//
// Created by muliana on 10/16/18.
//

#include "shiftmapping.h"

#ifndef MITSUBA_EXPLICIT_H
#define MITSUBA_EXPLICIT_H

MTS_NAMESPACE_BEGIN

class ShiftMappingExplicit : public ShiftMapping {
private:
    struct BaseState {
        Intersection its;
        Spectrum throughput;
        Float pdf;
        RayDifferential ray;
        const BSDF *bsdf;
        BSDFSampleResult sampled_bsdf;
    };
    struct BaseDirect {
        Point2 sample;
        DirectSamplingRecord dRec;
        bool visible;
        Spectrum emittedRadiance;
    };

    std::vector<BaseState> state = {};
    std::vector<BaseDirect> direct = {};
public:
    ShiftMappingExplicit(const GradientPathTracerConfig *config) :
            ShiftMapping(config) {}

    void evaluate(RayState &main, RayState *shiftedRays,
                  int secondaryCount, Spectrum &out_veryDirect) override {

        state.clear();
        direct.clear();

        // Nothing
        trace_base(main, out_veryDirect);
        Spectrum base_value = evalute_base(main);
        Float q = std::min(std::max(base_value.max(), m_config->rrShift), 1.0);
        if (q >= main.rRec.sampler->next1D()) {
            Float w = 1.0 / q;
            evaluate_gradients(main, shiftedRays, secondaryCount);
            for (int i = 0; i < secondaryCount; ++i) {
                RayState &shifted = shiftedRays[i];
                shifted.scaleValue(w);
            }
            main.scaleValue(w);
        }
    }

    void evaluateReuse(RayState *rays, int secondaryCount, Spectrum &out_veryDirect, int id_main,
                       std::vector<GradientInfo> &gradients) override {
        state.clear();
        direct.clear();

        trace_base(rays[id_main], out_veryDirect);
        Spectrum base_value = evalute_base(rays[id_main]);
        Float q = std::min(std::max(base_value.max(), m_config->rrShift), 1.0);
        if (q >= rays[id_main].rRec.sampler->next1D()) {
            Float w = 1.0 / q;
            // TODO: Clean this code
            evaluate_gradients_reuse(rays, secondaryCount, id_main, gradients, w);
            for (int i = 0; i < secondaryCount; ++i) {
                RayState &shifted = rays[i];
                shifted.scaleValue(w);
            }
        }
    }

    Spectrum evalute_base(RayState &mainTraced) {
        const Scene *scene = mainTraced.rRec.scene;

        Spectrum Li(0.f);
        int bounce = 1;
        while (bounce < m_config->m_maxDepth || m_config->m_maxDepth < 0) {
            if (state.size() < bounce) {
                return Li;
            }

            // Direct sampling
            {
                BaseState &mainState = state[bounce - 1];
                Intersection &mainIts = mainState.its;
                const BSDF *mainBSDF = mainState.bsdf;
                if (mainBSDF->getType() & BSDF::ESmooth) {
                    mitsuba::Point2 lightSample = mainTraced.rRec.nextSample2D();
                    DirectSamplingRecord dRec(mainState.its);

                    std::pair<Spectrum, bool> emitterTuple = scene->sampleEmitterDirectVisible(dRec, lightSample);
                    Spectrum mainEmitterRadiance = emitterTuple.first * dRec.pdf;
                    bool mainEmitterVisible = emitterTuple.second;

                    const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

                    // If the emitter sampler produces a non-emitter, that's a problem.
                    SAssert(emitter != nullptr);

                    // Add radiance and gradients to the base path and its offset path.
                    // Query the BSDF to the emitter's direction.
                    BSDFSamplingRecord mainBRec(mainIts, mainIts.toLocal(dRec.d), ERadiance);

                    // Evaluate BSDF * cos(theta).
                    Spectrum mainBSDFValue = mainBSDF->eval(mainBRec);

                    // Calculate the probability density of having generated the sampled path segment by BSDF sampling. Note that if the emitter is not visible, the probability density is zero.
                    // Even if the BSDF sampler has zero probability density, the light sampler can still sample it.
                    Float mainBsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle && mainEmitterVisible)
                                        ? mainBSDF->pdf(mainBRec) : 0;

                    // There values are probably needed soon for the Jacobians.
                    Float mainDistanceSquared = (mainIts.p - dRec.p).lengthSquared();
                    Float mainOpposingCosine = dot(dRec.n, (mainIts.p - dRec.p)) / sqrt(mainDistanceSquared);

                    // Power heuristic weights for the following strategies: light sample from base, BSDF sample from base.
                    Float mainWeightNumerator = mainState.pdf * dRec.pdf;
                    Float mainWeightDenominator =
                            (mainState.pdf * mainState.pdf) * ((dRec.pdf * dRec.pdf) + (mainBsdfPdf * mainBsdfPdf));

                    Float shiftedWeightDenominator = Float(0);
                    Float weight = mainWeightNumerator / (D_EPSILON + mainWeightDenominator);
                    Li += weight * mainState.throughput * (mainBSDFValue * mainEmitterRadiance);

                    direct.emplace_back(
                            BaseDirect{
                                    lightSample,
                                    dRec,
                                    mainEmitterVisible,
                                    mainEmitterRadiance,
                            }
                    );
                }
            }

            // BSDF Sampling
            if (state.size() > bounce) {
                // TODO: Change the function for getting the vertex type
                BaseState &prevState = state[bounce - 1];
                BaseState &nextState = state[bounce];

                // Reconstruct the MIS
                bool mainHitEmitter = false;
                Spectrum mainEmitterRadiance = Spectrum((Float) 0);
                DirectSamplingRecord mainDRec(prevState.its);
                const BSDF *mainBSDF = prevState.bsdf;

                // TODO: No env map shifting in this case
                // Intersected something - check if it was a luminaire.
                if (nextState.its.isEmitter()) {
                    mainEmitterRadiance = nextState.its.Le(-nextState.ray.d);

                    mainDRec.setQuery(nextState.ray, nextState.its);
                    mainHitEmitter = true;
                }

                // Sub-surface scattering.
                if (nextState.its.hasSubsurface()) {
                    mainEmitterRadiance += nextState.its.LoSub(scene, mainTraced.rRec.sampler, -nextState.ray.d,
                                                               bounce);
                }

                // Continue the shift.
                Float mainBsdfPdf = prevState.sampled_bsdf.pdf;
                Float mainPreviousPdf = prevState.pdf;

                // Compute the probability density of generating base path's direction using the implemented direct illumination sampling technique.
                const Float mainLumPdf = (mainHitEmitter && bounce + 1 >= m_config->m_minDepth &&
                                          !(prevState.sampled_bsdf.bRec.sampledType & BSDF::EDelta)) ?
                                         scene->pdfEmitterDirect(mainDRec) : 0;

                // Power heuristic weights for the following strategies: light sample from base, BSDF sample from base.
                Float mainWeightNumerator = mainPreviousPdf * prevState.sampled_bsdf.pdf;
                Float mainWeightDenominator =
                        (mainPreviousPdf * mainPreviousPdf) * ((mainLumPdf * mainLumPdf) + (mainBsdfPdf * mainBsdfPdf));


                Float weight = mainWeightNumerator / (D_EPSILON + mainWeightDenominator);
                Li += weight * nextState.throughput * mainEmitterRadiance;
            }

            if (bounce++ >= m_config->m_rrDepth) {
//                SLog(EError, "not implemented");
            }
        }

    }

    void evaluate_gradients_reuse(RayState *shiftedRays, int secondaryCount,
                                  int id_main, std::vector<GradientInfo> &gradients, Float rrShiftWeight) {
        const Scene *scene = shiftedRays[id_main].rRec.scene;

        // Perform the same first ray intersection for the offset paths.
        for (int i = 0; i < secondaryCount; ++i) {
            if (i == id_main) continue;
            RayState &shifted = shiftedRays[i];
            shifted.rRec.rayIntersect(shifted.ray);
            shifted.ray.mint = Epsilon;
        }

        // If no intersection of an offset ray could be found, its offset paths can not be generated.
        for (int i = 0; i < secondaryCount; ++i) {
            if (i == id_main) continue;
            RayState &shifted = shiftedRays[i];
            if (!shifted.rRec.its.isValid()) {
                shifted.alive = false;
            }
        }

        // Main path tracing loop.
        int bounce = 1;
        int smooth_surface = 1;
        while (bounce < m_config->m_maxDepth || m_config->m_maxDepth < 0) {
            if (state.size() < bounce) {
                return;
            }

            if (m_config->m_strictNormals) {
                for (int i = 0; i < secondaryCount; ++i) {
                    if (i == id_main) continue;
                    RayState &shifted = shiftedRays[i];
                    if (dot(shifted.ray.d, shifted.rRec.its.geoFrame.n) * Frame::cosTheta(shifted.rRec.its.wi) >= 0) {
                        // This is an impossible offset path.
                        shifted.alive = false;
                    }
                }
            }

            {
                BaseState &mainState = state[bounce - 1];
                Intersection &mainIts = mainState.its;
                const BSDF *mainBSDF = mainState.bsdf;
                if (mainBSDF->getType() & BSDF::ESmooth) {
                    const BaseDirect &baseDirect = direct[smooth_surface - 1];
                    smooth_surface += 1;
                    const Emitter *emitter = static_cast<const Emitter *>(baseDirect.dRec.object);
                    const DirectSamplingRecord &dRec = baseDirect.dRec;
                    SAssert(emitter != nullptr);

                    // Add radiance and gradients to the base path and its offset path.
                    // Query the BSDF to the emitter's direction.
                    BSDFSamplingRecord mainBRec(mainIts, mainIts.toLocal(dRec.d), ERadiance);

                    // Evaluate BSDF * cos(theta).
                    Spectrum mainBSDFValue = mainBSDF->eval(mainBRec);

                    // Calculate the probability density of having generated the sampled path segment by BSDF sampling. Note that if the emitter is not visible, the probability density is zero.
                    // Even if the BSDF sampler has zero probability density, the light sampler can still sample it.
                    Float mainBsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle && baseDirect.visible)
                                        ? mainBSDF->pdf(mainBRec) : 0;

                    // There values are probably needed soon for the Jacobians.
                    Float mainDistanceSquared = (mainIts.p - dRec.p).lengthSquared();
                    Float mainOpposingCosine = dot(dRec.n, (mainIts.p - dRec.p)) / sqrt(mainDistanceSquared);

                    // Power heuristic weights for the following strategies: light sample from base, BSDF sample from base.
                    Float weightNumerator = 1.0;
                    Float weightDenominator = mainState.pdf * (dRec.pdf + mainBsdfPdf);

                    std::vector<Spectrum> contribution(secondaryCount, Spectrum(0.f));
                    // Strict normals check to produce the same results as bidirectional methods when normal mapping is used.
                    if (!m_config->m_strictNormals ||
                        dot(mainIts.geoFrame.n, dRec.d) * Frame::cosTheta(mainBRec.wo) > 0) {
                        // The base path is good. Add radiance differences to offset paths.
                        for (int i = 0; i < secondaryCount; ++i) {
                            if (i == id_main) {
                                contribution[i] = mainState.throughput * (mainBSDFValue * baseDirect.emittedRadiance);
                                continue;
                            }

                            // Evaluate and apply the gradient.
                            RayState &shifted = shiftedRays[i];
                            Spectrum shiftedContribution(Float(0));
                            bool shiftSuccessful = shifted.alive;

                            // Construct the offset path.
                            if (shiftSuccessful) {
                                // Generate the offset path.
                                if (shifted.connection_status == RAY_CONNECTED) {
                                    // Follow the base path. All relevant vertices are shared.
                                    Float shiftedBsdfPdf = mainBsdfPdf;
                                    Float shiftedDRecPdf = dRec.pdf;
                                    Spectrum shiftedBsdfValue = mainBSDFValue;
                                    Spectrum shiftedEmitterRadiance = baseDirect.emittedRadiance;
                                    Float jacobian = (Float) 1;

                                    // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                                    weightDenominator += (jacobian * shifted.pdf) * (shiftedDRecPdf + shiftedBsdfPdf);
                                    shiftedContribution =
                                            jacobian * shifted.throughput * (shiftedBsdfValue * shiftedEmitterRadiance);

                                    // Note: The Jacobians were baked into shifted.pdf and shifted.throughput at connection phase.
                                } else if (shifted.connection_status == RAY_RECENTLY_CONNECTED) {
                                    // Follow the base path. The current vertex is shared, but the incoming directions differ.
                                    Vector3 incomingDirection = normalize(shifted.rRec.its.p - mainIts.p);

                                    BSDFSamplingRecord bRec(mainState.its, mainState.its.toLocal(incomingDirection),
                                                            mainState.its.toLocal(dRec.d), ERadiance);

                                    // Sample the BSDF.
                                    Float shiftedBsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle &&
                                                            baseDirect.visible) ? mainBSDF->pdf(bRec)
                                                                                : 0; // The BSDF sampler can not sample occluded path segments.
                                    Float shiftedDRecPdf = dRec.pdf;
                                    Spectrum shiftedBsdfValue = mainBSDF->eval(bRec);
                                    Spectrum shiftedEmitterRadiance = baseDirect.emittedRadiance;
                                    Float jacobian = (Float) 1;

                                    // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                                    weightDenominator += (jacobian * shifted.pdf) * (shiftedDRecPdf + shiftedBsdfPdf);
                                    shiftedContribution =
                                            jacobian * shifted.throughput * (shiftedBsdfValue * shiftedEmitterRadiance);

                                    // Note: The Jacobians were baked into shifted.pdf and shifted.throughput at connection phase.
                                } else {
                                    // Reconnect to the sampled light vertex. No shared vertices.
                                    SAssert(shifted.connection_status == RAY_NOT_CONNECTED);

                                    const BSDF *shiftedBSDF = shifted.rRec.its.getBSDF(shifted.ray);

                                    // This implementation uses light sampling only for the reconnect-shift.
                                    // When one of the BSDFs is very glossy, light sampling essentially reduces to a failed shift anyway.
                                    bool mainAtPointLight = (dRec.measure == EDiscrete);
                                    VertexType mainVertexType = getVertexType(mainBSDF, mainIts, *m_config,
                                                                              BSDF::ESmooth);
                                    VertexType shiftedVertexType = getVertexType(shifted, *m_config, BSDF::ESmooth);

                                    if (mainAtPointLight || (mainVertexType == VERTEX_TYPE_DIFFUSE &&
                                                             shiftedVertexType == VERTEX_TYPE_DIFFUSE)) {
                                        // Get emitter radiance.
                                        DirectSamplingRecord shiftedDRec(shifted.rRec.its);
                                        std::pair<Spectrum, bool> emitterTuple = scene->sampleEmitterDirectVisible(
                                                shiftedDRec, baseDirect.sample);
                                        bool shiftedEmitterVisible = emitterTuple.second;

                                        Spectrum shiftedEmitterRadiance = emitterTuple.first * shiftedDRec.pdf;
                                        Float shiftedDRecPdf = shiftedDRec.pdf;

                                        // Sample the BSDF.
                                        Float shiftedDistanceSquared = (dRec.p - shifted.rRec.its.p).lengthSquared();
                                        Vector emitterDirection =
                                                (dRec.p - shifted.rRec.its.p) / sqrt(shiftedDistanceSquared);
                                        Float shiftedOpposingCosine = -dot(dRec.n, emitterDirection);

                                        BSDFSamplingRecord bRec(shifted.rRec.its,
                                                                shifted.rRec.its.toLocal(emitterDirection), ERadiance);

                                        // Strict normals check, to make the output match with bidirectional methods when normal maps are present.
                                        if (m_config->m_strictNormals &&
                                            dot(shifted.rRec.its.geoFrame.n, emitterDirection) *
                                            Frame::cosTheta(bRec.wo) < 0) {
                                            // Invalid, non-samplable offset path.
                                            shiftSuccessful = false;
                                        } else {
                                            Spectrum shiftedBsdfValue = shiftedBSDF->eval(bRec);
                                            Float shiftedBsdfPdf = (emitter->isOnSurface() &&
                                                                    dRec.measure == ESolidAngle &&
                                                                    shiftedEmitterVisible) ? shiftedBSDF->pdf(bRec) : 0;
                                            Float jacobian = std::abs(shiftedOpposingCosine * mainDistanceSquared) /
                                                             (Epsilon +
                                                              std::abs(mainOpposingCosine * shiftedDistanceSquared));

                                            // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                                            weightDenominator +=
                                                    (jacobian * shifted.pdf) * (shiftedDRecPdf + shiftedBsdfPdf);
                                            shiftedContribution = jacobian * shifted.throughput *
                                                                  (shiftedBsdfValue * shiftedEmitterRadiance);
                                        }
                                    }
                                }
                            }

                            contribution[i] = shiftedContribution;

                        } // for(int i = 0; i < secondaryCount; ++i)
                    } // Strict normals

                    if (weightDenominator != 0) {
                        Float weight = weightNumerator / weightDenominator;
                        for (int i = 0; i < secondaryCount; ++i) {
                            shiftedRays[i].addRadiance(contribution[i], weight);
                        }
                        for (auto &g: gradients) {
                            g.value += (contribution[g.pixel_a_index] - contribution[g.pixel_b_index]) * weight *
                                       rrShiftWeight;
                        }
                    }
                }
            } // Sample incoming radiance from lights.


            if (state.size() > bounce) {
                // TODO: Change the function for getting the vertex type
                BaseState &prevState = state[bounce - 1];
                BaseState &nextState = state[bounce];

                // Reconstruct the MIS
                bool mainHitEmitter = false;
                Spectrum mainEmitterRadiance = Spectrum((Float) 0);
                DirectSamplingRecord mainDRec(prevState.its);
                const BSDF *mainBSDF = prevState.bsdf;


                // Update the vertex types.
                VertexType mainVertexType = getVertexType(mainBSDF, prevState.its, *m_config,
                                                          prevState.sampled_bsdf.bRec.sampledType);
                VertexType mainNextVertexType;


                // TODO: No env map shifting in this case
                // Intersected something - check if it was a luminaire.
                if (nextState.its.isEmitter()) {
                    mainEmitterRadiance = nextState.its.Le(-nextState.ray.d);

                    mainDRec.setQuery(nextState.ray, nextState.its);
                    mainHitEmitter = true;
                }

                // Sub-surface scattering.
                if (nextState.its.hasSubsurface()) {
                    mainEmitterRadiance += nextState.its.LoSub(scene,
                                                               shiftedRays[id_main].rRec.sampler,
                                                               -nextState.ray.d,
                                                               bounce);
                }

                // Update the vertex type.
                // FIXME: This part is wrong I guess
                mainNextVertexType = getVertexType(nextState.bsdf, nextState.its, *m_config,
                                                   prevState.sampled_bsdf.bRec.sampledType);


                // Continue the shift.
                Float mainBsdfPdf = prevState.sampled_bsdf.pdf;
                Float mainPreviousPdf = prevState.pdf;

                // Compute the probability density of generating base path's direction using the implemented direct illumination sampling technique.
                const Float mainLumPdf = (mainHitEmitter && bounce + 1 >= m_config->m_minDepth &&
                                          !(prevState.sampled_bsdf.bRec.sampledType & BSDF::EDelta)) ?
                                         scene->pdfEmitterDirect(mainDRec) : 0;

                // Power heuristic weights for the following strategies: light sample from base, BSDF sample from base.
                Float weightNumerator = 1.0;
                Float weightDenominator = mainPreviousPdf * (mainLumPdf + mainBsdfPdf);

                std::vector<Spectrum> contribution(secondaryCount, Spectrum(0.f));
                for (int i = 0; i < secondaryCount; ++i) {
                    if (i == id_main) {
                        contribution[i] = nextState.throughput * mainEmitterRadiance;
                        continue;
                    }
                    RayState &shifted = shiftedRays[i];
                    Spectrum shiftedContribution(Float(0));
                    Float shiftDenominatorWeight(0);

                    bool postponedShiftEnd = false; // Kills the shift after evaluating the current radiance.
                    if (shifted.alive) {
                        // The offset path is still good, so it makes sense to continue its construction.
                        Float shiftedPreviousPdf = shifted.pdf;

                        if (shifted.connection_status == RAY_CONNECTED) {
                            // The offset path keeps following the base path.
                            // As all relevant vertices are shared, we can just reuse the sampled values.
                            Spectrum shiftedBsdfValue = prevState.sampled_bsdf.weight * prevState.sampled_bsdf.pdf;
                            Float shiftedBsdfPdf = mainBsdfPdf;
                            Float shiftedLumPdf = mainLumPdf;
                            Spectrum shiftedEmitterRadiance = mainEmitterRadiance;

                            // Update throughput and pdf.
                            shifted.throughput *= shiftedBsdfValue;
                            shifted.pdf *= shiftedBsdfPdf;
                            shifted.eta *= prevState.sampled_bsdf.bRec.eta;

                            // Power heuristic between light sample from base, BSDF sample from base,
                            // light sample from offset, BSDF sample from offset.
                            shiftDenominatorWeight = shiftedPreviousPdf * (shiftedLumPdf + shiftedBsdfPdf);
                            shiftedContribution =
                                    shifted.throughput *
                                    shiftedEmitterRadiance; // Note: Jacobian baked into .throughput.
                        } else if (shifted.connection_status == RAY_RECENTLY_CONNECTED) {
                            // Recently connected - follow the base path but evaluate BSDF to the new direction.
                            Vector3 incomingDirection = normalize(shifted.rRec.its.p - nextState.ray.o);
                            BSDFSamplingRecord bRec(prevState.its, prevState.its.toLocal(incomingDirection),
                                                    prevState.its.toLocal(nextState.ray.d), ERadiance);

                            // Note: mainBSDF is the BSDF at previousMainIts, which is the current position of the offset path.
                            EMeasure measure = (prevState.sampled_bsdf.bRec.sampledType & BSDF::EDelta) ? EDiscrete
                                                                                                        : ESolidAngle;

                            Spectrum shiftedBsdfValue = mainBSDF->eval(bRec, measure);
                            Float shiftedBsdfPdf = mainBSDF->pdf(bRec, measure);

                            Float shiftedLumPdf = mainLumPdf;
                            Spectrum shiftedEmitterRadiance = mainEmitterRadiance;

                            // Update throughput and pdf.
                            shifted.throughput *= shiftedBsdfValue;
                            shifted.pdf *= shiftedBsdfPdf;
                            shifted.eta *= prevState.sampled_bsdf.bRec.eta;

                            shifted.connection_status = RAY_CONNECTED;

                            // Power heuristic between light sample from base, BSDF sample from base,
                            // light sample from offset, BSDF sample from offset.
                            shiftDenominatorWeight = shiftedPreviousPdf * (shiftedLumPdf + shiftedBsdfPdf);
                            shiftedContribution =
                                    shifted.throughput *
                                    shiftedEmitterRadiance; // Note: Jacobian baked into .throughput.
                        } else {
                            // Not connected - apply either reconnection or half-vector duplication shift.

                            const BSDF *shiftedBSDF = shifted.rRec.its.getBSDF(shifted.ray);

                            // Update the vertex type of the offset path.
                            VertexType shiftedVertexType = getVertexType(shifted, *m_config,
                                                                         prevState.sampled_bsdf.bRec.sampledType);

                            if (mainVertexType == VERTEX_TYPE_DIFFUSE && mainNextVertexType == VERTEX_TYPE_DIFFUSE &&
                                shiftedVertexType == VERTEX_TYPE_DIFFUSE) {
                                // Use reconnection shift.
                                bool lastSegment = (bounce + 1 == m_config->m_maxDepth);
                                // Optimization: Skip the last raycast and BSDF evaluation for the offset path when it won't contribute and isn't needed anymore.
                                if (!lastSegment || mainHitEmitter || nextState.its.hasSubsurface()) {
                                    ReconnectionShiftResult shiftResult;
                                    // FIXME: Environment map
//                                    bool environmentConnection = false;

                                    if (nextState.its.isValid()) {
                                        // This is an actual reconnection shift.
                                        shiftResult = reconnectShift(scene, nextState.ray.o, nextState.its.p,
                                                                     shifted.rRec.its.p, nextState.its.geoFrame.n,
                                                                     nextState.ray.time);
                                    } else {
                                        SLog(EError, "No support of environement map for now...");
                                    }

                                    if (!shiftResult.success) {
                                        // Failed to construct the offset path.
                                        shifted.alive = false;
                                        goto shift_failed;
                                    }

                                    Vector3 incomingDirection = -shifted.ray.d;
                                    Vector3 outgoingDirection = shiftResult.wo;

                                    BSDFSamplingRecord bRec(shifted.rRec.its,
                                                            shifted.rRec.its.toLocal(incomingDirection),
                                                            shifted.rRec.its.toLocal(outgoingDirection), ERadiance);

                                    // Strict normals check.
                                    if (m_config->m_strictNormals &&
                                        dot(outgoingDirection, shifted.rRec.its.geoFrame.n) *
                                        Frame::cosTheta(bRec.wo) <=
                                        0) {
                                        shifted.alive = false;
                                        goto shift_failed;
                                    }

                                    // Evaluate the BRDF to the new direction.
                                    Spectrum shiftedBsdfValue = shiftedBSDF->eval(bRec);
                                    Float shiftedBsdfPdf = shiftedBSDF->pdf(bRec);

                                    // Update throughput and pdf.
                                    shifted.throughput *= shiftedBsdfValue * shiftResult.jacobian;
                                    shifted.pdf *= shiftedBsdfPdf * shiftResult.jacobian;
                                    //shifted.eta *= bRec.eta;
                                    // ignore eta here

                                    shifted.connection_status = RAY_RECENTLY_CONNECTED;

                                    if (mainHitEmitter || nextState.its.hasSubsurface()) {
                                        // Also the offset path hit the emitter, as visibility was checked at reconnectShift or environmentShift.

                                        // Evaluate radiance to this direction.
                                        Spectrum shiftedEmitterRadiance(Float(0));
                                        Float shiftedLumPdf = Float(0);

                                        if (nextState.its.isValid()) {
                                            // Hit an object.
                                            if (mainHitEmitter) {
                                                shiftedEmitterRadiance = nextState.its.Le(-outgoingDirection);

                                                // Evaluate the light sampling PDF of the new segment.
                                                DirectSamplingRecord shiftedDRec;
                                                shiftedDRec.p = mainDRec.p;
                                                shiftedDRec.n = mainDRec.n;
                                                shiftedDRec.dist = (mainDRec.p - shifted.rRec.its.p).length();
                                                shiftedDRec.d = (mainDRec.p - shifted.rRec.its.p) / shiftedDRec.dist;
                                                shiftedDRec.ref = mainDRec.ref;
                                                shiftedDRec.refN = shifted.rRec.its.shFrame.n;
                                                shiftedDRec.object = mainDRec.object;

                                                shiftedLumPdf = scene->pdfEmitterDirect(shiftedDRec);
                                            }

                                            // Sub-surface scattering. Note: Should use the same random numbers as the base path!
                                            if (nextState.its.hasSubsurface()) {
                                                shiftedEmitterRadiance += nextState.its.LoSub(scene,
                                                                                              shifted.rRec.sampler,
                                                                                              -outgoingDirection,
                                                                                              bounce);
                                            }
                                        } else {
                                            // Hit the environment.
                                            shiftedEmitterRadiance = mainEmitterRadiance;
                                            shiftedLumPdf = mainLumPdf;
                                        }

                                        // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                                        shiftDenominatorWeight = shiftedPreviousPdf * (shiftedLumPdf + shiftedBsdfPdf);
                                        shiftedContribution = shifted.throughput *
                                                              shiftedEmitterRadiance; // Note: Jacobian baked into .throughput.
                                    }
                                }
                            } else {

                                // Use half-vector duplication shift. These paths could not have been sampled by light sampling (by our decision).
                                Vector3 tangentSpaceIncomingDirection = shifted.rRec.its.toLocal(-shifted.ray.d);
                                Vector3 tangentSpaceOutgoingDirection;
                                Spectrum shiftedEmitterRadiance(Float(0));

                                const BSDF *shiftedBSDF = shifted.rRec.its.getBSDF(shifted.ray);

                                HalfVectorShiftResult shiftResult;
                                EMeasure measure;
                                BSDFSamplingRecord bRec(shifted.rRec.its, tangentSpaceIncomingDirection,
                                                        tangentSpaceOutgoingDirection, ERadiance);
                                Vector3 outgoingDirection;
                                VertexType shiftedVertexType;

                                // Deny shifts between Dirac and non-Dirac BSDFs.
                                bool bothDelta = (prevState.sampled_bsdf.bRec.sampledType & BSDF::EDelta) &&
                                                 (shiftedBSDF->getType() & BSDF::EDelta);
                                bool bothSmooth = (prevState.sampled_bsdf.bRec.sampledType & BSDF::ESmooth) &&
                                                  (shiftedBSDF->getType() & BSDF::ESmooth);
                                if (!(bothDelta || bothSmooth)) {
                                    shifted.alive = false;
                                    goto half_vector_shift_failed;
                                }

                                SAssert(fabs(shifted.ray.d.lengthSquared() - 1) < 0.01);

                                // Apply the local shift.
                                shiftResult = halfVectorShift(prevState.sampled_bsdf.bRec.wi,
                                                              prevState.sampled_bsdf.bRec.wo,
                                                              shifted.rRec.its.toLocal(-shifted.ray.d),
                                                              mainBSDF->getEta(),
                                                              shiftedBSDF->getEta());
                                bRec.wo = shiftResult.wo;

                                if (prevState.sampled_bsdf.bRec.sampledType & BSDF::EDelta) {
                                    // Dirac delta integral is a point evaluation - no Jacobian determinant!
                                    shiftResult.jacobian = Float(1);
                                }

                                if (shiftResult.success) {
                                    // Invertible shift, success.
                                    shifted.throughput *= shiftResult.jacobian;
                                    shifted.pdf *= shiftResult.jacobian;
                                    tangentSpaceOutgoingDirection = shiftResult.wo;
                                } else {
                                    // The shift is non-invertible so kill it.
                                    shifted.alive = false;
                                    goto half_vector_shift_failed;
                                }


                                outgoingDirection = shifted.rRec.its.toWorld(tangentSpaceOutgoingDirection);

                                // Update throughput and pdf.
                                measure = (prevState.sampled_bsdf.bRec.sampledType & BSDF::EDelta) ? EDiscrete
                                                                                                   : ESolidAngle;

                                shifted.throughput *= shiftedBSDF->eval(bRec, measure);
                                shifted.pdf *= shiftedBSDF->pdf(bRec, measure);

                                if (shifted.pdf == Float(0)) {
                                    // Offset path is invalid!
                                    shifted.alive = false;
                                    goto half_vector_shift_failed;
                                }

                                // Strict normals check to produce the same results as bidirectional methods when normal mapping is used.
                                if (m_config->m_strictNormals &&
                                    dot(outgoingDirection, shifted.rRec.its.geoFrame.n) * Frame::cosTheta(bRec.wo) <=
                                    0) {
                                    shifted.alive = false;
                                    goto half_vector_shift_failed;
                                }


                                // Update the vertex type.
                                shiftedVertexType = getVertexType(shifted, *m_config,
                                                                  prevState.sampled_bsdf.bRec.sampledType);

                                // Trace the next hit point.
                                shifted.ray = Ray(shifted.rRec.its.p, outgoingDirection, nextState.ray.time);

                                if (!scene->rayIntersect(shifted.ray, shifted.rRec.its)) {
                                    // Hit nothing - Evaluate environment radiance.
                                    const Emitter *env = scene->getEnvironmentEmitter();
                                    if (!env) {
                                        // Since base paths that hit nothing are not shifted, we must be symmetric and kill shifts that hit nothing.
                                        shifted.alive = false;
                                        goto half_vector_shift_failed;
                                    }
                                    if (nextState.its.isValid()) {
                                        // Deny shifts between env and non-env.
                                        shifted.alive = false;
                                        goto half_vector_shift_failed;
                                    }

                                    if (mainVertexType == VERTEX_TYPE_DIFFUSE &&
                                        shiftedVertexType == VERTEX_TYPE_DIFFUSE) {
                                        // Environment reconnection shift would have been used for the reverse direction!
                                        shifted.alive = false;
                                        goto half_vector_shift_failed;
                                    }

                                    // The offset path is no longer valid after this path segment.
                                    shiftedEmitterRadiance = env->evalEnvironment(shifted.ray);
                                    postponedShiftEnd = true;
                                } else {
                                    // Hit something.
                                    if (!nextState.its.isValid()) {
                                        // Deny shifts between env and non-env.
                                        shifted.alive = false;
                                        goto half_vector_shift_failed;
                                    }

                                    // FIXME: Can be wrong as well
                                    VertexType shiftedNextVertexType = getVertexType(shifted, *m_config,
                                                                                     prevState.sampled_bsdf.bRec.sampledType);

                                    // Make sure that the reverse shift would use this same strategy!
                                    // ==============================================================

                                    if (mainVertexType == VERTEX_TYPE_DIFFUSE &&
                                        shiftedVertexType == VERTEX_TYPE_DIFFUSE &&
                                        shiftedNextVertexType == VERTEX_TYPE_DIFFUSE) {
                                        // Non-invertible shift: the reverse-shift would use another strategy!
                                        shifted.alive = false;
                                        goto half_vector_shift_failed;
                                    }

                                    if (shifted.rRec.its.isEmitter()) {
                                        // Hit emitter.
                                        shiftedEmitterRadiance = shifted.rRec.its.Le(-shifted.ray.d);
                                    }
                                    // Sub-surface scattering. Note: Should use the same random numbers as the base path!
                                    if (shifted.rRec.its.hasSubsurface()) {
                                        shiftedEmitterRadiance += shifted.rRec.its.LoSub(scene, shifted.rRec.sampler,
                                                                                         -shifted.ray.d,
                                                                                         bounce);
                                    }

                                }


                                half_vector_shift_failed:
                                if (shifted.alive) {
                                    // Evaluate radiance difference using power heuristic between BSDF samples from base and offset paths.
                                    // Note: No MIS with light sampling since we don't use it for this connection type.
                                    shiftDenominatorWeight = shifted.pdf;
                                    shiftedContribution = shifted.throughput *
                                                          shiftedEmitterRadiance; // Note: Jacobian baked into .throughput.
                                } else {
                                    // Handle the failure without taking MIS with light sampling, as we decided not to use it in the half-vector-duplication case.
                                    // Could have used it, but so far there has been no need. It doesn't seem to be very useful.
                                    shiftDenominatorWeight = 0.0;
                                    shiftedContribution = Spectrum(Float(0));

                                    // Disable the failure detection below since the failure was already handled.
                                    shifted.alive = true;
                                    postponedShiftEnd = true;

                                    // (TODO: Restructure into smaller functions and get rid of the gotos... Although this may mean having lots of small functions with a large number of parameters.)
                                }
                            }
                        }
                    }

                    shift_failed:
                    if (!shifted.alive) {
                        // The offset path cannot be generated; Set offset PDF and offset throughput to zero.
                        shiftDenominatorWeight = 0.0;
                        shiftedContribution = Spectrum(Float(0));
                    }

                    contribution[i] = shiftedContribution;
                    weightDenominator += shiftDenominatorWeight;

                    if (postponedShiftEnd) {
                        shifted.alive = false;
                    }
                }

                if (weightDenominator != 0) {
                    Float weight = weightNumerator / weightDenominator;
                    for (int i = 0; i < secondaryCount; ++i) {
                        shiftedRays[i].addRadiance(contribution[i], weight);
                    }
                    for (auto &g: gradients) {
                        g.value += (contribution[g.pixel_a_index] - contribution[g.pixel_b_index]) * weight *
                                   rrShiftWeight;
                    }
                }
            }

            if (bounce++ >= m_config->m_rrDepth) {
                SLog(EError, "not implemented");
            }
        }
    }

    void evaluate_gradients(RayState &mainTraced, RayState *shiftedRays, int secondaryCount) {
        const Scene *scene = mainTraced.rRec.scene;

        // Perform the same first ray intersection for the offset paths.
        for (int i = 0; i < secondaryCount; ++i) {
            RayState &shifted = shiftedRays[i];
            shifted.rRec.rayIntersect(shifted.ray);
            shifted.ray.mint = Epsilon;
        }

        // If no intersection of an offset ray could be found, its offset paths can not be generated.
        for (int i = 0; i < secondaryCount; ++i) {
            RayState &shifted = shiftedRays[i];
            if (!shifted.rRec.its.isValid()) {
                shifted.alive = false;
            }
        }

        // Main path tracing loop.
        size_t bounce = 1;
        int smooth_surface = 1;
        while (bounce < m_config->m_maxDepth || m_config->m_maxDepth < 0) {
            if (state.size() < bounce) {
                return;
            }

            if (m_config->m_strictNormals) {
                for (int i = 0; i < secondaryCount; ++i) {
                    RayState &shifted = shiftedRays[i];

                    if (dot(shifted.ray.d, shifted.rRec.its.geoFrame.n) * Frame::cosTheta(shifted.rRec.its.wi) >= 0) {
                        // This is an impossible offset path.
                        shifted.alive = false;
                    }
                }
            }

            {
                BaseState &mainState = state[bounce - 1];
                Intersection &mainIts = mainState.its;
                const BSDF *mainBSDF = mainState.bsdf;
                if (mainBSDF->getType() & BSDF::ESmooth) {
                    const BaseDirect &baseDirect = direct[smooth_surface - 1];
                    smooth_surface += 1;
                    const Emitter *emitter = static_cast<const Emitter *>(baseDirect.dRec.object);
                    const DirectSamplingRecord &dRec = baseDirect.dRec;
                    SAssert(emitter != nullptr);

                    // Add radiance and gradients to the base path and its offset path.
                    // Query the BSDF to the emitter's direction.
                    BSDFSamplingRecord mainBRec(mainIts, mainIts.toLocal(dRec.d), ERadiance);

                    // Evaluate BSDF * cos(theta).
                    Spectrum mainBSDFValue = mainBSDF->eval(mainBRec);

                    // Calculate the probability density of having generated the sampled path segment by BSDF sampling. Note that if the emitter is not visible, the probability density is zero.
                    // Even if the BSDF sampler has zero probability density, the light sampler can still sample it.
                    Float mainBsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle && baseDirect.visible)
                                        ? mainBSDF->pdf(mainBRec) : 0;

                    // There values are probably needed soon for the Jacobians.
                    Float mainDistanceSquared = (mainIts.p - dRec.p).lengthSquared();
                    Float mainOpposingCosine = dot(dRec.n, (mainIts.p - dRec.p)) / sqrt(mainDistanceSquared);

                    // Power heuristic weights for the following strategies: light sample from base, BSDF sample from base.
                    Float mainWeightNumerator = mainState.pdf * dRec.pdf;
                    Float mainWeightDenominator =
                            (mainState.pdf * mainState.pdf) * ((dRec.pdf * dRec.pdf) + (mainBsdfPdf * mainBsdfPdf));

                    // Strict normals check to produce the same results as bidirectional methods when normal mapping is used.
                    if (!m_config->m_strictNormals ||
                        dot(mainIts.geoFrame.n, dRec.d) * Frame::cosTheta(mainBRec.wo) > 0) {
                        // The base path is good. Add radiance differences to offset paths.
                        for (int i = 0; i < secondaryCount; ++i) {
                            // Evaluate and apply the gradient.
                            RayState &shifted = shiftedRays[i];

                            Spectrum mainContribution(Float(0));
                            Spectrum shiftedContribution(Float(0));
                            Float weight = Float(0);

                            bool shiftSuccessful = shifted.alive;

                            // Construct the offset path.
                            if (shiftSuccessful) {
                                // Generate the offset path.
                                if (shifted.connection_status == RAY_CONNECTED) {
                                    // Follow the base path. All relevant vertices are shared.
                                    Float shiftedBsdfPdf = mainBsdfPdf;
                                    Float shiftedDRecPdf = dRec.pdf;
                                    Spectrum shiftedBsdfValue = mainBSDFValue;
                                    Spectrum shiftedEmitterRadiance = baseDirect.emittedRadiance;
                                    Float jacobian = (Float) 1;

                                    // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                                    Float shiftedWeightDenominator =
                                            (jacobian * shifted.pdf) * (jacobian * shifted.pdf) *
                                            ((shiftedDRecPdf * shiftedDRecPdf) + (shiftedBsdfPdf * shiftedBsdfPdf));
                                    weight = mainWeightNumerator /
                                             (D_EPSILON + shiftedWeightDenominator + mainWeightDenominator);

                                    mainContribution =
                                            mainState.throughput * (mainBSDFValue * baseDirect.emittedRadiance);
                                    shiftedContribution =
                                            jacobian * shifted.throughput * (shiftedBsdfValue * shiftedEmitterRadiance);

                                    // Note: The Jacobians were baked into shifted.pdf and shifted.throughput at connection phase.
                                } else if (shifted.connection_status == RAY_RECENTLY_CONNECTED) {
                                    // Follow the base path. The current vertex is shared, but the incoming directions differ.
                                    Vector3 incomingDirection = normalize(shifted.rRec.its.p - mainIts.p);

                                    BSDFSamplingRecord bRec(mainState.its, mainState.its.toLocal(incomingDirection),
                                                            mainState.its.toLocal(dRec.d), ERadiance);

                                    // Sample the BSDF.
                                    Float shiftedBsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle &&
                                                            baseDirect.visible) ? mainBSDF->pdf(bRec)
                                                                                : 0; // The BSDF sampler can not sample occluded path segments.
                                    Float shiftedDRecPdf = dRec.pdf;
                                    Spectrum shiftedBsdfValue = mainBSDF->eval(bRec);
                                    Spectrum shiftedEmitterRadiance = baseDirect.emittedRadiance;
                                    Float jacobian = (Float) 1;

                                    // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                                    Float shiftedWeightDenominator =
                                            (jacobian * shifted.pdf) * (jacobian * shifted.pdf) *
                                            ((shiftedDRecPdf * shiftedDRecPdf) + (shiftedBsdfPdf * shiftedBsdfPdf));
                                    weight = mainWeightNumerator /
                                             (D_EPSILON + shiftedWeightDenominator + mainWeightDenominator);

                                    mainContribution =
                                            mainState.throughput * (mainBSDFValue * baseDirect.emittedRadiance);
                                    shiftedContribution =
                                            jacobian * shifted.throughput * (shiftedBsdfValue * shiftedEmitterRadiance);

                                    // Note: The Jacobians were baked into shifted.pdf and shifted.throughput at connection phase.
                                } else {
                                    // Reconnect to the sampled light vertex. No shared vertices.
                                    SAssert(shifted.connection_status == RAY_NOT_CONNECTED);

                                    const BSDF *shiftedBSDF = shifted.rRec.its.getBSDF(shifted.ray);

                                    // This implementation uses light sampling only for the reconnect-shift.
                                    // When one of the BSDFs is very glossy, light sampling essentially reduces to a failed shift anyway.
                                    bool mainAtPointLight = (dRec.measure == EDiscrete);
                                    VertexType mainVertexType = getVertexType(mainBSDF, mainIts, *m_config,
                                                                              BSDF::ESmooth);
                                    VertexType shiftedVertexType = getVertexType(shifted, *m_config, BSDF::ESmooth);

                                    if (mainAtPointLight || (mainVertexType == VERTEX_TYPE_DIFFUSE &&
                                                             shiftedVertexType == VERTEX_TYPE_DIFFUSE)) {
                                        // Get emitter radiance.
                                        DirectSamplingRecord shiftedDRec(shifted.rRec.its);
                                        std::pair<Spectrum, bool> emitterTuple = scene->sampleEmitterDirectVisible(
                                                shiftedDRec, baseDirect.sample);
                                        bool shiftedEmitterVisible = emitterTuple.second;

                                        Spectrum shiftedEmitterRadiance = emitterTuple.first * shiftedDRec.pdf;
                                        Float shiftedDRecPdf = shiftedDRec.pdf;

                                        // Sample the BSDF.
                                        Float shiftedDistanceSquared = (dRec.p - shifted.rRec.its.p).lengthSquared();
                                        Vector emitterDirection =
                                                (dRec.p - shifted.rRec.its.p) / sqrt(shiftedDistanceSquared);
                                        Float shiftedOpposingCosine = -dot(dRec.n, emitterDirection);

                                        BSDFSamplingRecord bRec(shifted.rRec.its,
                                                                shifted.rRec.its.toLocal(emitterDirection), ERadiance);

                                        // Strict normals check, to make the output match with bidirectional methods when normal maps are present.
                                        if (m_config->m_strictNormals &&
                                            dot(shifted.rRec.its.geoFrame.n, emitterDirection) *
                                            Frame::cosTheta(bRec.wo) < 0) {
                                            // Invalid, non-samplable offset path.
                                            shiftSuccessful = false;
                                        } else {
                                            Spectrum shiftedBsdfValue = shiftedBSDF->eval(bRec);
                                            Float shiftedBsdfPdf = (emitter->isOnSurface() &&
                                                                    dRec.measure == ESolidAngle &&
                                                                    shiftedEmitterVisible) ? shiftedBSDF->pdf(bRec) : 0;
                                            Float jacobian = std::abs(shiftedOpposingCosine * mainDistanceSquared) /
                                                             (Epsilon +
                                                              std::abs(mainOpposingCosine * shiftedDistanceSquared));

                                            // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                                            Float shiftedWeightDenominator =
                                                    (jacobian * shifted.pdf) * (jacobian * shifted.pdf) *
                                                    ((shiftedDRecPdf * shiftedDRecPdf) +
                                                     (shiftedBsdfPdf * shiftedBsdfPdf));
                                            weight = mainWeightNumerator /
                                                     (D_EPSILON + shiftedWeightDenominator + mainWeightDenominator);

                                            mainContribution =
                                                    mainState.throughput * (mainBSDFValue * baseDirect.emittedRadiance);
                                            shiftedContribution = jacobian * shifted.throughput *
                                                                  (shiftedBsdfValue * shiftedEmitterRadiance);
                                        }
                                    }
                                }
                            }

                            if (!shiftSuccessful) {
                                // The offset path cannot be generated; Set offset PDF and offset throughput to zero. This is what remains.

                                // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset. (Offset path has zero PDF)
                                Float shiftedWeightDenominator = Float(0);
                                weight = mainWeightNumerator / (D_EPSILON + mainWeightDenominator);

                                mainContribution = mainState.throughput * (mainBSDFValue * baseDirect.emittedRadiance);
                                shiftedContribution = Spectrum((Float) 0);
                            }

                            // Note: Using also the offset paths for the throughput estimate, like we do here, provides some advantage when a large reconstruction alpha is used,
                            // but using only throughputs of the base paths doesn't usually lose by much.
                            if(m_config->reusePrimal) {
                              mainTraced.addRadiance(mainContribution, weight);
                              shifted.addRadiance(shiftedContribution, weight);
                            }
                            shifted.addGradient(shiftedContribution - mainContribution, weight);
                        } // for(int i = 0; i < secondaryCount; ++i)
                    } // Strict normals
                }
            } // Sample incoming radiance from lights.


            if (state.size() > bounce) {
                // TODO: Change the function for getting the vertex type
                BaseState &prevState = state[bounce - 1];
                BaseState &nextState = state[bounce];

                // Reconstruct the MIS
                bool mainHitEmitter = false;
                Spectrum mainEmitterRadiance = Spectrum((Float) 0);
                DirectSamplingRecord mainDRec(prevState.its);
                const BSDF *mainBSDF = prevState.bsdf;


                // Update the vertex types.
                VertexType mainVertexType = getVertexType(mainBSDF, prevState.its, *m_config,
                                                          prevState.sampled_bsdf.bRec.sampledType);
                VertexType mainNextVertexType;


                // TODO: No env map shifting in this case
                // Intersected something - check if it was a luminaire.
                if (nextState.its.isEmitter()) {
                    mainEmitterRadiance = nextState.its.Le(-nextState.ray.d);

                    mainDRec.setQuery(nextState.ray, nextState.its);
                    mainHitEmitter = true;
                }

                // Sub-surface scattering.
                if (nextState.its.hasSubsurface()) {
                    mainEmitterRadiance += nextState.its.LoSub(scene, mainTraced.rRec.sampler, -nextState.ray.d,
                                                               bounce);
                }

                // Update the vertex type.
                // FIXME: This part is wrong I guess
                mainNextVertexType = getVertexType(nextState.bsdf, nextState.its, *m_config,
                                                   prevState.sampled_bsdf.bRec.sampledType);


                // Continue the shift.
                Float mainBsdfPdf = prevState.sampled_bsdf.pdf;
                Float mainPreviousPdf = prevState.pdf;

                // Compute the probability density of generating base path's direction using the implemented direct illumination sampling technique.
                const Float mainLumPdf = (mainHitEmitter && bounce + 1 >= m_config->m_minDepth &&
                                          !(prevState.sampled_bsdf.bRec.sampledType & BSDF::EDelta)) ?
                                         scene->pdfEmitterDirect(mainDRec) : 0;

                // Power heuristic weights for the following strategies: light sample from base, BSDF sample from base.
                Float mainWeightNumerator = mainPreviousPdf * prevState.sampled_bsdf.pdf;
                Float mainWeightDenominator =
                        (mainPreviousPdf * mainPreviousPdf) * ((mainLumPdf * mainLumPdf) + (mainBsdfPdf * mainBsdfPdf));

                for (int i = 0; i < secondaryCount; ++i) {
                    RayState &shifted = shiftedRays[i];
                    Spectrum mainContribution(Float(0));
                    Spectrum shiftedContribution(Float(0));
                    Float weight(0);

                    bool postponedShiftEnd = false; // Kills the shift after evaluating the current radiance.
                    if (shifted.alive) {
                        // The offset path is still good, so it makes sense to continue its construction.
                        Float shiftedPreviousPdf = shifted.pdf;

                        if (shifted.connection_status == RAY_CONNECTED) {
                            // The offset path keeps following the base path.
                            // As all relevant vertices are shared, we can just reuse the sampled values.
                            Spectrum shiftedBsdfValue = prevState.sampled_bsdf.weight * prevState.sampled_bsdf.pdf;
                            Float shiftedBsdfPdf = mainBsdfPdf;
                            Float shiftedLumPdf = mainLumPdf;
                            Spectrum shiftedEmitterRadiance = mainEmitterRadiance;

                            // Update throughput and pdf.
                            shifted.throughput *= shiftedBsdfValue;
                            shifted.pdf *= shiftedBsdfPdf;
                            shifted.eta *= prevState.sampled_bsdf.bRec.eta;

                            // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                            Float shiftedWeightDenominator = (shiftedPreviousPdf * shiftedPreviousPdf) *
                                                             ((shiftedLumPdf * shiftedLumPdf) +
                                                              (shiftedBsdfPdf * shiftedBsdfPdf));
                            weight = mainWeightNumerator /
                                     (D_EPSILON + shiftedWeightDenominator + mainWeightDenominator);

                            mainContribution = nextState.throughput * mainEmitterRadiance;
                            shiftedContribution =
                                    shifted.throughput *
                                    shiftedEmitterRadiance; // Note: Jacobian baked into .throughput.
                        } else if (shifted.connection_status == RAY_RECENTLY_CONNECTED) {
                            // Recently connected - follow the base path but evaluate BSDF to the new direction.
                            Vector3 incomingDirection = normalize(shifted.rRec.its.p - nextState.ray.o);
                            BSDFSamplingRecord bRec(prevState.its, prevState.its.toLocal(incomingDirection),
                                                    prevState.its.toLocal(nextState.ray.d), ERadiance);

                            // Note: mainBSDF is the BSDF at previousMainIts, which is the current position of the offset path.
                            EMeasure measure = (prevState.sampled_bsdf.bRec.sampledType & BSDF::EDelta) ? EDiscrete
                                                                                                        : ESolidAngle;

                            Spectrum shiftedBsdfValue = mainBSDF->eval(bRec, measure);
                            Float shiftedBsdfPdf = mainBSDF->pdf(bRec, measure);

                            Float shiftedLumPdf = mainLumPdf;
                            Spectrum shiftedEmitterRadiance = mainEmitterRadiance;

                            // Update throughput and pdf.
                            shifted.throughput *= shiftedBsdfValue;
                            shifted.pdf *= shiftedBsdfPdf;
                            shifted.eta *= prevState.sampled_bsdf.bRec.eta;

                            shifted.connection_status = RAY_CONNECTED;

                            // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                            Float shiftedWeightDenominator = (shiftedPreviousPdf * shiftedPreviousPdf) *
                                                             ((shiftedLumPdf * shiftedLumPdf) +
                                                              (shiftedBsdfPdf * shiftedBsdfPdf));
                            weight = mainWeightNumerator /
                                     (D_EPSILON + shiftedWeightDenominator + mainWeightDenominator);

                            mainContribution = nextState.throughput * mainEmitterRadiance;
                            shiftedContribution =
                                    shifted.throughput *
                                    shiftedEmitterRadiance; // Note: Jacobian baked into .throughput.
                        } else {
                            // Not connected - apply either reconnection or half-vector duplication shift.

                            const BSDF *shiftedBSDF = shifted.rRec.its.getBSDF(shifted.ray);

                            // Update the vertex type of the offset path.
                            VertexType shiftedVertexType = getVertexType(shifted, *m_config,
                                                                         prevState.sampled_bsdf.bRec.sampledType);

                            if (mainVertexType == VERTEX_TYPE_DIFFUSE && mainNextVertexType == VERTEX_TYPE_DIFFUSE &&
                                shiftedVertexType == VERTEX_TYPE_DIFFUSE) {
                                // Use reconnection shift.
                                bool lastSegment = (bounce + 1 == m_config->m_maxDepth);
                                // Optimization: Skip the last raycast and BSDF evaluation for the offset path when it won't contribute and isn't needed anymore.
                                if (!lastSegment || mainHitEmitter || nextState.its.hasSubsurface()) {
                                    ReconnectionShiftResult shiftResult;
                                    // FIXME: Environment map
//                                    bool environmentConnection = false;

                                    if (nextState.its.isValid()) {
                                        // This is an actual reconnection shift.
                                        shiftResult = reconnectShift(scene, nextState.ray.o, nextState.its.p,
                                                                     shifted.rRec.its.p, nextState.its.geoFrame.n,
                                                                     nextState.ray.time);
                                    } else {
                                        SLog(EError, "No support of environement map for now...");
                                    }

                                    if (!shiftResult.success) {
                                        // Failed to construct the offset path.
                                        shifted.alive = false;
                                        goto shift_failed;
                                    }

                                    Vector3 incomingDirection = -shifted.ray.d;
                                    Vector3 outgoingDirection = shiftResult.wo;

                                    BSDFSamplingRecord bRec(shifted.rRec.its,
                                                            shifted.rRec.its.toLocal(incomingDirection),
                                                            shifted.rRec.its.toLocal(outgoingDirection), ERadiance);

                                    // Strict normals check.
                                    if (m_config->m_strictNormals &&
                                        dot(outgoingDirection, shifted.rRec.its.geoFrame.n) *
                                        Frame::cosTheta(bRec.wo) <=
                                        0) {
                                        shifted.alive = false;
                                        goto shift_failed;
                                    }

                                    // Evaluate the BRDF to the new direction.
                                    Spectrum shiftedBsdfValue = shiftedBSDF->eval(bRec);
                                    Float shiftedBsdfPdf = shiftedBSDF->pdf(bRec);

                                    // Update throughput and pdf.
                                    shifted.throughput *= shiftedBsdfValue * shiftResult.jacobian;
                                    shifted.pdf *= shiftedBsdfPdf * shiftResult.jacobian;
                                    //shifted.eta *= bRec.eta;
                                    // ignore eta here

                                    shifted.connection_status = RAY_RECENTLY_CONNECTED;

                                    if (mainHitEmitter || nextState.its.hasSubsurface()) {
                                        // Also the offset path hit the emitter, as visibility was checked at reconnectShift or environmentShift.

                                        // Evaluate radiance to this direction.
                                        Spectrum shiftedEmitterRadiance(Float(0));
                                        Float shiftedLumPdf = Float(0);

                                        if (nextState.its.isValid()) {
                                            // Hit an object.
                                            if (mainHitEmitter) {
                                                shiftedEmitterRadiance = nextState.its.Le(-outgoingDirection);

                                                // Evaluate the light sampling PDF of the new segment.
                                                DirectSamplingRecord shiftedDRec;
                                                shiftedDRec.p = mainDRec.p;
                                                shiftedDRec.n = mainDRec.n;
                                                shiftedDRec.dist = (mainDRec.p - shifted.rRec.its.p).length();
                                                shiftedDRec.d = (mainDRec.p - shifted.rRec.its.p) / shiftedDRec.dist;
                                                shiftedDRec.ref = mainDRec.ref;
                                                shiftedDRec.refN = shifted.rRec.its.shFrame.n;
                                                shiftedDRec.object = mainDRec.object;

                                                shiftedLumPdf = scene->pdfEmitterDirect(shiftedDRec);
                                            }

                                            // Sub-surface scattering. Note: Should use the same random numbers as the base path!
                                            if (nextState.its.hasSubsurface()) {
                                                shiftedEmitterRadiance += nextState.its.LoSub(scene,
                                                                                              shifted.rRec.sampler,
                                                                                              -outgoingDirection,
                                                                                              bounce);
                                            }
                                        } else {
                                            // Hit the environment.
                                            shiftedEmitterRadiance = mainEmitterRadiance;
                                            shiftedLumPdf = mainLumPdf;
                                        }

                                        // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                                        Float shiftedWeightDenominator = (shiftedPreviousPdf * shiftedPreviousPdf) *
                                                                         ((shiftedLumPdf * shiftedLumPdf) +
                                                                          (shiftedBsdfPdf * shiftedBsdfPdf));
                                        weight = mainWeightNumerator /
                                                 (D_EPSILON + shiftedWeightDenominator + mainWeightDenominator);

                                        mainContribution = nextState.throughput * mainEmitterRadiance;
                                        shiftedContribution = shifted.throughput *
                                                              shiftedEmitterRadiance; // Note: Jacobian baked into .throughput.
                                    }
                                }
                            } else {

                                // Use half-vector duplication shift. These paths could not have been sampled by light sampling (by our decision).
                                Vector3 tangentSpaceIncomingDirection = shifted.rRec.its.toLocal(-shifted.ray.d);
                                Vector3 tangentSpaceOutgoingDirection;
                                Spectrum shiftedEmitterRadiance(Float(0));

                                const BSDF *shiftedBSDF = shifted.rRec.its.getBSDF(shifted.ray);

                                HalfVectorShiftResult shiftResult;
                                EMeasure measure;
                                BSDFSamplingRecord bRec(shifted.rRec.its, tangentSpaceIncomingDirection,
                                                        tangentSpaceOutgoingDirection, ERadiance);
                                Vector3 outgoingDirection;
                                VertexType shiftedVertexType;

                                // Deny shifts between Dirac and non-Dirac BSDFs.
                                bool bothDelta = (prevState.sampled_bsdf.bRec.sampledType & BSDF::EDelta) &&
                                                 (shiftedBSDF->getType() & BSDF::EDelta);
                                bool bothSmooth = (prevState.sampled_bsdf.bRec.sampledType & BSDF::ESmooth) &&
                                                  (shiftedBSDF->getType() & BSDF::ESmooth);
                                if (!(bothDelta || bothSmooth)) {
                                    shifted.alive = false;
                                    goto half_vector_shift_failed;
                                }

                                SAssert(fabs(shifted.ray.d.lengthSquared() - 1) < 0.01);

                                // Apply the local shift.
                                shiftResult = halfVectorShift(prevState.sampled_bsdf.bRec.wi,
                                                              prevState.sampled_bsdf.bRec.wo,
                                                              shifted.rRec.its.toLocal(-shifted.ray.d),
                                                              mainBSDF->getEta(),
                                                              shiftedBSDF->getEta());
                                bRec.wo = shiftResult.wo;

                                if (prevState.sampled_bsdf.bRec.sampledType & BSDF::EDelta) {
                                    // Dirac delta integral is a point evaluation - no Jacobian determinant!
                                    shiftResult.jacobian = Float(1);
                                }

                                if (shiftResult.success) {
                                    // Invertible shift, success.
                                    shifted.throughput *= shiftResult.jacobian;
                                    shifted.pdf *= shiftResult.jacobian;
                                    tangentSpaceOutgoingDirection = shiftResult.wo;
                                } else {
                                    // The shift is non-invertible so kill it.
                                    shifted.alive = false;
                                    goto half_vector_shift_failed;
                                }


                                outgoingDirection = shifted.rRec.its.toWorld(tangentSpaceOutgoingDirection);

                                // Update throughput and pdf.
                                measure = (prevState.sampled_bsdf.bRec.sampledType & BSDF::EDelta) ? EDiscrete
                                                                                                   : ESolidAngle;

                                shifted.throughput *= shiftedBSDF->eval(bRec, measure);
                                shifted.pdf *= shiftedBSDF->pdf(bRec, measure);

                                if (shifted.pdf == Float(0)) {
                                    // Offset path is invalid!
                                    shifted.alive = false;
                                    goto half_vector_shift_failed;
                                }

                                // Strict normals check to produce the same results as bidirectional methods when normal mapping is used.
                                if (m_config->m_strictNormals &&
                                    dot(outgoingDirection, shifted.rRec.its.geoFrame.n) * Frame::cosTheta(bRec.wo) <=
                                    0) {
                                    shifted.alive = false;
                                    goto half_vector_shift_failed;
                                }


                                // Update the vertex type.
                                shiftedVertexType = getVertexType(shifted, *m_config,
                                                                  prevState.sampled_bsdf.bRec.sampledType);

                                // Trace the next hit point.
                                shifted.ray = Ray(shifted.rRec.its.p, outgoingDirection, nextState.ray.time);

                                if (!scene->rayIntersect(shifted.ray, shifted.rRec.its)) {
                                    // Hit nothing - Evaluate environment radiance.
                                    const Emitter *env = scene->getEnvironmentEmitter();
                                    if (!env) {
                                        // Since base paths that hit nothing are not shifted, we must be symmetric and kill shifts that hit nothing.
                                        shifted.alive = false;
                                        goto half_vector_shift_failed;
                                    }
                                    if (nextState.its.isValid()) {
                                        // Deny shifts between env and non-env.
                                        shifted.alive = false;
                                        goto half_vector_shift_failed;
                                    }

                                    if (mainVertexType == VERTEX_TYPE_DIFFUSE &&
                                        shiftedVertexType == VERTEX_TYPE_DIFFUSE) {
                                        // Environment reconnection shift would have been used for the reverse direction!
                                        shifted.alive = false;
                                        goto half_vector_shift_failed;
                                    }

                                    // The offset path is no longer valid after this path segment.
                                    shiftedEmitterRadiance = env->evalEnvironment(shifted.ray);
                                    postponedShiftEnd = true;
                                } else {
                                    // Hit something.
                                    if (!nextState.its.isValid()) {
                                        // Deny shifts between env and non-env.
                                        shifted.alive = false;
                                        goto half_vector_shift_failed;
                                    }

                                    // FIXME: Can be wrong as well
                                    VertexType shiftedNextVertexType = getVertexType(shifted, *m_config,
                                                                                     prevState.sampled_bsdf.bRec.sampledType);

                                    // Make sure that the reverse shift would use this same strategy!
                                    // ==============================================================

                                    if (mainVertexType == VERTEX_TYPE_DIFFUSE &&
                                        shiftedVertexType == VERTEX_TYPE_DIFFUSE &&
                                        shiftedNextVertexType == VERTEX_TYPE_DIFFUSE) {
                                        // Non-invertible shift: the reverse-shift would use another strategy!
                                        shifted.alive = false;
                                        goto half_vector_shift_failed;
                                    }

                                    if (shifted.rRec.its.isEmitter()) {
                                        // Hit emitter.
                                        shiftedEmitterRadiance = shifted.rRec.its.Le(-shifted.ray.d);
                                    }
                                    // Sub-surface scattering. Note: Should use the same random numbers as the base path!
                                    if (shifted.rRec.its.hasSubsurface()) {
                                        shiftedEmitterRadiance += shifted.rRec.its.LoSub(scene, shifted.rRec.sampler,
                                                                                         -shifted.ray.d,
                                                                                         bounce);
                                    }

                                }


                                half_vector_shift_failed:
                                if (shifted.alive) {
                                    // Evaluate radiance difference using power heuristic between BSDF samples from base and offset paths.
                                    // Note: No MIS with light sampling since we don't use it for this connection type.
                                    weight = nextState.pdf /
                                             (D_EPSILON + shifted.pdf * shifted.pdf + nextState.pdf * nextState.pdf);
                                    mainContribution = nextState.throughput * mainEmitterRadiance;
                                    shiftedContribution = shifted.throughput *
                                                          shiftedEmitterRadiance; // Note: Jacobian baked into .throughput.
                                } else {
                                    // Handle the failure without taking MIS with light sampling, as we decided not to use it in the half-vector-duplication case.
                                    // Could have used it, but so far there has been no need. It doesn't seem to be very useful.
                                    weight = Float(1) / (D_EPSILON + nextState.pdf);
                                    mainContribution = nextState.throughput * mainEmitterRadiance;
                                    shiftedContribution = Spectrum(Float(0));

                                    // Disable the failure detection below since the failure was already handled.
                                    shifted.alive = true;
                                    postponedShiftEnd = true;

                                    // (TODO: Restructure into smaller functions and get rid of the gotos... Although this may mean having lots of small functions with a large number of parameters.)
                                }
                            }
                        }
                    }

                    shift_failed:
                    if (!shifted.alive) {
                        // The offset path cannot be generated; Set offset PDF and offset throughput to zero.
                        weight = mainWeightNumerator / (D_EPSILON + mainWeightDenominator);
                        mainContribution = nextState.throughput * mainEmitterRadiance;
                        shiftedContribution = Spectrum((Float) 0);
                    }

                    // Note: Using also the offset paths for the throughput estimate, like we do here, provides some advantage when a large reconstruction alpha is used,
                    // but using only throughputs of the base paths doesn't usually lose by much.
                    if (bounce + 1 >= m_config->m_minDepth) {
                      if(m_config->reusePrimal) {
                        mainTraced.addRadiance(mainContribution, weight);
                        shifted.addRadiance(shiftedContribution, weight);
                      }
                      shifted.addGradient(shiftedContribution - mainContribution, weight);
                    }

                    if (postponedShiftEnd) {
                        shifted.alive = false;
                    }
                }
            }

            if (bounce++ >= m_config->m_rrDepth) {
                SLog(EError, "not implemented");
            }
        }
    }

    void trace_base(RayState &main, Spectrum &out_veryDirect) {
        const Scene *scene = main.rRec.scene;

        // Perform the first ray intersection for the base path (or ignore if the intersection has already been provided).
        main.rRec.rayIntersect(main.ray);
        main.ray.mint = Epsilon;

        if (!main.rRec.its.isValid()) {
            // First hit is not in the scene so can't continue. Also there there are no paths to shift.
            // Add potential very direct light from the environment as gradients are not used for that.
            if (main.rRec.type & RadianceQueryRecord::EEmittedRadiance) {
                out_veryDirect += main.throughput * scene->evalEnvironment(main.ray);
            }
            return;
        }

        // Add very direct light from non-environment.
        // Include emitted radiance if requested.
        if (main.rRec.its.isEmitter() && (main.rRec.type & RadianceQueryRecord::EEmittedRadiance)) {
            out_veryDirect += main.throughput * main.rRec.its.Le(-main.ray.d);
        }

        if (m_config->m_strictNormals) {
            // If 'strictNormals'=true, when the geometric and shading normals classify the incident direction to the same side, then the main path is still good.
            if (dot(main.ray.d, main.rRec.its.geoFrame.n) * Frame::cosTheta(main.rRec.its.wi) >= 0) {
                // This is an impossible base path.
                return;
            }
        }

        // Main path tracing loop.
        main.rRec.depth = 1;
        while (main.rRec.depth < m_config->m_maxDepth || m_config->m_maxDepth < 0) {
            if (m_config->m_strictNormals) {
                if (dot(main.ray.d, main.rRec.its.geoFrame.n) * Frame::cosTheta(main.rRec.its.wi) >= 0) {
                    // This is an impossible main path, and there are no more paths to shift.
                    return;
                }
            }

            /* ==================================================================== */
            /*               BSDF sampling and emitter hits                         */
            /* ==================================================================== */
            // Sample a new direction from BSDF * cos(theta).
            {
                auto bsdfResult = sampleBSDF(main);
                {
                    const BSDF *mainBSDF = main.rRec.its.getBSDF(main.ray);
                    state.emplace_back(
                            BaseState{
                                    main.rRec.its,
                                    main.throughput,
                                    main.pdf,
                                    main.ray,
                                    mainBSDF,
                                    bsdfResult,
                            });
                }

                if (bsdfResult.pdf <= (Float) 0.0) {
                    // Impossible base path.
                    break;
                }

                const Vector mainWo = main.rRec.its.toWorld(bsdfResult.bRec.wo);

                // Prevent light leaks due to the use of shading normals.
                Float mainWoDotGeoN = dot(main.rRec.its.geoFrame.n, mainWo);
                if (m_config->m_strictNormals && mainWoDotGeoN * Frame::cosTheta(bsdfResult.bRec.wo) <= 0) {
                    break;
                }

                auto dRec = DirectSamplingRecord(main.rRec.its);
                const BSDF *mainBSDF = main.rRec.its.getBSDF(main.ray);

                // Update the vertex types.
                VertexType mainVertexType = getVertexType(main, *m_config, bsdfResult.bRec.sampledType);
                main.ray = Ray(main.rRec.its.p, mainWo, main.ray.time);
                if (!scene->rayIntersect(main.ray, main.rRec.its)) {
                    // TODO: The environment map are not handle
                    break;
                }

                main.throughput *= bsdfResult.weight * bsdfResult.pdf;
                main.pdf *= bsdfResult.pdf;
                main.eta *= bsdfResult.bRec.eta;
            }

            /* ==================================================================== */
            /*                    Russian roulette                                  */
            /* ==================================================================== */
            // Stop if the base path hit the environment.
            main.rRec.type = RadianceQueryRecord::ERadianceNoEmission;
            if (!main.rRec.its.isValid() || !(main.rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance)) {
                break;
            }

            if (main.rRec.depth + 1 == m_config->m_maxDepth) {
                // This is the last vertex, need to push it
                const BSDF *mainBSDF = main.rRec.its.getBSDF(main.ray);
                auto bsdfResult = sampleBSDF(main); // TODO: Select the one last
                state.emplace_back(
                        BaseState{
                                main.rRec.its,
                                main.throughput,
                                main.pdf,
                                main.ray,
                                mainBSDF,
                                bsdfResult,
                        });
            }

            if (main.rRec.depth++ >= m_config->m_rrDepth) {
                SLog(EError, "Not supported");
                /* Russian roulette: try to keep path weights equal to one,
                   while accounting for the solid angle compression at refractive
                   index boundaries. Stop with at least some probability to avoid
                   getting stuck (e.g. due to total internal reflection) */

                Float q = std::min((main.throughput / main.pdf).max() * main.eta * main.eta, (Float) 0.95f);
                if (main.rRec.nextSample1D() >= q)
                    break;

                main.pdf *= q;
            }
        }
        return;
    }
};

MTS_NAMESPACE_END

#endif //MITSUBA_EXPLICIT_H
