//
// Created by muliana on 10/30/18.
//

#include "shiftmapping.h"

#ifndef MITSUBA_ORIGINAL_VOLUME_H
#define MITSUBA_ORIGINAL_VOLUME_H

#define MEDIUM_DIRECT_TRANS 1
#define MEDIUM_INTER_TRANS 1
#define SURFACE_DIRECT_TRANS 1
#define SURFACE_INTER_TRANS 1

MTS_NAMESPACE_BEGIN

class ShiftMappingVolumeOriginal : public ShiftMapping {
public:
    ShiftMappingVolumeOriginal(const GradientPathTracerConfig *config) : ShiftMapping(config) {}

    void evaluateReuse(RayState *rays, int secondaryCount,
                       Spectrum &out_veryDirect, int id_main,
                       std::vector<GradientInfo> &gradients) override {
        SLog(EError, "Impossible to use reuse with the original shift mapping");
    }

    void evaluate(RayState &main, RayState *shiftedRays, int secondaryCount, Spectrum &out_veryDirect) override {
        const Scene *scene = main.rRec.scene;
        const Medium *med = main.rRec.medium;
        bool mediumScattered = true; // Allow to have surfaces

        // Perform the first ray intersection for the base path (or ignore if the intersection has already been provided).
        main.rRec.rayIntersect(main.ray);
        main.ray.mint = Epsilon;

        // Perform the same first ray intersection for the offset paths.
        for (int i = 0; i < secondaryCount; ++i) {
            RayState &shifted = shiftedRays[i];
            shifted.rRec.rayIntersect(shifted.ray);
            shifted.ray.mint = Epsilon;
        }

        if (!main.rRec.its.isValid()) {
            // First hit is not in the scene so can't continue. Also there there are no paths to shift.

            // Add potential very direct light from the environment as gradients are not used for that.
            if ((main.rRec.type & RadianceQueryRecord::EEmittedRadiance) && mediumScattered) {
                out_veryDirect += main.throughput * scene->evalEnvironment(main.ray);
            }

            //SLog(EInfo, "Main ray(%d): First hit not in scene.", rayCount);
            return;
        }

        // Add very direct light from non-environment.
        {
            // Include emitted radiance if requested.
            if (main.rRec.its.isEmitter() &&
                (main.rRec.type & RadianceQueryRecord::EEmittedRadiance) &&
                mediumScattered) {
                out_veryDirect += main.throughput * main.rRec.its.Le(-main.ray.d);
            }

            // Include radiance from a subsurface scattering model if requested. Note: Not tested!
            if (main.rRec.its.hasSubsurface() && (main.rRec.type & RadianceQueryRecord::ESubsurfaceRadiance) &&
                mediumScattered) {
                out_veryDirect += main.throughput * main.rRec.its.LoSub(scene, main.rRec.sampler, -main.ray.d, 0);
            }
        }

        // If no intersection of an offset ray could be found, its offset paths can not be generated.
        for (int i = 0; i < secondaryCount; ++i) {
            RayState &shifted = shiftedRays[i];
            if (!shifted.rRec.its.isValid()) {
                shifted.alive = false;
            }
        }

        // Strict normals check to produce the same results as bidirectional methods when normal mapping is used.
        if (m_config->m_strictNormals) {
            // If 'strictNormals'=true, when the geometric and shading normals classify the incident direction to the same side, then the main path is still good.
            if (dot(main.ray.d, main.rRec.its.geoFrame.n) * Frame::cosTheta(main.rRec.its.wi) >= 0) {
                // This is an impossible base path.
                return;
            }

            for (int i = 0; i < secondaryCount; ++i) {
                RayState &shifted = shiftedRays[i];

                if (dot(shifted.ray.d, shifted.rRec.its.geoFrame.n) * Frame::cosTheta(shifted.rRec.its.wi) >= 0) {
                    // This is an impossible offset path.
                    shifted.alive = false;
                }
            }
        }


        // Main path tracing loop.
        int nbSurfaceInteractions = 0;
        main.rRec.depth = 1;
        // To know if we have bounced inside the media or not
        bool sampledMediumEvent = false;
        // If we are dealing with a surface interaction on the main path
        // We need to store

        while (main.rRec.depth < m_config->m_maxDepth || m_config->m_maxDepth < 0) {

            // Strict normals check to produce the same results as bidirectional methods when normal mapping is used.
            // If 'strictNormals'=true, when the geometric and shading normals classify the incident direction to the same side, then the main path is still good.
            if (m_config->m_strictNormals) {
                SLog(EError, "Not possible to do strict normal here");
            }

            // Some optimizations can be made if this is the last traced segment.
            bool lastSegment = (main.rRec.depth + 1 == m_config->m_maxDepth);

            if (main.rRec.depth == 1) {
                // At the first loop we need to sample the medium from the mainPath
                if (med && med->sampleDistance(Ray(main.ray, 0, main.rRec.its.t), main.mRec, main.rRec.sampler)) {
                    sampledMediumEvent = true;
                } else {
                    sampledMediumEvent = false;
                }
            }

            if (sampledMediumEvent) {
                if (nbSurfaceInteractions < m_config->m_minCameraDepth) {
                    return; // Early kill the path
                }

                /* ====================================================================
                 *  We are inside the medium
                 */
                mediumScattered = true; // Mark that we can compute direct lighting

                // Update the transmittance for the main path
                const PhaseFunction *mainPhase = main.mRec.getPhaseFunction();
                main.throughput *= main.mRec.sigmaS * main.mRec.transmittance / main.mRec.pdfSuccess;
                main.pdf *= main.mRec.pdfSuccess;

                // Here more tricky, we need to update the transmittance for shifted path
                // There is two possibility:
                // 1) No reconnection is done yet, we check that we enough distance to shift
                // 2) We have done the reconnection, so we use main path transmittance sampling
                for (int i = 0; i < secondaryCount; ++i) {
                    // Evaluate and apply the gradient.
                    RayState &shifted = shiftedRays[i];

                    bool shiftSuccessful = shifted.alive;

                    // Construct the offset path.
                    if (shiftSuccessful) {
                        if (shifted.connection_status == RAY_CONNECTED) {
                            // We have done the reconnection, as the incomming direction does not matter
                            // Follow the base path. All relevant vertices are shared.
                            shifted.throughput *= main.mRec.sigmaS * main.mRec.transmittance / main.mRec.pdfSuccess;
                            shifted.pdf *= main.mRec.pdfSuccess;

                        } else if (shifted.connection_status == RAY_RECENTLY_CONNECTED) {
                            // We recently reconnected the path to the base path
                            // However, the travel distance inside the smoke can be different
                            // compare to the base path, so we need to re-evaluate
                            // based on the previous shifted path position
                            Vector d = main.mRec.p - shifted.prevPos();
                            Float distReco = d.length();
                            d /= distReco;

                            Ray shiftRay(shifted.prevPos(), d, 0, distReco, main.mRec.time);
                            MediumSamplingRecord mRecReco; // We use a another medium rec to not erase previous one
                            med->eval(shiftRay, mRecReco);

                            // Update shifted values
                            if (mRecReco.transmittance.isZero()) {
                                shifted.alive = false;
                            } else {
                                shifted.throughput *= mRecReco.sigmaS * mRecReco.transmittance / main.mRec.pdfSuccess;
                                shifted.pdf *= mRecReco.pdfSuccess;
                            }
                        } else {
                            // We did not have reconnected yet the path
                            // We just check that we can copy the distance
                            // If it is not the case, the shift is aborded
                            //SAssert(!shifted.lastVolume);
                            if (shifted.rRec.its.t < main.mRec.t) {
                                shifted.alive = false;
                            } else {
                                Ray shiftRay = shifted.ray;
                                shiftRay.maxt = main.mRec.t;
                                med->eval(shiftRay, shifted.mRec);

                                shifted.throughput *=
                                        shifted.mRec.sigmaS * shifted.mRec.transmittance / main.mRec.pdfSuccess;
                                shifted.pdf *= shifted.mRec.pdfSuccess;

                                // Weird but mitsuba does not fill mRec when we eval
                                shifted.mRec.t = shiftRay.maxt - shiftRay.mint;
                                shifted.mRec.p = shiftRay(shifted.mRec.t);
                                shifted.lastVolume = true;
                            }
                        }
                    }
                }

                /* ==================================================================== */
                /*                     Direct illumination sampling                     */
                /* ==================================================================== */
                if (mediumScattered) {
                    // TODO: Always assuming that the phase function is smooth
                    if (main.rRec.type & RadianceQueryRecord::EDirectMediumRadiance
                        && main.rRec.depth + 1 >= m_config->m_minDepth &&
                        nbSurfaceInteractions >= m_config->m_minCameraDepth) {
                        // Sample an emitter and evaluate f = f/p * p for it. */
                        DirectSamplingRecord dRec(main.mRec.p, main.mRec.time);

                        mitsuba::Point2 lightSample = main.rRec.nextSample2D();
                        int interactions = 0; //m_config->m_maxDepth - main.rRec.depth - 1; // -> no interaction...
                        std::pair<Spectrum, bool> emitterTupleValue = scene->sampleAttenuatedEmitterDirectVisible(dRec,
                                                                                                                  med,
                                                                                                                  interactions,
                                                                                                                  lightSample,
                                                                                                                  main.rRec.sampler,
                                                                                                                  false);
                        Spectrum mainEmitterRadiance = emitterTupleValue.first;
                        bool mainEmitterVisible = emitterTupleValue.second;
                        const auto *emitter = dynamic_cast<const Emitter *>(dRec.object);

                        // If the emitter sampler produces a non-emitter, that's a problem.
                        SAssert(emitter != nullptr);

                        // Add radiance and gradients to the base path and its offset path.
                        // Query the BSDF to the emitter's direction.
                        PhaseFunctionSamplingRecord mainPRec(main.mRec,
                                                             -main.ray.d,
                                                             dRec.d,
                                                             ERadiance);

                        // Evaluate BSDF * cos(theta).
                        Float mainPhaseValue = mainPhase->eval(mainPRec); //main->eval(mainBRec);

                        // Calculate the probability density of having generated the sampled path segment by BSDF sampling. Note that if the emitter is not visible, the probability density is zero.
                        // Even if the BSDF sampler has zero probability density, the light sampler can still sample it.
                        Float mainPhasePdf =
                                (emitter->isOnSurface() && dRec.measure == ESolidAngle && mainEmitterVisible)
                                ? mainPhase->pdf(mainPRec) : 0;

                        // There values are probably needed soon for the Jacobians.
                        Float mainDistanceSquared = (main.mRec.p - dRec.p).lengthSquared();
                        Float mainOpposingCosine = dot(dRec.n, (main.mRec.p - dRec.p)) / sqrt(mainDistanceSquared);

                        // Power heuristic weights for the following strategies: light sample from base, BSDF sample from base.
                        Float mainWeightDenominator = (main.pdf) * (dRec.pdf + mainPhasePdf);
                        Float mainWeightNumerator = main.pdf * dRec.pdf;
#ifdef CENTRAL_RADIANCE
                        if (main.rRec.depth + 1 >= m_config->m_minDepth && mainWeightDenominator != 0.f && nbSurfaceInteractions >= m_config->m_minCameraDepth) {
#if MEDIUM_DIRECT_TRANS
              main.addRadiance(main.throughput * (mainPhaseValue * mainEmitterRadiance),
                               mainWeightNumerator / (mainWeightDenominator));
#endif
            }
#endif

                        // The base path is good. Add radiance differences to offset paths.
                        for (int i = 0; i < secondaryCount; ++i) {
                            // Evaluate and apply the gradient.
                            RayState &shifted = shiftedRays[i];

                            Spectrum mainContribution(Float(0));
                            Spectrum shiftedContribution(Float(0));
                            auto weight = Float(0);

                            bool shiftSuccessful = shifted.alive;

                            // Construct the offset path.
                            if (shiftSuccessful) {
                                // Generate the offset path.
                                if (shifted.connection_status == RAY_CONNECTED) {
                                    // Follow the base path. All relevant vertices are shared.
                                    Float shiftedPhasePdf = mainPhasePdf;
                                    Float shiftedDRecPdf = dRec.pdf;
                                    Float shiftedPhaseValue = mainPhaseValue;
                                    Spectrum shiftedEmitterRadiance = mainEmitterRadiance;
                                    auto jacobian = (Float) 1;

                                    // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                                    Float shiftedWeightDenominator =
                                            (jacobian * shifted.pdf) * (shiftedDRecPdf + shiftedPhasePdf);
                                    weight = mainWeightNumerator /
                                             (D_EPSILON + shiftedWeightDenominator + mainWeightDenominator);

                                    mainContribution = main.throughput * (shiftedPhaseValue * mainEmitterRadiance);
                                    shiftedContribution = jacobian * shifted.throughput
                                                          * (shiftedPhaseValue * shiftedEmitterRadiance);

                                    // Note: The Jacobians were baked into shifted.pdf and shifted.throughput at connection phase.
                                } else if (shifted.connection_status == RAY_RECENTLY_CONNECTED) {
                                    // Follow the base path. The current vertex is shared, but the incoming directions differ.
                                    Vector3 incomingDirection = normalize(shifted.prevPos() - main.mRec.p);
                                    PhaseFunctionSamplingRecord PRec(main.mRec,
                                                                     incomingDirection,
                                                                     dRec.d,
                                                                     ERadiance);

                                    // Sample the BSDF.
                                    Float shiftedPhasePdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle
                                                             && mainEmitterVisible) ? mainPhase->pdf(PRec)
                                                                                    : 0; // The BSDF sampler can not sample occluded path segments.
                                    Float shiftedDRecPdf = dRec.pdf;
                                    Float shiftedPhaseValue = mainPhase->eval(PRec);
                                    Spectrum shiftedEmitterRadiance = mainEmitterRadiance;
                                    auto jacobian = (Float) 1;

                                    // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                                    Float shiftedWeightDenominator =
                                            (jacobian * shifted.pdf) * ((shiftedDRecPdf + shiftedPhasePdf));
                                    weight = mainWeightNumerator
                                             / (D_EPSILON + shiftedWeightDenominator + mainWeightDenominator);

                                    mainContribution = main.throughput * (mainPhaseValue * mainEmitterRadiance);
                                    shiftedContribution = jacobian * shifted.throughput
                                                          * (shiftedPhaseValue * shiftedEmitterRadiance);

                                    // Note: The Jacobians were baked into shifted.pdf and shifted.throughput at connection phase.
                                } else {
                                    // Reconnect to the sampled light vertex. No shared vertices.
                                    SAssert(shifted.connection_status == RAY_NOT_CONNECTED);

                                    // This implementation uses light sampling only for the reconnect-shift.
                                    // When one of the BSDFs is very glossy, light sampling essentially reduces to a failed shift anyway.
                                    bool mainAtPointLight = (dRec.measure == EDiscrete);

                                    // We force to treat phase function as always smooth
                                    // for the current implementation, it is assume that the phase function is diffuse
                                    VertexType mainVertexType = VertexType::VERTEX_TYPE_DIFFUSE;
                                    VertexType shiftedVertexType = VertexType::VERTEX_TYPE_DIFFUSE;

                                    if (mainAtPointLight || (mainVertexType == VERTEX_TYPE_DIFFUSE
                                                             && shiftedVertexType == VERTEX_TYPE_DIFFUSE)) {
                                        // Get emitter radiance.
                                        DirectSamplingRecord shiftedDRec(shifted.mRec.p, shifted.mRec.time);
                                        std::pair<Spectrum, bool>
                                                emitterTupleShiftValue = scene->sampleAttenuatedEmitterDirectVisible(
                                                shiftedDRec,
                                                med,
                                                interactions,
                                                lightSample,
                                                shifted.rRec.sampler, false);
                                        bool shiftedEmitterVisible = emitterTupleShiftValue.second;
                                        Spectrum shiftedEmitterRadiance =
                                                emitterTupleShiftValue.first * shiftedDRec.pdf;
                                        Float shiftedDRecPdf = shiftedDRec.pdf;

                                        // Sample the BSDF.
                                        Float shiftedDistanceSquared =
                                                (dRec.p - shifted.mRec.p).lengthSquared();
                                        Vector emitterDirection =
                                                (dRec.p - shifted.mRec.p) / sqrt(shiftedDistanceSquared);
                                        Float shiftedOpposingCosine = -dot(dRec.n, emitterDirection);

                                        PhaseFunctionSamplingRecord PRec(shifted.mRec,
                                                                         -shifted.ray.d,
                                                                         emitterDirection,
                                                                         ERadiance);

                                        // Strict normals check, to make the output match with bidirectional methods when normal maps are present.
                                        Float shiftedPhaseValue = mainPhase->eval(PRec);
                                        Float shiftedPhasePdf =
                                                (emitter->isOnSurface() && dRec.measure == ESolidAngle
                                                 && shiftedEmitterVisible) ? mainPhase->pdf(PRec) : 0;
                                        Float jacobian = std::abs(shiftedOpposingCosine * mainDistanceSquared)
                                                         / (Epsilon + std::abs(
                                                mainOpposingCosine * shiftedDistanceSquared));

                                        // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                                        Float shiftedWeightDenominator =
                                                (jacobian * shifted.pdf) * ((shiftedDRecPdf + shiftedPhasePdf));

                                        weight = mainWeightNumerator / (D_EPSILON + shiftedWeightDenominator
                                                                        + mainWeightDenominator);

                                        mainContribution =
                                                main.throughput * (shiftedPhaseValue * mainEmitterRadiance);
                                        if (weight != 0) {
                                            shiftedContribution = jacobian * shifted.throughput
                                                                  * (shiftedPhaseValue * shiftedEmitterRadiance) /
                                                                  dRec.pdf;
                                        }
                                    }

                                }
                            }

                            if (!shiftSuccessful) {
                                // The offset path cannot be generated; Set offset PDF and offset throughput to zero. This is what remains.

                                // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset. (Offset path has zero PDF)
                                weight = mainWeightNumerator / (D_EPSILON + mainWeightDenominator);

                                mainContribution = main.throughput * (mainPhaseValue * mainEmitterRadiance);
                                shiftedContribution = Spectrum((Float) 0);
                            }

                            // Note: Using also the offset paths for the throughput estimate, like we do here, provides some advantage when a large reconstruction alpha is used,
                            // but using only throughputs of the base paths doesn't usually lose by much.
                            if (main.rRec.depth + 1 >= m_config->m_minDepth &&
                                nbSurfaceInteractions >= m_config->m_minCameraDepth) {
#if MEDIUM_DIRECT_TRANS
#ifndef CENTRAL_RADIANCE
                                main.addRadiance(mainContribution, weight);
                                shifted.addRadiance(shiftedContribution, weight);
#endif
                                shifted.addGradient(shiftedContribution - mainContribution, weight);
#endif
                            }
                        } // for(int i = 0; i < secondaryCount; ++i)
                    }
                } // Sample incoming radiance from lights.

                /* ==================================================================== */
                /*                     Phase function sampling                          */
                /* ==================================================================== */

                // Sample a new direction from BSDF * cos(theta).
                PhaseSampleResult mainPhaseResult = samplePhase(main, main.mRec);
                if (mainPhaseResult.pdf <= (Float) 0.0) {
                    // Impossible base path.
                    break;
                }
                const Vector mainWo = mainPhaseResult.pRec.wo;
                MediumSamplingRecord previousMainMediumRec = main.mRec;

                // Trace a ray in the sampled direction.
                bool mainHitEmitter = false;
                Spectrum mainEmitterRadiance = Spectrum((Float) 0);

                // FIXME: We modify the record if we really hit the object
                DirectSamplingRecord mainDRec(main.mRec.p, main.mRec.time);

                // Update the vertex types.
                // We are inside the volume so the BSDF is diffuse
                VertexType mainVertexType = VertexType::VERTEX_TYPE_DIFFUSE;
                VertexType mainNextVertexType;

                // Use the position inside the volume
                main.ray = Ray(main.mRec.p, mainWo, main.ray.time);
                if (scene->rayIntersect(main.ray, main.rRec.its)) {
                    sampledMediumEvent = med->sampleDistance(Ray(main.ray, 0, main.rRec.its.t), main.mRec,
                                                             main.rRec.sampler);
                    mainNextVertexType = getVertexType(main, *m_config, BSDF::EGlossy);
#if ONLY_EXPLICIT_CONNECTION
                    mainHitEmitter = false;
#else
                    // Intersected something - check if it was a luminaire.
                    if (main.rRec.its.isEmitter()) {
                        mainEmitterRadiance = main.rRec.its.Le(-main.ray.d);
                        // Evaluate transmittance
                        mainEmitterRadiance *= med->evalTransmittance(Ray(main.ray, 0, main.rRec.its.t),
                                                                      main.rRec.sampler);
                        mainDRec.setQuery(main.ray, main.rRec.its);
                        mainHitEmitter = true;
                    }

                    // Sub-surface scattering.
                    if (main.rRec.its.hasSubsurface() && (main.rRec.type & RadianceQueryRecord::ESubsurfaceRadiance)) {
                        SLog(EError, "Not tested BSSRDF");
                        mainEmitterRadiance +=
                                main.rRec.its.LoSub(scene, main.rRec.sampler, -main.ray.d, main.rRec.depth);
                    }
#endif
                } else {
                    // Intersected nothing -- perhaps there is an environment map?
                    const Emitter *env = scene->getEnvironmentEmitter();

                    if (env) {
                        SLog(EError, "No support of env map");
                    } else {
                        // Nothing to do anymore.
                        break;
                    }
                    mainNextVertexType = VERTEX_TYPE_DIFFUSE;

                    // Quit
                    return;
                }

                // Continue the shift.
                Float mainPhasePdf = mainPhaseResult.pdf;
                Float mainPreviousPdf = main.pdf;

                main.throughput *= mainPhaseResult.weight;
                main.pdf *= mainPhaseResult.pdf;

                // Compute the probability density of generating base path's direction using the implemented direct illumination sampling technique.
                // FIXME: Change pdfEmitterDirect to attenuated (check other code)
                const Float mainLumPdf = (mainHitEmitter && main.rRec.depth + 1 >= m_config->m_minDepth &&
                                          nbSurfaceInteractions >= m_config->m_minCameraDepth) ?
                                         scene->pdfEmitterDirect(mainDRec) : 0;

                // Power heuristic weights for the following strategies: light sample from base, BSDF sample from base.
                // FIXME: Maybe the weight are wrong
                Float mainWeightNumerator = mainPreviousPdf * mainPhaseResult.pdf;
                Float mainWeightDenominator = (mainPreviousPdf * ((mainLumPdf + mainPhasePdf)));

#ifdef CENTRAL_RADIANCE
                if (main.rRec.depth + 1 >= m_config->m_minDepth && mainWeightDenominator != 0.f && nbSurfaceInteractions >= m_config->m_minCameraDepth) {
#if MEDIUM_INTER_TRANS
          main.addRadiance(main.throughput * mainEmitterRadiance,
                           mainWeightNumerator / (mainWeightDenominator));
#endif
        }
#endif

                // Construct the offset paths and evaluate emitter hits.
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
                            Float shiftedPhaseValue = mainPhaseResult.weight;
                            Float shiftedPhasePdf = mainPhasePdf;
                            Float shiftedLumPdf = mainLumPdf;
                            Spectrum shiftedEmitterRadiance = mainEmitterRadiance;

                            // Update throughput and pdf.
                            shifted.throughput *= shiftedPhaseValue;
                            shifted.pdf *= shiftedPhasePdf;

                            // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                            Float shiftedWeightDenominator = shiftedPreviousPdf * ((shiftedLumPdf + shiftedPhasePdf));
                            weight = mainWeightNumerator
                                     / (D_EPSILON + shiftedWeightDenominator + mainWeightDenominator);

                            mainContribution = main.throughput * mainEmitterRadiance;
                            shiftedContribution = shifted.throughput
                                                  * shiftedEmitterRadiance; // Note: Jacobian baked into .throughput.
                        } else if (shifted.connection_status == RAY_RECENTLY_CONNECTED) {
                            // Recently connected - follow the base path but evaluate BSDF to the new direction.
                            Vector3 incomingDirection = normalize(shifted.prevPos() - main.ray.o);
                            PhaseFunctionSamplingRecord pRec(previousMainMediumRec,
                                                             incomingDirection,
                                                             main.ray.d);

                            // Note: mainBSDF is the BSDF at previousMainIts, which is the current position of the offset path.
                            Float shiftedPhaseValue = mainPhase->eval(pRec);
                            Float shiftedPhasePdf = mainPhase->pdf(pRec);

                            Float shiftedLumPdf = mainLumPdf;
                            Spectrum shiftedEmitterRadiance = mainEmitterRadiance;

                            // Update throughput and pdf.
                            shifted.throughput *= shiftedPhaseValue / mainPhasePdf;
                            shifted.pdf *= shiftedPhasePdf;

                            shifted.connection_status = RAY_CONNECTED;

                            // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                            Float shiftedWeightDenominator = shiftedPreviousPdf * ((shiftedLumPdf + shiftedPhasePdf));
                            weight = mainWeightNumerator
                                     / (D_EPSILON + shiftedWeightDenominator + mainWeightDenominator);

                            mainContribution = main.throughput * mainEmitterRadiance;
                            shiftedContribution = shifted.throughput
                                                  * shiftedEmitterRadiance; // Note: Jacobian baked into .throughput.
                        } else {
                            // Not connected - apply either reconnection or half-vector duplication shift.
                            // Update the vertex type of the offset path.
                            // We have force that the shift path is inside the volume
                            // So we can classify as diffuse
                            VertexType shiftedVertexType = VertexType::VERTEX_TYPE_DIFFUSE;

                            if (mainVertexType == VERTEX_TYPE_DIFFUSE && mainNextVertexType == VERTEX_TYPE_DIFFUSE
                                && shiftedVertexType == VERTEX_TYPE_DIFFUSE) {
                                // Use reconnection shift.

                                // Optimization: Skip the last raycast and BSDF evaluation for the offset path when it won't contribute and isn't needed anymore.
                                if (!lastSegment || mainHitEmitter || main.rRec.its.hasSubsurface()) {
                                    ReconnectionShiftResult shiftResult;
                                    if (main.rRec.its.isValid()) {
                                        // This is an actual reconnection shift.
                                        if (!sampledMediumEvent) {
                                            shiftResult = reconnectShift(scene,
                                                                         main.ray.o, // Previous intersection
                                                                         main.rRec.its.p, // New intersection
                                                                         shifted.mRec.p, // Previous intersection shift
                                                                         main.rRec.its.geoFrame.n, // Normal (OK)
                                                                         main.ray.time);
                                        } else {
                                            shiftResult = reconnectShift(scene,
                                                                         main.ray.o, // Previous intersection
                                                                         main.mRec.p, // New intersection
                                                                         shifted.mRec.p, // Previous intersection shift
                                                                         Normal(0.f), // FIXME: Normal ignored here
                                                                         main.ray.time, true);
                                        }
                                    } else {
                                        // This is a reconnection at infinity in environment direction.
                                        SLog(EError, "Not implemented");
                                        shifted.alive = false;
                                        goto shift_failed_volume;
                                    }

                                    if (!shiftResult.success) {
                                        // Failed to construct the offset path.
                                        shifted.alive = false;
                                        goto shift_failed_volume;
                                    }

                                    // FIXME: FIXME FIXME FIXME FIXME
                                    // FIXME: The transmittance is not taken into account here
                                    Vector3 incomingDirection = -shifted.ray.d;
                                    Vector3 outgoingDirection = shiftResult.wo;
                                    PhaseFunctionSamplingRecord pRec(previousMainMediumRec,
                                                                     incomingDirection,
                                                                     outgoingDirection);

                                    // Evaluate the BRDF to the new direction.
                                    // FIXME: Assuming spatial constant phase function
                                    Float shiftedPhaseValue = mainPhase->eval(pRec);
                                    Float shiftedPhasePdf = mainPhase->pdf(pRec);

                                    // Update throughput and pdf.
                                    shifted.throughput *= shiftedPhaseValue * shiftResult.jacobian / mainPhasePdf;
                                    shifted.pdf *= shiftedPhasePdf * shiftResult.jacobian;

                                    shifted.connection_status = RAY_RECENTLY_CONNECTED;

                                    if (mainHitEmitter) {
                                        // Also the offset path hit the emitter, as visibility was checked at reconnectShift or environmentShift.

                                        // Evaluate radiance to this direction.
                                        Spectrum shiftedEmitterRadiance(Float(0));
                                        auto shiftedLumPdf = Float(0);

                                        if (main.rRec.its.isValid()) {
                                            // Hit an object.
                                            if (mainHitEmitter) {
                                                shiftedEmitterRadiance = main.rRec.its.Le(-outgoingDirection);
                                                // FIXME: Fix position of computation
                                                // Evaluate the light sampling PDF of the new segment.
                                                DirectSamplingRecord shiftedDRec;
                                                shiftedDRec.p = mainDRec.p;
                                                shiftedDRec.n = mainDRec.n; // OK as it is the light source
                                                shiftedDRec.dist = (mainDRec.p - shifted.mRec.p).length();
                                                shiftedDRec.d = (mainDRec.p - shifted.mRec.p) / shiftedDRec.dist;
                                                shiftedDRec.ref = mainDRec.ref;
                                                shiftedDRec.refN = mainDRec.refN;
                                                shiftedDRec.object = mainDRec.object;
                                                shiftedDRec.measure = ESolidAngle;

                                                shiftedEmitterRadiance *= med->evalTransmittance(Ray(shifted.mRec.p,
                                                                                                     shiftedDRec.d,
                                                                                                     0,
                                                                                                     shiftedDRec.dist,
                                                                                                     shifted.mRec.time),
                                                                                                 shifted.rRec.sampler);
                                                shiftedLumPdf = scene->pdfEmitterDirect(shiftedDRec);
                                            }

                                            // Sub-surface scattering. Note: Should use the same random numbers as the base path!
                                            if (main.rRec.its.hasSubsurface()
                                                && (main.rRec.type & RadianceQueryRecord::ESubsurfaceRadiance)) {
                                                SLog(EError, "No BSSRDF support");
                                                shiftedEmitterRadiance += main.rRec.its.LoSub(scene,
                                                                                              shifted.rRec.sampler,
                                                                                              -outgoingDirection,
                                                                                              main.rRec.depth);
                                            }
                                        } else {
                                            SLog(EError, "Not implemented!");
                                        }

                                        // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                                        Float shiftedWeightDenominator =
                                                shiftedPreviousPdf * ((shiftedLumPdf + shiftedPhasePdf));
                                        weight = mainWeightNumerator
                                                 / (D_EPSILON + shiftedWeightDenominator + mainWeightDenominator);

                                        mainContribution = main.throughput * mainEmitterRadiance;
                                        shiftedContribution = shifted.throughput
                                                              *
                                                              shiftedEmitterRadiance; // Note: Jacobian baked into .throughput.
                                    }
                                }
                            } else {
                                // Use half-vector duplication shift. These paths could not have been sampled by light sampling (by our decision).
                                Spectrum shiftedEmitterRadiance(Float(0));
#if ORIGINAL_TESTS
                                SAssert(fabs(shifted.ray.d.lengthSquared() - 1) < 0.000001);
#else
                                if (fabs(shifted.ray.d.lengthSquared() - 1) >= 0.000001) {
                                    shifted.alive = false;
                                    goto half_vector_shift_failed_volume;
                                }
#endif
                                {
                                    // Apply the local shift.
                                    HalfVectorShiftResult shiftResult = halfVectorShiftVolume(mainPhaseResult.pRec.wi,
                                                                                              mainPhaseResult.pRec.wo,
                                                                                              -shifted.ray.d);

                                    if (shiftResult.success) {
                                        // Invertible shift, success.
                                        shifted.throughput *= shiftResult.jacobian;
                                        shifted.pdf *= shiftResult.jacobian;
                                    } else {
                                        // The shift is non-invertible so kill it.
                                        shifted.alive = false;
                                        goto half_vector_shift_failed_volume;
                                    }

                                    {
                                        // FIXME: Maybe the phase function is assumed to be constant
                                        PhaseFunctionSamplingRecord pRec(previousMainMediumRec,
                                                                         -shifted.ray.d,
                                                                         shiftResult.wo);

                                        // Evaluate the phase function
                                        shifted.throughput *= mainPhase->eval(pRec) / mainPhasePdf;
                                        shifted.pdf *= mainPhase->pdf(pRec);

                                        if (shifted.pdf == Float(0)) {
                                            // Offset path is invalid!
                                            shifted.alive = false;
                                            goto half_vector_shift_failed_volume;
                                        }

                                        // Update the vertex type.
                                        // as we are inside the volume, we always assume diffuse vertex type
                                        shiftedVertexType = VERTEX_TYPE_DIFFUSE;

                                        // Trace the next hit point.
                                        shifted.ray = Ray(shifted.mRec.p, shiftResult.wo, main.ray.time);
                                        if (!scene->rayIntersect(shifted.ray, shifted.rRec.its)) {
                                            // Hit nothing - Evaluate environment radiance.
                                            const Emitter *env = scene->getEnvironmentEmitter();
                                            if (!env) {
                                                // Since base paths that hit nothing are not shifted, we must be symmetric and kill shifts that hit nothing.
                                                shifted.alive = false;
                                                goto half_vector_shift_failed_volume;
                                            }
                                            if (main.rRec.its.isValid()) {
                                                // Deny shifts between env and non-env.
                                                shifted.alive = false;
                                                goto half_vector_shift_failed_volume;
                                            }

                                            SLog(EError, "No env map support");
                                            postponedShiftEnd = true;
                                        } else {
                                            // Hit something.

                                            if (!main.rRec.its.isValid()) {
                                                // Deny shifts between env and non-env.
                                                shifted.alive = false;
                                                goto half_vector_shift_failed_volume;
                                            }

                                            // FIXME: Seems ok
                                            VertexType shiftedNextVertexType = getVertexType(shifted, *m_config,
                                                                                             BSDF::EGlossy);

                                            // Make sure that the reverse shift would use this same strategy!
                                            // ==============================================================
                                            if (mainVertexType == VERTEX_TYPE_DIFFUSE
                                                && shiftedVertexType == VERTEX_TYPE_DIFFUSE
                                                && shiftedNextVertexType == VERTEX_TYPE_DIFFUSE) {
                                                // Non-invertible shift: the reverse-shift would use another strategy!
                                                shifted.alive = false;
                                                goto half_vector_shift_failed_volume;
                                            }

                                            if (mediumScattered) {
                                                if (shifted.rRec.its.isEmitter()) {
                                                    // Hit emitter.
                                                    shiftedEmitterRadiance = shifted.rRec.its.Le(-shifted.ray.d);
                                                }
                                                // Sub-surface scattering. Note: Should use the same random numbers as the base path!
                                                if (shifted.rRec.its.hasSubsurface()
                                                    && (shifted.rRec.type & RadianceQueryRecord::ESubsurfaceRadiance)) {
                                                    shiftedEmitterRadiance += shifted.rRec.its.LoSub(scene,
                                                                                                     shifted.rRec.sampler,
                                                                                                     -shifted.ray.d,
                                                                                                     main.rRec.depth);
                                                }
                                            }
                                        }
                                    }
                                }

                                half_vector_shift_failed_volume:
                                if (shifted.alive) {
                                    // Evaluate radiance difference using power heuristic between BSDF samples from base and offset paths.
                                    // Note: No MIS with light sampling since we don't use it for this connection type.
                                    weight = main.pdf / (shifted.pdf + main.pdf);
                                    mainContribution = main.throughput * mainEmitterRadiance;
                                    shiftedContribution = shifted.throughput
                                                          *
                                                          shiftedEmitterRadiance; // Note: Jacobian baked into .throughput.
                                } else {
                                    // Handle the failure without taking MIS with light sampling, as we decided not to use it in the half-vector-duplication case.
                                    // Could have used it, but so far there has been no need. It doesn't seem to be very useful.
                                    weight = Float(1);
                                    mainContribution = main.throughput * mainEmitterRadiance;
                                    shiftedContribution = Spectrum(Float(0));

                                    // Disable the failure detection below since the failure was already handled.
                                    shifted.alive = true;
                                    postponedShiftEnd = true;
                                }
                            }
                        }
                    }

                    shift_failed_volume:
                    if (!shifted.alive) {
                        // The offset path cannot be generated; Set offset PDF and offset throughput to zero.
                        weight = mainWeightNumerator / (D_EPSILON + mainWeightDenominator);
                        mainContribution = main.throughput * mainEmitterRadiance;
                        shiftedContribution = Spectrum((Float) 0);
                    }

                    // Note: Using also the offset paths for the throughput estimate, like we do here, provides some advantage when a large reconstruction alpha is used,
                    // but using only throughputs of the base paths doesn't usually lose by much.
                    if (main.rRec.depth + 1 >= m_config->m_minDepth &&
                        nbSurfaceInteractions >= m_config->m_minCameraDepth) {
#if MEDIUM_INTER_TRANS
#ifndef CENTRAL_RADIANCE
                        main.addRadiance(mainContribution, weight);
                        shifted.addRadiance(shiftedContribution, weight);
#endif
                        shifted.addGradient(shiftedContribution - mainContribution, weight);
#endif
                    }

                    if (postponedShiftEnd) {
                        shifted.alive = false;
                    }
                }

                // Stop if the base path hit the environment.
                main.rRec.type = RadianceQueryRecord::ERadianceNoEmission;
                if ((!sampledMediumEvent && !main.rRec.its.isValid()) ||
                    !(main.rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance)) {
                    break;
                }

            } else {
                /* ==================================================================== */
                /*                     Sampled point on surface                         */
                /* ==================================================================== */

                // Update with the medium weights
                if (med) {
                    if (main.mRec.transmittance.isZero()) {
                        return; // Dead main
                    }
                    main.throughput *= main.mRec.transmittance / main.mRec.pdfFailure;
                    main.pdf *= main.mRec.pdfFailure;

                    for (int i = 0; i < secondaryCount; ++i) {
                        // Evaluate and apply the gradient.
                        RayState &shifted = shiftedRays[i];

                        if (shifted.alive) {
                            if (shifted.connection_status == RAY_CONNECTED) {
                                // We have done the reconnection, as the incomming direction does not matter
                                // We treat CONNECTED and RECENTLY_CONNECTED as the same way
                                // Follow the base path. All relevant vertices are shared.
                                shifted.throughput *= main.mRec.transmittance / main.mRec.pdfFailure;
                                shifted.pdf *= main.mRec.pdfFailure;
                            } else if (shifted.connection_status == RAY_RECENTLY_CONNECTED) {
                                // We recently reconnected the path to the base path
                                // However, the travel distance inside the smoke can be different
                                // compare to the base path, so we need to re-evaluate
                                // based on the previous shifted path position
                                Vector d = main.rRec.its.p - shifted.prevPos();
                                Float distReco = d.length();
                                d /= distReco;

                                Ray shiftRay(shifted.prevPos(), d, 0, distReco, main.mRec.time);
                                MediumSamplingRecord mRecReco; // We use a another medium rec to not erase previous one
                                med->eval(shiftRay, mRecReco);

                                if (mRecReco.transmittance.isZero()) {
                                    // The path is dead
                                    shifted.alive = false;
                                } else {
                                    // Update shifted values
                                    shifted.throughput *= mRecReco.transmittance / main.mRec.pdfFailure;
                                    shifted.pdf *= mRecReco.pdfFailure;
                                }
                            } else {
                                // We did not have reconnected yet the path
                                // We just check that we can copy the distance
                                shifted.ray.maxt = shifted.rRec.its.t;
                                med->eval(shifted.ray, shifted.mRec);
                                if (shifted.mRec.transmittance.isZero()) {
                                    shifted.alive = false;
                                } else {
                                    // FIXME: do we need to use main pdfFailure or shifted one???
                                    // FIXME: Need to derive the equations
                                    shifted.throughput *= shifted.mRec.transmittance / main.mRec.pdfFailure;
                                    shifted.pdf *= shifted.mRec.pdfFailure;
                                    shifted.lastVolume = false;
                                }
                            }
                        } // isAlive
                    }
                } // If med

                if (!mediumScattered) {
                    nbSurfaceInteractions += 1;
                }

                /* ==================================================================== */
                /*                     Direct illumination sampling                     */
                /* ==================================================================== */

                // Sample incoming radiance from lights (next event estimation).
                if (mediumScattered) {
                    const BSDF *mainBSDF = main.rRec.its.getBSDF(main.ray);

                    if (main.rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance
                        && mainBSDF->getType() & BSDF::ESmooth && main.rRec.depth + 1 >= m_config->m_minDepth &&
                        nbSurfaceInteractions >= m_config->m_minCameraDepth) {
                        // Sample an emitter and evaluate f = f/p * p for it. */
                        DirectSamplingRecord dRec(main.rRec.its);

                        mitsuba::Point2 lightSample = main.rRec.nextSample2D();
                        int interactions = 0; //m_config->m_maxDepth - main.rRec.depth - 1;
                        std::pair<Spectrum, bool> emitterTupleValue = scene->sampleAttenuatedEmitterDirectVisible(dRec,
                                                                                                                  med,
                                                                                                                  interactions,
                                                                                                                  lightSample,
                                                                                                                  main.rRec.sampler);

                        Spectrum mainEmitterRadiance = emitterTupleValue.first;
                        bool mainEmitterVisible = emitterTupleValue.second;
                        const auto *emitter = dynamic_cast<const Emitter *>(dRec.object);

                        // If the emitter sampler produces a non-emitter, that's a problem.
                        SAssert(emitter != nullptr);

                        // Add radiance and gradients to the base path and its offset path.
                        // Query the BSDF to the emitter's direction.
                        BSDFSamplingRecord mainBRec(main.rRec.its, main.rRec.its.toLocal(dRec.d), ERadiance);

                        // Evaluate BSDF * cos(theta).
                        Spectrum mainBSDFValue = mainBSDF->eval(mainBRec);

                        // Calculate the probability density of having generated the sampled path segment by BSDF sampling. Note that if the emitter is not visible, the probability density is zero.
                        // Even if the BSDF sampler has zero probability density, the light sampler can still sample it.
                        Float mainBsdfPdf =
                                (emitter->isOnSurface() && dRec.measure == ESolidAngle && mainEmitterVisible)
                                ? mainBSDF->pdf(mainBRec) : 0;

                        // There values are probably needed soon for the Jacobians.
                        Float mainDistanceSquared = (main.rRec.its.p - dRec.p).lengthSquared();
                        Float mainOpposingCosine = dot(dRec.n, (main.rRec.its.p - dRec.p)) / sqrt(mainDistanceSquared);

                        // Power heuristic weights for the following strategies: light sample from base, BSDF sample from base.
                        Float mainWeightNumerator = main.pdf * dRec.pdf;
                        Float mainWeightDenominator = main.pdf * ((dRec.pdf + mainBsdfPdf));

#ifdef CENTRAL_RADIANCE
                        if (main.rRec.depth + 1 >= m_config->m_minDepth && mainWeightDenominator != 0.f && nbSurfaceInteractions >= m_config->m_minCameraDepth) {
#if SURFACE_DIRECT_TRANS
              main.addRadiance(main.throughput * (mainBSDFValue * mainEmitterRadiance),
                               mainWeightNumerator / (mainWeightDenominator));
#endif
            }
#endif

                        // Strict normals check to produce the same results as bidirectional methods when normal mapping is used.
                        if (!m_config->m_strictNormals
                            || dot(main.rRec.its.geoFrame.n, dRec.d) * Frame::cosTheta(mainBRec.wo) > 0) {
                            // The base path is good. Add radiance differences to offset paths.
                            for (int i = 0; i < secondaryCount; ++i) {
                                // Evaluate and apply the gradient.
                                RayState &shifted = shiftedRays[i];

                                Spectrum mainContribution(Float(0));
                                Spectrum shiftedContribution(Float(0));
                                auto weight = Float(0);

                                bool shiftSuccessful = shifted.alive;

                                // Construct the offset path.
                                if (shiftSuccessful) {
                                    // Generate the offset path.
                                    if (shifted.connection_status == RAY_CONNECTED) {
                                        // Follow the base path. All relevant vertices are shared.
                                        Float shiftedBsdfPdf = mainBsdfPdf;
                                        Float shiftedDRecPdf = dRec.pdf;
                                        Spectrum shiftedBsdfValue = mainBSDFValue;
                                        Spectrum shiftedEmitterRadiance = mainEmitterRadiance;
                                        auto jacobian = (Float) 1;

                                        // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                                        Float shiftedWeightDenominator =
                                                (jacobian * shifted.pdf) * ((shiftedDRecPdf + shiftedBsdfPdf));
                                        weight = mainWeightNumerator
                                                 / (D_EPSILON + shiftedWeightDenominator + mainWeightDenominator);

                                        mainContribution = main.throughput * (mainBSDFValue * mainEmitterRadiance);
                                        shiftedContribution = jacobian * shifted.throughput
                                                              * (shiftedBsdfValue * shiftedEmitterRadiance);

                                        // Note: The Jacobians were baked into shifted.pdf and shifted.throughput at connection phase.
                                    } else if (shifted.connection_status == RAY_RECENTLY_CONNECTED) {
                                        // Follow the base path. The current vertex is shared, but the incoming directions differ.
                                        Vector3 incomingDirection = normalize(shifted.prevPos() - main.rRec.its.p);

                                        BSDFSamplingRecord bRec(main.rRec.its,
                                                                main.rRec.its.toLocal(incomingDirection),
                                                                main.rRec.its.toLocal(dRec.d),
                                                                ERadiance);

                                        // Sample the BSDF.
                                        Float shiftedBsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle
                                                                && mainEmitterVisible) ? mainBSDF->pdf(bRec)
                                                                                       : 0; // The BSDF sampler can not sample occluded path segments.
                                        Float shiftedDRecPdf = dRec.pdf;
                                        Spectrum shiftedBsdfValue = mainBSDF->eval(bRec);
                                        Spectrum shiftedEmitterRadiance = mainEmitterRadiance;
                                        auto jacobian = (Float) 1;

                                        // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                                        Float shiftedWeightDenominator =
                                                (jacobian * shifted.pdf) * ((shiftedDRecPdf + shiftedBsdfPdf));
                                        weight = mainWeightNumerator
                                                 / (D_EPSILON + shiftedWeightDenominator + mainWeightDenominator);

                                        mainContribution = main.throughput * (mainBSDFValue * mainEmitterRadiance);
                                        shiftedContribution = jacobian * shifted.throughput
                                                              * (shiftedBsdfValue * shiftedEmitterRadiance);

                                        // Note: The Jacobians were baked into shifted.pdf and shifted.throughput at connection phase.
                                    } else {
                                        // Reconnect to the sampled light vertex. No shared vertices.
                                        SAssert(shifted.connection_status == RAY_NOT_CONNECTED);

                                        const BSDF *shiftedBSDF = shifted.rRec.its.getBSDF(shifted.ray);

                                        // This implementation uses light sampling only for the reconnect-shift.
                                        // When one of the BSDFs is very glossy, light sampling essentially reduces to a failed shift anyway.
                                        bool mainAtPointLight = (dRec.measure == EDiscrete);

                                        // The component have been already choosen before end.
                                        VertexType mainVertexType = getVertexType(main, *m_config, BSDF::ESmooth);
                                        VertexType shiftedVertexType = getVertexType(shifted, *m_config, BSDF::ESmooth);
                                        if (mainAtPointLight || (mainVertexType == VERTEX_TYPE_DIFFUSE
                                                                 && shiftedVertexType == VERTEX_TYPE_DIFFUSE)) {
                                            // Get emitter radiance.
                                            DirectSamplingRecord shiftedDRec(shifted.rRec.its);
                                            std::pair<Spectrum, bool>
                                                    emitterTupleShiftValue = scene->sampleAttenuatedEmitterDirectVisible(
                                                    shiftedDRec,
                                                    med, interactions,
                                                    lightSample,
                                                    shifted.rRec.sampler);
                                            //sampleEmitterDirectVisible(shiftedDRec, lightSample);
                                            bool shiftedEmitterVisible = emitterTupleShiftValue.second;
                                            Spectrum shiftedEmitterRadiance =
                                                    emitterTupleShiftValue.first * shiftedDRec.pdf;
                                            Float shiftedDRecPdf = shiftedDRec.pdf;

                                            // Sample the BSDF.
                                            Float shiftedDistanceSquared =
                                                    (dRec.p - shifted.rRec.its.p).lengthSquared();
                                            Vector emitterDirection =
                                                    (dRec.p - shifted.rRec.its.p) / sqrt(shiftedDistanceSquared);
                                            Float shiftedOpposingCosine = -dot(dRec.n, emitterDirection);

                                            BSDFSamplingRecord bRec(shifted.rRec.its,
                                                                    shifted.rRec.its.toLocal(emitterDirection),
                                                                    ERadiance);

                                            // Strict normals check, to make the output match with bidirectional methods when normal maps are present.
                                            if (m_config->m_strictNormals
                                                && dot(shifted.rRec.its.geoFrame.n, emitterDirection)
                                                   * Frame::cosTheta(bRec.wo) < 0) {
                                                // Invalid, non-samplable offset path.
                                                shiftSuccessful = false;
                                            } else {
                                                Spectrum shiftedBsdfValue = shiftedBSDF->eval(bRec);
                                                Float shiftedBsdfPdf =
                                                        (emitter->isOnSurface() && dRec.measure == ESolidAngle
                                                         && shiftedEmitterVisible) ? shiftedBSDF->pdf(bRec) : 0;
                                                Float jacobian = std::abs(shiftedOpposingCosine * mainDistanceSquared)
                                                                 / (Epsilon + std::abs(
                                                        mainOpposingCosine * shiftedDistanceSquared));

                                                // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                                                Float shiftedWeightDenominator =
                                                        (jacobian * shifted.pdf) * ((shiftedDRecPdf + shiftedBsdfPdf));
                                                weight = mainWeightNumerator / (D_EPSILON + shiftedWeightDenominator
                                                                                + mainWeightDenominator);

                                                mainContribution =
                                                        main.throughput * (mainBSDFValue * mainEmitterRadiance);
                                                if (weight != 0.f) {
                                                    shiftedContribution = jacobian * shifted.throughput
                                                                          *
                                                                          (shiftedBsdfValue * shiftedEmitterRadiance) /
                                                                          dRec.pdf;
                                                }
                                            }
                                        }
                                    }
                                }

                                if (!shiftSuccessful) {
                                    // The offset path cannot be generated; Set offset PDF and offset throughput to zero. This is what remains.

                                    // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset. (Offset path has zero PDF)
                                    weight = mainWeightNumerator / (D_EPSILON + mainWeightDenominator);

                                    mainContribution = main.throughput * (mainBSDFValue * mainEmitterRadiance);
                                    shiftedContribution = Spectrum((Float) 0);
                                }

                                // Note: Using also the offset paths for the throughput estimate, like we do here, provides some advantage when a large reconstruction alpha is used,
                                // but using only throughputs of the base paths doesn't usually lose by much.
                                if (main.rRec.depth + 1 >= m_config->m_minDepth &&
                                    nbSurfaceInteractions >= m_config->m_minCameraDepth) {
#if SURFACE_DIRECT_TRANS
#ifndef CENTRAL_RADIANCE
                                    main.addRadiance(mainContribution, weight);
                                    shifted.addRadiance(shiftedContribution, weight);
#endif
                                    shifted.addGradient(shiftedContribution - mainContribution, weight);
#endif
                                }
                            } // for(int i = 0; i < secondaryCount; ++i)
                        } // Strict normals
                    }
                } // Sample incoming radiance from lights.

                /* ==================================================================== */
                /*               BSDF sampling and emitter hits                         */
                /* ==================================================================== */

                // Sample a new direction from BSDF * cos(theta).
                BSDFSampleResult mainBsdfResult = sampleBSDF(main);
                if (mainBsdfResult.pdf <= (Float) 0.0) {
                    // Impossible base path.
                    break;
                }

                const Vector mainWo = main.rRec.its.toWorld(mainBsdfResult.bRec.wo);

                // Prevent light leaks due to the use of shading normals.
                Float mainWoDotGeoN = dot(main.rRec.its.geoFrame.n, mainWo);
                if (m_config->m_strictNormals && mainWoDotGeoN * Frame::cosTheta(mainBsdfResult.bRec.wo) <= 0) {
                    break;
                }

                // The old intersection structure is still needed after main.rRec.its gets updated.
                Intersection previousMainIts = main.rRec.its;

                // Trace a ray in the sampled direction.
                bool mainHitEmitter = false;
                Spectrum mainEmitterRadiance = Spectrum((Float) 0);

                DirectSamplingRecord mainDRec(main.rRec.its);
                const BSDF *mainBSDF = main.rRec.its.getBSDF(main.ray);

                // Check if we are in good med, and update it
                if (main.rRec.its.isMediumTransition()) {
                    med = main.rRec.its.getTargetMedium(mainWo);
                }

                // Update the vertex types.
                VertexType mainVertexType = getVertexType(main, *m_config, mainBsdfResult.bRec.sampledType);
                VertexType mainNextVertexType;

                // Handle environment connection as diffuse (that's ~infinitely far away).
                main.ray = Ray(main.rRec.its.p, mainWo, main.ray.time);
                if (scene->rayIntersect(main.ray, main.rRec.its)) {
                    if (med) {
                        sampledMediumEvent = med->sampleDistance(Ray(main.ray, 0, main.rRec.its.t), main.mRec,
                                                                 main.rRec.sampler);
                    } else {
                        sampledMediumEvent = false;
                    }

                    // Intersected something - check if it was a luminaire.
                    if (mediumScattered) {
                        if (main.rRec.its.isEmitter()) {
                            mainEmitterRadiance = main.rRec.its.Le(-main.ray.d);
                            // Count for the transmittance
                            if (med) {
                                mainEmitterRadiance *= med->evalTransmittance(Ray(main.ray, 0, main.rRec.its.t),
                                                                              main.rRec.sampler);
                            }
                            mainDRec.setQuery(main.ray, main.rRec.its);
                            mainHitEmitter = true;
                        }

                        // Sub-surface scattering.
                        if (main.rRec.its.hasSubsurface() &&
                            (main.rRec.type & RadianceQueryRecord::ESubsurfaceRadiance)) {
                            mainEmitterRadiance +=
                                    main.rRec.its.LoSub(scene, main.rRec.sampler, -main.ray.d, main.rRec.depth);
                        }
                    }

                    // Update the vertex type for the next surface vertex
                    mainNextVertexType = VERTEX_TYPE_DIFFUSE;
                    if (!sampledMediumEvent) {
                        // We have reach a surface, we need to check that we do not have any
                        // component where we need to bounce on. If it is the case,
                        // we cannot connect to it ...
                        if (main.rRec.its.isValid()) {
                            mainNextVertexType = getVertexType(main, *m_config, mainBsdfResult.bRec.sampledType);
                        }
                    } else {
                        // We sample an event inside the volume
                        mainNextVertexType = VERTEX_TYPE_DIFFUSE;
                    }
                } else {
                    sampledMediumEvent = false;

                    // Intersected nothing -- perhaps there is an environment map?
                    const Emitter *env = scene->getEnvironmentEmitter();

                    if (env) {
                        SLog(EError, "No support of env map");
                        // Hit the environment map.
                        mainEmitterRadiance = env->evalEnvironment(main.ray);
                        if (!env->fillDirectSamplingRecord(mainDRec, main.ray))
                            break;
                        mainHitEmitter = true;
                        return;
                    } else {
                        // Nothing to do anymore.
                        break;
                    }
                }

                // Continue the shift.
                Float mainBsdfPdf = mainBsdfResult.pdf;
                Float mainPreviousPdf = main.pdf;

                main.throughput *= mainBsdfResult.weight;
                main.pdf *= mainBsdfResult.pdf;
                main.eta *= mainBsdfResult.bRec.eta;

                // Compute the probability density of generating base path's direction using the implemented direct illumination sampling technique.
                const Float mainLumPdf = (mainHitEmitter && main.rRec.depth + 1 >= m_config->m_minDepth &&
                                          nbSurfaceInteractions >= m_config->m_minCameraDepth
                                          && !(mainBsdfResult.bRec.sampledType & BSDF::EDelta)) ?
                                         scene->pdfEmitterDirect(mainDRec) : 0;

                // Power heuristic weights for the following strategies: light sample from base, BSDF sample from base.
                Float mainWeightNumerator = mainPreviousPdf * mainBsdfResult.pdf;
                Float mainWeightDenominator = mainPreviousPdf * ((mainLumPdf + mainBsdfPdf));

#ifdef CENTRAL_RADIANCE
                if (main.rRec.depth + 1 >= m_config->m_minDepth && mainWeightDenominator != 0.f && nbSurfaceInteractions >= m_config->m_minCameraDepth) {
#if SURFACE_INTER_TRANS
          main.addRadiance(main.throughput * mainEmitterRadiance, mainWeightNumerator / (mainWeightDenominator));
#endif
        }
#endif

                // Construct the offset paths and evaluate emitter hits.

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
                            Spectrum shiftedBsdfValue = mainBsdfResult.weight;
                            Float shiftedBsdfPdf = mainBsdfPdf;
                            Float shiftedLumPdf = mainLumPdf;
                            Spectrum shiftedEmitterRadiance = mainEmitterRadiance;

                            // Update throughput and pdf.
                            shifted.throughput *= shiftedBsdfValue;
                            shifted.pdf *= shiftedBsdfPdf;

                            // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                            Float shiftedWeightDenominator = shiftedPreviousPdf * ((shiftedLumPdf + shiftedBsdfPdf));
                            weight = mainWeightNumerator
                                     / (D_EPSILON + shiftedWeightDenominator + mainWeightDenominator);

                            mainContribution = main.throughput * mainEmitterRadiance;
                            shiftedContribution = shifted.throughput
                                                  * shiftedEmitterRadiance; // Note: Jacobian baked into .throughput.
                        } else if (shifted.connection_status == RAY_RECENTLY_CONNECTED) {
                            // Recently connected - follow the base path but evaluate BSDF to the new direction.
                            Vector3 incomingDirection = normalize(shifted.prevPos() - main.ray.o);
                            BSDFSamplingRecord bRec(previousMainIts,
                                                    previousMainIts.toLocal(incomingDirection),
                                                    previousMainIts.toLocal(main.ray.d),
                                                    ERadiance);

                            // Note: mainBSDF is the BSDF at previousMainIts, which is the current position of the offset path.
                            EMeasure measure =
                                    (mainBsdfResult.bRec.sampledType & BSDF::EDelta) ? EDiscrete : ESolidAngle;

                            Spectrum shiftedBsdfValue = mainBSDF->eval(bRec, measure);
                            Float shiftedBsdfPdf = mainBSDF->pdf(bRec, measure);
                            shiftedBsdfPdf *= mainBSDF->pdfComponent(bRec);

                            Float shiftedLumPdf = mainLumPdf;
                            Spectrum shiftedEmitterRadiance = mainEmitterRadiance;

                            // Update throughput and pdf.
                            shifted.throughput *= shiftedBsdfValue / mainBsdfPdf;
                            shifted.pdf *= shiftedBsdfPdf;

                            shifted.connection_status = RAY_CONNECTED;

                            // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                            Float shiftedWeightDenominator = shiftedPreviousPdf * ((shiftedLumPdf + shiftedBsdfPdf));
                            weight = mainWeightNumerator
                                     / (D_EPSILON + shiftedWeightDenominator + mainWeightDenominator);

                            mainContribution = main.throughput * mainEmitterRadiance;
                            shiftedContribution = shifted.throughput
                                                  * shiftedEmitterRadiance; // Note: Jacobian baked into .throughput.
                        } else {
                            // Not connected - apply either reconnection or half-vector duplication shift.
                            const BSDF *shiftedBSDF = shifted.rRec.its.getBSDF(shifted.ray);

                            // Check that the two path (main and shifted) are inside the same volume
                            if (shifted.rRec.its.isMediumTransition()) {
                                const Medium *shiftMed = shifted.rRec.its.getTargetMedium(shifted.ray.d);
                                if (med != shiftMed) {
                                    // The ray is dead due to medium inconsistency.
                                    shifted.alive = false;
                                    goto shift_failed;
                                }
                            }

                            VertexType shiftedVertexType = getVertexType(shifted, *m_config,
                                                                         mainBsdfResult.bRec.sampledType);

                            if (mainVertexType == VERTEX_TYPE_DIFFUSE && mainNextVertexType == VERTEX_TYPE_DIFFUSE
                                && shiftedVertexType == VERTEX_TYPE_DIFFUSE) {
                                // Use reconnection shift.

                                // Optimization: Skip the last raycast and BSDF evaluation for the offset path when it won't contribute and isn't needed anymore.
                                if (!lastSegment || mainHitEmitter || main.rRec.its.hasSubsurface()) {
                                    ReconnectionShiftResult shiftResult;
                                    if (main.rRec.its.isValid()) {
                                        // This is an actual reconnection shift.
                                        if (!sampledMediumEvent) {
                                            shiftResult = reconnectShift(scene,
                                                                         main.ray.o,
                                                                         main.rRec.its.p,
                                                                         shifted.rRec.its.p,
                                                                         main.rRec.its.geoFrame.n,
                                                                         main.ray.time);
                                        } else {
                                            shiftResult = reconnectShift(scene,
                                                                         main.ray.o, // Previous intersection
                                                                         main.mRec.p, // New intersection
                                                                         shifted.rRec.its.p, // Previous intersection shift
                                                                         main.ray.d, // Normal ignored here
                                                                         main.ray.time, true);
                                        }
                                    } else {
                                        SLog(EError, "Do not support env map");
                                        shifted.alive = false;
                                        goto shift_failed;
                                        // This is a reconnection at infinity in environment direction.
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
                                                            shifted.rRec.its.toLocal(outgoingDirection),
                                                            ERadiance);

                                    // Strict normals check.
                                    if (m_config->m_strictNormals && dot(outgoingDirection, shifted.rRec.its.geoFrame.n)
                                                                     * Frame::cosTheta(bRec.wo) <= 0) {
                                        shifted.alive = false;
                                        goto shift_failed;
                                    }

                                    // Evaluate the BRDF to the new direction.
                                    Spectrum shiftedBsdfValue = shiftedBSDF->eval(bRec);
                                    Float shiftedBsdfPdf = shiftedBSDF->pdf(bRec);
                                    shiftedBsdfPdf *= shiftedBSDF->pdfComponent(bRec);

                                    // Update throughput and pdf.
                                    shifted.throughput *= shiftedBsdfValue * shiftResult.jacobian / mainBsdfPdf;
                                    shifted.pdf *= shiftedBsdfPdf * shiftResult.jacobian;
                                    shifted.connection_status = RAY_RECENTLY_CONNECTED;

                                    if (mainHitEmitter || main.rRec.its.hasSubsurface()) {
                                        // Also the offset path hit the emitter, as visibility was checked at reconnectShift or environmentShift.

                                        // Evaluate radiance to this direction.
                                        Spectrum shiftedEmitterRadiance(Float(0));
                                        auto shiftedLumPdf = Float(0);

                                        if (main.rRec.its.isValid()) {
                                            // Hit an object.
                                            if (mainHitEmitter) {
                                                shiftedEmitterRadiance = main.rRec.its.Le(-outgoingDirection);

                                                // Evaluate the light sampling PDF of the new segment.
                                                DirectSamplingRecord shiftedDRec;
                                                shiftedDRec.p = mainDRec.p;
                                                shiftedDRec.n = mainDRec.n;
                                                shiftedDRec.dist = (mainDRec.p - shifted.rRec.its.p).length();
                                                shiftedDRec.d = (mainDRec.p - shifted.rRec.its.p) / shiftedDRec.dist;
                                                shiftedDRec.ref = mainDRec.ref;
                                                shiftedDRec.refN = shifted.rRec.its.shFrame.n;
                                                shiftedDRec.object = mainDRec.object;
                                                shiftedDRec.measure = ESolidAngle;

                                                // Add the transmittance
                                                if (med) {
                                                    shiftedEmitterRadiance *= med->evalTransmittance(
                                                            Ray(shifted.rRec.its.p,
                                                                shiftedDRec.d,
                                                                0,
                                                                shiftedDRec.dist,
                                                                shifted.mRec.time),
                                                            shifted.rRec.sampler);
                                                }
                                                shiftedLumPdf = scene->pdfEmitterDirect(shiftedDRec);
                                            }

                                            // Sub-surface scattering. Note: Should use the same random numbers as the base path!
                                            if (main.rRec.its.hasSubsurface()
                                                && (main.rRec.type & RadianceQueryRecord::ESubsurfaceRadiance)) {
                                                shiftedEmitterRadiance += main.rRec.its.LoSub(scene,
                                                                                              shifted.rRec.sampler,
                                                                                              -outgoingDirection,
                                                                                              main.rRec.depth);
                                            }
                                        } else {
                                            // Hit the environment.
                                            shiftedEmitterRadiance = mainEmitterRadiance;
                                            shiftedLumPdf = mainLumPdf;
                                        }

                                        // Power heuristic between light sample from base, BSDF sample from base, light sample from offset, BSDF sample from offset.
                                        Float shiftedWeightDenominator =
                                                shiftedPreviousPdf * ((shiftedLumPdf + shiftedBsdfPdf));
                                        weight = mainWeightNumerator
                                                 / (D_EPSILON + shiftedWeightDenominator + mainWeightDenominator);

                                        mainContribution = main.throughput * mainEmitterRadiance;
                                        shiftedContribution = shifted.throughput
                                                              *
                                                              shiftedEmitterRadiance; // Note: Jacobian baked into .throughput.
                                    }
                                }
                            } else {
                                // Use half-vector duplication shift. These paths could not have been sampled by light sampling (by our decision).
                                Vector3 tangentSpaceIncomingDirection = shifted.rRec.its.toLocal(-shifted.ray.d);
                                Vector3 tangentSpaceOutgoingDirection;
                                Spectrum shiftedEmitterRadiance(Float(0));

                                // Deny shifts between Dirac and non-Dirac BSDFs.
                                bool bothDelta = (mainBsdfResult.bRec.sampledType & BSDF::EDelta)
                                                 && (shiftedBSDF->getType() & BSDF::EDelta);
                                bool bothSmooth = (mainBsdfResult.bRec.sampledType & BSDF::ESmooth)
                                                  && (shiftedBSDF->getType() & BSDF::ESmooth);
                                if (!(bothDelta || bothSmooth)) {
                                    shifted.alive = false;
                                    goto half_vector_shift_failed;
                                }

#if ORIGINAL_TESTS
                                SAssert(fabs(shifted.ray.d.lengthSquared() - 1) < 0.000001);
#else
                                if (fabs(shifted.ray.d.lengthSquared() - 1) >= 0.000001) {
                                    shifted.alive = false;
                                    goto half_vector_shift_failed;
                                }
#endif
                                {
                                    // Apply the local shift.
                                    HalfVectorShiftResult shiftResult = halfVectorShift(mainBsdfResult.bRec.wi,
                                                                                        mainBsdfResult.bRec.wo,
                                                                                        shifted.rRec.its.toLocal(
                                                                                                -shifted.ray.d),
                                                                                        mainBSDF->getEta(),
                                                                                        shiftedBSDF->getEta());

                                    if (mainBsdfResult.bRec.sampledType & BSDF::EDelta) {
                                        // Dirac delta integral is a point evaluation - no Jacobian determinant!
                                        shiftResult.jacobian = Float(1);
                                    }

                                    if (shiftResult.success) {
                                        // Invertible shift, success.
                                        shifted.throughput *= shiftResult.jacobian;
                                        shifted.pdf *= shiftResult.jacobian;
                                        tangentSpaceOutgoingDirection = shiftResult.wo;
                                    } else {
                                        // The shift is non-invertible so kill it
                                        shifted.alive = false;
                                        goto half_vector_shift_failed;
                                    }

                                    Vector3 outgoingDirection = shifted.rRec.its.toWorld(tangentSpaceOutgoingDirection);

                                    // Update throughput and pdf.
                                    BSDFSamplingRecord bRec(shifted.rRec.its,
                                                            tangentSpaceIncomingDirection,
                                                            tangentSpaceOutgoingDirection,
                                                            ERadiance);
                                    EMeasure measure = (mainBsdfResult.bRec.sampledType & BSDF::EDelta) ? EDiscrete
                                                                                                        : ESolidAngle;

                                    shifted.throughput *= shiftedBSDF->eval(bRec, measure) / mainBsdfPdf;
                                    shifted.pdf *= shiftedBSDF->pdf(bRec, measure);
                                    shifted.pdf *= shiftedBSDF->pdfComponent(bRec);

                                    if (shifted.pdf == Float(0)) {
                                        // Offset path is invalid!
                                        shifted.alive = false;
                                        goto half_vector_shift_failed;
                                    }

                                    // Strict normals check to produce the same results as bidirectional methods when normal mapping is used.
                                    if (m_config->m_strictNormals && dot(outgoingDirection, shifted.rRec.its.geoFrame.n)
                                                                     * Frame::cosTheta(bRec.wo) <= 0) {
                                        shifted.alive = false;
                                        goto half_vector_shift_failed;
                                    }

                                    // Update the vertex type.
                                    shiftedVertexType = getVertexType(shifted, *m_config,
                                                                      mainBsdfResult.bRec.sampledType);

                                    // Trace the next hit point.
                                    shifted.ray = Ray(shifted.rRec.its.p, outgoingDirection, main.ray.time);
                                    if (!scene->rayIntersect(shifted.ray, shifted.rRec.its)) {
                                        // Hit nothing - Evaluate environment radiance.
                                        const Emitter *env = scene->getEnvironmentEmitter();
                                        if (!env) {
                                            // Since base paths that hit nothing are not shifted, we must be symmetric and kill shifts that hit nothing.
                                            shifted.alive = false;
                                            goto half_vector_shift_failed;
                                        }
                                        if (main.rRec.its.isValid()) {
                                            // Deny shifts between env and non-env.
                                            shifted.alive = false;
                                            goto half_vector_shift_failed;
                                        }

                                        if (mainVertexType == VERTEX_TYPE_DIFFUSE
                                            && shiftedVertexType == VERTEX_TYPE_DIFFUSE) {
                                            // Environment reconnection shift would have been used for the reverse direction!
                                            shifted.alive = false;
                                            goto half_vector_shift_failed;
                                        }

                                        SLog(EError, "No env support");
                                        // The offset path is no longer valid after this path segment.
                                        shiftedEmitterRadiance = env->evalEnvironment(shifted.ray);
                                        postponedShiftEnd = true;
                                    } else {
                                        // Hit something.

                                        if (!main.rRec.its.isValid()) {
                                            // Deny shifts between env and non-env.
                                            shifted.alive = false;
                                            goto half_vector_shift_failed;
                                        }

                                        // FIXME: Need to check which shift we currently dealing
                                        // FIXME: If it is a surface shift, we need to be sure we have the same BSDF
                                        // FIXME: If it is the case, we can copy the base path classification
                                        // FIXME: This might be not optimal but good enough
                                        VertexType shiftedNextVertexType;
                                        if (!sampledMediumEvent) {
                                            shiftedNextVertexType = getVertexType(shifted, *m_config,
                                                                                  mainBsdfResult.bRec.sampledType);
                                            if (!shifted.rRec.its.isEmitter()) {
                                                // We need to check if we have the same BSDF
                                                // If not, we reject the shift
                                                const BSDF *mainNextBSDF = main.rRec.its.getBSDF(main.ray);
                                                const BSDF *shiftNextBSDF = shifted.rRec.its.getBSDF(shifted.ray);
                                                if (mainNextBSDF != shiftNextBSDF) {
                                                    shifted.alive = false;
                                                    goto half_vector_shift_failed;
                                                }
                                            }
                                        } else {
                                            shiftedNextVertexType = VERTEX_TYPE_DIFFUSE;
                                        }
                                        //    getVertexType(shifted, *m_config, mainBsdfResult.bRec.sampledType);

                                        // Make sure that the reverse shift would use this same strategy!
                                        // ==============================================================

                                        if (mainVertexType == VERTEX_TYPE_DIFFUSE
                                            && shiftedVertexType == VERTEX_TYPE_DIFFUSE
                                            && shiftedNextVertexType == VERTEX_TYPE_DIFFUSE) {
                                            // Non-invertible shift: the reverse-shift would use another strategy!
                                            shifted.alive = false;
                                            goto half_vector_shift_failed;
                                        }

                                        if (mediumScattered) {
                                            if (shifted.rRec.its.isEmitter()) {
                                                // Hit emitter.
                                                shiftedEmitterRadiance = shifted.rRec.its.Le(-shifted.ray.d);
                                            }


                                            // Sub-surface scattering. Note: Should use the same random numbers as the base path!
                                            if (shifted.rRec.its.hasSubsurface()
                                                && (shifted.rRec.type & RadianceQueryRecord::ESubsurfaceRadiance)) {
                                                shiftedEmitterRadiance += shifted.rRec.its.LoSub(scene,
                                                                                                 shifted.rRec.sampler,
                                                                                                 -shifted.ray.d,
                                                                                                 main.rRec.depth);
                                            }
                                        }
                                    }
                                }

                                half_vector_shift_failed:
                                if (shifted.alive) {
                                    // Evaluate radiance difference using power heuristic between BSDF samples from base and offset paths.
                                    // Note: No MIS with light sampling since we don't use it for this connection type.
                                    weight = main.pdf / (shifted.pdf + main.pdf);
                                    mainContribution = main.throughput * mainEmitterRadiance;
                                    if (!shiftedEmitterRadiance.isZero()) {
                                        shiftedContribution = shifted.throughput
                                                              *
                                                              shiftedEmitterRadiance; // Note: Jacobian baked into .throughput.
                                    } else {
                                        shiftedContribution = Spectrum(0.f);
                                    }
                                } else {
                                    // Handle the failure without taking MIS with light sampling, as we decided not to use it in the half-vector-duplication case.
                                    // Could have used it, but so far there has been no need. It doesn't seem to be very useful.
                                    weight = Float(1); // main.pdf;
                                    mainContribution = main.throughput * mainEmitterRadiance;
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
                        mainContribution = main.throughput * mainEmitterRadiance;
                        shiftedContribution = Spectrum((Float) 0);
                    }

                    // Note: Using also the offset paths for the throughput estimate, like we do here, provides some advantage when a large reconstruction alpha is used,
                    // but using only throughputs of the base paths doesn't usually lose by much.
                    if (main.rRec.depth + 1 >= m_config->m_minDepth &&
                        nbSurfaceInteractions >= m_config->m_minCameraDepth) {
#if SURFACE_INTER_TRANS
#ifndef CENTRAL_RADIANCE
                        main.addRadiance(mainContribution, weight);
                        shifted.addRadiance(shiftedContribution, weight);
#endif
                        shifted.addGradient(shiftedContribution - mainContribution, weight);
#endif
                    }

                    if (postponedShiftEnd) {
                        shifted.alive = false;
                    }
                }

                // Stop if the base path hit the environment.
                main.rRec.type = RadianceQueryRecord::ERadianceNoEmission;
                if (!main.rRec.its.isValid() || !(main.rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance)) {
                    break;
                }
            }

            if (main.rRec.depth++ >= m_config->m_rrDepth) {
                /* Russian roulette: try to keep path weights equal to one,
                   while accounting for the solid angle compression at refractive
                   index boundaries. Stop with at least some probability to avoid
                   getting stuck (e.g. due to total internal reflection) */

                Float q = std::min((main.throughput / main.pdf).max() * main.eta * main.eta, (Float) 0.95f);
                if (main.rRec.nextSample1D() >= q)
                    break;

                main.pdf *= q;
                for (int i = 0; i < secondaryCount; ++i) {
                    RayState &shifted = shiftedRays[i];
                    shifted.pdf *= q;
                }
            }
        }
    }
};
MTS_NAMESPACE_END

#endif //MITSUBA_ORIGINAL_VOLUME_H
