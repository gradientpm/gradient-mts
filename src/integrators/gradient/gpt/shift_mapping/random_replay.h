//
// Created by muliana on 10/18/18.
//

#include "shiftmapping.h"

#ifndef MITSUBA_RANDOM_REPLAY_H
#define MITSUBA_RANDOM_REPLAY_H

MTS_NAMESPACE_BEGIN

class ShiftMappingRandom: public ShiftMapping {
private:
    Sampler* m_sampler = nullptr;
    std::vector<Float> m_random = {};
    size_t m_index = 0;

    Float next1D() {
        auto v = [&]() -> Float {
            if(m_index < m_random.size()) {
                return m_random[m_index];
            } else {
                auto v = m_sampler->next1D();
                m_random.push_back(v);
                return v;
            }
        }();
        m_index++;
        return v;
    }

    Point2 next2D() {
        return Point2(next1D(), next1D());
    }

    void reset() {
        m_index = 0;
    }
public:
    ShiftMappingRandom(const GradientPathTracerConfig* config): ShiftMapping(config) {}
    void evaluateReuse(RayState *rays, int secondaryCount,
                       Spectrum &out_veryDirect, int id_main,
                       std::vector<GradientInfo>& gradients) override {
        m_sampler = rays[id_main].rRec.sampler;
        m_random.clear();

        Spectrum base_value = Li(rays[id_main].ray, rays[id_main].rRec);
        Float q = std::min(std::max(base_value.max(), m_config->rrShift), 1.0);
        if (q >= rays[id_main].rRec.sampler->next1D()) {
            std::vector<Spectrum> contributions(secondaryCount, Spectrum(0.f));
            Float w = 1.0 / q;
            w *= 1.0 / secondaryCount;
            rays[id_main].radiance += base_value * w;
            contributions[id_main] = base_value;
            for (int i = 0; i < secondaryCount; ++i) {
                if(i == id_main) continue;
                reset();
                RayState &shifted = rays[i];
                contributions[i] = Li(shifted.ray, shifted.rRec);
                shifted.radiance += contributions[i] * w;
            }

            for (auto &g: gradients) {
                g.value += (contributions[g.pixel_a_index] - contributions[g.pixel_b_index]) * w;
            }

        }
    }
    void evaluate(RayState &main, RayState *shiftedRays, int secondaryCount, Spectrum &out_veryDirect) override {
        m_sampler = main.rRec.sampler;
        m_random.clear();

        Spectrum base_value = Li(main.ray, main.rRec);
        Float q = std::min(std::max(base_value.max(), m_config->rrShift), 1.0);
        if (q >= main.rRec.sampler->next1D()) {
            Float w = 1.0 / q;
            for (int i = 0; i < secondaryCount; ++i) {
                reset();
                RayState &shifted = shiftedRays[i];
                Spectrum shift_value = Li(shifted.ray, shifted.rRec);
                if(m_config->reusePrimal) {
                  shifted.radiance += shift_value * 0.5 * w;
                  main.radiance += base_value * 0.5 * w;
                }
                shifted.gradient += (shift_value - base_value) * 0.5 * w;
            }
        }
        // TODO: Need to check if we put this code to RR or not.
        if(!m_config->reusePrimal) {
          main.radiance += base_value;
        }
    }

    inline Float miWeight(Float pdfA, Float pdfB) {
        pdfA *= pdfA;
        pdfB *= pdfB;
        return pdfA / (pdfA + pdfB);
    }

    // FIXME: Does support very direct
    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) {
        // Duplicate of MIPathTracer::Li to support sub-surface scattering initialization.

        /* Some aliases and local variables */
        const Scene *scene = rRec.scene;
        Intersection &its = rRec.its;
        RayDifferential ray(r);
        Spectrum Li(0.0f);
        bool scattered = false;

        /* Perform the first ray intersection (or ignore if the
            intersection has already been provided). */
        rRec.rayIntersect(ray);
        ray.mint = Epsilon;

        Spectrum throughput(1.0f);
        Float eta = 1.0f;
        bool m_hideEmitters = false; // FIXME
        while (rRec.depth <= m_config->m_maxDepth || m_config->m_maxDepth < 0) {
            if (!its.isValid()) {
                /* If no intersection could be found, potentially return
                    radiance from a environment luminaire if it exists */
                if(m_config->m_minDepth <= 0) {
                  if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
                      && (!m_hideEmitters || scattered))
                    Li += throughput * scene->evalEnvironment(ray);
                }
                break;
            }

            const BSDF *bsdf = its.getBSDF(ray);

            /* Possibly include emitted radiance if requested */
            if(m_config->m_minDepth <= 0) {
              if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
                  && (!m_hideEmitters || scattered))
                Li += throughput * its.Le(-ray.d);

              /* Include radiance from a subsurface scattering model if requested */
              if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance)) {
                SLog(EError, "Subsurface is not supported");
                Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);
              }
            }

            if ((rRec.depth >= m_config->m_maxDepth && m_config->m_maxDepth > 0)
                || (m_config->m_strictNormals && dot(ray.d, its.geoFrame.n)
                                       * Frame::cosTheta(its.wi) >= 0)) {

                /* Only continue if:
                    1. The current path length is below the specifed maximum
                    2. If 'strictNormals'=true, when the geometric and shading
                        normals classify the incident direction to the same side */
                break;
            }

            /* ==================================================================== */
            /*                     Direct illumination sampling                     */
            /* ==================================================================== */

            /* Estimate the direct illumination if this is requested */
            DirectSamplingRecord dRec(its);
            Point2 randSampleLight = next2D();
            if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
                (bsdf->getType() & BSDF::ESmooth)) {
                Spectrum value = scene->sampleEmitterDirect(dRec, randSampleLight);
                if (!value.isZero()) {
                    const auto *emitter = dynamic_cast<const Emitter *>(dRec.object);

                    /* Allocate a record for querying the BSDF */
                    BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);

                    /* Evaluate BSDF * cos(theta) */
                    const Spectrum bsdfVal = bsdf->eval(bRec);

                    /* Prevent light leaks due to the use of shading normals */
                    if (!bsdfVal.isZero() && (!m_config->m_strictNormals
                                              || dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0)) {

                        /* Calculate prob. of having generated that direction
                            using BSDF sampling */
                        Float bsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle)
                                        ? bsdf->pdf(bRec) : 0;

                        /* Weight using the power heuristic */
                        Float weight = miWeight(dRec.pdf, bsdfPdf);
                        if(rRec.depth + 1 > m_config->m_minDepth) {
#if DO_DIRECT_COMPUTATION
                          Li += throughput * value * bsdfVal * weight;
#endif
                        }
                    }
                }
            }

            /* ==================================================================== */
            /*                            BSDF sampling                             */
            /* ==================================================================== */

            /* Sample BSDF * cos(theta) */
            Float bsdfPdf;
            BSDFSamplingRecord bRec(its, nullptr, ERadiance);
            Point2 randomBSDF = next2D();
            Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, randomBSDF);
            if (bsdfWeight.isZero())
                break;

            scattered |= bRec.sampledType != BSDF::ENull;

            /* Prevent light leaks due to the use of shading normals */
            const Vector wo = its.toWorld(bRec.wo);
            Float woDotGeoN = dot(its.geoFrame.n, wo);
            if (m_config->m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
                break;

            bool hitEmitter = false;
            Spectrum value;

            /* Trace a ray in this direction */
            ray = Ray(its.p, wo, ray.time);
            if (scene->rayIntersect(ray, its)) {
                /* Intersected something - check if it was a luminaire */
                if (its.isEmitter()) {
                    value = its.Le(-ray.d);
                    dRec.setQuery(ray, its);
                    hitEmitter = true;
                }
            } else {
                /* Intersected nothing -- perhaps there is an environment map? */
                const Emitter *env = scene->getEnvironmentEmitter();

                if (env) {
                    if (m_hideEmitters && !scattered)
                        break;

                    value = env->evalEnvironment(ray);
                    if (!env->fillDirectSamplingRecord(dRec, ray))
                        break;
                    hitEmitter = true;
                } else {
                    break;
                }
            }

            /* Keep track of the throughput and relative
                refractive index along the path */
            throughput *= bsdfWeight;
            eta *= bRec.eta;

            /* If a luminaire was hit, estimate the local illumination and
                weight using the power heuristic */
            if (hitEmitter &&
                (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
                /* Compute the prob. of generating that direction using the
                    implemented direct illumination sampling technique */
                const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
                                     scene->pdfEmitterDirect(dRec) : 0;
                if(rRec.depth + 1 > m_config->m_minDepth) {
#if DO_BSDF_COMPUTATION
                  Li += throughput * value * miWeight(bsdfPdf, lumPdf);
#endif
                }
            }

            /* ==================================================================== */
            /*                         Indirect illumination                        */
            /* ==================================================================== */

            /* Set the recursive query type. Stop if no surface was hit by the
                BSDF sample or if indirect illumination was not requested */
            if (!its.isValid() || !(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
                break;
            rRec.type = RadianceQueryRecord::ERadianceNoEmission;

            if (rRec.depth++ >= m_config->m_rrDepth) {
                /* Russian roulette: try to keep path weights equal to one,
                    while accounting for the solid angle compression at refractive
                    index boundaries. Stop with at least some probability to avoid
                    getting stuck (e.g. due to total internal reflection) */

                Float q = std::min(throughput.max() * eta * eta, (Float) 0.95f);
                if (next1D() >= q)
                    break;
                throughput /= q;
            }
        }

        return Li;
    }
};

MTS_NAMESPACE_END

#endif //MITSUBA_RANDOM_REPLAY_H
