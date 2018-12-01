//
// Created by beltegeuse on 10/20/18.
//

#include "shiftmapping.h"

#ifndef MITSUBA_RANDOM_REPLAY_RECONNECT_H
#define MITSUBA_RANDOM_REPLAY_RECONNECT_H

MTS_NAMESPACE_BEGIN

class ShiftMappingRandomReconnection : public ShiftMapping {
private:

public:
  ShiftMappingRandomReconnection(const GradientPathTracerConfig *config) : ShiftMapping(config) {}
  void evaluateReuse(RayState *rays, int secondaryCount,
                     Spectrum &out_veryDirect, int id_main,
                     std::vector<GradientInfo> &gradients) override {
    SLog(EError, "Impossible to pathReuse");
  }

  void evaluate(RayState &main, RayState *shiftedRays, int secondaryCount, Spectrum &out_veryDirect) override {
    /***
     * FIXME: psLight is when we use path reuse. Better to debug this approach using classical GPT
     * FIXME: No MIS for the light sampling / bsdf sampling. Need to check this
     * FIXME: Does the global shift connection is good? Better to base on the material (I guess), no?
     * FIXME: If it is based on the material, need to compute the other probability to keep shifting (as it is a conditional PDF)
     */
    const Scene *scene = main.rRec.scene;

    /*******
     * First hit computation
     */
    // Perform the first ray intersection for the base path (or ignore if the intersection has already been provided).
    main.rRec.rayIntersect(main.ray);
    main.ray.mint = Epsilon;

    // Perform the same first ray intersection for the offset paths.
    for (int i = 0; i < secondaryCount; ++i) {
      RayState &shifted = shiftedRays[i];
      shifted.rRec.rayIntersect(shifted.ray);
      shifted.ray.mint = Epsilon;
    }

    if(m_config->m_minDepth <= 0) {
      if (!main.rRec.its.isValid()) {
        // First hit is not in the scene so can't continue. Also there there are no paths to shift.
        // Add potential very direct light from the environment as gradients are not used for that.
        if (main.rRec.type & RadianceQueryRecord::EEmittedRadiance) {
          out_veryDirect += main.throughput * scene->evalEnvironment(main.ray);
        }
        return;
      } else {
        // Include emitted radiance if requested.
        if (main.rRec.its.isEmitter() && (main.rRec.type & RadianceQueryRecord::EEmittedRadiance)) {
          out_veryDirect += main.throughput * main.rRec.its.Le(-main.ray.d);
        }
        // Include radiance from a subsurface scattering model if requested. Note: Not tested!
        if (main.rRec.its.hasSubsurface() && (main.rRec.type & RadianceQueryRecord::ESubsurfaceRadiance)) {
          out_veryDirect += main.throughput * main.rRec.its.LoSub(scene, main.rRec.sampler, -main.ray.d, 0);
        }
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

    int depth = 1;
    while (depth < m_config->m_maxDepth || m_config->m_maxDepth < 0) {
      // Only direct rendering
      const BSDF *mainBSDF = main.rRec.its.getBSDF(main.ray);
      if (m_config->m_strictNormals) {
        if (dot(main.ray.d, main.rRec.its.geoFrame.n) * Frame::cosTheta(main.rRec.its.wi) >= 0) {
          return;
        }
        for (int i = 0; i < secondaryCount; ++i) {
          RayState &shifted = shiftedRays[i];
          if (dot(shifted.ray.d, shifted.rRec.its.geoFrame.n) * Frame::cosTheta(shifted.rRec.its.wi) >= 0) {
            shifted.alive = false;
          }
        }
      }

#if COHERENT_FRAME
    main.rRec.its.wi = main.rRec.its.toWorld(main.rRec.its.wi);
    {
      main.rRec.its.shFrame = Frame(main.rRec.its.shFrame.n);
      main.rRec.its.wi = main.rRec.its.toLocal(main.rRec.its.wi);
    }
    for (size_t i = 0; i < secondaryCount; i++) {
      RayState &shifted = shiftedRays[i];
      shifted.rRec.its.wi = shifted.rRec.its.toWorld(shifted.rRec.its.wi);
      {
        shifted.rRec.its.shFrame = Frame(shifted.rRec.its.shFrame.n);
        shifted.rRec.its.wi = shifted.rRec.its.toLocal(shifted.rRec.its.wi);
      }
    }
#endif

      /// This code check the status of the shift
      /// To get some error
      if (!main.alive) {
        SLog(EError, "Invalid main (d: %i)", depth);
        return;
      }

      // Compute the shift only if the base path
      // if not dead
      // recMainLight.pdf != 0.f &&
      if(mainBSDF->getType() & BSDF::ESmooth && depth + 1 >= m_config->m_minDepth) {
        // Sample the list
        MainLightSamplingInfo mainInfo(main.rRec.its);
        mainInfo.sample = main.rRec.nextSample2D();
        std::pair<Spectrum, bool> emitterTuple = scene->sampleEmitterDirectVisible(mainInfo.dRec, mainInfo.sample);
        mainInfo.value = emitterTuple.first * mainInfo.dRec.pdf;
        mainInfo.visible = emitterTuple.second;
        const Emitter *emitter = static_cast<const Emitter *>(mainInfo.dRec.object);
        SAssert(emitter != nullptr);

        // Compute BSDF component for direct (and MIS)
        BSDFSamplingRecord mainBRec(main.rRec.its, main.rRec.its.toLocal(mainInfo.dRec.d), ERadiance);
        Spectrum mainBSDFValue = mainBSDF->eval(mainBRec);
        Float mainBsdfPdf = (emitter->isOnSurface() && mainInfo.dRec.measure == ESolidAngle && mainInfo.visible)
                            ? mainBSDF->pdf(mainBRec) : 0;

        // Precomputed values for the Jacobian
        mainInfo.distSquare = (main.rRec.its.p - mainInfo.dRec.p).lengthSquared();
        mainInfo.opCos = dot(mainInfo.dRec.n, (main.rRec.its.p - mainInfo.dRec.p)) / sqrt(mainInfo.distSquare);

        Float mainWeightDenominator = (main.pdf) * (mainInfo.dRec.pdf + mainBsdfPdf);
        Spectrum mainContribution = main.throughput * (mainBSDFValue * mainInfo.value);
        if(depth + 1 >= m_config->m_minDepth && (!m_config->reusePrimal)) {
#if DO_DIRECT_COMPUTATION
          main.addRadiance(mainContribution, 1.0 / (D_EPSILON + mainWeightDenominator));
#endif
        }
        for (size_t i = 0; i < secondaryCount; i++) {
          RayState &shifted = shiftedRays[i];
          auto shift_dRec = shiftLightSampling(shifted, main, mainInfo, mainBSDFValue, mainBsdfPdf, mainBSDF);

          if(depth + 1 >= m_config->m_minDepth) {
#if DO_DIRECT_COMPUTATION
            Float weight = 1.0 / (D_EPSILON + mainWeightDenominator + shift_dRec.weightDenominator);
            if(m_config->reusePrimal) {
              main.addRadiance(mainContribution, weight);
              shifted.addRadiance(shift_dRec.value, weight);
            }
            shifted.addGradient(shift_dRec.value - mainContribution, weight);
#endif
          }
        }
      }

      BSDFSampleResult mainBsdfResult = sampleBSDF(main);
      if (mainBsdfResult.pdf <= (Float) 0.0) {
        // Impossible base path.
        break;
      }

      const Vector mainWo = main.rRec.its.toWorld(mainBsdfResult.bRec.wo);
      Float mainWoDotGeoN = dot(main.rRec.its.geoFrame.n, mainWo);
      if (m_config->m_strictNormals && mainWoDotGeoN * Frame::cosTheta(mainBsdfResult.bRec.wo) <= 0) {
        break;
      }

      MainBSDFSamplingInfo mainInfo(main.rRec.its);
      mainInfo.bsdf_pdf = mainBsdfResult.pdf;
      mainInfo.bsdf_value = mainBsdfResult.weight * mainBsdfResult.pdf;
      mainInfo.sampledType = mainBsdfResult.bRec.sampledType;
      mainInfo.eta = mainBsdfResult.bRec.eta;

      // Update the vertex types and do the next intersection
      mainInfo.current = getVertexType(main, *m_config, mainBsdfResult.bRec.sampledType);
      main.ray = Ray(main.rRec.its.p, mainWo, main.ray.time);
      if (scene->rayIntersect(main.ray, main.rRec.its)) {
        // Intersected something - check if it was a luminaire.
        if (main.rRec.its.isEmitter()) {
          mainInfo.light_value = main.rRec.its.Le(-main.ray.d);
          mainInfo.dRec.setQuery(main.ray, main.rRec.its);
          mainInfo.hitEmitter = true;
          mainInfo.next = VERTEX_TYPE_DIFFUSE;
        } else {
          // FIXME: This code seems very strange ...
          mainInfo.next = getVertexType(main, *m_config, mainBsdfResult.bRec.sampledType);
        }

        // Sub-surface scattering.
        if (main.rRec.its.hasSubsurface() && (main.rRec.type & RadianceQueryRecord::ESubsurfaceRadiance)) {
          mainInfo.light_value += main.rRec.its.LoSub(scene, main.rRec.sampler, -main.ray.d, depth);
        }
      } else {
        // Intersected nothing -- perhaps there is an environment map?
        const Emitter *env = scene->getEnvironmentEmitter();

        if (env) {
          // Hit the environment map.
          mainInfo.light_value = env->evalEnvironment(main.ray);
          if (!env->fillDirectSamplingRecord(mainInfo.dRec, main.ray))
            break;
          mainInfo.hitEmitter = true;

          // Handle environment connection as diffuse (that's ~infinitely far away).
          // Update the vertex type.
          mainInfo.next = VERTEX_TYPE_DIFFUSE;
        } else {
          // Nothing to do anymore.
          break;
        }
      }

      // Continue the shift.
      Float mainBsdfPdf = mainBsdfResult.pdf;
      Float mainPreviousPdf = main.pdf;
      main.throughput *= mainBsdfResult.weight * mainBsdfResult.pdf;
      main.pdf *= mainBsdfResult.pdf;
      main.eta *= mainBsdfResult.bRec.eta;

      // Compute the weight for later
      mainInfo.light_pdf = (mainInfo.hitEmitter && depth + 1 >= m_config->m_minDepth &&
          !(mainBsdfResult.bRec.sampledType & BSDF::EDelta)) ? scene->pdfEmitterDirect(mainInfo.dRec) : 0;
      Spectrum mainContribution = main.throughput * mainInfo.light_value;
      Float mainWeightDenominator = (mainPreviousPdf) * (mainInfo.light_pdf + mainBsdfPdf);
#if DO_BSDF_COMPUTATION
      if(depth + 1 >= m_config->m_minDepth && (!m_config->reusePrimal)) {
        main.addRadiance(mainContribution, 1.0 / (D_EPSILON + mainWeightDenominator));
      }
#endif

      for (int i = 0; i < secondaryCount; ++i) {
        RayState &shifted = shiftedRays[i];
        auto shift_bsdf = shiftBSDFSampling(shifted, main, mainInfo, mainBSDF, mainBsdfResult.sample);

        if(depth + 1 >= m_config->m_minDepth) {
#if DO_BSDF_COMPUTATION
          auto weight = [&]() -> Float {
            if(shift_bsdf.exclude_light_sampling) {
              return 1.0 / (D_EPSILON + main.pdf + shift_bsdf.weightDenominator);
            } else {
              return 1.0 / (D_EPSILON + mainWeightDenominator + shift_bsdf.weightDenominator);
            }
          }();
          if(m_config->reusePrimal) {
            main.addRadiance(mainContribution, weight);
            shifted.addRadiance(shift_bsdf.value, weight);
          }
          shifted.addGradient(shift_bsdf.value - mainContribution, weight);
#endif
        }
      }


      // Stop if the base path hit the environment.
      main.rRec.type = RadianceQueryRecord::ERadianceNoEmission;
      if (!main.rRec.its.isValid() || !(main.rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance)) {
        break;
      }

      if (depth++ >= m_config->m_rrDepth) {
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

#endif //MITSUBA_RANDOM_REPLAY_RECONNECT_H
