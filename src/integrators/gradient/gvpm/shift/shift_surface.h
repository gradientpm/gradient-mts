#pragma once

#include "shift_utilities.h"
#include "shift_cameraPath.h"

#include "../gvpm_accel.h"

#ifndef MITSUBA_SHIFT_SURFACE_H
#define MITSUBA_SHIFT_SURFACE_H

MTS_NAMESPACE_BEGIN

/**
 * These structure is to store all the gather point data
 * useful for the gradient computation
 */
struct SurfaceGradientRecord {
  Scene *scene;
  GatherPoint *baseGather;
  std::vector<ShiftGatherPoint> &shiftGPs;

  Spectrum baseFlux;                      // map to total flux defined in GatherPoint struct
  Spectrum shiftedFlux[4];
  Spectrum weightedBaseFlux[4];

  const GPMConfig &config;
  GPMThreadData &thdata;

  SurfaceGradientRecord(Scene *scene, GatherPoint *gp,
                        const GPMConfig &config,
                        GPMThreadData &thdata, std::vector<ShiftGatherPoint> &_shiftGPs) :
      scene(scene), baseGather(gp), shiftGPs(_shiftGPs),
      baseFlux(0.f), config(config), thdata(thdata) {
    for (auto i = 0; i < shiftGPs.size(); ++i) {
      shiftedFlux[i] = Spectrum(0.0f);
      weightedBaseFlux[i] = Spectrum(0.0f);
    }
  }

  bool shiftPhoton(const Path *lt, int currVertex, const GatherPoint *shiftGather,
                   GradientSamplingResult &result);
  bool shiftPhotonDiffuse(GradientSamplingResult &result,
                          const Path *lt, int currVertex,
                          const Intersection &itsProj,
                          const GatherPoint *shiftGather);
  bool shiftPhotonManifold(GradientSamplingResult &result,
                           const Path *lt, int c, int b,
                           const Intersection &itsProj,
                           const GatherPoint *shiftGather);

  void estimateGradEmitter() {
    // Shift
    const Point2 basePixel = baseGather->path.vertex(1)->getSamplePosition();
    const std::array<Point2, 4> pixels = generateOffsetPos(basePixel);

    Vector2i filmSize = scene->getFilm()->getSize();

    for (auto i = 0; i < shiftGPs.size(); ++i) {
      Float miWeight = 1.f;
      shiftGPs[i].generate(scene, thdata.pool, *baseGather, pixels[i], true);
      if (shiftGPs[i].surfaceValid) {
        if (shiftGPs[i].lastIts().isEmitter()) {
          miWeight = 1.f /
              (1.f +
                  (shiftGPs[i].surfInfo().pdf / baseGather->surfInfo().pdf)
                      * shiftGPs[i].surfInfo().jacobian);
          if (miWeight > 1.f) {
            SLog(EError, "MIS problem");
          }

          if ((i == ERight && (int) basePixel.x == filmSize.x - 1) ||
              (i == ETop && (int) basePixel.y == filmSize.y - 1)) {
            miWeight = 1.0f;
          }

          baseGather->shiftedEmitterFlux[i] += shiftGPs[i].currEmission * miWeight;
        }
      } // Shift success

      baseGather->weightedEmitterFlux[i] += baseGather->currEmission * miWeight;
    }

  }

  std::tuple<bool, Spectrum> getSurfacePhotonContrib(const GPhoton &photon,
                                                     const GatherPoint *gp, Float otherRadii = -1) {
    if (config.bsdfInteractionMode != BSDF::EAll &&
        photon.prevComponentType > 0 && !(photon.prevComponentType & config.bsdfInteractionMode)) {
      return std::make_tuple(false, Spectrum(0.f));
    }

    // Test the radius
    // if an other radii is provided, use this one
    Float lengthSqr = (gp->lastIts().p - photon.its.p).lengthSquared();
    if (otherRadii > 0) {
      if ((otherRadii * otherRadii - lengthSqr) < 0)
        return std::make_tuple(false, Spectrum(0.f));
    } else {
      // No other radii is provided, just use GP radii
      if ((gp->radius * gp->radius - lengthSqr) < 0)
        return std::make_tuple(false, Spectrum(0.f));
    }

    Vector photonWi = photon.its.toWorld(photon.its.wi);
    Normal photonNormal = photon.its.geoFrame.n;

#ifndef MTS_NOSHADINGNORMAL
    Float wiDotGeoN = absDot(photonNormal, photonWi);
#endif

    if (dot(photonNormal, gp->lastIts().shFrame.n) < 1e-1f
#ifndef MTS_NOSHADINGNORMAL
        || wiDotGeoN < 1e-2f
#endif
        ) {
      return std::make_tuple(false, Spectrum(0.f));
    }

    // Accumulate the contribution of the photon
    BSDFSamplingRecord bRec(gp->lastIts(), gp->lastIts().toLocal(photonWi), gp->lastIts().wi, EImportance);
    bRec.component = gp->sampledComponent;

#ifdef MTS_NOSHADINGNORMAL
    Spectrum value = photon.weight * (gp->its.getBSDF()->eval(bRec)) / std::abs(Frame::cosTheta(bRec.wo));
#else
    Spectrum value = photon.weight * gp->lastIts().getBSDF()->eval(bRec);
#endif

    value /= gp->pdfComponent;

    if (value.isZero()) {
      // Return true because this is still a sampleable path by photon mapping
      return std::make_tuple(true, Spectrum(0.f));
    }

    // Account for non-symmetry due to shading normals
    // In photon mapping case, woDotGeoN term is cancelled because photon gathering
    // does not involve a cosine term.
#ifndef MTS_NOSHADINGNORMAL
    value *= std::abs(Frame::cosTheta(bRec.wi) / (wiDotGeoN * Frame::cosTheta(bRec.wo)));
#endif

    // Accumulate the results
    return std::make_tuple(true, value);
  }

  inline void operator()(const GPhotonNodeKD &nodePhoton) {
    GPhoton photon;
    nodePhoton.getData().getPhoton(photon);

    int pathLength = baseGather->depth + photon.depth;
    if (pathLength < config.minDepth ||
        (config.maxDepth > 0 && pathLength > config.maxDepth))
      return;

    // Check the type of path
    if ((config.lightingInteractionMode & EMedia2Surf) && (config.lightingInteractionMode & ESurf2Surf)) {
      // In this case, all the lighting effect for the surfaces are enable
      // So nothing to do in this case, just continue the computation
    } else {
      // In this case, only one lighting mode is available (from surface or media only)
      // So we need to check the previous vertex and may cancel the computation in this case
      const PathVertex *vPrev = nodePhoton.getData().lightPath->vertex(nodePhoton.getData().vertexId - 1);
      if (vPrev->isMediumInteraction() && !(config.lightingInteractionMode & EMedia2Surf)) {
        return; // Ignore this path
      }
      if ((vPrev->isSurfaceInteraction() || vPrev->isEmitterSample())
          && !(config.lightingInteractionMode & ESurf2Surf)) {
        return; // Ignore this path
      }
    }

    // Debug one shift by cancel it
    if (config.debugShift != EAllShift && config.debugShift != ENullShift) {
      int b;
      ELightShiftType currShift = getTypeShift(nodePhoton.getData().lightPath,
                                               nodePhoton.getData().vertexId, b);
      if (config.debugShift != currShift) {
        return; // Do not compute the photon contribution
      }
    }

    // Compute photon contribution and shift photons
    std::tuple<bool, Spectrum> rBase = getSurfacePhotonContrib(photon, baseGather);
    if (!std::get<0>(rBase)) {
      // Not a valid base path
      return;
    }

    Spectrum surfaceBaseFlux = baseGather->surfInfo().weight * std::get<1>(rBase);
    baseFlux += surfaceBaseFlux;

    // Shift
    const Point2 basePixel = baseGather->path.vertex(1)->getSamplePosition();
    const std::array<Point2, 4> pixels = generateOffsetPos(basePixel);
    Vector2i filmSize = scene->getFilm()->getSize();

    for (int i = 0; i < shiftGPs.size(); ++i) {
      GradientSamplingResult result;

      // If the gather point are not traced
      // Generate them
      shiftGPs[i].generate(scene, thdata.pool, *baseGather, pixels[i], false);
      if (shiftGPs[i].surfaceValid) {
        shiftPhoton(nodePhoton.getData().lightPath,
                    nodePhoton.getData().vertexId,
                    &shiftGPs[i], result);
      }

      // no reverse shift at right and top corners
      if ((i == ERight && (int) basePixel.x == filmSize.x - 1) ||
          (i == ETop && (int) basePixel.y == filmSize.y - 1)) {
        result.weight = 1.0f;
      }

      // Weight checking
      if (result.weight > 1.0f || result.weight < 0.f) {
        SLog(EError, "Weight invalid: %f", result.weight);
      }

      shiftedFlux[i] += result.weight * result.shiftedFlux;
      weightedBaseFlux[i] += result.weight * surfaceBaseFlux;
    }
  }
};

MTS_NAMESPACE_END

#endif //MITSUBA_SHIFT_SURFACE_H
