#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/scene.h>

#include "shift_cameraPath.h"
#include "shift_utilities.h"

#include "../gvpm_beams.h"

// For diffuse shift
#include "operation/shift_diffuse.h"

// For ME shift
#include "operation/shift_ME.h"
#include <mitsuba/bidir/manifold.h>
#include <mitsuba/bidir/mut_manifold.h>

#ifndef MITSUBA_SHIFT_VOLUME_BEAMS_H
#define MITSUBA_SHIFT_VOLUME_BEAMS_H

MTS_NAMESPACE_BEGIN

struct BeamKernelRecord {
  Float radius;
  Float v;            // the sampled distance on photon beam
  Float w;            // the sampled distance on camera ray
  Float pdfKernel;
  Float pdfEdgeFailure;   // short beam probability during photon beam sampling

  // 1D kernel
  Float u; // TODO: Do we need it?

  EVolumeTechnique volTechnique;
  Spectrum beamTrans;
  Spectrum contrib;
  Float weightKernel;

  BeamKernelRecord(const BeamKernelRecord &ori,
                   const Medium *medium,
                   const LTPhotonBeam *beam,
                   const Ray &cameraRay) : radius(ori.radius), u(0.f), volTechnique(ori.volTechnique),
                                           contrib(0.f) {
    if (ori.volTechnique == EBeamBeam1D) {
      eval(medium, nullptr, beam, cameraRay, 0, beam->getLength());
    } else if (volTechnique == EBeamBeam3D_Optimized) {
      // Pre-compute ray form for camera and beam
      // in order to call the intersection procedure
      Ray _cam = Ray(cameraRay(cameraRay.mint), cameraRay.d,
                     0, cameraRay.maxt - cameraRay.mint, 0.f);

      // In this sampling strategy, we consider the beam entirely.
      // However, as a beam can be represented by several subbeam
      // We discard the computation of non valid subbeams (to select only one subbeam at the end)
      Ray _beam = Ray(beam->getOri(), beam->getDir(),
                      0.f, beam->getLength(), 0.f); //TODO: Make a function

      // To compute the PDF
      double tNearBeam, tFarBeam;
      if (!cylinderIntersection(_cam, _beam, radius, tNearBeam, tFarBeam)) {
        return; // No intersection
      }

      // We didn't move the light path sampling point
      // So we can reuse it
      v = ori.v;
      pdfKernel = 1.0 / std::max(tFarBeam - tNearBeam, 0.0001);
      if (v < 0 || v > beam->getLength()) {
        return;
      }

      // Optimized: sample the point over the camera beam
      // on the overlapping region
      Point kernelCentroid = beam->getPos(v);

      // Find the projection of the photon position to the camera ray
      Float distToProj = dot(kernelCentroid - cameraRay.o, cameraRay.d);
      Float distSqr = (cameraRay(distToProj) - kernelCentroid).lengthSquared();
      Float radSqr = radius * radius;
      if (distSqr >= radSqr) {
        SLog(EError, "Instersection problem");
        return;  // ray does not intersect the kernel
      }

      // Sample camera segment [tc-, tc+] (uniform sampling)
      Float deltaT = math::safe_sqrt(radSqr - distSqr);    // equal to half of [tc-, tc+] range
      w = ori.w;
      pdfKernel *= 1.0 / std::max(2.0 * deltaT, 0.0001);

      // Check if it is in a acceptable distance for the camera
      if (w < cameraRay.mint || w > cameraRay.maxt) {
        SLog(EError, "Intersection problem (camera beam)");
        return;
      }

      // Photon beam evaluation
      // Note that in the null shift strategy, we keep w and v constant
      // So we can just copy the value of the ori record
      // This only works for homogenous media
#if 0
      Ray rayTrans(beam->getOri(), beam->getDir(), 0.f);
      rayTrans.mint = 0.f;
      rayTrans.maxt = v;
      MediumSamplingRecord mRecBeam;
      medium->eval(rayTrans, mRecBeam);

      // Camera evaluation
      Ray cameraRayEval(cameraRay);
      cameraRayEval.mint = 0.f;                     // contribute to the entire camera segment
      cameraRayEval.maxt = w;
      MediumSamplingRecord mRecCamera;
      medium->eval(cameraRayEval, mRecCamera);

      PhaseFunctionSamplingRecord pRec(mRecCamera, -beam->getDir(), -cameraRay.d, EImportance);
      Float phaseTerm = mRecCamera.getPhaseFunction()->eval(pRec);

      Float kernelVol = ((4.0 / 3.0) * M_PI * std::pow(radius, 3));

      contrib = beam->flux * mRecBeam.transmittance *
        mRecCamera.sigmaS * mRecCamera.transmittance * phaseTerm / pdfKernel;  // not yet include kernel weight
      weightKernel = 1.0 / kernelVol;

      if (!beam->longBeams) {
          contrib /= mRecBeam.pdfFailure;
          pdfEdgeFailure = mRecBeam.pdfFailure;
      } else {
          pdfEdgeFailure = 1.0f;
      }
#else
      contrib = ori.contrib * (ori.pdfKernel / pdfKernel);
      weightKernel = ori.weightKernel;
      beamTrans = ori.beamTrans;

      //FIXME: ASSERT if it is a heterogenous media (as we reneed to evaluate the transmittance)

      if (!beam->longBeams) {
        pdfEdgeFailure = ori.pdfEdgeFailure;
      } else {
        pdfEdgeFailure = 1.0f;
      }
#endif
    }
  }

  BeamKernelRecord(EVolumeTechnique _volTechnique, const Medium *medium, Sampler *sampler,
                   const LTPhotonBeam *beam, Float _radius,
                   const Ray &cameraRay,
                   Float tmin, Float tmax) : radius(_radius), v(0.f), w(0.f),
                                             pdfKernel(0.f), pdfEdgeFailure(0.f), u(0.f), volTechnique(_volTechnique),
                                             contrib(0.f), weightKernel(0.f) {
    eval(medium, sampler, beam, cameraRay, tmin, tmax);
  }

  // Do the evaluation
  void eval(const Medium *medium, Sampler *sampler,
            const LTPhotonBeam *beam,
            const Ray &cameraRay,
            Float tmin, Float tmax) {
    if (beam->isInvalid()) {
      return;
    }

    // TODO: See this correction
    // We might want to not restrict to maxDepth
    if (tmax > beam->getLength()) {
      tmax = beam->getLength();
    }

    if (volTechnique == EBeamBeam1D) {
      // Also update kRec when checking intersection
      if (!beam->rayIntersect1D(cameraRay,
                                tmin,
                                tmax,
                                u,
                                v,
                                w,
                                pdfKernel)) {
        return;
      }

      // Evaluate the medium transmittance along the beam ray
      Ray cameraRayEval(cameraRay);
      cameraRayEval.mint = 0.f;
      cameraRayEval.maxt = w;
      MediumSamplingRecord mRecCamera;
      medium->eval(cameraRayEval, mRecCamera);

      weightKernel = 0.5f / radius;

      pdfEdgeFailure = 0.f;
      contrib = beam->getContrib(v, mRecCamera, cameraRayEval.d, beamTrans, pdfEdgeFailure);
      if (!contrib.isZero()) {
        contrib /= pdfKernel;
      }
    } else if (volTechnique == EBeamBeam3D_Optimized) {
      // Pre-compute ray form for camera and beam
      // in order to call the intersection procedure
      Ray _cam = Ray(cameraRay(cameraRay.mint), cameraRay.d,
                     0, cameraRay.maxt - cameraRay.mint, 0.f);

      // In this sampling strategy, we consider the beam entirely.
      // However, as a beam can be represented by several subbeam
      // We discard the computation of non valid subbeams (to select only one subbeam at the end)
      Ray _beam = Ray(beam->getOri(), beam->getDir(),
                      0.f, beam->getLength(), 0.f); //TODO: Make a function

      // 1) Compute tNear and tFar of the beam to camera "beam" (Cylinder intersection)
      double tNearBeam, tFarBeam;
      if (!cylinderIntersection(_cam, _beam, radius, tNearBeam, tFarBeam)) {
        return; // No intersection
      }

      // 2) Check if it is the correct sub-beam (to avoid to integrate twice)
      if (tNearBeam < 0 && tmin <= Epsilon) {
        // Acceptable, it is the first subbeam
      } else if (tNearBeam > tmin && tNearBeam < tmax) {
        // The tNearBeam need to be inside the subbeam
      } else {
        return; // Another subbeam will take care of it
      }


      // 3) Sample a point on the photon beam
      v = tNearBeam + (tFarBeam - tNearBeam) * sampler->next1D();
      pdfKernel = 1.0 / std::max(tFarBeam - tNearBeam, 0.0001);
      if (v < 0 || v > beam->getLength()) {
        return;
      }

      // 4) Determine the camera position
      // Optimized: sample the point over the camera beam
      // on the overlapping region
      Point kernelCentroid = beam->getPos(v);

      // Find the projection of the photon position to the camera ray
      Float distToProj = dot(kernelCentroid - cameraRay.o, cameraRay.d);
      Float distSqr = (cameraRay(distToProj) - kernelCentroid).lengthSquared();
      Float radSqr = radius * radius;
      if (distSqr >= radSqr) {
        return;  // ray does not intersect the kernel
      }

      // Sample camera segment [tc-, tc+] (uniform sampling)
      Float deltaT = math::safe_sqrt(radSqr - distSqr);    // equal to half of [tc-, tc+] range
      w = distToProj - deltaT + 2 * deltaT * sampler->next1D(); // the random point is always inside the kernel
      pdfKernel *= 1.0 / std::max(2.0 * deltaT, 0.0001);


      // Check if it is in a acceptable distance for the camera
      if (w < cameraRay.mint || w > cameraRay.maxt) {
        return;
      }

      // Photon beam evaluation
      // Need to do it after rejection
      Ray rayTrans(beam->getOri(), beam->getDir(), 0.f); //TODO: Make a function
      rayTrans.mint = 0.f;
      rayTrans.maxt = v;
      MediumSamplingRecord mRecBeam;
      medium->eval(rayTrans, mRecBeam);

      // Camera evaluation
      Ray cameraRayEval(cameraRay);
      cameraRayEval.mint = 0.f;                     // contribute to the entire camera segment
      cameraRayEval.maxt = w;
      MediumSamplingRecord mRecCamera;
      medium->eval(cameraRayEval, mRecCamera);

      PhaseFunctionSamplingRecord pRec(mRecCamera, -beam->getDir(), -cameraRay.d, EImportance);
      Float phaseTerm = mRecCamera.getPhaseFunction()->eval(pRec);

      Float kernelVol = ((4.0 / 3.0) * M_PI * std::pow(radius, 3));

      contrib = beam->flux * mRecBeam.transmittance *
          mRecCamera.sigmaS * mRecCamera.transmittance * phaseTerm / pdfKernel;  // not yet include kernel weight
      weightKernel = 1.0 / kernelVol;
      beamTrans = mRecBeam.transmittance;

      if (!beam->longBeams) {
        contrib /= mRecBeam.pdfFailure;
        pdfEdgeFailure = mRecBeam.pdfFailure;
      } else {
        pdfEdgeFailure = 1.0f;
      }
    } else {
      SAssert(false);
    }
  }

  bool isValid() const {
    return !contrib.isZero();
  }

  Float pdf() const {
    return pdfEdgeFailure * pdfKernel;
  }

  Float kernelPDF(const Ray &cameraRay, const Point &orgBeam, const Vector &dBeam, const Float newDLength) const {
    if (volTechnique == EBeamBeam1D) {
      return std::sqrt(cross(cameraRay.d, dBeam).lengthSquared());
    } else if (volTechnique == EBeamBeam3D_Naive) {
      // SH: pdf of naive case is not tested
      SAssert(false); // Not tested properly
      // See commit 07/16 to see this code
      return 0.f;
    } else if (volTechnique == EBeamBeam3D_Optimized) {
      Ray photonRay(orgBeam, dBeam, 0.f);
      Float _pdfKernel = 0;

      Ray _beam(photonRay.o, photonRay.d, 0.f, photonRay.maxt, 0.f);
      Ray _cam(cameraRay.o, cameraRay.d, 0.f, cameraRay.maxt, 0.f);

      double tNearBeam, tFarBeam;
      if (cylinderIntersection(_cam, _beam, radius, tNearBeam, tFarBeam)) {
        // Uniform sampling on a line
        _pdfKernel = 1.0 / std::max(tFarBeam - tNearBeam, 0.0001);

        // Pdf of the sample point on camera beam
        Point kernelCentroid = photonRay.o + photonRay.d * newDLength;
        Float distToProj = dot(kernelCentroid - cameraRay.o, cameraRay.d);
        Float distSqr = (cameraRay(distToProj) - kernelCentroid).lengthSquared();
        Float radSqr = radius * radius;
        if (distSqr < radSqr) {
          Float deltaT = math::safe_sqrt(radSqr - distSqr); // equal to half of [tc-, tc+] range
          _pdfKernel *= 1.0 / std::max(2.0 * deltaT, 0.0001);
          return _pdfKernel;
        } else {
          return 0.f;
        }
      }
      return 0.f;
    } else {
      SLog(EError, "Unsupported pdf for shifted ray for volume technique %d", volTechnique);
      return 0.f;
    }
  }

};

class BeamGradRadianceQuery {
public:
  // Scene and camera path informations
  Scene *scene;
  GatherPoint *baseGather;
  std::vector<ShiftGatherPoint> &shiftGPs;

  // Extra information
  const GPMConfig &config;
  GPMThreadData &thdata;

  // Accumulate spectrum value
  Spectrum extraFlux = Spectrum(0.f);
  Spectrum mediumFlux;
  Spectrum shiftedMediumFlux[4];
  Spectrum weightedMediumFlux[4];
  const Ray &baseCameraRay;               // The query (camera) ray
  const Medium *medium;

  // The current camera edgeID
  size_t currCameraEdge;

  // Cache values
  Path cachePath;
  Float cacheDetSource;

  // Sampler
  Sampler *sampler;

#if HAVE_ADDITIONAL_STATS
  // To have some statistics about the shifts, per pixels
  ShiftStats shiftStats[4];
#endif

  BeamGradRadianceQuery(Scene *scene_,
                        GatherPoint *baseGP_,
                        std::vector<ShiftGatherPoint> &shiftGPs_,
                        const Ray &baseCameraRay_,
                        const Medium *medium_,
                        const GPMConfig &config_, GPMThreadData &thdata_,
                        size_t _currCameraEdge, Sampler *_sampler) :
      scene(scene_), baseGather(baseGP_), shiftGPs(shiftGPs_),
      config(config_), thdata(thdata_),
      baseCameraRay(baseCameraRay_), medium(medium_), currCameraEdge(_currCameraEdge),
      cacheDetSource(-1), sampler(_sampler) {
    for (int i = 0; i < 4; i++) {
      shiftedMediumFlux[i] = weightedMediumFlux[i] = Spectrum(0.f);
    }
    mediumFlux = Spectrum(0.f);
  }

  void clearCache() {
    //cachePath.release(thdata.pool);
    if (cachePath.length() != 0) {
      // Deallocate only the path
      thdata.pool.release(cachePath.edge(cachePath.edgeCount() - 1));
      thdata.pool.release(cachePath.vertex(cachePath.vertexCount() - 1));
      thdata.pool.release(cachePath.vertex(cachePath.vertexCount() - 2));
      cachePath.clear();
    }
    cacheDetSource = -1.f;
  }

  /**
   Entry point of the beam beam gradient implementation 
   
   beam: the photon beam found by some nearest neighbor searches.
   */
  bool operator()(const LTPhotonBeam *beam,
                  Float tmin = 0,
                  Float tmax = std::numeric_limits<Float>::infinity());

  Point getShiftPos(const Ray &bRay,
                    const Ray &sRay,
                    Float w, const Vector &u, Float radius,
                    Float newW, bool coherent = true);

  Point getShiftPos1D(const Ray &bRay,
                      const Ray &sRay,
                      const Point &a,
                      const Vector &bBeamDir,
                      Float w, Float u);

  /**
   * beam       : photon beam on the light subpath. 
   * shiftRay   : the camera path edge that we will seek intersection with the beam. Maxt is set up to the intersection only.
   * v          : distance to intersection on photon beam.
   * w          : distance to intersection on camera edge. 
   * u          : distance between two "intersections" on photon beam and camera edge. 
   */
  bool shiftBeam(const Point &newPos,
                 const LTPhotonBeam *beam,
                 const ShiftGatherPoint &shiftGP,
                 const Ray &shiftRay,
                 Float shiftW,
                 GradientSamplingResult &result,
                 const BeamKernelRecord &kRec,
                 int b,
                 ELightShiftType currShift);

  bool shiftBeamME(int b, int c,
                   const LTPhotonBeam *beam, const ShiftGatherPoint &shiftGP,
                   const Ray &shiftRay, Float shiftW, GradientSamplingResult &result,
                   const BeamKernelRecord &kRec, const Point &newPos);

  bool shiftBeamDiffuse(const LTPhotonBeam *beam, const ShiftGatherPoint &shiftGP,
                        const Ray &shiftRay, Float shiftW, GradientSamplingResult &result,
                        const BeamKernelRecord &kRec,
                        const Point &newPos);

  void cacheSourcePath(int b, int c,
                       const LTPhotonBeam *beam, const BeamKernelRecord &kRec);

  /**
   * When sampling is false, the kernel record has to be filled with 
   * v, w, radius, pdfV, pdfW, and pdfEdgeFailure. 
   */
  bool shiftNull3D(const LTPhotonBeam *beam, const ShiftGatherPoint &shiftGP,
                   const Ray &shiftRay, Float shiftW, GradientSamplingResult &result,
                   const BeamKernelRecord &kRec, BeamKernelRecord &kRecShift);

};

MTS_NAMESPACE_END

#endif //MITSUBA_SHIFT_VOLUME_BEAMS_H
