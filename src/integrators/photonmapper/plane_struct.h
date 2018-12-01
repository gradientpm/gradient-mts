//
// Created by beltegeuse on 8/29/17.
//



#ifndef MITSUBA_PLANE_STRUCT_H
#define MITSUBA_PLANE_STRUCT_H

#include <mitsuba/render/photonmap.h>
#include <mitsuba/render/scene.h>

#include "beams_struct.h"
#include "../volume_utils.h"

MTS_NAMESPACE_BEGIN

struct PhotonPlane {
public:
  struct IntersectionRecord {
    Float tCam, t0, t1, invDet;
  };

protected:
  // Photon beam
  Point _ori;
  Vector _w0;
  Float _length0;

  // The extension of the beam
  Vector _w1;
  Float _length1;
public:
  const Medium *medium;
  Spectrum _flux;
  int depth;

  PhotonPlane() {
    // FIXME: Intialize the values ...
  }

  PhotonPlane(const Point &ori_, const Vector &w0_, const Float length0_,
              const Vector &w1_, const Float length1_,
              const Medium *m, const Spectrum &f, int depth_) :
      _ori(ori_), _w0(w0_), _length0(length0_),
      _w1(w1_), _length1(length1_),
      medium(m), _flux(f), depth(depth_) {
    SAssert(_length0 != 0 && std::isfinite(_length0));
    SAssert(_length1 != 0 && std::isfinite(_length1));
  }

  /**
   * @return center point of the photon plane
   * @note this function is usefull for pointKD tree construction
   */
  Point getCenter() const {
    return _ori + _w0 * _length0 * 0.5 + _w1 * _length1 * 0.5;
  }

  AABB getAABB() const {
    AABB aabb(_ori);
    aabb.expandBy(_ori + _w0 * _length0);
    aabb.expandBy(_ori + _w1 * _length1);
    aabb.expandBy(_ori + _w1 * _length1 + _w0 * _length0);
    return aabb;
  }

  /**
   * Helper to transform photon beams into plane
   * @param beam
   * @return
   */
  static PhotonPlane transformBeam(const PhotonBeam &beam, Sampler *sampler) {
    // TODO: Assume homogenous media
    // So the direction and the location of sample distance
    // Doesn't really matter here.

    // 1) Sample a distance (without taking into account the visibility)
    // The transmittance value of the first edge
    Ray rayNew(beam.getOri(), beam.getDir(), 0.f);
    MediumSamplingRecord mRecNew;
    beam.medium->sampleDistance(rayNew, mRecNew, sampler); // This sample the distance without vis tests.

    // 2) Sample additional scattering direction
    PhaseFunctionSamplingRecord pRec(mRecNew, -beam.getDir(), EImportance);
    do {
      beam.medium->getPhaseFunction()->sample(pRec, sampler);
    } while (absDot(beam.getDir(), pRec.wo) == 1);

    // 3) Create the photon plane and return it
    return PhotonPlane(beam.getOri(), beam.getDir(), beam.getLength(),
                       pRec.wo, mRecNew.t, beam.medium, beam.flux, beam.depth);
  }

  /**
   * Möller–Trumbore intersection algorithm variant for planes
   * @param ray_ the camera ray
   * @param tCam the intersection distance from the camera
   * @param t0 the distance on the first direction
   * @param t1 the distance of the second direction
   * @param invDet the non-normalize inverse determinant (w0 . (w1 x ray_.d))
   * @return
   */
  bool intersectPlane0D(const Ray &ray_,
                        Float &tCam, Float &t0, Float &t1,
                        Float &invDet) const {
    Vector e0 = _w0 * _length0;
    Vector e1 = _w1 * _length1;
    Vector P = cross(ray_.d, e1);
    float det = dot(e0, P);
    if (std::abs(det) < 1e-5f)
      return false;

    invDet = 1.0f / det;
    Vector T = ray_.o - _ori;
    t0 = dot(T, P) * invDet;
    if (t0 < 0.0f || t0 > 1.0f)
      return false;

    Vector Q = cross(T, e0);
    t1 = dot(ray_.d, Q) * invDet;
    if (t1 < 0.0f || t1 > 1.0f)
      return false;

    tCam = dot(e1, Q) * invDet;
    if (tCam <= ray_.mint || tCam >= ray_.maxt)
      return false;

    // Scale to the correct distance
    // In order to use correctly transmittance sampling
    t1 *= _length1;
    t0 *= _length0;

    return true;
  }
  bool intersectPlane0D(const Ray &ray,
                        IntersectionRecord &rec) const {
    return intersectPlane0D(ray, rec.tCam, rec.t0, rec.t1, rec.invDet);
  }

  /**
   * Compute contribution of the beam 0D kernel
   * @param t0 the distance on the first edge of the plane
   * @param t1 the distance on the second edge of the plane
   * @param mRecCamera the camera values
   * @param d the direction of the camera
   * @param invDet the non-normalized inv determinant
   * @return
   */
  Spectrum getContrib0D(Float t0,
                        Float t1,
                        const MediumSamplingRecord &mRecCamera,
                        const Vector &d,
                        const Float &invDet) const {
    //WARN: The visibility is assumed to be correct

    PhaseFunctionSamplingRecord pRec(mRecCamera, -_w1, -d, EImportance);
    Float phaseTerm = mRecCamera.getPhaseFunction()->eval(pRec);

    // The transmittance value of the first edge
    Ray rayTrans0(_ori, _w0, 0.f);
    rayTrans0.mint = 0.f;
    rayTrans0.maxt = t0;
    MediumSamplingRecord mRec0;
    medium->eval(rayTrans0, mRec0);

    // The transmittance value of the second edge
    Ray rayTrans1(rayTrans0(t0), _w1, 0.f);
    rayTrans1.mint = 0.f;
    rayTrans1.maxt = t1;
    MediumSamplingRecord mRec1;
    medium->eval(rayTrans1, mRec1);

    Spectrum contrib = mRecCamera.transmittance
        * (mRecCamera.sigmaS) * (mRec0.sigmaS) // The bounce on mRec1
        * _flux * phaseTerm; // The _flux

    // Transmittance of the current path
    // This terms cancel out to gives only sigmaT
    contrib *= mRec1.transmittance * mRec0.transmittance;
    contrib /= mRec0.pdfFailure;
    contrib /= mRec1.pdfFailure;

    // The jacobian expression
    // as it is not normalized, we normalize it.
    //contrib *= std::abs(invDet) * _length0 * _length1;

    // Use the same jacobian inside the paper
    contrib *= invJacobian(d);

    return contrib;
  }
  Spectrum getContrib0D(const IntersectionRecord &rec,
                        const MediumSamplingRecord &mRecCamera,
                        const Vector &d) const {
    return getContrib0D(rec.t0, rec.t1, mRecCamera, d, rec.invDet);
  }

  const Float invJacobian(const Vector &k) const {
    return 1.0 / absDot(_w0, cross(_w1, k));
  }
  const Point ori() const {
    return _ori;
  }
  const Vector w1() const {
    return _w1;
  }
  const Vector w0() const {
    return _w0;
  }
  const Float area() const {
    return _length0 * _length1;
  }
  const Float length0() const {
    return _length0;
  }
  const Float length1() const {
    return _length1;
  }
  const Spectrum &flux() const {
    return _flux;
  }

};

// Query to compute photon plane contribution
class PhotonPlaneQuery {
public:
  PhotonPlaneQuery(const Ray &ray_,
                   const Medium *medium_,
                   int maxDepth_,
                   int minDepth_,
                   ref <Sampler> sampler_,
                   EVolumeTechnique volTechnique_) :
      baseCameraRay(ray_), medium(medium_),
      maxDepth(maxDepth_), minDepth(minDepth_), Li(0.f), sampler(sampler_), volTechnique(volTechnique_) {}

  bool operator()(const PhotonPlane *plane) {
    Float tCam, t0, t1, invDet;
    if (plane->intersectPlane0D(baseCameraRay, tCam, t0, t1, invDet)) {
      // Compute the transmittance of the camera
      Ray rayCamToPlane = baseCameraRay;
      rayCamToPlane.mint = 0.f;
      rayCamToPlane.maxt = tCam;

      MediumSamplingRecord mRecCam;
      medium->eval(rayCamToPlane, mRecCam);

      // Compute the total contrib of the plane
      Li += plane->getContrib0D(t0, t1, mRecCam, baseCameraRay.d, invDet);

      return true;
    }

    return false;
  }

  // The initial ray
  const Ray &baseCameraRay;
  const Medium *medium;
  int maxDepth;
  int minDepth;

  // Accumulate spectrum value
  Spectrum Li;

  ref <Sampler> sampler;
  EVolumeTechnique volTechnique;
};

MTS_NAMESPACE_END

#endif //MITSUBA_PLANE_STRUCT_H
