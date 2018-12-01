#ifndef MITSUBA_SHIFT_MEDIUM_H
#define MITSUBA_SHIFT_MEDIUM_H

#include <mitsuba/mitsuba.h>

#include "../shift_utilities.h"
#include "../../gvpm_beams.h"

MTS_NAMESPACE_BEGIN

bool intersectionMediumShift(const Ray &ray_, const Point &ori, const Vector &w0, const Vector &w1,
                             Float &t0, Float &t1) {
  // Sanity checks
  SAssert(fabs(w1.length() - 1.0) < 10e-5);
  SAssert(fabs(w0.length() - 1.0) < 10e-5);

  Vector P = cross(ray_.d, w1);
  Float det = dot(w0, P);
  if (std::abs(det) < 1e-8f) {
    SLog(EWarn, "Is it possible?");
    return false;
  }

  Float invDet = 1.0f / det;
  Vector T = ray_.o - ori;

  Vector Q = cross(T, w0);
  t1 = dot(ray_.d, Q) * invDet;
  if (t1 < 0.0f) {
    SLog(EError, "Impossible case?");
    return false;
  }

  t0 = dot(T, P) * invDet;
  return t0 >= 0.0f;

}

bool mediumRotationShift(ShiftRecord &sRec,
                         const Vector &baseSensorDir,
                         const Vector &offsetSensorDir,
                         const PathVertex *aPred,
                         const PathVertex *a0,
                         const PathEdge *e0,
                         const PathVertex *a1,
                         const PathEdge *e1,
                         const Point a2p,
                         Vector &newW1) {
  // In this shfit, we create a "virtual photon plane"
  // Where a0, a1 and a2 create this plane.
  // The intersection with the sensor path is mimic at the position a2
  Ray shiftRay(a2p - offsetSensorDir,
               offsetSensorDir, 0.f, 2, 0.f);

  // Construct orthogonal planes
  Vector orthNewW1 = a2p - (a0->getPosition() + e0->d * dot(a2p - a0->getPosition(), e0->d));
  orthNewW1 /= orthNewW1.length();

  // From W1Ortho' compute W1'
  // by copying the local dot product of the orignal W0 and W1 vectors
  Float w0Dot = dot(e0->d, e1->d);
  // sqrt(1 - (w0Dot*w0Dot)) -> pythagore
  newW1 = sqrt(1 - (w0Dot * w0Dot)) * orthNewW1 + e0->d * w0Dot;

  // Check that the new plane does is still valid. For that we need to check
  // - t0New and t1New are valid (> 0)
  // - visibility if t0New > t0
  // - visibility between a2p and a1p
  Float t0New, t1New;
  Vector d0New = e0->d;
  if (!intersectionMediumShift(shiftRay, a0->getPosition(), e0->d, newW1,
                               t0New, t1New)) {
    if (t0New < 0) {
      // Flip back the first vector
      d0New = -d0New;
      t0New = std::abs(t0New);
    } else {
      return false;
    }
  }

  // FIXME: Test VISIBILITY !

  // Prepare variable for the accumulations
  sRec.pdf = 1.f;
  sRec.throughtput = Spectrum(1.f);

  //////////////////
  // First vertex interaction
  /////////////////
  bool diracInteraction = false;
  if (a0->getType() == PathVertex::ESurfaceInteraction) {
    // parentIts could be on an emitter (even though type is surface interaction).
    // We allow this because photon tracing allows path to bounce on a light source surface.
    const Intersection &parentIts = a0->getIntersection();
    SAssert(parentIts.p == a0->getPosition());

    Vector3 pWo = parentIts.toLocal(d0New);
    Vector3 pWi = parentIts.wi;
    const BSDF *parentBSDF = parentIts.getBSDF();
    BSDFSamplingRecord bRec(parentIts,
                            pWi,
                            pWo,
                            EImportance);
    bRec.component = a0->sampledComponentIndex;
    if (bRec.component >= parentBSDF->getComponentCount()) {
      SLog(EWarn, "Invalid component request %d", bRec.component);
    }

    EMeasure measure = a0->isDegenerate() ? EDiscrete : ESolidAngle;
    if (a0->isDegenerate()) {
      diracInteraction = true;
    }

    // Evaluate the BRDF from the parent vertex to the offset vertex
    sRec.throughtput *= parentBSDF->eval(bRec, measure);
    sRec.pdf *= parentBSDF->pdf(bRec, measure);
    sRec.pdf *= parentBSDF->pdfComponent(bRec);

    // Adjoint BSDF for shading normals (see vertex.cpp for example)
    Float wiDotGeoN = dot(parentIts.geoFrame.n, parentIts.toWorld(pWi));
    Float woDotGeoN = dot(parentIts.geoFrame.n, e0->d);
    if (wiDotGeoN * Frame::cosTheta(bRec.wi) <= 0 ||
        woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
      return false;

    // FIXME: Missing adjoint compensation
  } else if (a0->getType() == PathVertex::EMediumInteraction) {
    const MediumSamplingRecord &pMRec = a0->getMediumSamplingRecord();
    Vector3 pWo = d0New;
    Vector3 pWi = normalize(aPred->getPosition() - pMRec.p);
    const PhaseFunction *parentPhase = pMRec.getPhaseFunction();
    PhaseFunctionSamplingRecord pRec(pMRec,
                                     pWi,
                                     pWo,
                                     EImportance);

    // No component inside phase function, no pdf for selecting it
    // Phase function are always expressed in solid angle
    sRec.throughtput *= pMRec.sigmaS * parentPhase->eval(pRec);
    sRec.pdf *= parentPhase->pdf(pRec);
  } else {
    SAssert(a0->getType() == PathVertex::EEmitterSample);
    SLog(EError, "Impossible for now");
  }

  // Put the first interaction in the area domain
  if (!diracInteraction) {
    Float GOp = 1 / (t0New * t0New);
    sRec.pdf *= GOp;
    sRec.throughtput *= GOp;
  }
  sRec.throughtput /= a0->pdf[EImportance];
  sRec.throughtput *= a0->rrWeight;

  ////////////////////
  // Second vertex interaction
  ////////////////////
  // We know that this second vertex is a medium interaction
  SAssert(a1->getType() == PathVertex::EMediumInteraction);
  {
    const MediumSamplingRecord &pMRec = a1->getMediumSamplingRecord();
    Vector3 pWo = newW1;
    Vector3 pWi = -d0New;
    const PhaseFunction *parentPhase = pMRec.getPhaseFunction();
    PhaseFunctionSamplingRecord pRec(pMRec,
                                     pWi,
                                     pWo,
                                     EImportance);

    // No component inside phase function, no pdf for selecting it
    // Phase function are always expressed in solid angle
    sRec.throughtput *= pMRec.sigmaS * parentPhase->eval(pRec);
    sRec.pdf *= parentPhase->pdf(pRec);
  }

  // Put the second interaction in the area domain
  {
    Float GOp = 1 / (t1New * t1New);
    sRec.pdf *= GOp;
    sRec.throughtput *= GOp;
  }
  sRec.throughtput /= a1->pdf[EImportance];
  sRec.throughtput *= a1->rrWeight;

  ///////////////////////
  // Evaluate transmittance changes
  //////////////////////
  // Evaluate the transmittance for the new
  // shifted light path
  Ray rayTrans0Shift(a0->getPosition(), d0New, 0.f);
  rayTrans0Shift.mint = 0.f;
  rayTrans0Shift.maxt = t0New;
  MediumSamplingRecord mRec0Shift;
  e0->medium->eval(rayTrans0Shift, mRec0Shift);

  Ray rayTrans1Shift(a0->getPosition(), newW1, 0.f);
  rayTrans1Shift.mint = 0.f;
  rayTrans1Shift.maxt = t1New;
  MediumSamplingRecord mRec1Shift;
  e0->medium->eval(rayTrans1Shift, mRec1Shift);

  // Compute the shift throughput
  // Update the throughput based on transmittance
  sRec.throughtput *= mRec0Shift.transmittance / e0->pdf[EImportance];
  sRec.throughtput *= mRec1Shift.transmittance / e1->pdf[EImportance];
  sRec.pdf *= mRec0Shift.pdfSuccess;
  sRec.pdf *= mRec1Shift.pdfSuccess;

  // FIXME: Do we need this Jacobian ratio? I do not think so...
  //throughputShift /= plane->invJacobian(baseCameraRay.d);
  //throughputShift *= 1.0 /  absDot(plane->w0(), cross(newW1, shiftRay.d));

  // The jacobian is to:
  // - decoupling the variables
  // - recoupling the variables
  sRec.jacobian = 1.f;
  Float ratioJ = absDot(d0New, cross(newW1, offsetSensorDir)) / absDot(e0->d, cross(e1->d, baseSensorDir));
  sRec.jacobian *= ratioJ;
  sRec.jacobian /= t1New / e1->length;
  if (diracInteraction) { // If parent from parent was a specular vertex
    sRec.jacobian /= t0New / e0->length;
  }

//  Float pdfBase = a0->pdf[EImportance] * a1->pdf[EImportance] * e0->pdf[EImportance] * e1->pdf[EImportance];

  return true;
}

MTS_NAMESPACE_END

#endif //MITSUBA_SHIFT_MEDIUM_H
