//
// Created by beltegeuse on 7/12/17.
//

#include "shift_diffuse.h"

#include "../shift_volume_beams.h"

MTS_NAMESPACE_BEGIN

bool diffuseReconnection(ShiftRecord &sRec, const Intersection &newIts,
                         const Vector &newD, const Float newDLength,
                         const PathVertex *parentVertex,
                         const PathEdge *edge, bool isVolumeBase,
                         const Point &predPos, bool adjointCorr) {

  Float pdfValue;
  if (parentVertex->getType() == PathVertex::ESurfaceInteraction) {
    // parentIts could be on an emitter (even though type is surface interaction).
    // We allow this because photon tracing allows path to bounce on a light source surface.
    const Intersection &parentIts = parentVertex->getIntersection();
    SAssert(parentIts.p == parentVertex->getPosition());

    Vector3 pWo = parentIts.toLocal(newD);
    Vector3 pWi = parentIts.wi;
    const BSDF *parentBSDF = parentIts.getBSDF();
    BSDFSamplingRecord bRec(parentIts,
                            pWi,
                            pWo,
                            EImportance);
    bRec.component = parentVertex->sampledComponentIndex;
    if (bRec.component >= parentBSDF->getComponentCount()) {
      SLog(EWarn, "Invalid component request %d", bRec.component);
    }
    EMeasure measure = ESolidAngle;

    // Evaluate the BRDF from the parent vertex to the offset vertex
    sRec.throughtput *= parentBSDF->eval(bRec, measure);
    pdfValue = parentBSDF->pdf(bRec, measure);
    pdfValue *= bRec.component == -1? 1.0 : parentBSDF->pdfComponent(bRec);

    // Adjoint BSDF for shading normals (see vertex.cpp for example)
    Float wiDotGeoN = dot(parentIts.geoFrame.n, parentIts.toWorld(pWi));
    Float woDotGeoN = dot(parentIts.geoFrame.n, newD);
    if (wiDotGeoN * Frame::cosTheta(bRec.wi) <= 0 ||
        woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
      return false;

    if (adjointCorr) {
      sRec.throughtput *= std::abs(
          (Frame::cosTheta(bRec.wi) * woDotGeoN) /
              (Frame::cosTheta(bRec.wo) * wiDotGeoN));
    }
  } else if (parentVertex->getType() == PathVertex::EMediumInteraction) {
    // Said that the parent is not on surface
    // Usefull for the Jacobian and pdf computation
    const MediumSamplingRecord &pMRec = parentVertex->getMediumSamplingRecord();
    Vector3 pWo = newD;
    Vector3 pWi = normalize(predPos - pMRec.p);
    const PhaseFunction *parentPhase = pMRec.getPhaseFunction();
    PhaseFunctionSamplingRecord pRec(pMRec,
                                     pWi,
                                     pWo,
                                     EImportance);

    // No component inside phase function, no pdf for selecting it
    // Phase function are always expressed in solid angle
    sRec.throughtput *= pMRec.sigmaS * parentPhase->eval(pRec);
    pdfValue = parentPhase->pdf(pRec);

    // No adjoint correction required because no surface
  } else {
    // Parent vertex samples an out-going direction from the light.
    // The position sampling on the light is handled by the emitter super node.
    SAssert(parentVertex->getType() == PathVertex::EEmitterSample);

    const AbstractEmitter *emitter = parentVertex->getAbstractEmitter();
    const PositionSamplingRecord &pRec = parentVertex->getPositionSamplingRecord();
    DirectionSamplingRecord dRec;
    dRec.d = newD;       // world
    dRec.measure = ESolidAngle;
    sRec.throughtput *= emitter->evalDirection(dRec, pRec);
    pdfValue = emitter->pdfDirection(dRec, pRec);
  }



  // All the values are in solid angle, convert them
  if (isVolumeBase) {
    Float GOp = 1 / (newDLength * newDLength);
    sRec.pdf = pdfValue * GOp;
    sRec.throughtput *= GOp;
  } else {
    Float GOp = geometryOpposingTerm(parentVertex->getPosition(),
                                     newIts.p,
                                     newIts.geoFrame.n);
    sRec.pdf = pdfValue * GOp;
    sRec.throughtput *= GOp;
  }

  if (parentVertex->pdf[EImportance] == 0.f) {
    SLog(EWarn, "Do not do diffuse reconnection as parent pdf is 0");
    sRec.pdf = 0.f;
    return false;
  }

  // Jacobian deduction
  sRec.jacobian = 1.0f;

  // In Area domain
  sRec.throughtput /= parentVertex->pdf[EImportance];
  sRec.throughtput *= parentVertex->rrWeight;

  if (edge->medium != nullptr) {
    MediumSamplingRecord mRecShift;
    Ray mRay(parentVertex->getPosition(), newD, 0.f, newDLength, 0.f);
    edge->medium->eval(mRay, mRecShift);
    // There is a medium between the two points,
    // Need to take into account the pdf change of traversing the media
    if (isVolumeBase) {
      sRec.pdf *= mRecShift.pdfSuccess; // add the prob to stop here
    } else {
      sRec.pdf *= mRecShift.pdfFailure; // add the prob of pass through the medium
    }

    // Update throughtput
    sRec.throughtput *= mRecShift.transmittance / edge->pdf[EImportance]; // mRecShift.pdfFailure;

    // No Jacobian is have to be taking into account
    // Because the shift doesn't change the decision (measure).
  }

  return true;
}

bool diffuseReconnectionPhotonBeam(ShiftRecord &sRec, const Point &newPos, const Point &basePos,
                                   const Vector &newD, const Float newDLength,
                                   const PathVertex *baseVertex,
                                   const PathVertex *parentVertex,
                                   const PathEdge *parentEdge,
                                   Float pdfEdgeAndKernel, bool longBeam, const Point &parentParentPos,
                                   bool adjointCorr) {
  /**
   * parentVertex --> basePos (intersection point with other beam) --> baseVertex (end point of this beam during photon sampling)
   *
   * We evaluate throughput and pdf up to basePos only.
   */
  Float pdfValueSA;
  if (parentVertex->getType() == PathVertex::ESurfaceInteraction) {
    // parentIts could be on an emitter (even though type is surface interaction).
    // We allow this because photon tracing allows path to bounce on a light source surface.
    const Intersection &parentIts = parentVertex->getIntersection();
    SAssert(parentIts.p == parentVertex->getPosition());

    Vector3 pWo = parentIts.toLocal(newD);
    Vector3 pWi = parentIts.wi;
    const BSDF *parentBSDF = parentIts.getBSDF();
    BSDFSamplingRecord bRec(parentIts,
                            pWi,
                            pWo,
                            EImportance);
    bRec.component = parentVertex->sampledComponentIndex;
    if (bRec.component >= parentBSDF->getComponentCount()) {
      SLog(EWarn, "Invalid component request %d", bRec.component);
    }
    EMeasure measure = ESolidAngle;

    // Evaluate the BRDF from the parent vertex to the offset vertex
    sRec.throughtput *= parentBSDF->eval(bRec, measure);
    pdfValueSA = parentBSDF->pdf(bRec, measure);
    pdfValueSA *= parentBSDF->pdfComponent(bRec);

    // Adjoint BSDF for shading normals (see vertex.cpp for example)
    Float wiDotGeoN = dot(parentIts.geoFrame.n, parentIts.toWorld(pWi));
    Float woDotGeoN = dot(parentIts.geoFrame.n, newD);
    if (wiDotGeoN * Frame::cosTheta(bRec.wi) <= 0 ||
        woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
      return false;

    if (adjointCorr) {
      sRec.throughtput *= std::abs(
          (Frame::cosTheta(bRec.wi) * woDotGeoN) /
              (Frame::cosTheta(bRec.wo) * wiDotGeoN));
    }
  } else if (parentVertex->getType() == PathVertex::EMediumInteraction) {
    // Parent vertex is in the medium
    // Useful for the Jacobian and pdf computation
    const MediumSamplingRecord &pMRec = parentVertex->getMediumSamplingRecord();
    Vector3 pWo = newD;
    Vector3 pWi = normalize(parentParentPos - pMRec.p);
    const PhaseFunction *parentPhase = pMRec.getPhaseFunction();
    PhaseFunctionSamplingRecord pRec(pMRec,
                                     pWi,
                                     pWo,
                                     EImportance);

    // No component inside phase function, no pdf for selecting it
    // Phase function are always expressed in solid angle
    sRec.throughtput *= pMRec.sigmaS * parentPhase->eval(pRec);
    pdfValueSA = parentPhase->pdf(pRec);

    // No adjoint correction required because no surface
  } else {
    // Parent vertex samples an out-going direction from the light.
    // The position sampling on the light is handled by the emitter super node.
    SAssert(parentVertex->getType() == PathVertex::EEmitterSample);

    const AbstractEmitter *emitter = parentVertex->getAbstractEmitter();
    const PositionSamplingRecord &pRec = parentVertex->getPositionSamplingRecord();
    DirectionSamplingRecord dRec;
    dRec.d = newD;       // world
    dRec.measure = ESolidAngle;
    sRec.throughtput *= emitter->evalDirection(dRec, pRec);
    pdfValueSA = emitter->pdfDirection(dRec, pRec);
  }

  // Convert solid angle measure to area
  // In volume shift, the new vertex is always not on surface, so no cosine angle involved for GOp.
  Float GOpNew = 1 / (newDLength * newDLength);
  sRec.pdf = pdfValueSA * GOpNew; // Area domain
  sRec.throughtput *= GOpNew; // Area domain

  // Jacobian in area measure is 1
  sRec.jacobian = 1.0f;

  // Go from area to solid angle measure (to remove old geometry factor)
  Float pdfBasePos =
      parentVertex->pdf[EImportance] * (parentVertex->getPosition() - baseVertex->getPosition()).lengthSquared();
  if (baseVertex->isOnSurface())
    pdfBasePos /= absDot(baseVertex->getGeometricNormal(), parentEdge->d);

  // Add back the "new geometry factor"
  Float GOpBase = 1.0f / (parentVertex->getPosition() - basePos).lengthSquared();
  pdfBasePos *= GOpBase; // Area domain

  if (pdfBasePos == 0.f) {
    SLog(EWarn, "pdfBase is 0... do not do diffuse reconnection");
    sRec.pdf = 0.f;
    return false;
  }

  // For a fast implementation, we use solid angle measure to compute the throughput
  sRec.throughtput /= pdfBasePos; // Area domain
  sRec.throughtput *= parentVertex->rrWeight;

  if (parentEdge->medium != nullptr) {

    // We have to "manually" compute here because we didn't create an PathEdge on the shifted path
    MediumSamplingRecord mRecShift;
    Ray mRay(parentVertex->getPosition(), newD, 0.f, newDLength, 0.f);
    parentEdge->medium->eval(mRay, mRecShift);
    // There is a medium between the two points,
    // Need to take into account the pdf change of traversing the media


    if (!longBeam) {
      // Long Beams, so the PDF is equal to 1
      sRec.pdf *= mRecShift.pdfFailure; // add the prob to continue
    } else {
      // during long beam tracing, no distance sampling is done so pdf is 1
    }

    // Update throughtput
    sRec.throughtput *= mRecShift.transmittance / pdfEdgeAndKernel;
  }

  return true;
}

MTS_NAMESPACE_END
