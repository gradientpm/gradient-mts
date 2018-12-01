#pragma once

#ifndef MITSUBA_SHIFT_VOLUME_PLANES_H
#define MITSUBA_SHIFT_VOLUME_PLANES_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/scene.h>

#include "shift_cameraPath.h"
#include "shift_utilities.h"

#include "../gvpm_beams.h"

// For diffuse shift
#include "operation/shift_diffuse.h"

MTS_NAMESPACE_BEGIN

class PlaneGradRadianceQuery {
public:
  // Scene and camera path informations
  Scene *scene;
  GatherPoint *baseGather;
  std::vector<ShiftGatherPoint> &shiftGPs;

  // Extra information
  const GPMConfig &config;
  GPMThreadData &thdata;

  // Accumulate spectrum value
  Spectrum mediumFlux;
  Spectrum shiftedMediumFlux[4];
  Spectrum weightedMediumFlux[4];

  const Ray &baseCameraRay;               // The query (camera) ray
  const Medium *medium;

  // The current camera edgeID
  int currCameraEdge;

  PlaneGradRadianceQuery(Scene *scene_,
                         GatherPoint *baseGP_,
                         std::vector<ShiftGatherPoint> &shiftGPs_,
                         const Ray &baseCameraRay_,
                         const Medium *medium_,
                         const GPMConfig &config_, GPMThreadData &thdata_,
                         int _currCameraEdge) :
      scene(scene_), baseGather(baseGP_), shiftGPs(shiftGPs_),
      config(config_), thdata(thdata_),
      baseCameraRay(baseCameraRay_), medium(medium_), currCameraEdge(_currCameraEdge) {
    for (int i = 0; i < 4; i++) {
      shiftedMediumFlux[i] = weightedMediumFlux[i] = Spectrum(0.f);
    }
    mediumFlux = Spectrum(0.f);
  }

  bool operator()(const LTPhotonPlane *plane) {
    PhotonPlane::IntersectionRecord bRec;
    if (!plane->intersectPlane0D(baseCameraRay, bRec)) {
      return false;
    }

    // Compute the transmittance of the camera
    Ray rayCamToPlane = baseCameraRay;
    rayCamToPlane.mint = 0.f;
    rayCamToPlane.maxt = bRec.tCam;

    MediumSamplingRecord mRecCam;
    medium->eval(rayCamToPlane, mRecCam);

    // Compute the total contrib of the plane
    Spectrum baseContrib = plane->getContrib0D(bRec, mRecCam, baseCameraRay.d);
    mediumFlux += baseContrib;

    // Prepare the shift informations
    const Point2 basePixel = baseGather->path.vertex(1)->getSamplePosition();
    const std::array<Point2, 4> pixels = generateOffsetPos(basePixel);
//      Vector2i filmSize = scene->getFilm()->getSize();

    for (int i = 0; i < 4; ++i) {
      shiftGPs[i].generate(scene, thdata.pool, *baseGather,
                           pixels[i], false);

      GradientSamplingResult result;
      if (shiftGPs[i].validVolumeEdge(currCameraEdge, plane->medium)) {
        // Retieve relevent information for the shift path
//              if(plane->edgeID > 2) {
//                  // Shift the origin of the plane to get
//                  // the same intersection parametrisation.
//                  diffuseShift(shiftGPs[i], bRec, plane, result);
//              } else {
        // Origin and w0 is lock
        specularShift(shiftGPs[i], bRec, plane, baseContrib, result);
//              }
      }

      shiftedMediumFlux[i] += result.weight * result.shiftedFlux;
      weightedMediumFlux[i] += result.weight * baseContrib;
    }
    return false;
  }

  /**
   * This shift is a basic shift that keep the local parametrisation
   * and only move the original photon plane to achieve the same
   * parametrisation
   */
  bool diffuseShift(const ShiftGatherPoint &shiftGP,
                    const PhotonPlane::IntersectionRecord &bRec,
                    const LTPhotonPlane *plane,
                    GradientSamplingResult &result) {
    const Path &shiftPath = shiftGP.path;
    Vector shiftDir = -shiftPath.edge(currCameraEdge)->d;
    Float shiftDist = shiftPath.edge(currCameraEdge)->length;
    Ray shiftRay(shiftPath.vertex(currCameraEdge)->getPosition(),
                 shiftDir, Epsilon, shiftDist, 0.f);

    // TODO: Only work if the two rays are parallel, see how to change the formulation otherwise
    // Determine the local vector v
    Vector v = plane->w0() * bRec.t0 + plane->w1() * bRec.t1;

    // Found the new origin
    Point newPlaneOri = shiftRay(bRec.tCam) - v;

    // Maintain the photon distribution at the shift gather point to be the same as the base
    int currVertex = plane->edgeID;
    const Path *source = plane->path;
    const PathVertex *basePhoton = source->vertex(currVertex);
    const PathVertex *parentPhoton = source->vertex(currVertex - 1);
    const PathEdge *parentEdge = source->edge(currVertex - 1);

    // TODO: Do the cosine checking of the parents if it is on a surface
    // TODO: Check the visibility of the new place vs the parent photon

    // Recompute the associated _flux of the unchanged path part
    Spectrum photonWeight(1.f);
    for (int i = 0; i < currVertex - 1; ++i) {
      photonWeight *= source->vertex(i)->weight[EImportance] *
          source->vertex(i)->rrWeight *
          source->edge(i)->weight[EImportance];
    }

    // Query the diffuse shift computation
    Vector dProj = newPlaneOri - parentPhoton->getPosition();
    Float lProj = dProj.length();
    dProj /= lProj;
    ShiftRecord sRec;
    {
      Intersection fakeInter;
      fakeInter.p = newPlaneOri;
      // SH: is this correct because there might still have vertex that is on surface?
      // In this case, give the parent parent location
      diffuseReconnection(sRec, fakeInter, dProj, lProj,
                          parentPhoton, parentEdge,
                          true,
                          currVertex >= 3 ? source->vertex(currVertex - 2)->getPosition() : Point(1.f));

      if (sRec.pdf == Float(0)) {
        result.weight = 1.0f;
        return false;
      }
    }

    ////////////////////////
    // 1) Reconnection to the parent vertex
    ////////////////////////
    // Update the jacobian/_flux with the diffuse shift
    result.jacobian *= sRec.jacobian;
    photonWeight *= sRec.throughtput;

    // Compute the photon plane contrib
    // Compute the transmittance of the camera
    // But also used for phase function eval
    Ray rayCamToPlane = shiftRay;
    rayCamToPlane.mint = 0.f;
    rayCamToPlane.maxt = bRec.tCam;
    MediumSamplingRecord mRecCam;
    medium->eval(rayCamToPlane, mRecCam);

    ////////////////////////
    // 2) Re-evaluate phase function new origin
    ////////////////////////
    // Compute the shift throughput on the new plane
    // origin vertex (rrWeights)
    PhaseFunctionSamplingRecord pRec(mRecCam,
                                     normalize(parentPhoton->getPosition() - basePhoton->getPosition()),
                                     plane->w0());
    PhaseFunctionSamplingRecord pRecNew(mRecCam,
                                        -dProj,
                                        plane->w0());
    photonWeight *= mRecCam.getPhaseFunction()->eval(pRecNew) * (mRecCam.sigmaS);
    photonWeight /= mRecCam.getPhaseFunction()->pdf(pRec);
    photonWeight *= basePhoton->rrWeight;

    ////////////////////////
    // 3) Re-compute the photon plane contrib
    ////////////////////////
    // Note: This is a hack to remove the flux here
    PhaseFunctionSamplingRecord pRecCam(mRecCam, -plane->w1(), -shiftRay.d, EImportance);
    Float phaseTerm = mRecCam.getPhaseFunction()->eval(pRecCam);

    // The transmittance value of the first edge
    Ray rayTrans0(newPlaneOri, plane->w0(), 0.f);
    rayTrans0.mint = 0.f;
    rayTrans0.maxt = bRec.t0;
    MediumSamplingRecord mRec0;
    medium->eval(rayTrans0, mRec0);

    // The transmittance value of the second edge
    Ray rayTrans1(rayTrans0(bRec.t0), plane->w1(), 0.f);
    rayTrans1.mint = 0.f;
    rayTrans1.maxt = bRec.t1;
    MediumSamplingRecord mRec1;
    medium->eval(rayTrans1, mRec1);

    photonWeight *= mRecCam.transmittance
        * (mRecCam.sigmaS) * (mRec0.sigmaS)
        * phaseTerm;

    // Transmittance of the current path
    // This terms cancel out to gives only sigmaT
    photonWeight *= mRec1.transmittance * mRec0.transmittance;
    photonWeight /= mRec0.pdfFailure;
    photonWeight /= mRec1.pdfFailure;

    // Use the same jacobian inside the paper
    photonWeight *= std::abs(bRec.invDet) * plane->length0() * plane->length1();

    result.shiftedFlux = photonWeight;
    result.weight = 0.5f;

    if (config.useMIS) {
      SAssert(currVertex - 1 > 0);

      Float basePdf = source->vertex(currVertex - 1)->pdf[EImportance];
      basePdf *= source->edge(currVertex - 1)->pdf[EImportance];
      basePdf *= mRecCam.getPhaseFunction()->pdf(pRec);

      Float offsetPdf = sRec.pdf;
      offsetPdf *= mRecCam.getPhaseFunction()->pdf(pRecNew);

      // TODO: Here the jacobian have not been taking into account
      // As we are in a easy case

      if (offsetPdf == Float(0) || basePdf == Float(0)) {
        //SLog(EWarn, "Invalid path");
        result.weight = 1.0f;
        return false;
      }

      const Float sensorPart = shiftGP.sensorMIS(currCameraEdge, *baseGather,
                                                 bRec.tCam, bRec.tCam);
      result.weight = 1.0f / (1.0f + sensorPart * result.jacobian * (offsetPdf / basePdf));
    }
  }

#define BETTERSHIFT 0
  /**
   * This shift lock w0 plane basis vector
   * and its origin as these edge is produced
   * by a specular interaction
   */
  bool specularShift(const ShiftGatherPoint &shiftGP,
                     const PhotonPlane::IntersectionRecord &bRec,
                     const LTPhotonPlane *plane,
                     const Spectrum &baseContrib,
                     GradientSamplingResult &result) {
    const Path &shiftPath = shiftGP.path;
    Vector shiftDir = -shiftPath.edge(currCameraEdge)->d;
    Float shiftDist = shiftPath.edge(currCameraEdge)->length;
    Ray shiftRay(shiftPath.vertex(currCameraEdge)->getPosition(),
                 shiftDir, Epsilon, shiftDist, 0.f);

    // Compute the base vector W1Ortho' that is orthogonal to W0
    // Note that W1Othro' and W0 form a plane that contains
    // the new intersection point
    Point newIntersection = shiftRay(bRec.tCam);

    // Construct orthogonal planes
    Vector
        orthNewW1 = newIntersection - (plane->ori() + plane->w0() * dot(newIntersection - plane->ori(), plane->w0()));
    orthNewW1 /= orthNewW1.length();

#if BETTERSHIFT
    // Try to copy a valid plane here
    // By deducing the new local parametrisation.
    // Note that this technique doesn't keep the local dot product constant
    // So we need to re-evaluate it later if this shift is succesful.
    Point oldIntersection = baseCameraRay(bRec.tCam);
    Vector orthW1 = oldIntersection - (plane->ori() + plane->w0()*dot(oldIntersection - plane->ori(), plane->w0()));
    orthW1 /= orthW1.length();
    Float oldAngle = dot(plane->w1(), orthW1);
    Float w0Dot = dot(plane->w0(), plane->w1());

    Vector newW1 = orthNewW1*oldAngle + w0Dot*plane->w0();
    newW1 /= newW1.length();
#else
    // From W1Ortho' compute W1'
    // by copying the local dot product of the orignal W0 and W1 vectors
    Float w0Dot = dot(plane->w0(), plane->w1());
    // sqrt(1 - (w0Dot*w0Dot)) -> pythagore
    Vector newW1 = sqrt(1 - (w0Dot * w0Dot)) * orthNewW1 + plane->w0() * w0Dot;
#endif

    // Recheck the intersection to be sure that
    // t0New and t1New are valid (> 0)
    Float t0New, t1New, tCamNew, invDetNew;
    if (!intersection(shiftRay, plane->ori(), plane->w0(), newW1,
                      tCamNew, t0New, t1New, invDetNew)) {
      result.weight = 1.0;
      return false;
    }

    // Remark: This code is a little bit complex and different to other shifts
    // We reuse the baseContrib to avoid to do some re-computation.

    // Do the computation
    // Original PDF
    Ray rayTrans1(plane->ori(), plane->w1(), 0.f);
    rayTrans1.mint = 0.f;
    rayTrans1.maxt = bRec.t1;
    MediumSamplingRecord mRec1;
    medium->eval(rayTrans1, mRec1);

    Ray rayTrans0(plane->ori(), plane->w0(), 0.f);
    rayTrans0.mint = 0.f;
    rayTrans0.maxt = bRec.t0;
    MediumSamplingRecord mRec0;
    medium->eval(rayTrans0, mRec0);

    // Shift PDF
    Ray rayTrans1Shift(plane->ori(), newW1, 0.f);
    rayTrans1Shift.mint = 0.f;
    rayTrans1Shift.maxt = t1New;
    MediumSamplingRecord mRec1Shift;
    medium->eval(rayTrans1Shift, mRec1Shift);

    Ray rayTrans0Shift(plane->ori(), plane->w0(), 0.f);
    rayTrans0Shift.mint = 0.f;
    rayTrans0Shift.maxt = t0New;
    MediumSamplingRecord mRec0Shift;
    medium->eval(rayTrans0Shift, mRec0Shift);

    // Compute the shift throughput
    Spectrum throughputShift = baseContrib;
    throughputShift *= mRec0Shift.transmittance / mRec0.transmittance;
    throughputShift *= mRec1Shift.transmittance / mRec1.transmittance;
    throughputShift /= plane->invJacobian(baseCameraRay.d);
    throughputShift *= 1.0 / absDot(plane->w0(), cross(newW1, shiftRay.d));

    // The jacobian is to:
    // - decoupling the variables
    // - recoupling the variables
    result.jacobian = plane->invJacobian(baseCameraRay.d);
    result.jacobian *= absDot(plane->w0(), cross(newW1, shiftRay.d));
    result.jacobian /= t1New / bRec.t1;

    if (plane->edgeID != 1) {
      result.jacobian /= t0New / bRec.t0;
    }

    PhaseFunctionSamplingRecord pRec(mRec0Shift, // Not used
                                     -plane->w1(),
                                     -baseCameraRay.d);
    PhaseFunctionSamplingRecord pRecNew(mRec0Shift,// Not used
                                        -newW1,
                                        -shiftRay.d);
    throughputShift *= mRec0Shift.getPhaseFunction()->eval(pRecNew); // Sigma S already inside
    throughputShift /= mRec0Shift.getPhaseFunction()->eval(pRec); // Remove the base contrib here

#if BETTERSHIFT
    PhaseFunctionSamplingRecord pRecW1(mRec0Shift, // Not used
                                     plane->w0(), plane->w1());
    PhaseFunctionSamplingRecord pRecW1New(mRec0Shift,// Not used
                                       plane->w0(), newW1);
    throughputShift *= mRec0Shift.getPhaseFunction()->eval(pRecW1New);
    throughputShift /= mRec0Shift.getPhaseFunction()->pdf(pRecW1New);
#endif

    result.weight = 0.5f;
    result.shiftedFlux = result.jacobian * throughputShift;

    if (config.useMIS) {
      ///////////////////////// BASE PDF
      // * Eye subpath
      Float basePdf = mRec0.pdfSuccess;
      basePdf *= mRec1.pdfSuccess;
      basePdf *= mRec0Shift.getPhaseFunction()->pdf(pRec);
      //basePdf /= plane->invJacobian(baseCameraRay.d); // Jacobian

      ///////////////////////// SHIFT PDF
      // * Eye subpath
      //offsetPdf *= absDot(plane->w0(), cross(newW1, shiftRay.d));
      Float offsetPdf = mRec0Shift.pdfSuccess;// * GOpNew0;
      offsetPdf *= mRec1Shift.pdfSuccess;// * GOpNew1;
      offsetPdf *= mRec0Shift.getPhaseFunction()->pdf(pRecNew);

#if BETTERSHIFT
      basePdf *= mRec0Shift.getPhaseFunction()->pdf(pRecW1);
      offsetPdf *= mRec0Shift.getPhaseFunction()->pdf(pRecW1New);
#endif

      if (offsetPdf == Float(0) || basePdf == Float(0)) {
        SLog(EWarn, "Invalid path: offset = %f | base = %f", offsetPdf, basePdf);
        result.weight = 1.0f;
        return false;
      }

      const Float sensorPart = shiftGP.sensorMIS(currCameraEdge, *baseGather, bRec.tCam, bRec.tCam);
      result.weight = 1.0f / (1.0f + sensorPart * result.jacobian * offsetPdf / basePdf);

    }

    return true;

  }

  /**
   * This is the same code as the intersection code inside photon plane structure
   * Two things changes here:
   *  - w1 and w0 are normalized, so no need to renormalize invDet or t0/t1
   *  - remove the condition on t0 > 1 or t1 > 1:
   *     * as these values are not expressed in the local space anymore
   *     * we don't want to miss the shift due to a too short photon plane
   */
  bool intersection(const Ray &ray_, const Point &ori, const Vector &w0, const Vector &w1,
                    Float &tCam, Float &t0, Float &t1,
                    Float &invDet) const {
    // Sanity checks
    // TODO: Remove these check when the shift will be correct
    SAssert(fabs(w1.length() - 1.0) < 10e-5);
    SAssert(fabs(w0.length() - 1.0) < 10e-5);

    Vector P = cross(ray_.d, w1);
    Float det = dot(w0, P);
    if (std::abs(det) < 1e-8f)
      return false;

    invDet = 1.0f / det;
    Vector T = ray_.o - ori;
    t0 = dot(T, P) * invDet;
    if (t0 < 0.0f)
      return false;

    Vector Q = cross(T, w0);
    t1 = dot(ray_.d, Q) * invDet;
    if (t1 < 0.0f)
      return false;

    tCam = dot(w1, Q) * invDet;
    return !(tCam <= ray_.mint || tCam >= ray_.maxt);

  }
};

MTS_NAMESPACE_END

#endif //MITSUBA_SHIFT_VOLUME_PLANES_H
