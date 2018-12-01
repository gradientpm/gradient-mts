#pragma once

#include <mitsuba/render/scene.h>
#include <mitsuba/bidir/mempool.h>

#include "../gvpm_struct.h"
#include "shift_utilities.h"

#ifndef MITSUBA_SHIFT_CAMERAPATH_H
#define MITSUBA_SHIFT_CAMERAPATH_H

MTS_NAMESPACE_BEGIN

class ShiftGatherPoint : public GatherPoint {
private:
  struct SVertexMeasure {
    EMeasure measure;
    size_t vertex_id;
  };

public:
  bool generated;
  bool surfaceValid;

  ShiftGatherPoint(): generated(false), surfaceValid(false) {
  }
  virtual ~ShiftGatherPoint() {}

  bool generate(Scene *scene, MemoryPool &pool,
                const GatherPoint &baseGP,
                const Point2 offsetPixel,
                bool emissionShift = false) {
    // Nothing to do
    if (generated)
      return false;
    generated = true;
    pixel = Point2i((int) offsetPixel.x, (int) offsetPixel.y);

    // Generate the shift path
    const Path &basePath = baseGP.path;
    surfaceValid = trace(scene, offsetPixel, pool, emissionShift, basePath);

    Float accumulatedGShift = 1.f;
    Float accumulatedJacobianCorrection = 1.f;
    // Compute the Jacobian and update the pdf in the surface space.
    for (size_t i = 2; i < path.vertexCount(); i++) {

      // FIXME: Becarefull if the code below is uncommented.
      // FIXME: If might produce wrong results...
      // We skip null intersection vertices
      // However, for the last vertex, we always do this computation
//      if(i != (path.vertexCount() - 1) && path.vertex(i)->isNullInteraction()) {
//        continue; // Skip this vertex as we need to ignore it
//      }

      SVertexPDF &currInfo = getVertexInfo(i - 1);

      auto previous_vertex = [&scene, i, &basePath]() -> SVertexMeasure {
        // Search the correct vertex (first)
        size_t t = i-1;
//        while(basePath.vertex(t)->isNullInteraction() && t > 1) { t--; };
        // Depending on the type of vertex selected, gather the measure information
        if(t == 1) {
          auto sensor = scene->getSensor();
          return {(sensor->getType() & Sensor::EDeltaDirection ? EDiscrete : ESolidAngle), t};
        } else {
          auto previous_valid_vertex = basePath.vertex(t);
          if(previous_valid_vertex->isNullInteraction()) {
            return {EDiscrete, t};
          } else {
            return {((previous_valid_vertex->getComponentType() & BSDF::EDelta) ? EDiscrete : ESolidAngle), t};
          }
        }
      }();

      switch (previous_vertex.measure) {
      case ESolidAngle: {
        // SH: GOpNew here has to account only to the beam intersection point?

        const Float GOpNew = geometryOpposingTerm(path.vertex(previous_vertex.vertex_id)->getPosition(),
                                                  path.vertex(i)->getPosition(),
                                                  path.vertex(i)->getGeometricNormal());
        const Float GOpBase = geometryOpposingTerm(basePath.vertex(previous_vertex.vertex_id)->getPosition(),
                                                   basePath.vertex(i)->getPosition(),
                                                   basePath.vertex(i)->getGeometricNormal());
        accumulatedGShift *= GOpNew;

        // Strictly speaking, we need
        // currInfo->jacobian *= GOpBase / GOpNew;     // convert Jacobian from solid angle measure to area measure
        // currInfo->weight *= GOpNew / GOpBase;       // evaluation in area measure
        // then weight *= jacobian
        // which results in the ratio of geometry term canels out
        //
        // Therefore, we simply do not update the weight and vertexWeight here.
        //
        // We still update the Jacobian as it is needed for MIS
        // (with some correction for the last segment pending during density estimation).
        accumulatedJacobianCorrection *= GOpBase / GOpNew;
        break;
      }
      case EDiscrete: {
        // Nothing for now as the measure
        // if ill defined in this case
        break;
      }
      default: {
        SLog(EError, "Not possible measure supported");
        break;
      }
      }

      // FIXME: This make creates some float inaccuracy
      // Apply the correction
      currInfo.pdf *= accumulatedGShift;
      currInfo.jacobian *= accumulatedJacobianCorrection;
    }

    radius = baseGP.radius;
    pixel = Point2i((int) offsetPixel.x, (int) offsetPixel.y);
    sampledComponent = baseGP.sampledComponent;
    pdfComponent = baseGP.pdfComponent;
    depth = (int) basePath.vertexCount() - 2;

    if (emissionShift) {
      const Intersection &shiftIts = lastIts();
      if (shiftIts.isEmitter()) {
        if (!baseGP.pureSpecular) {
          emission += surfInfo().weight * shiftIts.Le(shiftIts.toWorld(shiftIts.wi));
          currEmission += surfInfo().weight * shiftIts.Le(shiftIts.toWorld(shiftIts.wi));
        }
      }
    }

    return true;

  }

  inline bool validVolumeEdge(size_t idEdge, const Medium *med) const {
    if (path.edgeCount() <= idEdge)
      return false;
    const PathEdge *e = path.edge(idEdge);
    return e->medium == med;
  }

private:
  /// This function trace the shifted path
  /// by looking to the base path
  /// Note that Jacobian and PDF are not properly updated yet
  bool trace(Scene *scene, const Point2 offsetPixel,
             MemoryPool &pool, bool emissionShift,
             const Path &basePath) {
    Path &shiftPath = path;
    clear();

    // Find primary hit point
    Sensor *sensor = scene->getSensor();
    shiftPath.initialize(scene, 0.f, ERadiance, pool);

    // Sample the primary ray from the camera
    Point2i pixelID((int) offsetPixel.x, (int) offsetPixel.y);
    PathVertex *v1 = pool.allocVertex(), *v2 = pool.allocVertex();
    PathEdge *e0 = pool.allocEdge(), *e1 = pool.allocEdge();
    int t = shiftPath.vertex(0)->sampleSensor(scene, nullptr, pixelID, e0, v1, e1, v2, true, offsetPixel, true);
    if (t < 1 || t < 2) {
      pool.release(e0);
      pool.release(v1);
      pool.release(e1);
      pool.release(v2);
      shiftPath.release(pool);
      return false;
    }
    // Here we assume that e1 always have the same parametrisation than the base path
    shiftPath.append(e0, v1);
    shiftPath.append(e1, v2);


    // Find the measure depending of the sensor type
    EMeasure measure = sensor->getType() & Sensor::EDeltaDirection ? EDiscrete : ESolidAngle;
    // Note: If is solid angle, the pdf will be converted into area during the 1st loop iteration

    // Cache vertex info for shifted gather point
    {
      SVertexPDF *currInfo = createInfo();        // vertex 0
      currInfo->weight = Spectrum(1.f);
      currInfo->vertexWeight = Spectrum(1.f);
      currInfo->pdf = 1.0f;
      currInfo->jacobian = 1.0f;
    }
    SVertexPDF *currInfo = createInfo();                    // vertex 1

    // Now handle the shift of the primary hit point (vertex 2)
    // The Jacobian of the shift is stored at vertex 1.

    // We use pdf ratio to calculate the Jacobian term of the first hit point in area measure
    Float pdf1, pdf2;
    Point org = basePath.vertex(1)->getPosition();
    {
      Point first = basePath.vertex(2)->getPosition(); // first hit point;
      PositionSamplingRecord pRecSrc = basePath.vertex(1)->getPositionSamplingRecord();
      DirectionSamplingRecord dRecSrc;
      dRecSrc.d = normalize(first - org);
      dRecSrc.measure = measure;
      pdf1 = sensor->pdfDirection(dRecSrc, pRecSrc);
    }

    {
      PositionSamplingRecord pRecDst;
      pRecDst.measure = EArea;
      pRecDst.time = 0.0f;
      DirectionSamplingRecord dRecDst;
      dRecDst.d = -e1->d;
      dRecDst.measure = measure;
      pdf2 = sensor->pdfDirection(dRecDst, pRecDst);

      currInfo->pdf = pdf2;
      currInfo->weight = sensor->evalDirection(dRecDst, pRecDst);
    }

    // At boundary, pdf of the shifted point is zero.
    // We assume the Jacobian and weight to be 1 in such cases.
    if (pdf2 == 0.0f) {
      // FIXME: is this assumption correct?
      currInfo->jacobian = 1.0f;
      currInfo->weight = Spectrum(1.0f);
      currInfo->vertexWeight = Spectrum(1.0f);

      if (pdf1 == 0.f)
        SLog(EError, "Pdf1 0 prob. on a sampled ray");
      currInfo->pdf = pdf1;

      // SH: should we just return an empty camera path here?

    } else {
      // Strictly speaking, we must evaluate
      // currInfo->weight /= pdf1;
      // currInfo->weight *= currInfo->jacobian;
      // which eventually boils down to
      currInfo->weight /= pdf2;
      currInfo->vertexWeight = currInfo->weight;
      // which is equivalent to (assume perspective camera and importance sampling on sensor)
      //currInfo->weight = Spectrum(1.f);
      //currInfo->vertexWeight = Spectrum(1.f);

      currInfo->jacobian = pdf1 / pdf2;   // solid angle measure, will be converted into area measure in the first loop
    }

    const Medium *currMed = e1->medium;
    for (size_t i = 2; i < basePath.vertexCount(); ++i) {
      const PathVertex *baseVertex = basePath.vertex(i);

      // Update the weight based on the edge
      // FIXME: This is not necessary, no?
      if (currMed != basePath.edge(i - 1)->medium) {
        surfaceValid = false;
        return false;
      }
      if (basePath.edge(i - 1)->medium != nullptr) {
        // because camera path uses long beam, there is no sampling on the edge,
        // so the weight is equal to the transmittance
        currInfo->weight *= shiftPath.edge(i - 1)->weight[ERadiance];
      }

      if(baseVertex->isMediumInteraction() || shiftPath.vertex(i)->isMediumInteraction()) {
        // This can happens when we get a too thick participating media.
        // In this case we need to stop the shift operation as everything
        // is already done
        surfaceValid = false;
        return false; // TODO: Do we need to recompute some information on the shift edge???
        // TODO: Maybe we need to check if we was able to sample or not.
      }

      const Intersection &baseIts = baseVertex->getIntersection();
      const BSDF *baseBSDF = baseIts.getBSDF();
      const Intersection &shiftIts = shiftPath.vertex(i)->getIntersection();
      const BSDF *shiftedBSDF = shiftIts.getBSDF();
      // This check is not applicable when shifting the last vertex on a complete path
      if (!(emissionShift && i == basePath.vertexCount() - 1)) {
        // Deny shifts between Dirac and non-Dirac BSDFs
        bool bothDelta, bothSmooth;
        if (baseVertex->componentType != 0) {
          bothDelta = (baseVertex->componentType & BSDF::EDelta) && (shiftedBSDF->getType() & BSDF::EDelta);
          bothSmooth = (baseVertex->componentType & BSDF::ESmooth) && (shiftedBSDF->getType() & BSDF::ESmooth);
        } else {
          bothDelta = (baseBSDF->getType() & BSDF::EDelta) && (shiftedBSDF->getType() & BSDF::EDelta);
          bothSmooth = (baseBSDF->getType() & BSDF::ESmooth) && (shiftedBSDF->getType() & BSDF::ESmooth);
        }
        if (!(bothDelta || bothSmooth)) {
          return false;
        }
      }

      if (shiftedBSDF->getType() != baseBSDF->getType()) {
        return false;
      }

      if (i == (int) basePath.vertexCount() - 1)
        break;

      // Perform half vector copy
      currInfo = createInfo();

      // --- Compute the direction incident direction
      PathVertex *predVertex = basePath.vertex(i - 1);
      PathVertex *succVertex = basePath.vertex(i + 1);
      bool almostSameNormal = (shiftIts.shFrame.n - baseIts.shFrame.n).lengthSquared() < 0.01;
      Vector3 shiftWi;
      if (almostSameNormal) {
        shiftWi = shiftIts.toLocalCoherent(shiftPath.edge(i - 1)->d);
      } else {
        shiftWi = shiftIts.toLocal(shiftPath.edge(i - 1)->d);
      }

      // Compute Jacobian and outgoing direction
      Vector3 shiftWo;
      Float shiftJacobian = 1.0;
      if (baseVertex->isNullInteraction()) {
        // The case of null interaction is straight forward
        // We just need to keep the same direction
        // The rest of the computation is degenerated.
        shiftWo = -shiftWi;
      } else {
        Vector baseWi = normalize(predVertex->getPosition() - baseVertex->getPosition());
        Vector baseWo = normalize(succVertex->getPosition() - baseVertex->getPosition());

        HalfVectorShiftResult shiftResult;
        if (almostSameNormal) {
          shiftResult = halfVectorShift(baseIts.toLocalCoherent(baseWi),
                                        baseIts.toLocalCoherent(baseWo),
                                        shiftWi,
                                        baseIts.getBSDF()->getEta(), shiftedBSDF->getEta());
        } else {
          shiftResult = halfVectorShift(baseIts.toLocal(baseWi),
                                        baseIts.toLocal(baseWo),
                                        shiftWi,
                                        baseIts.getBSDF()->getEta(), shiftedBSDF->getEta());
        }

        if (baseVertex->componentType & BSDF::EDelta)
          shiftResult.jacobian = Float(1);
        if (!shiftResult.success) {
          return false;
        }
        shiftJacobian = shiftResult.jacobian;
        shiftWo = shiftResult.wo;
      }

      // Evaluate BSDF at the new vertex
      BSDFSamplingRecord bRec(shiftIts, shiftWi, shiftWo);
      bRec.component = baseVertex->sampledComponentIndex;
      EMeasure measure = (baseVertex->getComponentType() & BSDF::EDelta) ? EDiscrete : ESolidAngle;
      if (bRec.component >= shiftedBSDF->getComponentCount())
        SLog(EWarn,
             "Invalid component request %d",
             bRec.component);

      Spectrum bsdfValue = shiftedBSDF->eval(bRec, measure);

      Float shiftPdfWo = shiftedBSDF->pdf(bRec, measure);
      shiftPdfWo *= shiftedBSDF->pdfComponent(bRec);
      currInfo->pdf *= shiftPdfWo;
      if (currInfo->pdf == Float(0)) {
        return false;
      }

      // Account for the probability of the next vertex on base path
      // Turn existing area density into solid angle
      // Pdf of the component already accounted in pdf[2] array.
      Float baseWoPdf;
      if (baseVertex->measure != EDiscrete) {
        baseWoPdf = baseVertex->pdf[ERadiance] / geometryOpposingTerm(baseVertex->getPosition(),
                                                                      succVertex->getPosition(),
                                                                      succVertex->getGeometricNormal());
      } else {
        baseWoPdf = baseVertex->pdf[ERadiance];
      }
      currInfo->weight *= baseVertex->rrWeight;
      currInfo->weight *= bsdfValue / baseWoPdf;
      currInfo->vertexWeight = bsdfValue / baseWoPdf;
      currInfo->vertexWeight *= baseVertex->rrWeight;

      // Account Jacobian
      currInfo->weight *= shiftJacobian;
      currInfo->vertexWeight *= shiftJacobian;
      currInfo->jacobian *= shiftJacobian;

      // Prepare for next tracing step
      Vector3 outgoingDirection;
      if (almostSameNormal) {
        outgoingDirection = shiftIts.toWorldCoherent(shiftWo);
      } else {
        outgoingDirection = shiftIts.toWorld(shiftWo);
      }

      // Check the medium, if needed, update the current medium/
      const Medium *nextMed = basePath.edge(i)->medium;
      if (shiftIts.isMediumTransition())
        currMed = shiftIts.getTargetMedium(outgoingDirection);
      if (currMed != nextMed) {
        return false;
      }

      // Generate the new ray and compute the intersection
      Ray ray(shiftIts.p, outgoingDirection, 0.f);
      PathEdge *ei = pool.allocEdge();
      PathVertex *vi = pool.allocVertex();
      ei->medium = currMed; // assign the medium to the next edge that will be generated
      if (!ei->sampleNext(scene, nullptr, shiftPath.vertex(i), ray, vi, ERadiance, true)) {
        // camera path is long beam
        pool.release(ei);
        pool.release(vi);
        return false;
      }
      shiftPath.append(ei, vi);
    }
    return true;
  }

};

MTS_NAMESPACE_END

#endif //MITSUBA_SHIFT_CAMERAPATH_H
