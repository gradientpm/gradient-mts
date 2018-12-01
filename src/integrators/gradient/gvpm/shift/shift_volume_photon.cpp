#include "shift_volume_photon.h"

// For the statistics
#include <mitsuba/core/statistics.h>

// Diffuse shift
#include "operation/shift_diffuse.h"

// For ME shift
#include "operation/shift_ME.h"
#include <mitsuba/bidir/manifold.h>

// For Medium shift
#include "operation/shift_medium.h"

MTS_NAMESPACE_BEGIN

StatsCounter MEMisSuccessVol("GPM",
                             "MIS computation succesfull", EPercentage);

StatsCounter MEShiftVol("GVPM",
                        "Percentage ME Vol Shift    : ", EPercentage);
StatsCounter HVShiftVol("GVPM",
                        "Percentage HVCopy Vol Shift: ", EPercentage);
StatsCounter DiffShiftVol("GVPM",
                          "Percentage Diff. Vol Shift : ", EPercentage);
StatsCounter DiffMediumVol("GVPM",
                           "Percentage Medium Vol Shift : ", EPercentage);
StatsCounter InvalidShiftVol("GVPM",
                             "Percentage Inval Vol Shift : ", EPercentage);
StatsCounter noVisShiftVol("GVPM",
                           "Percentage non visible Vol Shift : ", EPercentage);
StatsCounter noPDFVol("GVPM",
                           "Percentage 0 PDF vol : ", EPercentage);

StatsCounter shiftNullVol("GVPM",
                          "Percentage of null shift : ", EPercentage);
StatsCounter additionalShiftVol("GVPM",
                                "Percentage of additional shift: ", EPercentage);
StatsCounter cancelShortDistance("GVPM",
                                 "Deny due to shorter shift camera beam: ", EPercentage);

#define DEBUG_MIS 0

#if DEBUG_MIS
StatsCounter avgMISWeightME("GVPM", "Avg MIS weight ME", EAverage);
#endif

bool AbstractVolumeGradientRecord::shiftPhoton(const Point &offsetPos, const Path *source, size_t currVertex,
                                               const ShiftGatherPoint &shiftGP,
                                               const Ray &shiftRay, const MediumSamplingRecord& shiftMRec,
                                               GradientSamplingResult &result, Float photonRadius,
                                               Float pdfBaseRay, Float pdfShiftRay,
                                               Float additionalJacobian) {
  shiftNullVol.incrementBase();

  // Deduce the shift type
  int b;
  ELightShiftType type = getTypeShift(source, currVertex, b);

  MEShiftVol.incrementBase();
  HVShiftVol.incrementBase();
  DiffShiftVol.incrementBase();
  DiffMediumVol.incrementBase();
  InvalidShiftVol.incrementBase();

  if (type == EInvalidShift) {
    // Error?
    ++InvalidShiftVol;
    return false;
  } else if (type == EDiffuseShift) {
    // Just use diffuse reconnection here
    ++DiffShiftVol;
    return shiftPhotonDiffuse(offsetPos, source, currVertex,
                              shiftGP, shiftRay, shiftMRec, result, pdfBaseRay, pdfShiftRay,
                              additionalJacobian);
  } else if (type == EMediumShift) {
    if (config.noMediumShift) {
      ++DiffShiftVol;
      return shiftPhotonDiffuse(offsetPos, source, currVertex,
                                shiftGP, shiftRay, shiftMRec, result, pdfBaseRay, pdfShiftRay,
                                additionalJacobian);
    } else {
      ++DiffMediumVol;
      return shiftPhotonMedium(offsetPos, source, currVertex,
                               shiftGP, shiftRay, shiftMRec, result, pdfBaseRay, pdfShiftRay,
                               additionalJacobian);
    }
  } else if (type == EManifoldShift) {
    if (!config.useManifold) {
      ++InvalidShiftVol;
      return false; // Impossible to handle this path
    } else {
      ++MEShiftVol;
      return shiftPhotonManifold(b, currVertex, offsetPos, source,
                                 shiftGP, shiftRay, shiftMRec, result, photonRadius,
                                 pdfBaseRay, pdfShiftRay,
                                 additionalJacobian);
    }
  }
  SLog(EError, "Invalid type");
  ++InvalidShiftVol;
  return false;

}

bool AbstractVolumeGradientRecord::shiftNull(const Spectrum &photonFlux,
                                             const Vector3 &photonWi,
                                             const ShiftGatherPoint &shiftGP,
                                             const Ray &shiftRay,
                                             const MediumSamplingRecord& shiftMRec,
                                             GradientSamplingResult &resultNull,
                                             Float pdfBaseRay,
                                             Float pdfShiftRay,
                                             Float additionalJacobian) {
  shiftNullVol.incrementBase();
  ++shiftNullVol;

  // In the volume, the shift is always correct
  Spectrum contrib = getVolumePhotonContrib(photonFlux,shiftMRec,
                                            photonWi, -shiftRay.d, EImportance);
  Spectrum eyeShiftContrib = shiftGP.getWeightBeam(currEdge - 1) * // cummulative vertex(i-1).weight * edge(i-1).weight
      shiftGP.getWeightVertex(currEdge);    // vertex(i).weight
  resultNull.jacobian *= additionalJacobian;
  resultNull.shiftedFlux = shiftMRec.transmittance * contrib * eyeShiftContrib * resultNull.jacobian;
  resultNull.weight = 0.5f;

  // MIS in solid angle measure
  // FIXME: if parentIndex == 0, should not convert to solid angle
  if (config.useMIS) {

    Float basePdf = pdfBaseRay;
    Float offsetPdf = pdfShiftRay;
    if (offsetPdf == Float(0) || basePdf == Float(0)) {
      //SLog(EWarn, "Invalid path");
      resultNull.weight = 1.0f;
      return false;
    }

    const Float sensorPart = shiftGP.sensorMIS(currEdge, *baseGather,
                                               shiftRay.maxt, baseRay.maxt);
    resultNull.weight = 1.0f / (1.0f + sensorPart * offsetPdf * resultNull.jacobian / basePdf);
  }

  return true;
}

bool AbstractVolumeGradientRecord::shiftPhotonManifold(int b, int c,
                                                       const Point &offsetPos, const Path *lt,
                                                       const ShiftGatherPoint &shiftGP,
                                                       const Ray &shiftRay, const MediumSamplingRecord& shiftMRec,
                                                       GradientSamplingResult &result, Float photonRadius,
                                                       Float pdfBaseRay, Float pdfShiftRay,
                                                       Float additionalJacobian) {


  // Perform a walk through a specular chain to a non-specular vertex
  const Path &source = (Path &) *lt;
  Path proposal;
  {
    // Construct the new vertex
    PathVertex shiftVertex;
    memset(&shiftVertex, 0, sizeof(PathVertex));
    MediumSamplingRecord &cacheMRec = shiftVertex.getMediumSamplingRecord();
    cacheMRec.t = 0.f;
    cacheMRec.p = offsetPos; // Get the original pos.
    cacheMRec.medium = medium;
    shiftVertex.type = PathVertex::EMediumInteraction;
    shiftVertex.measure = EArea;
    shiftVertex.sampledComponentIndex = -1; // Because in media no component sampling

    if (!generateShiftPathME(source, proposal, b, c, thdata.pool, thdata.offsetGenerator.get(),
                             shiftVertex, photonRadius * config.relaxME,
                             baseRay(baseRay.maxt),
                             shiftRay(shiftRay.maxt))) {
      for (int i = b; i <= c; ++i) {
        thdata.pool.release(proposal.edge(i - 1));
        thdata.pool.release(proposal.vertex(i));
      }
      result.weight = 1.0f;
      return false;
    }
  }

  ShiftRecord sRecME;
  if (!ShiftME(sRecME, source, proposal, b, c, false)) {
    for (int i = b; i <= c; ++i) {
      thdata.pool.release(proposal.edge(i - 1));
      thdata.pool.release(proposal.vertex(i));
    }
    result.weight = 1.0f;
    return false;
  }

  // Evaluate Jacobian (area measure)
  result.jacobian *= sRecME.jacobian * additionalJacobian;  // middle term | dy_s' / dy_s |

  // Jacobian in area integral in original space
  SpecularManifold *manifold = thdata.offsetGenerator->getSpecularManifold();

  // Jacobian computed using det
  const Float detProposed = manifold->det(proposal, b, c);
  const Float detSource = manifold->det(source, b, c);
  result.jacobian *= detProposed / detSource;

  if (result.jacobian <= 0.0 || !std::isfinite(result.jacobian)) {
    SLog(EWarn, "Invalid jacobian %g %d %d", result.jacobian, b, c);
    result.weight = 1.0f;
    return false;
  }

  Spectrum photonWeight = sRecME.throughtput;

  // In the volume, the shift is always correct
  Spectrum contrib = getVolumePhotonContrib(photonWeight,
                                            shiftMRec,
                                            normalize(proposal.vertex(c - 1)->getPosition() - offsetPos),
                                            -shiftRay.d,
                                            EImportance);
  Spectrum eyeShiftContrib = shiftGP.getWeightBeam(currEdge - 1) * // cummulative vertex(i-1).weight * edge(i-1).weight
      shiftGP.getWeightVertex(currEdge);    // vertex(i).weight
  // As every vertex on the camera path is simply surface vertex, the weight

  // Default no MIS weight
  result.weight = 0.5f;
  result.shiftedFlux = shiftMRec.transmittance * contrib * eyeShiftContrib * result.jacobian;

  // Pdf check for offset path
  // FIXME: Double check this code
  Float offsetPdf = sRecME.pdf * pdfShiftRay;
  if (offsetPdf == Float(0)) {
    // This path has zero density due to manifold walk.
    // It cannot be sampled by particle tracing, and hence no reversible.
    result.weight = 1.0f;
    if (!result.shiftedFlux.isZero()) {
      SLog(EWarn, "0 PDF path but with a throughput: %s\n Set to 0", result.shiftedFlux.toString().c_str());
      result.shiftedFlux = Spectrum(0.f);
    }
    // When this assert fails, it could be due to pdfComponent is not correctly implemented.
  }

  if (config.useMIS) {
    MEMisSuccessVol.incrementBase();

    Float basePdf = pdfBaseRay;
    // No need to consider pdf of the unchanged segments
    for (int i = b; i < c; ++i) {
      basePdf *= source.vertex(i)->pdf[EImportance];
      basePdf *= source.edge(i)->pdf[EImportance];
    }

    if (basePdf == Float(0)) {
      SLog(EWarn, "Invalid base path. This case should not happen.");
      result.weight = 0.0f;
    } else {
      const Float sensorPart = shiftGP.sensorMIS(currEdge, *baseGather,
                                                 shiftRay.maxt, baseRay.maxt);
      if (config.powerHeuristic) {
        result.weight = 1.0f / (1.0f + powerOfTwo(sensorPart * result.jacobian * (offsetPdf / basePdf)));
      } else {
        result.weight = 1.0f / (1.0f + sensorPart * offsetPdf * result.jacobian / basePdf);
      }
#if DEBUG_MIS
      avgMISWeightME += result.weight*100;
      avgMISWeightME.incrementBase();
#endif
    }

  } else {
#if DEBUG_MIS
    avgMISWeightME += 50;
    avgMISWeightME.incrementBase();
#endif
  }

  // Clean up
  for (int i = b; i <= c; ++i) {
    thdata.pool.release(proposal.edge(i - 1));
    thdata.pool.release(proposal.vertex(i));
  }

  return true;
}

bool AbstractVolumeGradientRecord::shiftPhotonMedium(const Point &offsetPos,
                                                     const Path *source,
                                                     size_t currVertex,
                                                     const ShiftGatherPoint &shiftGP,
                                                     const Ray &shiftRay,
                                                     const MediumSamplingRecord& shiftMRec,
                                                     GradientSamplingResult &result,
                                                     Float pdfBaseRay,
                                                     Float pdfShiftRay,
                                                     Float additionalJacobian) {
  SAssert(false);
  // Gather all necessary information for the shift
  const PathVertex *a1 = source->vertex(currVertex - 1);
  const PathVertex *a0 = source->vertex(currVertex - 2);
  const PathEdge *e0 = source->edge(currVertex - 2);
  const PathEdge *e1 = source->edge(currVertex - 1);

  // Recompute the associated flux of the unchanged path part
  Spectrum photonWeight = source->vertex(0)->weight[EImportance] *
      source->vertex(0)->rrWeight *
      source->edge(0)->weight[EImportance];
  for (size_t i = 1; i < currVertex - 2; ++i) {
    photonWeight *= source->vertex(i)->weight[EImportance] *
        source->vertex(i)->rrWeight *
        source->edge(i)->weight[EImportance];
  }

  // Query the diffuse shift computation
  ShiftRecord sRec;
  Vector newW1;
  {
    PathVertex *aPred = nullptr;
    if (currVertex - 3 >= 1) {
      aPred = source->vertex(currVertex - 3);
    }
    bool success = mediumRotationShift(sRec, baseRay.d, shiftRay.d, aPred, a0, e0, a1, e1, offsetPos, newW1);
    if (!success || sRec.pdf == 0.f) {
      // Irreversible shift
      result.weight = 1.0f;
      return false;
    }
  }

  result.jacobian *= sRec.jacobian * additionalJacobian;
  photonWeight *= sRec.throughtput;
  Spectrum contrib = getVolumePhotonContrib(photonWeight, shiftMRec,
                                            -newW1, -shiftRay.d, EImportance);

  // We have to use the cache values here because the offset eye path is not explicitly constructed
  // currEdge = i means we are on the way from vertex i to vertex i + 1.
  Spectrum
      eyeShiftContrib = shiftGP.getWeightBeam(currEdge - 1) *    // cummulative vertex(i-1).weight * edge(i-1).weight
      shiftGP.getWeightVertex(currEdge);          // vertex(i).weight

  result.shiftedFlux = shiftMRec.transmittance * contrib * eyeShiftContrib * result.jacobian;
  result.weight = 0.5f;

  if (config.useMIS) {
    SAssert(currVertex - 1 > 0);

    Float basePdf = pdfBaseRay;
    basePdf *= a0->pdf[EImportance];
    basePdf *= e0->pdf[EImportance];
    basePdf *= a1->pdf[EImportance];
    basePdf *= e1->pdf[EImportance];

    Float offsetPdf = sRec.pdf * pdfShiftRay;
    if (offsetPdf == Float(0) || basePdf == Float(0)) {
      result.weight = 1.0f;
      return false;
    }

    const Float sensorPart = shiftGP.sensorMIS(currEdge, *baseGather,
                                               shiftRay.maxt, baseRay.maxt);
    if (config.powerHeuristic) {
      result.weight = 1.0f / (1.0f + powerOfTwo(sensorPart * result.jacobian * (offsetPdf / basePdf)));
    } else {
      result.weight = 1.0f / (1.0f + sensorPart * result.jacobian * (offsetPdf / basePdf));
    }

  }
  return true;

}

bool AbstractVolumeGradientRecord::shiftPhotonDiffuse(const Point &offsetPos, const Path *source, size_t currVertex,
                                                      const ShiftGatherPoint &shiftGP,
                                                      const Ray &shiftRay, const MediumSamplingRecord& shiftMRec,
                                                      GradientSamplingResult &result,
                                                      Float pdfBaseRay, Float pdfShiftRay,
                                                      Float additionalJacobian) {
  // Maintain the photon distribution at the shift gather point to be the same as the base
  const PathVertex *parentPhoton = source->vertex(currVertex - 1);
  const PathEdge *parentEdge = source->edge(currVertex - 1);

  // Check visibility of the current photon shift
  Vector dProj = offsetPos - parentPhoton->getPosition();
  Float lProj = dProj.length();
  dProj /= lProj;
  Ray projRay(parentPhoton->getPosition(), dProj, Epsilon, lProj * ShadowEpsilon, 0.f);
  noVisShiftVol.incrementBase();
  if (scene->rayIntersect(projRay)) {
    ++noVisShiftVol;
    // If there is an intersection, the current shift is impossible
    return false;
  }

  if (!parentPhoton->isMediumInteraction()) {
    const Float signDot = dot(parentPhoton->getGeometricNormal(), dProj) /
        dot(parentPhoton->getGeometricNormal(), source->edge(currVertex - 1)->d);
    if (signDot < 0.f) {
      // This mean a change of sign between the direction
      ++noVisShiftVol;
      return false;
    }
  }

  // Recompute the associated flux of the unchanged path part
  Spectrum photonWeight = source->vertex(0)->weight[EImportance] *
      source->vertex(0)->rrWeight *
      source->edge(0)->weight[EImportance];
  for (size_t i = 1; i < currVertex - 1; ++i) {
    photonWeight *= source->vertex(i)->weight[EImportance] *
        source->vertex(i)->rrWeight *
        source->edge(i)->weight[EImportance];
  }

  // Query the diffuse shift computation
  ShiftRecord sRec;
  {
    Intersection fakeInter;
    fakeInter.p = offsetPos;
    // SH: is this correct because there might still have vertex that is on surface?
    // In this case, give the parent parent location
    diffuseReconnection(sRec, fakeInter, dProj, lProj,
                        parentPhoton, parentEdge,
                        true,
                        currVertex >= 3 ? source->vertex(currVertex - 2)->getPosition() : Point(1.f));
    noPDFVol.incrementBase();
    if (sRec.pdf == Float(0)) {
      ++noPDFVol;
      result.weight = 1.0f;
      return false;
    }
  }

  // Jacobian for camera path contribution is already accounted in the camera vertex weight
  // Jacobian for light path contribution
  result.jacobian *= sRec.jacobian * additionalJacobian;
  photonWeight *= sRec.throughtput;

  Spectrum contrib = getVolumePhotonContrib(photonWeight, shiftMRec,
                                            -dProj, -shiftRay.d, EImportance);

  // We have to use the cache values here because the offset eye path is not explicitly constructed
  // currEdge = i means we are on the way from vertex i to vertex i + 1.
  Spectrum
      eyeShiftContrib = shiftGP.getWeightBeam(currEdge - 1) *    // cummulative vertex(i-1).weight * edge(i-1).weight
      shiftGP.getWeightVertex(currEdge);          // vertex(i).weight
  // FIXME: make sure BSDF value / area pdf is correct in shift_cameraPath.cpp for glossy material.

  result.shiftedFlux = shiftMRec.transmittance * contrib * eyeShiftContrib * result.jacobian;
  result.weight = 0.5f;

  // MIS for area measure is supported for now

  if (config.useMIS) {
    SAssert(currVertex > 1);

    Float basePdf = pdfBaseRay;
    basePdf *= source->vertex(currVertex - 1)->pdf[EImportance];
    basePdf *= source->edge(currVertex - 1)->pdf[EImportance];

    Float offsetPdf = sRec.pdf * pdfShiftRay;
    if (offsetPdf == Float(0) || basePdf == Float(0)) {
      //SLog(EWarn, "Invalid path");
      result.weight = 1.0f;
      return false;
    }

    const Float sensorPart = shiftGP.sensorMIS(currEdge, *baseGather,
                                               shiftRay.maxt, baseRay.maxt);
    if (config.powerHeuristic) {
      result.weight = 1.0f / (1.0f + powerOfTwo(sensorPart * result.jacobian * (offsetPdf / basePdf)));
    } else {
      result.weight = 1.0f / (1.0f + sensorPart * result.jacobian * (offsetPdf / basePdf));
    }
  }
  return true;
}

/////////////////////// VOLUME PHOTON
void VolumeGradientPositionQuery::operator()(const GPhotonNodeKD &nodePhoton) {
  if (baseGather->path.vertexCount() < 2)
    SLog(EError, "Problem path");

  GPhoton photon;
  nodePhoton.getData().getPhoton(photon);
  Point pos = baseRay(baseRay.maxt);

  // Test the radius
  Float lengthSqr = (pos - photon.its.p).lengthSquared();
  if ((searchRadius * searchRadius - lengthSqr) < 0)
    return;

  // Test the depth of the photon
  size_t pathLength = currEdge + photon.depth;
  if ((config.maxDepth > 0 && pathLength > config.maxDepth))
    return;

  // Test interaction modes
  const PathVertex *vPrev = nodePhoton.getData().lightPath->vertex(nodePhoton.getData().vertexId - 1);
  if (!computeVolumeContribution(config, vPrev))
    return;

  // Debug one shift by cancel it
  if (config.debugShift != EAllShift && config.debugShift != ENullShift) {
    int b;
    ELightShiftType currShift = getTypeShift(nodePhoton.getData().lightPath,
                                             nodePhoton.getData().vertexId, b);
    if (config.debugShift != currShift) {
      return; // Do not compute the photon contribution
    }
  }

  Spectrum photonContrib = getVolumePhotonContrib(photon.weight, *baseMRec,
                                                  photon.its.wi, -baseRay.d, EImportance);

  Spectrum eyeContrib =
      baseGather->getWeightBeam(currEdge - 1) *   // cummulative vertex(i-1).weight * edge(i-1).weight
          baseGather->getWeightVertex(currEdge);  // * vertex(i).weight

  Spectrum baseContrib = eyeContrib * baseMRec->transmittance * photonContrib;
  Float kernelVol = (4.0 / 3.0) * M_PI * std::pow(searchRadius, 3);
  mediumFlux += baseContrib / (kernelVol * pdfBaseRay());

  // Shift
  const Point2 basePixel = baseGather->path.vertex(1)->getSamplePosition();
  const std::array<Point2, 4> pixels = generateOffsetPos(basePixel);
  Vector2i filmSize = scene->getFilm()->getSize();

  for (int i = 0; i < shiftGPs.size(); ++i) {
    GradientSamplingResult result;

    // If the gather point are not traced
    // Generate them
    shiftGPs[i].generate(scene, thdata.pool, *baseGather,
                         pixels[i], false);

    // The distance shift for the camera
    Float additionalJacobian = 1.f;
    if (shiftGPs[i].validVolumeEdge(currEdge, baseMRec->medium) && !shiftMRecIntialized) {
      const GatherPoint &shiftGather = shiftGPs[i];
      // Generate the ray that represent the beam segment
      Vector shiftDir = -shiftGather.path.edge(currEdge)->d;
      Float shiftDistMax = shiftGather.path.edge(currEdge)->length;

      // Just copy the base camera distance
      if (shiftDistMax >= baseRay.maxt) {
        validShiftDist[i] = true;

        // constant distance
        shiftDistCamera[i] = baseRay.maxt;

        // Just evaluate the transmittance
        Ray shiftRay(shiftGather.path.vertex(currEdge)->getPosition(),
                     shiftDir, Epsilon, shiftDistMax, 0.f);
        shiftMRec[i].t = shiftDistCamera[i];
        photon.medium->eval(shiftRay, shiftMRec[i], EDistanceAlwaysValid);
      }
    }

    cancelShortDistance.incrementBase();
    if (validShiftDist[i]) {
#if HAVE_ADDITIONAL_STATS
      shiftStats[i].nbLightShifts += 1;
#endif
      if (shiftDistCamera[i] == -1)
        SLog(EError, "Distance not initialized");

      ShiftGatherPoint &shiftGather = shiftGPs[i];
      Ray shiftRay(shiftGather.path.vertex(currEdge)->getPosition(),
                   -shiftGather.path.edge(currEdge)->d,
                   Epsilon, shiftDistCamera[i], 0.f);

      bool alreadyShifted = false;
      // Take care about the simple shift here
      if (config.useShiftNull) {
        // Compute the projection (if we keep the same position on the shift path)

        Float distSqr = (nodePhoton.getPosition() - shiftRay(shiftRay.maxt)).lengthSquared();
        if (distSqr < searchRadius * searchRadius) {
          alreadyShifted = true;
          shiftNull(photon.weight, photon.its.wi,
                    shiftGather, shiftRay,
                    shiftMRec[i], result,
                    pdfBaseRay(), pdfShiftRay(i, photon.medium, shiftDistCamera[i]),
                    additionalJacobian);
#if HAVE_ADDITIONAL_STATS
          shiftStats[i].nullShifts += 1;
          shiftStats[i].nbSuccessLightShifts += 1;
#endif
        }
      }

      if (!alreadyShifted) {
        // Compute the offset shift
        Point offsetPos = getShiftPos(shiftRay, searchRadius, nodePhoton.getPosition());

        // Repair the path
        if (config.debugShift != ENullShift) {
          bool successLightShift = shiftPhoton(offsetPos,
                                               nodePhoton.getData().lightPath,
                                               nodePhoton.getData().vertexId,
                                               shiftGather,
                                               shiftRay,
                                               shiftMRec[i],
                                               result,
                                               searchRadius,
                                               pdfBaseRay(),
                                               pdfShiftRay(i, photon.medium, shiftDistCamera[i]),
                                               additionalJacobian);
#if HAVE_ADDITIONAL_STATS
          int b;
          ELightShiftType currShift = getTypeShift(nodePhoton.getData().lightPath,
                                                   nodePhoton.getData().vertexId, b);
          if (currShift == EDiffuseShift) {
              shiftStats[i].DiffuseShifts += 1;
          } else {
              shiftStats[i].MEShifts += 1;
          }
          if (successLightShift) {
              shiftStats[i].nbSuccessLightShifts += 1;
          }
#endif

        }
      }
    } else {
      ++cancelShortDistance;
      result.weight = 1.f;
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

    shiftedMediumFlux[i] += result.weight * result.shiftedFlux / (kernelVol * pdfBaseRay());
    weightedMediumFlux[i] += result.weight * baseContrib / (kernelVol * pdfBaseRay());
  }
  shiftMRecIntialized = true;
}

/////////////////////// BRE
void VolumeGradientBREQuery::operator()(const GPhotonNodeKD &nodePhoton, Float photonRadius, Float randValue) {
  SAssert(nodePhoton.getData().vertexId - 1 >= 0);
  const size_t currDepth = nodePhoton.getData().vertexId - 1;
  const Spectrum &currFlux = nodePhoton.getData().weight;
  const Vector &currWi = -nodePhoton.getData().lightPath->edge(currDepth)->d; // currVertex - 1;
  const Medium *currMed = nodePhoton.getData().lightPath->edge(currDepth)->medium;

  // Contrary to regular photon mapping,
  // Here we do not need to test the photon radius
  // Because the photon radius os kept same during the shift

  // Test the depth of the photon
  if (config.maxDepth > 0 && int(currDepth + currEdge) > config.maxDepth)
    return;
  if (config.minDepth != 0 && int(currDepth + currEdge) < config.minDepth)
    return;

  const PathVertex *vPrev = nodePhoton.getData().lightPath->vertex(nodePhoton.getData().vertexId - 1);
  if (!computeVolumeContribution(config, vPrev))
    return;

  // Debug one shift by cancel it
  if (config.debugShift != EAllShift && config.debugShift != ENullShift) {
    int b;
    ELightShiftType currShift = getTypeShift(nodePhoton.getData().lightPath,
                                             nodePhoton.getData().vertexId, b);
    if (config.debugShift != currShift) {
      return; // Do not compute the photon contribution
    }
  }

  Float rrGlobalWeight = 1;
  if (config.pathSet) {
    const Point2 posPix = baseGather->path.vertex(1)->getSamplePosition();
    size_t currentGroup = (int(posPix.x) + int(posPix.y)) % 2;
    if (nodePhoton.getData().pathID % 2 != currentGroup) {
      return; // Reject it
    }
    rrGlobalWeight = 2;
  }

  // We intersect disks aligned perpendicular to the ray direction
  // so the density kernel is 2D
  Float kernelVol = M_PI * std::pow(photonRadius, 2);
  Float pdfCameraPos = 1.f;

  // Keep the photon projection for the next computations
  bool validBaseDistance = true;
  const Float baseProjDist = baseRay.maxt;
  if (EVolumeTechniqueHelper::use3DKernel(config.volTechnique)) {
    kernelVol = ((4.0 / 3.0) * M_PI * std::pow(photonRadius, 3));

    Float distSqr = (baseRay(baseRay.maxt) - nodePhoton.getPosition()).lengthSquared();

    // Determine delta T
    Float deltaT = math::safe_sqrt(photonRadius * photonRadius - distSqr);

    // Compute the random new distance
    Float tminKernel = baseRay.maxt - deltaT;
    Float diskDistanceRand = tminKernel + (deltaT * 2) * randValue;

    if (diskDistanceRand < baseRay.mint || diskDistanceRand > baseGather->path.edge(currEdge)->length) {
      validBaseDistance = false;
    }

    // Replace the distance by the sampled one and compute the PDF
    baseRay.maxt = diskDistanceRand;
    pdfCameraPos = 1.f / std::max(deltaT * 2.0, 0.0001);
  } else {
    if (baseProjDist > baseGather->path.edge(currEdge)->length) {
      // Not possible

    }
  }

  // Compute the photon contribution.
  Spectrum baseContrib(0.f);
  if (validBaseDistance) {
    // Evaluate transmittance
    MediumSamplingRecord mRecBase;
    medium->eval(baseRay, mRecBase);

    // Compute the base contrib
    Spectrum contrib = getVolumePhotonContrib(currFlux, mRecBase, currWi, -baseRay.d, EImportance);
    Spectrum eyeContrib =
        baseGather->getWeightBeam(currEdge - 1) *   // cummulative vertex(i-1).weight * edge(i-1).weight
            baseGather->getWeightVertex(currEdge);  // * vertex(i).weight
    baseContrib = mRecBase.transmittance * contrib * eyeContrib;

    // Only do here
    mediumFlux += (baseContrib / (kernelVol * pdfCameraPos)) * rrGlobalWeight;
  } else {
    return; // Impossible to gather, we just skip the path
  }

  const Point2 &basePixel = baseGather->path.vertex(1)->getSamplePosition();
  const std::array<Point2, 4> pixels = generateOffsetPos(basePixel);
  const Vector2i &filmSize = scene->getFilm()->getSize();

  for (int i = 0; i < shiftGPs.size(); ++i) {
    GradientSamplingResult result;
    // If the gather point are not traced
    // Generate them
    shiftGPs[i].generate(scene, thdata.pool, *baseGather,
                         pixels[i], false);

    if (shiftGPs[i].validVolumeEdge(currEdge, currMed)) {
      const ShiftGatherPoint &shiftGather = shiftGPs[i];
      // Generate the ray that represent the beam segement
      const Vector &shiftDir = -shiftGather.path.edge(currEdge)->d;
      const Float shiftDistTotal = shiftGather.path.edge(currEdge)->length;
      Ray shiftRay(shiftGather.path.vertex(currEdge)->getPosition(),
                   shiftDir, Epsilon, baseRay.maxt, 0.f); // Copy the new location on the beam ray

      bool alreadyShift = false;
      const Path *source = nodePhoton.getData().lightPath;
      const size_t currVertex = nodePhoton.getData().vertexId;
      // Take care about the simple shift here
      if (config.useShiftNull) {
        if (!EVolumeTechniqueHelper::use3DKernel(config.volTechnique)) {
          SLog(EError, "Not possible :)");
        }

        // || z'_* - y_* ||
        const Float ZPtoY = (shiftRay(shiftRay.maxt) - nodePhoton.getPosition()).lengthSquared();
        if (ZPtoY < photonRadius * photonRadius &&
            shiftRay.maxt < shiftDistTotal) {
          Vector originToCenter = nodePhoton.getPosition() - shiftRay.o;
          Float diskDistance = dot(originToCenter, shiftRay.d);

          // Compute the extra probability for the shifted
          // This probability can be different due to some occlu
          const Float distSqr = (shiftRay(diskDistance) - nodePhoton.getPosition()).lengthSquared();
          const Float deltaT = math::safe_sqrt(photonRadius * photonRadius - distSqr);
          const Float pdfShiftPos = 1.f / std::max(2.0 * deltaT, 0.0001);

          // The photon is inside, pick it
          MediumSamplingRecord mRecShift;
          currMed->eval(shiftRay, mRecShift);
          shiftNull(currFlux, currWi, shiftGather, shiftRay, mRecShift,
                    result, pdfCameraPos, pdfShiftPos);

          // Make already shifted to avoid multiple shifts
          alreadyShift = true;
        }

      }


      // Check if the distance is ok for the shifted gatherpoint
      // if not, mark the gather as invalid
      if (!alreadyShift) {
        if (shiftDistTotal >= shiftRay.maxt) {
          Point offsetPos = getShiftPos(shiftRay,
                                        photonRadius,
                                        nodePhoton.getPosition(),
                                        !EVolumeTechniqueHelper::use3DKernel(config.volTechnique));

          Float pdfShiftPos = 1.f;
          if (EVolumeTechniqueHelper::use3DKernel(config.volTechnique)) {
            // Compute the possible PDF for the offset photon
            Vector originToCenter = offsetPos - shiftRay.o;
            Float diskDistance = dot(originToCenter, shiftRay.d);
            Float distSqr = (shiftRay(diskDistance) - offsetPos).lengthSquared();

            // Determine delta T for the shift path
            Float deltaT = math::safe_sqrt(photonRadius * photonRadius - distSqr);
            pdfShiftPos = 1.f / std::max(2.0 * deltaT, 0.0001); // We can avoid this computation
          }

          // Repair the photon path
          if (config.debugShift != ENullShift) {
            MediumSamplingRecord mRecShift;
            currMed->eval(shiftRay, mRecShift);
            shiftPhoton(offsetPos, nodePhoton.getData().lightPath,
                        nodePhoton.getData().vertexId,
                        shiftGather, shiftRay,
                        mRecShift, result, photonRadius,
                        pdfCameraPos, pdfShiftPos);
          }
        }
      }
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

    weightedMediumFlux[i] += rrGlobalWeight * result.weight * baseContrib / (kernelVol * pdfCameraPos);
    shiftedMediumFlux[i] += rrGlobalWeight * result.weight * result.shiftedFlux / (kernelVol * pdfCameraPos);
  }
}

Point AbstractVolumeGradientRecord::getShiftPos(const Ray &shiftRay, Float radius, Point basePhotonPos, bool coherent) {

  // Just preserve distance
  Float baseW = baseRay.maxt;
  Float shiftW = shiftRay.maxt;
  Point offsetPos = shiftRay(shiftW) + (basePhotonPos - baseRay(baseW));

  // Compute the position of the candidate
  // Use the original projection on the base ray
  if (coherent) {
    Frame bL(baseRay.d);
    coordinateSystemCoherent(bL.n, bL.s, bL.t);
    Frame nL(shiftRay.d);
    coordinateSystemCoherent(nL.n, nL.s, nL.t);
    const Vector localD = bL.toLocal(basePhotonPos - baseRay(baseW));
    offsetPos = shiftRay(shiftW) + nL.toWorld(localD);
  }

  if (config.useShiftNull) {
    additionalShiftVol.incrementBase();
    Float offDistSqr = (baseRay(baseW) - offsetPos).lengthSquared();
    if (offDistSqr < radius * radius) {
      // The current photon is inside the 3D kernel
      // Make the additional shift
      ++additionalShiftVol;
#if 0
      // Naive additional shift
      offsetPos = offsetPos + (shiftRay(shiftRay.maxt) - offsetPos)*2.f;
#else
      // More coherent additional shift
      Vector dShift = shiftRay(shiftRay.maxt) - baseRay(baseRay.maxt);
      dShift /= dShift.length();
      const Float cosD = dot(dShift, -(offsetPos - shiftRay(shiftRay.maxt)));
      offsetPos += dShift * cosD * 2;
#endif
    }
  }
  return offsetPos;
}

MTS_NAMESPACE_END
