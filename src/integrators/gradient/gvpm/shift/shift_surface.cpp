//
// Created by beltegeuse on 12/1/16.
//

#include "shift_surface.h"

// For the statistics
#include <mitsuba/core/statistics.h>
// For ME shift
#include <mitsuba/bidir/manifold.h>
// For the diffuse shift
#include "operation/shift_diffuse.h"
// For the ME shift
#include "operation/shift_ME.h"

MTS_NAMESPACE_BEGIN

StatsCounter intersectionSucessRatio("GPM",
                                     "Interesected surface projection", EPercentage);
StatsCounter contributiveShiftedPRatio("GPM",
                                       "Have contribution photon", EPercentage);

StatsCounter MEMisSuccess("GPM",
                          "MIS computation succesfull", EPercentage);

StatsCounter MEShift("GPM",
                     "Percentage ME Shift    : ", EPercentage);
StatsCounter HVShift("GPM",
                     "Percentage HVCopy Shift: ", EPercentage);
StatsCounter DiffShift("GPM",
                       "Percentage Diff. Shift : ", EPercentage);
StatsCounter InvalidShift("GPM",
                          "Percentage Inval Shift : ", EPercentage);

bool SurfaceGradientRecord::shiftPhoton(const Path *source, int currVertex, const GatherPoint *shiftGather,
                                        GradientSamplingResult &result) {
  if (shiftGather == nullptr || shiftGather->depth != baseGather->depth) {
    // Invalid gather point
    return false;
  }

  // Maintain the photon distribution at the shift gather point to be the same as the base
  const PathVertex *basePhoton = source->vertex(currVertex);

  Vector localOffset = baseGather->lastIts().toLocalCoherent(basePhoton->getPosition() - baseGather->lastIts().p);
  Point newVLoc = shiftGather->lastIts().p + shiftGather->lastIts().toWorldCoherent(localOffset);

  // Just test the visibility aka. projection from previous point
  // So compute the new direction of photon emission from the previous vertex
  const PathVertex *prevVertex = source->vertex(currVertex - 1);
  Vector dProj = normalize(newVLoc - prevVertex->getPosition());

  // Project the photon to the shifted photon
  Intersection itsProj;
  Ray projRay(prevVertex->getPosition(), dProj, 0.f);      // FIXME: zero time
  if (!intersectNULL(scene, projRay, itsProj)) {
    return false;
  }

  // Minor: Add a term to account for tracing event from parent photon
  // Inspired from Fig. 8 in HSLT paper
  // Due to tracing, dx'_k / dx_k is not exactly 1.
  int c = currVertex;
  Vector imgNormal = shiftGather->lastIts().geoFrame.n;
  result.jacobian = geometryOpposingTerm(source->vertex(c - 1)->getPosition(), newVLoc, imgNormal)
      / geometryOpposingTerm(source->vertex(c - 1)->getPosition(), itsProj.p, itsProj.geoFrame.n);
  if (result.jacobian == 0.0f) {
    // Sometimes the direction to the candidate offset photon is perpendicular to imgNormal,
    // causing jacobian to be zero.
    // We ignore such cases by resetting jacobian to 1,
    // as eventually we only need the on-surface offset photon
    result.jacobian = 1.0f;
  }

  // Search next diffuse vertex
  int b;
  ELightShiftType type = getTypeShift(source, currVertex, b);
  RayDifferential projRayDiff(projRay);

  MEShift.incrementBase();
  HVShift.incrementBase();
  DiffShift.incrementBase();
  InvalidShift.incrementBase();

  if (type == EInvalidShift) {
    // Error?
    ++InvalidShift;
    return false;
  } else if (type == EDiffuseShift) {
    // Just use diffuse reconnection here
    ++DiffShift;
    return shiftPhotonDiffuse(result, source, currVertex, itsProj, shiftGather);
  } else if (type == EManifoldShift) {
    if (!config.useManifold) {
      ++InvalidShift;
      return false; // Impossible to handle this path
    } else {
      ++MEShift;
      return shiftPhotonManifold(result, source, currVertex, b, itsProj, shiftGather);
    }
  }
  SLog(EError, "Invalid type");
  ++InvalidShift;
  return false;
}

bool SurfaceGradientRecord::shiftPhotonManifold(GradientSamplingResult &result,
                                                const Path *lt, int c, int b,
                                                const Intersection &itsProj,
                                                const GatherPoint *shiftGather) {
  // Perform a walk through a specular chain to a non-specular vertex
  const Path &source = (Path &) *lt;
  Path proposal;
  {
    // Construct the new vertex
    PathVertex shiftVertex;
    memset(&shiftVertex, 0, sizeof(PathVertex));
    shiftVertex.getIntersection() = itsProj;
    shiftVertex.type = PathVertex::ESurfaceInteraction;
    shiftVertex.measure = EArea;
    shiftVertex.sampledComponentIndex = source.vertex(c)->sampledComponentIndex;

    if (!generateShiftPathME(source, proposal, b, c, thdata.pool, thdata.offsetGenerator.get(),
                             shiftVertex, 0.f, //photonRadius * config.relaxME
                             baseGather->lastIts().p,
                             shiftGather->lastIts().p)) {
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
  result.jacobian *= sRecME.jacobian;  // middle term | dy_s' / dy_s |

  // Jacobian in area integral in original space
  SpecularManifold *manifold = thdata.offsetGenerator->getSpecularManifold();

  // Jacobian computed using det
  const Float detProposed = manifold->det(proposal, b, c);
  const Float detSource = manifold->det(source, b, c);
  result.jacobian *= detProposed / detSource;

  if (result.jacobian  <= 0.0 || !std::isfinite(result.jacobian )) {
    SLog(EWarn, "Invalid jacobian %g %d %d", result.jacobian , b, c);
	for (int i = b; i <= c; ++i) {
	  thdata.pool.release(proposal.edge(i - 1));
	  thdata.pool.release(proposal.vertex(i));
	}
    result.weight = 1.0f;
    return false;
  }

  Spectrum photonWeight = sRecME.throughtput;

  // Perform photon mapping
  GPhoton offsetPhoton(itsProj, nullptr, -1, photonWeight, c - 1, GPhoton::getPrevComponentType(source.vertex(c - 1)));
  offsetPhoton.its.wi =
      offsetPhoton.its.toLocal(normalize(proposal.vertex(c - 1)->getPosition() - proposal.vertex(c)->getPosition()));

  std::tuple<bool, Spectrum>
      photonContribution = getSurfacePhotonContrib(offsetPhoton, shiftGather, baseGather->radius);
  if (!std::get<0>(photonContribution)) {
    // The projected photon doesn't have any contribution:
    // 1) Not inside the radius of the shifted gather point.
    // 2) Share not the same normal as the gather points.
    // And hence, the shift is not reversible.
	for (int i = b; i <= c; ++i) {
	  thdata.pool.release(proposal.edge(i - 1));
	  thdata.pool.release(proposal.vertex(i));
	}
    result.weight = 1.0f;
    return false;
  }

  // Default no MIS weight
  result.weight = 0.5f;
  result.shiftedFlux = std::get<1>(photonContribution) * shiftGather->surfInfo().weight * result.jacobian;

  // Pdf check for offset path
  // So we do not need to recompute it for the MIS computation
  // and it is a sanity check
  Float offsetPdf = shiftGather->surfInfo().pdf * sRecME.pdf;
  if (offsetPdf == Float(0)) {
    //SLog(EWarn, "Invalid offset path");
    // This path has zero density due to manifold walk.
    // It cannot be sampled by particle tracing, and hence no reversanywayible.
    result.weight = 1.0f;
    if (!result.shiftedFlux.isZero()) {
      SLog(EWarn, "0 PDF path but with a throughput: %s\n Set to 0", result.shiftedFlux.toString().c_str());
      result.shiftedFlux = Spectrum(0.f);
    }
    // When this assert fails, it could be due to pdfComponent is not correctly implemented.
  }

  if (config.useMIS) {
    MEMisSuccess.incrementBase();

    Float basePdf = baseGather->surfInfo().pdf;
    // No need to consider pdf of the unchanged segments
    for (int i = b; i < c; ++i) {
      basePdf *= source.vertex(i)->pdf[EImportance];
    }

    Float allJacobian = shiftGather->surfInfo().jacobian * result.jacobian;

    if (basePdf == Float(0)) {
      SLog(EWarn, "Invalid base path. This case should not happen.");
      result.weight = 0.0f;
    } else {
      result.weight = 1.0f / (1.0f + offsetPdf * allJacobian / basePdf);
    }

  }

  // Clean up
  for (int i = b; i <= c; ++i) {
	thdata.pool.release(proposal.edge(i - 1));
	thdata.pool.release(proposal.vertex(i));
  }

  return true;
}

bool SurfaceGradientRecord::shiftPhotonDiffuse(GradientSamplingResult &result,
                                               const Path *lt, int currVertex,
                                               const Intersection &itsProj,
                                               const GatherPoint *shiftGather) {

  // Maintain the photon distribution at the shift gather point to be the same as the base
  const PathVertex *parentVertex = lt->vertex(currVertex - 1);
  const PathEdge *edge = lt->edge(currVertex - 1); // Get the edge between the two surface points

  Vector newD = itsProj.p - parentVertex->getPosition();
  Float newDLength = newD.length();
  newD /= newDLength;

  // Check visibility of the current photon shift
  Ray projRay(parentVertex->getPosition(), newD, Epsilon, newDLength * ShadowEpsilon, 0.f);
  if (scene->rayIntersect(projRay)) {
    // If there is an intersection, the current shift is impossible
    return false;
  }

  ShiftRecord sRec;
  if (currVertex >= 3) {
    diffuseReconnection(sRec, itsProj, newD, newDLength,
                        parentVertex, edge, false, lt->vertex(currVertex - 2)->getPosition());
  } else {
    diffuseReconnection(sRec, itsProj, newD, newDLength,
                        parentVertex, edge, false, Point(1.f));
  }
  if (sRec.pdf == Float(0)) {
    result.weight = 1.0f;
    return false;
  }

  // To be strict, we must re-evaluate the BSDF value at the parent
  // Compute the flux of parent photon
  Spectrum photonWeight = lt->vertex(0)->weight[EImportance] *
      lt->vertex(0)->rrWeight *
      lt->edge(0)->weight[EImportance];

  for (int i = 1; i < currVertex - 1; ++i) {
    photonWeight *= lt->vertex(i)->weight[EImportance] *
        lt->vertex(i)->rrWeight *
        lt->edge(i)->weight[EImportance];
  }

  // Jacobian for solid angle integral does not need to account for photon projection
  photonWeight *= sRec.throughtput * sRec.jacobian;

  // Assign the first intersection, remember to assign the incoming direction for the offset photon
  GPhoton offsetPhoton(itsProj, nullptr, -1, photonWeight, currVertex - 1,
                       GPhoton::getPrevComponentType(lt->vertex(currVertex - 1)));
  offsetPhoton.its.wi = offsetPhoton.its.toLocal(normalize(parentVertex->getPosition() - itsProj.p));

  // Evaluate BSDF at the photon
  std::tuple<bool, Spectrum>
      photonContribution = getSurfacePhotonContrib(offsetPhoton, shiftGather, baseGather->radius);

  contributiveShiftedPRatio.incrementBase(1);
  if (!std::get<0>(photonContribution)) {
    // The projected photon doesn't any contribution:
    // 1) Not inside the radius of the shifted gp.
    // 2) Share not the same normal as the gather points
    result.weight = 1.0f;
    return false;
  }
  ++contributiveShiftedPRatio;

  // Finish
  result.shiftedFlux = std::get<1>(photonContribution) * shiftGather->surfInfo().weight;
  result.weight = 0.5f;

  if (config.useMIS) {
    Float basePdf = baseGather->surfInfo().pdf;
    basePdf *= lt->vertex(currVertex - 1)->pdf[EImportance];
    basePdf *= lt->edge(currVertex - 1)->pdf[EImportance];

    Float offsetPdf = shiftGather->surfInfo().pdf * sRec.pdf;
    if (offsetPdf == Float(0) || basePdf == Float(0)) {
      //SLog(EWarn, "Invalid path");
      result.weight = 1.0f;
      return false;
    }

    Float allJacobian = shiftGather->surfInfo().jacobian * sRec.jacobian;

    result.weight = 1.0f / (1.0f + offsetPdf * allJacobian / basePdf);
  }
  return true;
}

MTS_NAMESPACE_END
