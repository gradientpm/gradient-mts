// For the statistics
#include <mitsuba/core/statistics.h>

#include "shift_volume_beams.h"

// For diffuse shift

// For ME shift

MTS_NAMESPACE_BEGIN

StatsCounter MEShiftVolBeam("GVPM",
                            "Percentage ME Vol Shift    : ", EPercentage);
StatsCounter HVShiftVolBeam("GVPM",
                            "Percentage HVCopy Vol Shift: ", EPercentage);
StatsCounter DiffShiftVolBeam("GVPM",
                              "Percentage Diff. Vol Shift : ", EPercentage);
StatsCounter InvalidShiftVolBeam("GVPM",
                                 "Percentage Inval Vol Shift : ", EPercentage);
StatsCounter NullShiftVolBeam("GVPM",
                              "Percentage Null Vol Shift : ", EPercentage);
StatsCounter additionalShiftVolBeam("GVPM",
                                    "Percentage of additional shift: ", EPercentage);

StatsCounter MEShiftVolBeamShiftFailed("GVPM",
                                       "Percentage of shift failed in ME: ", EPercentage);
StatsCounter MEShiftVolBeamContribFailed("GVPM",
                                         "Percentage of failed contribution in ME: ", EPercentage);
StatsCounter MEShiftVolBeamContribZero("GVPM",
                                       "Percentage of zero contribution in ME: ", EPercentage);
StatsCounter SingleScatteringShiftSkip("GVPM",
                                       "Percentage of skipped single scattering path: ", EPercentage);

/*
 * Build a basis with r.d and a projection vector as basis
 */
Frame localMatrix(const Ray &r, const Point &a) {
  Float d = dot(a - r.o, r.d);
  Vector s = normalize(a - r(d));
  Vector t = cross(r.d, s);
  return {r.d, s, t};
}

/**
 * Generate a new point such that it stays on a special plane the contains vertex a.
 * The distance between the camera ray to this plane is still equal to u.
 * @param r
 * @param a
 * @param u
 * @param w
 * @param flip Only two points can be generated. This parameter selects one of them.
 * @return
 */
Point shift(const Ray &r, const Point &a, Float u, Float w, bool flip = false) {
  // Construct the new basis
  Frame lT = localMatrix(r, a);

  // Compute the translation for that a is inside the plane
  // defined by the basis lT
  Float d = dot(a - r.o, r.d);
  Point tD = r(d);

  // Express a in the locally coordinates
  Vector localA = lT.toLocal(a - tD);

  // In this local space, founding the intersection point
  // is straightforward
#ifdef WIN32
  const Float M_PI_2 = 1.57079632679489661923f;
#endif
  Float phi = M_PI_2 - math::safe_asin(u / fabs(localA.y));
  if (flip)
    phi = -phi;
  Vector localW = Vector(0, u * cos(phi), u * sin(phi));

  // Re-express the new direction in world space
  Point worldU = Point(lT.toWorld(localW));
  return r(w) + worldU;
}

Point BeamGradRadianceQuery::getShiftPos1D(const Ray &bRay,
                                           const Ray &sRay,
                                           const Point &a,
                                           const Vector &bBeamDir,
                                           Float w, Float u) {
  // First compute the shift on the base path
  // the idea is to found if phi angle need to be positive
  // or negative
  Vector baseShitedBack = normalize(shift(bRay, a, u, w) - a);
  bool flipAngle = false;
  if ((baseShitedBack - bBeamDir).lengthSquared() > 0.001) {
    flipAngle = true;
  }

  return shift(sRay, a, u, w, flipAngle);
}

Point BeamGradRadianceQuery::getShiftPos(const Ray &bRay,
                                         const Ray &sRay,
                                         Float w, const Vector &u,
                                         Float radius,
                                         Float newW, bool coherent) {
  Point newPos = sRay(newW) + u;

  // Force vector u on base and offset path to be in the same direction
  // Vector u goes from the eye segment to light beam segment
  // In beam 3D, we define u direction from the eye segment in the kernel to the kernel center
  if (coherent) {
    Frame bL(bRay.d);
    coordinateSystemCoherent(bL.n, bL.s, bL.t);
    Frame nL(sRay.d);
    coordinateSystemCoherent(nL.n, nL.s, nL.t);
    newPos = sRay(newW) + nL.toWorld(bL.toLocal(u));
  }


  // Change the shift pos if we use shift
  if (config.useShiftNull) {
    const Point bCamW = bRay(w);
    additionalShiftVolBeam.incrementBase();

    // If shifted photon can cover the base camera ray, we still use null shift
    Float offDistSqr = (bCamW - newPos).lengthSquared();
    if (offDistSqr < radius * radius) {
      // The current photon is inside the 3D kernel
      Vector dShift = sRay(newW) - bCamW;
      dShift /= dShift.length();
      const Float cosD = dot(dShift, -(newPos - sRay(newW)));
      newPos += dShift * cosD * 2;

      ++additionalShiftVolBeam;
    }
  }
  // Note that in shift mapping, only w on the camera ray is preserved, and then vector u.
  // Value v on the shifted photon ray and the base photon ray could be different.
  return newPos;
}

bool BeamGradRadianceQuery::operator()(const LTPhotonBeam *beam, Float tmin, Float tmax) {
  Float rrGlobalWeight = 1;

  // Test the depth of the beam
  int pathLength = currCameraEdge + beam->depth;
  if ((config.maxDepth > 0 && pathLength > config.maxDepth))
    return false;

  // Discard if the contribution if needed
  const PathVertex *vPrev = beam->path->vertex(beam->edgeID);
  if (!computeVolumeContribution(config, vPrev))
    return false;

  // Before defining the shift type
  // Let check if we need to shift the path
  // If it was a L(S*) light path
  // This is a simple approximation only on the vertex type (pure specular)
  // Improvement will be needed.
  bool excludeLater = false;
  if (false) {
    const Path *source = beam->path;
    const int currVertex = beam->edgeID + 1;
    SingleScatteringShiftSkip.incrementBase();
    if (currVertex == 2 || currVertex == 1) {
      // If we have direct photon from the light source,
      // we still allow to do the shift
    } else {
      bool specularOnly = true;
      for (int i = currVertex - 1; i > 1 && specularOnly; i--) {
        EVertexType vertexType = VertexClassifier::type(*source->vertex(i), source->vertex(i)->sampledComponentIndex);
        specularOnly = vertexType != VERTEX_TYPE_DIFFUSE;
      }
      // We have found only specular vertex, so we mark this shift as failed.
      if (specularOnly) {
        ++SingleScatteringShiftSkip;
        excludeLater = true;
      }
    }
  }

  // If we exclude later, we do not split the path
  if (config.pathSet && (!excludeLater)) {
    const Point2 posPix = baseGather->path.vertex(1)->getSamplePosition();
    size_t currentGroup = (int(posPix.x) + int(posPix.y)) % 2;
    if (beam->pathID % 2 != currentGroup) {
      return false; // Reject it
    }
    rrGlobalWeight = 2;
  }

  const Float radius = beam->getRadius();
  BeamKernelRecord kRec = BeamKernelRecord(config.volTechnique,
                                           medium,
                                           sampler,
                                           beam, radius,
                                           baseCameraRay,
                                           tmin, tmax);
  if (!kRec.isValid()) {
    return false;
  }

  // cummulative vertex(i-1).weight * edge(i-1).weight
  // * vertex(i).weight = bsdf / p(w)
  Spectrum eyeContrib = baseGather->getWeightBeam(currCameraEdge - 1) *
          baseGather->getWeightVertex(currCameraEdge);

  // Phase function evaluation at the intersection is already in beam->getContrib()

  Spectrum baseContrib = eyeContrib * kRec.contrib * kRec.weightKernel;
  if (excludeLater) {
    extraFlux += baseContrib * rrGlobalWeight;
    return true; // Stop here
  } else {
    mediumFlux += baseContrib * rrGlobalWeight;
  }

  const Point2 basePixel = baseGather->path.vertex(1)->getSamplePosition();
  const std::array<Point2, 4> pixels = generateOffsetPos(basePixel);
  Vector2i filmSize = scene->getFilm()->getSize();

  int bID;
  ELightShiftType currShift;
  // Debug and count the number of ME in this area
  {
    currShift = getTypeShift(beam->path, beam->edgeID + 1, bID);
    if (config.debugShift != EAllShift && config.debugShift != ENullShift && config.debugShift != currShift) {
      return false; // Do not compute the photon contribution
    }
  }

  clearCache();

  for (int i = 0; i < 4; ++i) {
#if HAVE_ADDITIONAL_STATS
    shiftStats[i].nbLightShifts += 1;
#endif
    shiftGPs[i].generate(scene, thdata.pool, *baseGather,
                         pixels[i], false);

    GradientSamplingResult result;
    if (shiftGPs[i].validVolumeEdge(currCameraEdge, beam->medium)) {
      const Path &shiftPath = shiftGPs[i].path;
      Float shiftDistMAX = shiftPath.edge(currCameraEdge)->length;
      Ray shiftRay(shiftPath.vertex(currCameraEdge)->getPosition(),
                   -shiftPath.edge(currCameraEdge)->d,
                   Epsilon, shiftDistMAX, 0.f);

      // Keep the same camera distance
      Float shiftW = kRec.w;

      // Keep track of the shifts statistics
      NullShiftVolBeam.incrementBase();
      MEShiftVolBeam.incrementBase();
      HVShiftVolBeam.incrementBase();
      DiffShiftVolBeam.incrementBase();
      InvalidShiftVolBeam.incrementBase();

      bool alreadyShift = false;

      // Take care about the simple shift here
      if (config.useShiftNull) {
        // Old way to do the null shift:
        // We make the null shift only if the shift path intersect
        // the current kernel position over the base beam.
        Point kernelPos = beam->getPos(kRec.v);
        const Float ZPtoY = (shiftRay(shiftW) - kernelPos).lengthSquared();

        if (config.volTechnique == EBeamBeam1D) {
          // Ignored in case of Beam 1D kernel
        } else if (config.volTechnique == EBeamBeam3D_Naive ||
            config.volTechnique == EBeamBeam3D_EGSR ||
            config.volTechnique == EBeamBeam3D_Optimized) {
          if (ZPtoY < radius * radius &&
              kRec.w <= shiftDistMAX) {
            BeamKernelRecord kRecShift = BeamKernelRecord(kRec, medium,
                                                          beam,
                                                          shiftRay);
            if (kRecShift.isValid()) {
              ++NullShiftVolBeam;
              // Could be true or false depending on beam contrib evaluation
              shiftNull3D(beam, shiftGPs[i], shiftRay, shiftW, result, kRec, kRecShift);
              alreadyShift = true;
#if HAVE_ADDITIONAL_STATS
              shiftStats[i].nullShifts += 1;
              shiftStats[i].nbSuccessLightShifts += 1;
#endif
            }
          }

        } else {
          SLog(EError, "Unsupported shifting volume technique %d", config.volTechnique);
        }
      }

      // Check if the distance is ok for the shifted gatherpoint
      // if not, mark the gather as invalid
      if (!alreadyShift && kRec.w <= shiftDistMAX) {
        if (!config.newShiftBeam) {
          Float minDistSqr = (beam->getOri() - shiftRay(dot(beam->getOri() - shiftRay.o, shiftRay.d))).lengthSquared();
          if (minDistSqr > kRec.u * kRec.u) {
            Point offsetPos = getShiftPos(baseCameraRay, shiftRay,
                                          kRec.w, beam->getPos(kRec.v) - baseCameraRay(kRec.w),
                                          radius, shiftW);
            shiftBeam(offsetPos, beam, shiftGPs[i], shiftRay, shiftW, result, kRec, bID, currShift);

          } else {
            result.weight = 1.f;
          }
        } else {
          Point offsetPos = getShiftPos1D(baseCameraRay,
                                          shiftRay, beam->getOri(),
                                          beam->getDir(), kRec.w, kRec.u);
          shiftBeam(offsetPos, beam, shiftGPs[i], shiftRay, shiftW, result, kRec, bID, currShift);
        }

        // Update the stats
#if HAVE_ADDITIONAL_STATS
        int b;
        ELightShiftType currShift = getTypeShift(beam->path,
                                                 beam->edgeID + 1, b);
        if (currShift == EDiffuseShift) {
          shiftStats[i].DiffuseShifts += 1;
        } else {
          shiftStats[i].MEShifts += 1;
        }
        if (succesfullShift) {
          shiftStats[i].nbSuccessLightShifts += 1;
        }
#endif
      }
    } else {
      result.weight = 1.f;
    }

    // Finally, include gather kernel
    // For beam1D, sinTheta is already included inside "shiftedFlux"
    result.shiftedFlux *= kRec.weightKernel;

    if ((i == ERight && (int) basePixel.x == filmSize.x - 1) ||
        (i == ETop && (int) basePixel.y == filmSize.y - 1)) {
      result.weight = 1.0f;
    }

    shiftedMediumFlux[i] += result.weight * result.shiftedFlux * rrGlobalWeight;
    weightedMediumFlux[i] += result.weight * baseContrib * rrGlobalWeight;
  }

  // Release memory if any
  clearCache();

  return true;
}

bool BeamGradRadianceQuery::shiftBeam(const Point &newPos,
                                      const LTPhotonBeam *beam, const ShiftGatherPoint &shiftGP,
                                      const Ray &shiftRay, Float shiftW, GradientSamplingResult &result,
                                      const BeamKernelRecord &kRec, int bID, ELightShiftType currShift) {
  if (config.debugShift == ENullShift) {
    result.weight = 1.0f;
    return false;
  }

  if (shiftW > shiftRay.maxt) {
    result.weight = 1.0f;
    return false;
  }

  const int currVertex = beam->edgeID + 1;
  bool shiftedEnoughRough = true; // FIXME: We consider that we are always enough rough inside the PM

  if (currShift == EInvalidShift) {
    // Error?
    ++InvalidShiftVolBeam;
    return false;
  } else if (currShift == EDiffuseShift && shiftedEnoughRough) {
    // Just use diffuse reconnection here
    ++DiffShiftVolBeam;
    return shiftBeamDiffuse(beam, shiftGP, shiftRay, shiftW, result, kRec, newPos);
  } else if (currShift == EMediumShift) {
    if (config.noMediumShift) {
      ++DiffShiftVolBeam;
      return shiftBeamDiffuse(beam, shiftGP, shiftRay, shiftW, result, kRec, newPos);
    } else {
      SLog(EError, "Not implemented");
    }
  } else if (currShift == EDiffuseShift && !shiftedEnoughRough) {
    if (!config.useManifold) {
      ++InvalidShiftVolBeam;
      return false; // Impossible to handle this path
    } else {
      ++MEShiftVolBeam;
      return shiftBeamME(bID, currVertex, beam, shiftGP, shiftRay, shiftW, result, kRec, newPos);
    }
  } else if (currShift == EManifoldShift) {
    if (!config.useManifold) {
      ++InvalidShiftVolBeam;
      return false; // Impossible to handle this path
    } else {
      ++MEShiftVolBeam;
      return shiftBeamME(bID, currVertex, beam, shiftGP, shiftRay, shiftW, result, kRec, newPos);
    }
  } else {
    SLog(EError, "Invalid shift type");
  }
  ++InvalidShiftVolBeam;
  return false;
}

bool BeamGradRadianceQuery::shiftBeamDiffuse(const LTPhotonBeam *beam, const ShiftGatherPoint &shiftGP,
                                             const Ray &shiftRay, Float shiftW, GradientSamplingResult &result,
                                             const BeamKernelRecord &kRec,
                                             const Point &newPos) {
  // Compute the new beam direction that passing throught the shifted position
  Vector newPBDir = (newPos - beam->getOri());
  Float newPBDist = newPBDir.length();
  newPBDir /= newPBDist;

  // Visibility test
  Ray newRayPB(beam->getOri(), newPBDir, Epsilon, newPBDist, 0.f);
  if (scene->rayIntersect(newRayPB)) {
    // The beam is not visible, do not continue
    // to compute the shift information
    result.weight = 1.0f;
    return false;
  }

  // Get the information of the original photon beam path
  const Path *source = beam->path;
  const int currVertex = beam->edgeID;
  const PathVertex *parentPhoton = source->vertex(currVertex);
  const PathEdge *parentEdge = source->edge(currVertex);
  const PathVertex *baseVertex = source->vertex(currVertex + 1);

  // Recompute the position of intersection of the base photon beam
  SAssert(!beam->isInvalid());
  Point basePos = beam->getPos(kRec.v);

  // Compute the weight for the unchanged part
  Spectrum shiftPhotonWeight = source->vertex(0)->weight[EImportance] *
      source->vertex(0)->rrWeight *
      source->edge(0)->weight[EImportance];
  for (int i = 1; i <= currVertex - 1; ++i) {
    shiftPhotonWeight *= source->vertex(i)->weight[EImportance] *
        source->vertex(i)->rrWeight *
        source->edge(i)->weight[EImportance];
  }

  // Compute the weight towards the new photon
  Float pdfKernelAndDist = kRec.pdf();
  ShiftRecord sRec;
  {
    // In this case, give the parent parent location
    diffuseReconnectionPhotonBeam(sRec, newPos, basePos,
                                  newPBDir, newPBDist,
                                  baseVertex, parentPhoton, parentEdge, pdfKernelAndDist, beam->longBeams,
                                  currVertex >= 2 ? source->vertex(currVertex - 1)->getPosition() : Point(1.f));

    // pdfV and pdfW of the shifted path are already multiplied to sRec.pdf
    if (sRec.pdf == Float(0)) {
      result.weight = 1.0f;
      return false;
    }
  }

  // Jacobian for reconnection
  result.jacobian *= sRec.jacobian;

  Float shiftKernelPDF = kRec.kernelPDF(shiftRay,
                                        parentPhoton->getPosition(),
                                        newPBDir,
                                        newPBDist);
  if (shiftKernelPDF == 0) {
    result.weight = 1.0f;
    return false;
  }

  shiftPhotonWeight *= sRec.throughtput;

  // The shift camera throughput
  Spectrum
      eyeShiftContrib = shiftGP.getWeightBeam(currCameraEdge - 1) * // cummulative vertex(i-1).weight * edge(i-1).weight
      shiftGP.getWeightVertex(currCameraEdge);    // vertex(i).weight

  // Adding the term on the shift camera path (transmittance, scattering coef, phase function)
  MediumSamplingRecord mRecShift;
  Ray shiftRayEval(shiftRay.o, shiftRay.d, 0.f, shiftW, 0.f);            // transmittance up to w only
  medium->eval(shiftRayEval, mRecShift);

  PhaseFunctionSamplingRecord pRec(mRecShift, -newPBDir, -shiftRay.d, EImportance);
  Float phaseTerm = mRecShift.getPhaseFunction()->eval(pRec);
  shiftPhotonWeight *= mRecShift.transmittance * mRecShift.sigmaS * phaseTerm;

  // Finish
  result.shiftedFlux = shiftPhotonWeight * eyeShiftContrib * result.jacobian;
  result.weight = 0.5f;

  if (config.useMIS) {
    ///////////////////////// BASE PDF
    // * Eye subpath
    // * Light subpath
    Float basePdf = source->vertex(currVertex)->pdf[EImportance];            // parent vertex
    // Only account area measure up to the beam intersection
    // by first removing the previous geometry term
    basePdf *= (parentPhoton->getPosition() - baseVertex->getPosition()).lengthSquared();
    if (baseVertex->isOnSurface())
      basePdf /= absDot(baseVertex->getGeometricNormal(), parentEdge->d);
    // and compute the geometry factor up to the beam intersection
    // (as the intersection is in volume, so no cosine at the intersection is needed)
    basePdf /= (parentPhoton->getPosition() - basePos).lengthSquared();
    // Account for sampling due to
    // and then account for the probability of sampling the intersection
    basePdf *= pdfKernelAndDist;

    ///////////////////////// SHIFT PDF
    // * Eye subpath
    Float offsetPdf = shiftKernelPDF;
    // * Light subpath
    offsetPdf *= sRec.pdf;
    if (offsetPdf == Float(0) || basePdf == Float(0)) {
      //SLog(EWarn, "Invalid path");
      result.weight = 1.0f;
      return false;
    }

    const Float sensorPart = shiftGP.sensorMIS(currCameraEdge, *baseGather,
                                               shiftW, kRec.w);

    if (config.powerHeuristic) {
      result.weight = 1.0f / (1.0f + powerOfTwo(sensorPart * offsetPdf * result.jacobian / basePdf));
    } else {
      result.weight = 1.0f / (1.0f + sensorPart * offsetPdf * result.jacobian / basePdf);
    }

  }

  return true;

}

struct VertexRecord {
  Point p;
  Vector wi, wo;          // wi points to light, wo points to camera (EImportance mode)
  Intersection its;       // only for surface classification use
  const Medium *medium;
  PathVertex::EVertexType type;
  VertexRecord() : medium(nullptr), type(PathVertex::EInvalid) {

  }
};

void BeamGradRadianceQuery::cacheSourcePath(int b, int c,
                                            const LTPhotonBeam *beam, const BeamKernelRecord &kRec) {
  const auto &oriSource = (const Path &) *beam->path;
  if (cachePath.length() == 0) {
    // We have to move the original vertex c to the actual beam-beam intersection
    // So we modify:
    // * vertex(c - 1) due to geometry term change
    // * edge(c - 1) due to transmittance change
    // * vertex(c) due to position change
    cachePath.append(oriSource, 0, c - 1);    // Up to vertex c - 2 only
    cachePath.append(oriSource.edge(c - 2));  // Add the edge to c - 1

    PathEdge *cacheNewEdge;
    PathVertex *cacheVertexVolume;
    PathVertex *cacheVertexParent;

    // --- Vertex c - 1
    cacheVertexParent = oriSource.vertex(c - 1)->clone(thdata.pool);
    cachePath.append(cacheVertexParent);

    // --- Copy and modify the last edge on the base light path
    cacheNewEdge = oriSource.edge(c - 1)->clone(thdata.pool);
    cacheNewEdge->length = kRec.v;
    cacheNewEdge->pdf[EImportance] = kRec.pdf();        // Included sinTheta in Beam1D case
    cacheNewEdge->weight[EImportance] = kRec.beamTrans / cacheNewEdge->pdf[EImportance];
    cachePath.append(cacheNewEdge);

    // --- Create the vertex (inside the media)
    cacheVertexVolume = thdata.pool.allocVertex();
    memset(cacheVertexVolume, 0, sizeof(PathVertex));
    MediumSamplingRecord &cacheMRec = cacheVertexVolume->getMediumSamplingRecord();
    cacheMRec.t = kRec.v;
    cacheMRec.time = 0.f;
    cacheMRec.p = beam->getPos(kRec.v); // Get the original pos.
    cacheMRec.medium = beam->medium;
    cacheVertexVolume->type = PathVertex::EMediumInteraction;
    cacheVertexVolume->measure = EArea;
    cacheVertexVolume->sampledComponentIndex = -1;
    cachePath.append(cacheVertexVolume);

    // Adjust the pdf of the parent vertex to sample up to this intersection only
    if (oriSource.vertex(c - 1)->measure != EDiscrete) {
      Float oldG = fastGOp(oriSource, c - 1, c);
      Float newG = fastGOp(cachePath, c - 1, c);
      cacheVertexParent->pdf[EImportance] *= newG / oldG;
    }
  }
}

bool BeamGradRadianceQuery::shiftBeamME(int b, int c,
                                        const LTPhotonBeam *beam, const ShiftGatherPoint &shiftGP,
                                        const Ray &shiftRay, Float shiftW, GradientSamplingResult &result,
                                        const BeamKernelRecord &kRec, const Point &newPos) {

  // Perform a walk through a specular chain to a non-specular vertex
  this->cacheSourcePath(b, c, beam, kRec);

  // Generates the proposal path
  Path proposal;
  {
    // Construct the new vertex
    PathVertex shiftVertex;
    memset(&shiftVertex, 0, sizeof(PathVertex));
    MediumSamplingRecord &cacheMRec = shiftVertex.getMediumSamplingRecord();
    cacheMRec.t = 0.f;
    cacheMRec.p = newPos; // Get the original pos.
    cacheMRec.medium = beam->medium;
    shiftVertex.type = PathVertex::EMediumInteraction;
    shiftVertex.measure = EArea;
    shiftVertex.sampledComponentIndex = -1;

    MEShiftVolBeamShiftFailed.incrementBase();
    Float radius = beam->getRadius();
    if (!generateShiftPathME(cachePath, proposal, b, c, thdata.pool,
                             thdata.offsetGenerator.get(), shiftVertex,
                             radius * config.relaxME,
                             baseCameraRay(shiftW - baseCameraRay.mint),
                             shiftRay(shiftW - shiftRay.mint))) {
      for (int i = b; i <= c; ++i) {
        thdata.pool.release(proposal.edge(i - 1));
        thdata.pool.release(proposal.vertex(i));
      }
      result.weight = 1.0f;
      ++MEShiftVolBeamShiftFailed;
      return false;
    }
  }

  MEShiftVolBeamContribFailed.incrementBase();
  ShiftRecord sRecME;
  if (!ShiftME(sRecME, cachePath, proposal, b, c, true)) {
    for (int i = b; i <= c; ++i) {
      thdata.pool.release(proposal.edge(i - 1));
      thdata.pool.release(proposal.vertex(i));
    }
    result.weight = 1.0f;
    ++MEShiftVolBeamContribFailed;
    return false;
  }

  MEShiftVolBeamContribZero.incrementBase();
  if (sRecME.throughtput.isZero() || sRecME.pdf == 0.0f) {
    ++MEShiftVolBeamContribZero;
  }

  Float shiftKernelPDF = kRec.kernelPDF(shiftRay,
                                        proposal.vertex(c - 1)->getPosition(),
                                        proposal.edge(c - 1)->d,
                                        proposal.edge(c - 1)->length);
  if (shiftKernelPDF == 0) {
    for (int i = b; i <= c; ++i) {
      thdata.pool.release(proposal.edge(i - 1));
      thdata.pool.release(proposal.vertex(i));
    }
    result.weight = 1.0f;
    return false;
  }

  // Evaluate Jacobian (area measure)
  result.jacobian *= sRecME.jacobian;

  // Jacobian in area integral in original space
  SpecularManifold *manifold = thdata.offsetGenerator->getSpecularManifold();

  // Jacobian computed using det
  const Float detProposed = manifold->det(proposal, b, c);
  if (cacheDetSource == -1) {
    cacheDetSource = manifold->det(cachePath, b, c);
  }
  result.jacobian *= detProposed / cacheDetSource;

  // TEST:
  //result.jacobian = 1.0f;

  if (result.jacobian <= 0.0 || !std::isfinite(result.jacobian)) {
    SLog(EWarn, "Invalid jacobian %g %d %d", result.jacobian, b, c);
    result.weight = 1.0f;
    return false;
  }

  Spectrum shiftPhotonWeight = sRecME.throughtput;

  // Adding the term on the shift camera path (transmittance, scattering coef, phase function)
  MediumSamplingRecord mRecShift;
  Ray shiftRayEval(shiftRay.o, shiftRay.d, 0.f, shiftW, 0.f);
  medium->eval(shiftRayEval, mRecShift);

  PhaseFunctionSamplingRecord pRec(mRecShift, -proposal.edge(c - 1)->d, -shiftRay.d, EImportance);
  Float phaseTerm = mRecShift.getPhaseFunction()->eval(pRec);
  shiftPhotonWeight *= mRecShift.transmittance * mRecShift.sigmaS * phaseTerm;

  // The shift camera throughput (with base path pdf multiplied)
  Spectrum
      eyeShiftContrib = shiftGP.getWeightBeam(currCameraEdge - 1) * // cummulative vertex(i-1).weight * edge(i-1).weight
      shiftGP.getWeightVertex(currCameraEdge);    // vertex(i).weight

  // Finish
  result.shiftedFlux = shiftPhotonWeight * eyeShiftContrib * result.jacobian;
  result.weight = 0.5f;
  // result.jacobian assigned above

  if (config.useMIS) {
    // No need to multiply pdf of the common part in light subpath
    Float offsetPdf = sRecME.pdf                                                         // light subpath
        * shiftKernelPDF;

    Float basePdf = 1.f; // camera subpath, up to intersection
    for (int i = b; i < c; ++i) {   // light subpath
      basePdf *= cachePath.vertex(i)->pdf[EImportance];
      basePdf *= cachePath.edge(i)->pdf[EImportance];
    }

    if (basePdf == Float(0)) {
      SLog(EWarn, "Invalid base path. This case should not happen.");
      result.weight = 0.0f;
    } else if (offsetPdf == Float(0)) {
      SLog(EWarn, "Invalid offset path. This case should not happen: %f, %f", sRecME.pdf, shiftKernelPDF);
      result.weight = 1.0f;
    } else {
      const Float sensorPart = shiftGP.sensorMIS(currCameraEdge, *baseGather,
                                                 shiftW, kRec.w);
      if (config.powerHeuristic) {
        result.weight = 1.0f / (1.0f + powerOfTwo(sensorPart * result.jacobian * (offsetPdf / basePdf)));
      } else {
        result.weight = 1.0f / (1.0f + sensorPart * result.jacobian * (offsetPdf / basePdf));
      }
    }
  }

  for (int i = b; i <= c; ++i) {
    thdata.pool.release(proposal.edge(i - 1));
    thdata.pool.release(proposal.vertex(i));
  }
  return true;
}

bool BeamGradRadianceQuery::shiftNull3D(const LTPhotonBeam *beam, const ShiftGatherPoint &shiftGP,
                                        const Ray &shiftRay, Float shiftW, GradientSamplingResult &result,
                                        const BeamKernelRecord &kRec, BeamKernelRecord &kRecShift) {
  if (!kRecShift.isValid()) {
    result.weight = 1.0f;
    return false;
  }

  Spectrum
      eyeShiftContrib = shiftGP.getWeightBeam(currCameraEdge - 1) * // cummulative vertex(i-1).weight * edge(i-1).weight
      shiftGP.getWeightVertex(currCameraEdge);    // vertex(i).weight


  kRecShift.contrib *= kRecShift.pdf() / kRec.pdf(); // Use the right PDF
  result.jacobian = 1.0; // 3D changes
  result.shiftedFlux = kRecShift.contrib * eyeShiftContrib * result.jacobian;
  result.weight = 0.5f;

  if (config.useMIS) {
    Float basePdf = kRec.pdf();
    Float offsetPdf = kRecShift.pdf();

    if (offsetPdf == Float(0) || basePdf == Float(0)) {
      result.weight = 1.0f;
      return false;
    }

    const Float sensorPart = shiftGP.sensorMIS(currCameraEdge, *baseGather,
                                               kRecShift.w, kRec.w);
    if (config.powerHeuristic) {
      result.weight = 1.0f / (1.0f + powerOfTwo(sensorPart * result.jacobian * (offsetPdf / basePdf)));
    } else {
      result.weight = 1.0f / (1.0f + sensorPart * result.jacobian * (offsetPdf / basePdf));
    }

  }

  return true;
}

MTS_NAMESPACE_END
