#pragma once

#include "../gvpm_geoOps.h"
#include "../gvpm_struct.h"

#include <array>

#ifndef MITSUBA_SHIFT_UTILITIES_H
#define MITSUBA_SHIFT_UTILITIES_H

MTS_NAMESPACE_BEGIN

const Float D_EPSILON = (Float)(1e-14);

/// The structure to get the results after a shift
struct GradientSamplingResult {
  Spectrum shiftedFlux;           // contribution of the shifted path (weighted by density and Jacobian)
  Float weight;                   // MIS weight
  Float jacobian;              // Jacobian
  GradientSamplingResult() : shiftedFlux(0.0f), weight(1.0f), jacobian(1.f) {
  }
};

/// General shift record
struct ShiftRecord {
  Spectrum throughtput;
  Float pdf;
  Float jacobian;

  ShiftRecord() : throughtput(1.0f), pdf(0.0), jacobian(1.f) {
  }
};

/// Result of a half-vector duplication shift.
struct HalfVectorShiftResult {
  bool success;   ///< Whether the shift succeeded.
  Float jacobian; ///< Local Jacobian determinant of the shift.
  Vector3 wo;     ///< Tangent space outgoing vector for the shift.
};

/// Calculates the outgoing direction of a shift by duplicating the local half-vector.
inline HalfVectorShiftResult halfVectorShift(Vector3 tangentSpaceMainWi, Vector3 tangentSpaceMainWo,
                                             Vector3 tangentSpaceShiftedWi, Float mainEta, Float shiftedEta) {
  HalfVectorShiftResult result;

  if (Frame::cosTheta(tangentSpaceMainWi) * Frame::cosTheta(tangentSpaceMainWo) < (Float) 0) {
    // Refraction.

    // Refuse to shift if one of the Etas is exactly 1. This causes degenerate half-vectors.
    if (mainEta == (Float) 1 || shiftedEta == (Float) 1) {
      // This could be trivially handled as a special case if ever needed.
      result.success = false;
      return result;
    }

    // Get the non-normalized half vector.
    Vector3 tangentSpaceHalfVectorNonNormalizedMain;
    if (Frame::cosTheta(tangentSpaceMainWi) < (Float) 0) {
      tangentSpaceHalfVectorNonNormalizedMain = -(tangentSpaceMainWi * mainEta + tangentSpaceMainWo);
    } else {
      tangentSpaceHalfVectorNonNormalizedMain = -(tangentSpaceMainWi + tangentSpaceMainWo * mainEta);
    }

    // Get the normalized half vector.
    Vector3 tangentSpaceHalfVector = normalize(tangentSpaceHalfVectorNonNormalizedMain);

    // Refract to get the outgoing direction.
    Vector3 tangentSpaceShiftedWo = refract(tangentSpaceShiftedWi, tangentSpaceHalfVector, shiftedEta);

    // Refuse to shift between transmission and full internal reflection.
    // This shift would not be invertible: reflections always shift to other reflections.
    if (tangentSpaceShiftedWo.isZero()) {
      result.success = false;
      return result;
    }

    // Calculate the Jacobian.
    Vector3 tangentSpaceHalfVectorNonNormalizedShifted;
    if (Frame::cosTheta(tangentSpaceShiftedWi) < (Float) 0) {
      tangentSpaceHalfVectorNonNormalizedShifted = -(tangentSpaceShiftedWi * shiftedEta + tangentSpaceShiftedWo);
    } else {
      tangentSpaceHalfVectorNonNormalizedShifted = -(tangentSpaceShiftedWi + tangentSpaceShiftedWo * shiftedEta);
    }

    Float hLengthSquared = tangentSpaceHalfVectorNonNormalizedShifted.lengthSquared()
        / (D_EPSILON + tangentSpaceHalfVectorNonNormalizedMain.lengthSquared());
    Float WoDotH = std::abs(dot(tangentSpaceMainWo, tangentSpaceHalfVector))
        / (D_EPSILON + std::abs(dot(tangentSpaceShiftedWo, tangentSpaceHalfVector)));

    // Output results.
    result.success = true;
    result.wo = tangentSpaceShiftedWo;
    result.jacobian = hLengthSquared * WoDotH;
  } else {
    // Reflection.
    Vector3 tangentSpaceHalfVector = normalize(tangentSpaceMainWi + tangentSpaceMainWo);
    Vector3 tangentSpaceShiftedWo = reflect(tangentSpaceShiftedWi, tangentSpaceHalfVector);

    // FIXME: according to the paper, it must be the constant direction here, not the sampled direction.
    // In this case it does not matter because the cosine is the same.
    Float WoDotH = dot(tangentSpaceShiftedWo, tangentSpaceHalfVector) / dot(tangentSpaceMainWo, tangentSpaceHalfVector);
    Float jacobian = std::abs(WoDotH);

    result.success = true;
    result.wo = tangentSpaceShiftedWo;
    result.jacobian = jacobian;
  }

  return result;
}

inline ELightShiftType getTypeShift(const Path *lt, size_t currVertex, int &b) {
  SAssert(currVertex > 0);

  // Search next diffuse vertex
  b = -1;
  for (size_t i = currVertex - 1; i > 0 && b == -1; --i) {
    EVertexType vertexType = VertexClassifier::type(*lt->vertex(i), lt->vertex(i)->sampledComponentIndex);
    b = vertexType == VERTEX_TYPE_DIFFUSE ? i : -1;
  }

  if (b == -1) {
    // Error?
    return EInvalidShift;
  } else if (b + 1 == currVertex) {
    // Just use diffuse reconnection here
    return EDiffuseShift;
  } else if (lt->vertex(currVertex - 1)->getType() == PathVertex::EMediumInteraction) {
    // If the parent vertex is inside the medium,
    // we can use medium shift
    return EMediumShift;
  }

  // In this case, we use our more robust shift
  return EManifoldShift;
}

inline EVertexType vertexType(Intersection &its, int comp) {
  const BSDF *bsdf = its.getBSDF();
  return VertexClassifier::type(bsdf, its, comp);
}

inline EVertexType vertexType(Intersection &its, RayDifferential &rayDiff, int comp) {
  const BSDF *bsdf = its.getBSDF(rayDiff);
  return VertexClassifier::type(bsdf, its, comp);
}

/// Result of a reconnection shift.
struct ReconnectionShiftResult {
  bool success;   ///< Whether the shift succeeded.
  Float jacobian; ///< Local Jacobian determinant of the shift.
  Vector3 wo;     ///< World space outgoing vector for the shift.
};

/// Tries to connect the offset path to a specific vertex of the main path.
inline ReconnectionShiftResult reconnectShift(const Scene *scene,
                                              Point3 mainSourceVertex,
                                              Point3 targetVertex,
                                              Point3 shiftSourceVertex,
                                              Vector3 targetNormal,
                                              Float time) {
  ReconnectionShiftResult result;

  // Check visibility of the connection.
  Ray shadowRay(shiftSourceVertex,
                targetVertex - shiftSourceVertex,
                Epsilon, 1.0 - ShadowEpsilon, time);
  if (scene->rayIntersect(shadowRay)) {
    // Since this is not a light sample, we cannot allow shifts through occlusion.
    result.success = false;
    return result;
  }

  // Calculate the Jacobian.
  Vector3 mainEdge = mainSourceVertex - targetVertex;
  Vector3 shiftedEdge = shiftSourceVertex - targetVertex;

  Float mainEdgeLengthSquared = mainEdge.lengthSquared();
  Float shiftedEdgeLengthSquared = shiftedEdge.lengthSquared();

  Vector3 shiftedWo = -shiftedEdge / math::safe_sqrt(shiftedEdgeLengthSquared);

  Float mainOpposingCosine = dot(mainEdge, targetNormal) / math::safe_sqrt(mainEdgeLengthSquared);
  Float shiftedOpposingCosine = dot(shiftedWo, targetNormal);

  Float jacobian = std::abs(shiftedOpposingCosine * mainEdgeLengthSquared)
      / (D_EPSILON + std::abs(mainOpposingCosine * shiftedEdgeLengthSquared));

  // Return the results.
  result.success = true;
  result.jacobian = jacobian;
  result.wo = shiftedWo;
  return result;
}

/// Describes the state of a ray that is being traced in the scene.
struct RayState {
  RayState() :
      throughput(Spectrum(0.0f)),
      eta(1.0f),
      alive(true) {
    pdf[EPdfSolidAngle] = pdf[EPdfArea] = 0.0f;
    jacobian[EPdfSolidAngle] = jacobian[EPdfArea] = 1.0f;
  }

  RayDifferential ray;             ///< Current ray.

  Spectrum throughput;             ///< Current throughput of the path.
  Float pdf[2];                       ///< PDF of the path. When the path is an offset path, this stores
  /// the PDF as if the path is sampled with the same events as of the base path.
  Float jacobian[2];

  Intersection prevIts;
  Intersection its;

  Float eta;                       ///< Current refractive index of the ray.
  ///< For R.R use only.
  bool
      alive;                      ///< Whether the path matching to the ray is still good. Otherwise it's an invalid offset path with zero PDF and throughput.
};

static int getVertexComponentType(const PathVertex *prev) {
  if (prev->isEmitterSample() ||
      prev->isMediumInteraction())  // FIXME: to support "glossy" phase function
  {
    // As emitter component type is 0, treat it as diffuse
    return BSDF::EDiffuseReflection;
  } else {
    return prev->componentType;
  }
}

inline bool computeVolumeContribution(const GPMConfig &config, const PathVertex *vPrev) {
  // Check the type of path
  if ((config.lightingInteractionMode & ESurf2Media) && (config.lightingInteractionMode & EMedia2Media)) {
    // In this case, all the lighting effect for the volume are enable
    // So nothing to do in this case, just continue the computation
  } else {
    // In this case, only one lighting mode is available (from surface or media only)
    // So we need to check the previous vertex and may cancel the computation in this case
    if (vPrev->isMediumInteraction() && !(config.lightingInteractionMode & EMedia2Media)) {
      return false; // Ignore this path
    }
    if ((vPrev->isSurfaceInteraction() || vPrev->isEmitterSample())
        && !(config.lightingInteractionMode & ESurf2Media)) {
      return false; // Ignore this path
    }
  }

  return !(config.bsdfInteractionMode != BSDF::EAll &&
      getVertexComponentType(vPrev) > 0 && !(getVertexComponentType(vPrev) & config.bsdfInteractionMode));

}

inline std::array<Point2, 4> generateOffsetPos(const Point2 &basePixel) {
  const Point2 rightPixel = basePixel + Point2(1, 0);
  const Point2 leftPixel = basePixel + Point2(-1, 0);
  const Point2 bottomPixel = basePixel + Point2(0, -1);
  const Point2 topPixel = basePixel + Point2(0, 1);
  return {leftPixel, rightPixel, topPixel, bottomPixel};
}

inline bool checkVisibility(const Scene *sc, const Point &p1, const Point &p2) {
  Ray r(p1, p2 - p1, Epsilon, 1 - Epsilon, 0.f);
  return !sc->rayIntersectAll(r);
}

MTS_NAMESPACE_END

#endif //MITSUBA_SHIFT_UTILITIES_H
