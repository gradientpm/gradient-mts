#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/bidir/path.h>

#ifndef MITSUBA_GVPM_GEOOPS_H
#define MITSUBA_GVPM_GEOOPS_H

MTS_NAMESPACE_BEGIN

inline Float geometryOpposingTerm(const Point &p, const Point &q, const Vector &qn) {
  Vector pq = q - p;
  return absDot(qn, normalize(pq)) / pq.lengthSquared();
}

inline Float geometryOpposingTerm(const Path &path, int b, int c) {
  const Point p = path.vertex(b)->getPosition();
  const Point q = path.vertex(c)->getPosition();
  const Vector pq = q - p;
  if (path.vertex(c)->isOnSurface()) {
    const Vector qn = path.vertex(c)->getGeometricNormal();
    return absDot(qn, normalize(pq)) / (pq.lengthSquared());
  } else {
    return 1.0 / pq.lengthSquared();
  }
}

inline Float fastGOp(const Path &path, int b, int c) {
  SAssert(b == c - 1);
  Float v = 1.f;
  if (path.vertex(c)->isOnSurface()) {
    v = absDot(path.edge(b)->d, path.vertex(c)->getGeometricNormal());
  }
  v /= (path.edge(b)->length * path.edge(b)->length);
#if 0
  Float vOri = geometryOpposingTerm(path, b, c);
  if(std::abs(1.f - (v / vOri)) > Epsilon) {
      SLog(EWarn, "Different values ");
  }
#endif
  return v;
}

inline Float cosine(const Path &path, int b, int c) {
  // TODO: can use edge->d
  Point p = path.vertex(b)->getPosition();
  Point q = path.vertex(c)->getPosition();
  Vector pn = path.vertex(b)->getGeometricNormal();
  return absDot(pn, normalize(q - p));
}

static inline Float cosineRatio(const Point &source, const Vector &sourceNormal,
                                const Point &target, const Vector &targetNormal) {
  Vector edge = normalize(target - source);
  return std::abs(dot(edge, targetNormal) / dot(edge, sourceNormal));
}

/// Returns whether point1 sees point2.
inline bool testVisibility(const Scene *scene, const Point3 &point1, const Point3 &point2, Float time) {
  Ray shadowRay;
  shadowRay.setTime(time);
  shadowRay.setOrigin(point1);
  shadowRay.setDirection(point2 - point1);
  shadowRay.mint = Epsilon;
  shadowRay.maxt = (Float) 1.0 - ShadowEpsilon;

  return !scene->rayIntersect(shadowRay);
}

/** This function go through null surface (medium transition)
 * @param scene
 * @param ray
 * @param its
 * @return return if a valid (not null) intersection is found
 */
inline bool intersectNULL(const Scene *scene, Ray &ray, Intersection &its) {
  int maxIntersection = 100;
  while (scene->rayIntersect(ray, its) && maxIntersection > 0) {
    maxIntersection--; // Protect to looping due to bad intersection

    // Check if it is a null intersection
    // If it is, just update the ray mint
    // FIXME: For now, assume that null BSDF is always medium transition
    if (its.isMediumTransition()) {
      ray.mint = its.t + Epsilon;
    } else {
      return true;
    }
  }
  if (maxIntersection < 0) {
    SLog(EWarn, "Problem loop intersection");
  }
  return false;
}

inline Float computeStepSize(const Float dist) {
  return std::max((Float) 0.1, dist / (Float) 10.0);
}

MTS_NAMESPACE_END

#endif //MITSUBA_GVPM_GEOOPS_H
