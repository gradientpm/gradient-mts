#pragma once

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

inline Float getRadiusRayDifferential(const RayDifferential &ray, Float totalDist) {

  if (ray.hasDifferentials) {  // nbComponentExtra == 0 &&
    Point posProj = ray.o + ray.d * totalDist;
    Point rX = ray.rxOrigin + ray.rxDirection * totalDist;
    Point rY = ray.ryOrigin + ray.ryDirection * totalDist;
    Float dX = (rX - posProj).length();
    Float dY = (rY - posProj).length();

    Float r = std::max(dX, dY);
    return r;
  } else {
    SLog(EError, "No ray differential");
    return 0.f;
  }
}

MTS_NAMESPACE_END
