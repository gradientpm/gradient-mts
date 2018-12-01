#if !defined(__BEAMS_3D_INTERSECT_H)
#define __BEAMS__BEAMS_3D_INTERSECT_H_H

#include <mitsuba/render/photonmap.h>
#include <mitsuba/render/gatherproc.h>

MTS_NAMESPACE_BEGIN

/**
 * Cylinder intersection with buggy cap
 * @param rCylinder The ray that describe the cylinder
 * @param view The view ray which we want to compute the intersection
 * @param radius The radius of the cylinder
 * @param tNear out tNear
 * @param tFar out tFar
 * @return return true if there is a valid overlap with the cylinder.
 */
inline bool cylinderIntersectionBug(const Ray &rCylinder, const Ray &view, Float radius, double &tNear, double &tFar) {
  // This function come use 1D intersection to compute
  // 3D kernel intersection.
  // Compare to the slow intersection, a "bias" can be made
  // On the beam sphere cup (that is not there).

  const Vector d1d2c = cross(view.d, rCylinder.d);
  const float
      sinThetaSqr = dot(d1d2c, d1d2c); // Square of the sine between the two lines (||cross(d1, d2)|| = sinTheta).

  // Slower code to test if the lines are too far apart.
  // oDistance = absDot((O2 - O1), d1d2c) / d1d2c.size();
  // if(oDistance*oDistance >= maxDistSqr) return false;

  const float ad = dot((rCylinder.o - view.o), d1d2c);

  // Lines too far apart.
  if (ad * ad >= (radius * radius) * sinThetaSqr)
    return false;

  // Cosine between the two lines.
  const float d1d2 = dot(view.d, rCylinder.d);
  const float d1d2Sqr = d1d2 * d1d2;
  const float d1d2SqrMinus1 = d1d2Sqr - 1.0f;

  // Parallel lines?
  if (d1d2SqrMinus1 < 1e-5f && d1d2SqrMinus1 > -1e-5f)
    return false;

  const float d1O1 = dot(view.d, Vector(view.o));
  const float d1O2 = dot(view.d, Vector(rCylinder.o));

  Float w =
      (d1O1 - d1O2 - d1d2 * (dot(rCylinder.d, Vector(view.o)) - dot(rCylinder.d, Vector(rCylinder.o)))) / d1d2SqrMinus1;

  // Out of range on ray 1.
  if (w <= view.mint || w >= view.maxt)
    return false;

  Float v = (w + d1O1 - d1O2) / d1d2;

  // Out of range on ray 2.
  if (v <= 0.0 || v >= rCylinder.maxt || std::isnan(v))
    return false;

  // Out of range due to beam section
  if (rCylinder.mint >= v || rCylinder.maxt < v)
    return false;

  const float sinThetaConst = std::sqrt(sinThetaSqr);
  Float u = std::abs(ad) / sinThetaConst;

  Float rad = math::safe_sqrt(radius * radius - u * u) / sinThetaConst;
  tNear = w - rad;
  tFar = w + rad;

  return true;
}

inline bool cylinderIntersection(const Ray &rCylinder, const Ray &view, Float radius, double &tNear, double &tFar) {
  // SH: is Epsilon acceptable?
  SAssert(rCylinder.mint == 0.f);

  // Early intersection test (from optimized beam-beam 1D intersection)
  // TODO: Verify that is test is correct or not.
  const Vector d1d2c = cross(view.d, rCylinder.d);
  const float sinThetaSqr = dot(d1d2c, d1d2c);
  const float ad = dot((rCylinder.o - view.o), d1d2c);
  if (ad * ad >= (radius * radius) * sinThetaSqr)
    return false; // No intersection

  // TODO: See if we can get rid of this function
  // and compute the intersection directly in the world domain.
  Transform objToWorld = Transform::translate(Vector(rCylinder.o)) *
      Transform::fromFrame(Frame(rCylinder.d));
  Transform worldToObject = objToWorld.inverse();
  const Float lMax = rCylinder.maxt;

  ///// Test the intersection
  Ray ray = view;
  worldToObject(view, ray); // Put in local coordinate

  const double
      ox = ray.o.x,
      oy = ray.o.y,
      dx = ray.d.x,
      dy = ray.d.y;

  const double A = dx * dx + dy * dy;
  const double B = 2 * (dx * ox + dy * oy);
  const double C = ox * ox + oy * oy - radius * radius;

  if (!solveQuadraticDouble(A, B, C, tNear, tFar))
    return false; // No intersection

  if (tNear > view.maxt || tFar < 0) {
    return false; // Non valid intersection
  }

  const double zPosNear = ray.o.z + ray.d.z * tNear;
  const double zPosFar = ray.o.z + ray.d.z * tFar;

  // http://woo4.me/wootracer/cylinder-intersection/
  // beam intersection with flat cap.
  if (zPosNear < 0) {
    if (zPosFar < 0)
      return false;
    // Hit the c
    float th = tNear + (tFar - tNear) * (zPosNear) / (zPosNear - zPosFar);
    tNear = th;
    return true;
  } else if (zPosNear >= 0 && zPosNear < lMax) {
    // Hit the cylinder, keep all values
    return true;
  } else if (zPosNear > lMax) {
    if (zPosFar > lMax)
      return false;
    float th = tNear + (tFar - tNear) * (zPosNear - lMax) / (zPosNear - zPosFar);
    tNear = th;
    return true;
  }
  return false;
}

inline bool cylinderIntersectionDebug(const Ray &rCylinder,
                                      const Ray &view,
                                      Float radius,
                                      double &tNear,
                                      double &tFar) {
  // SH: is Epsilon acceptable?
  SAssert(rCylinder.mint == 0.f);

  // Early intersection test (from optimized beam-beam 1D intersection)
  // TODO: Verify that is test is correct or not.
  /*
  const Vector d1d2c = cross(view.d, rCylinder.d);
  const float sinThetaSqr = dot(d1d2c, d1d2c);
  const float ad = dot((rCylinder.o - view.o), d1d2c);
  if (ad * ad >= (radius * radius) * sinThetaSqr) {
      SLog(EInfo, "Failed at sintheta");
      return false; // No intersection
  }*/

  // TODO: See if we can get rid of this function
  // and compute the intersection directly in the world domain.
  Transform objToWorld = Transform::translate(Vector(rCylinder.o)) *
      Transform::fromFrame(Frame(rCylinder.d));
  Transform worldToObject = objToWorld.inverse();
  const Float lMax = rCylinder.maxt;

  ///// Test the intersection
  Ray ray;
  worldToObject(view, ray); // Put in local coordinate

  Ray rLocalCylinder;
  worldToObject(rCylinder, rLocalCylinder);
  SLog(EInfo,
       "rLocalCylinder %f %f %f %f %f %f",
       rLocalCylinder.o.x,
       rLocalCylinder.o.y,
       rLocalCylinder.o.z,
       rLocalCylinder.d.x,
       rLocalCylinder.d.y,
       rLocalCylinder.d.z);
  SLog(EInfo, "ray %f %f %f %f %f %f", ray.o.x, ray.o.y, ray.o.z, ray.d.x, ray.d.y, ray.d.z);

  const double
      ox = ray.o.x,
      oy = ray.o.y,
      dx = ray.d.x,
      dy = ray.d.y;

  const double A = dx * dx + dy * dy;
  const double B = 2 * (dx * ox + dy * oy);
  const double C = ox * ox + oy * oy - radius * radius;
  SLog(EInfo, "A B C %f %f %f", A, B, C);

  if (!solveQuadraticDouble(A, B, C, tNear, tFar)) {
    SLog(EInfo, "Failed at quadratic");
    return false; // No intersection
  }

  if (tNear > view.maxt || tFar < 0) {
    SLog(EInfo, "Failed at no valid intersection");
    return false; // No valid intersection
  }

  const double zPosNear = ray.o.z + ray.d.z * tNear;
  const double zPosFar = ray.o.z + ray.d.z * tFar;

  // http://woo4.me/wootracer/cylinder-intersection/
  // beam intersection with flat cap.
  if (zPosNear < 0) {
    if (zPosFar < 0) {
      SLog(EInfo, "Failed at zPosFar");
      return false;
    }
    // Hit the c
    float th = tNear + (tFar - tNear) * (zPosNear) / (zPosNear - zPosFar);
    tNear = th;
    return true;
  } else if (zPosNear >= 0 && zPosNear < lMax) {
    // Hit the cylinder, keep all values
    return true;
  } else if (zPosNear > lMax) {
    if (zPosFar > lMax) {
      SLog(EInfo, "Failed at zPos lMax");
      return false;
    }
    float th = tNear + (tFar - tNear) * (zPosNear - lMax) / (zPosNear - zPosFar);
    tNear = th;
    return true;
  }

  SLog(EInfo, "Failed at the end");
  return false;
}

MTS_NAMESPACE_END

#endif