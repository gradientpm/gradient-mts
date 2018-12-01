//
// Created by beltegeuse on 11/30/16.
//

#ifndef MITSUBA_BEAMS_STRUCT_H
#define MITSUBA_BEAMS_STRUCT_H

#include <mitsuba/render/photonmap.h>
#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

// PhotonBeam struct defines a beam from vertex i to vertex i + 1
// (see LTPhotonBeam struct in gvpm_beams.h for implementation).
// 
// If we evaluate contribution on a point on the beam, 
// flux has to be accounted only up to that point. There are two factors to adjust: geometry term, and transmittance. 
//
// For primal domain, we do not need to adjust geometry term because it cancels out in Monte Carlo estimation. 
// 
// For gradient domain, we have to take care of both factors. 

extern Float getRadiusRayDifferential(const RayDifferential &ray, Float totalDist);

struct PhotonBeam {
private:
  Point p1;
  Point p2;
  Vector dir;
  Float length;

public:
  const Medium *medium;
  Spectrum flux;
  int depth;
  Float radius;   //for differential beam this is the top caps radius

  // Just a protection if p2 is not set...
  bool invalid;
  bool longBeams;

  //differential data
  RayDifferential ray;

  PhotonBeam() :
      medium(0), flux(0.f), depth(-1), invalid(true) {
    p1 = Point(0.f);
    p2 = Point(0.f);
    dir = Vector(0.f);
    length = 0.f;
    longBeams = false;
    radius = 0.f;
  }

  PhotonBeam(const Point &p,
             const Medium *m,
             const Spectrum &f,
             int curDepth,
             Float beamRadius = 0.f) :
      p1(p), p2(0.f), medium(m), flux(f), depth(curDepth), radius(beamRadius) {
    invalid = true;
    longBeams = false;
    ray.hasDifferentials = false;
    length = 0.f;
  }

  void setDifferential(const RayDifferential &rayDiff) {
    ray = rayDiff;
    ray.o = p1;
  }

  void setEndPoint(const Point &p) {
    invalid = false;
    p2 = p;

    // Setup normalize dir and length of the beam
    dir = p2 - p1;
    length = dir.length();
    dir /= length;
  }

  bool isInvalid() const {
    return invalid;
  }

  void convertToLong(const Scene *scene) {
    Ray r(p1, normalize(p2 - p1), 0.f);
    Intersection its;
    if (!scene->rayIntersect(r, its))
      SLog(EError, "Impossible to convert to long beam...");
    p2 = its.p; // Change the last points position
    longBeams = true;
    //Fixme : @nicolas : must call setupDifferential
  }

  Float getProbSuccess(Float v) const {
    if (longBeams) {
      return 1.f;
    } else {
      Ray rayTrans(p1, getDir(), 0.f);
      rayTrans.mint = 0.f;
      rayTrans.maxt = v;

      MediumSamplingRecord mRec;
      medium->eval(rayTrans, mRec);

      return mRec.pdfSuccess;
    }
  }

  const Point &getOri() const {
    return p1;
  }
  const Point getEndPoint() const {
    return p1 + dir * length;
  }

  Point getPos(Float v) const {
    return p1 + getDir() * v;
  }

  const Vector &getDir() const {
    return dir;
  }

  Float getLength() const {
    return length;
  }

  Spectrum getContrib(Float v, const MediumSamplingRecord &mRecCamera, const Vector &d) const {
    Spectrum beamTransmittance;
    Float pdfFailure;
    return getContrib(v, mRecCamera, d, beamTransmittance, pdfFailure);
  }

  Spectrum getContrib(Float v, const MediumSamplingRecord &mRecCamera, const Vector &d,
                      Spectrum &beamTransmittance, Float &pdfFailure) const {
    // Compute the contribution of this beam to a receiving beam.
    // The gather point could be in the middle of the beam edge.
    //
    // v          : distance to a point on the current beam
    // mRecCamera : for transmittance evaluation at the receiving beam
    // d          : camera ray direction (wo = -d)


    Ray rayTrans(p1, getDir(), 0.f);
    rayTrans.mint = 0.f;
    rayTrans.maxt = v;

    MediumSamplingRecord mRec;
    medium->eval(rayTrans, mRec);
    beamTransmittance = mRec.transmittance;

    PhaseFunctionSamplingRecord pRec(mRecCamera, -getDir(), -d, EImportance);
    Float phaseTerm = mRecCamera.getPhaseFunction()->eval(pRec);

    Spectrum beamContrib = mRec.transmittance * mRecCamera.transmittance * mRec.sigmaS *
        flux * phaseTerm;

    if (!longBeams) {
      // If we didn't sample long beam during photon tracing,
      // need to include the russian roulette value to account for short beam
      if (mRec.pdfFailure == 0 && (!mRec.transmittance.isZero())) {
        SLog(EWarn, "Long beam creation is not possible: %s\n depth: %i", mRec.toString().c_str(), depth);
        pdfFailure = mRec.pdfFailure;
        return Spectrum(0.f);
      }

      // Here mRec.pdfFailure is the probability for sampling a beam with length > v
      // so the gather point can exist
      // See "Unifying Points, Beams, and Paths in Volumetric Light Transport Simulation" [Sec. 7.3]
      beamContrib /= mRec.pdfFailure;
      pdfFailure = mRec.pdfFailure;
    } else {
      pdfFailure = 1.f;
    }

    if (!beamContrib.isValid()) {
      SLog(EError, "Invalid beam contrib: %s [%f, %f, %f, %f, %f, %f]", beamContrib.toString().c_str(),
           mRec.transmittance.getLuminance(), mRecCamera.transmittance.getLuminance(),
           mRec.sigmaS, flux.getLuminance(), phaseTerm);
    }

    return beamContrib;
  }

  bool rayIntersect1D(const Ray &_ray,
                      Float tminBeam,
                      Float tmaxBeam,
                      Float &u,
                      Float &v,
                      Float &w,
                      Float &sinTheta) const {
    return rayIntersectInternal1D(getRadius(), _ray, tminBeam, tmaxBeam,
                                  u, v, w,
                                  sinTheta);
  }

  bool rayIntersectInternalShift(const Ray &_ray, const Vector3 &newBeamDir,
                                 Float &u, Float &v, Float &w,
                                 Float &sinTheta) const {

    const Vector d1d2c = cross(_ray.d, newBeamDir);
    const float
        sinThetaSqr = dot(d1d2c, d1d2c); // Square of the sine between the two lines (||cross(d1, d2)|| = sinTheta).

    // Slower code to test if the lines are too far apart.
    // oDistance = absDot((O2 - O1), d1d2c) / d1d2c.size();
    // if(oDistance*oDistance >= maxDistSqr) return false;

    const float ad = dot((p1 - _ray.o), d1d2c);

    // Lines too far apart.
    // We don't care with the shift
    //if (ad * ad >= (radius * radius) * sinThetaSqr)
    //    return false;

    // Cosine between the two lines.
    const float d1d2 = dot(_ray.d, newBeamDir);
    const float d1d2Sqr = d1d2 * d1d2;
    const float d1d2SqrMinus1 = d1d2Sqr - 1.0f;

    // Parallel lines?
    if (d1d2SqrMinus1 < 1e-5f && d1d2SqrMinus1 > -1e-5f)
      return false;

    const float d1O1 = dot(_ray.d, Vector(_ray.o));
    const float d1O2 = dot(_ray.d, Vector(p1));

    w = (d1O1 - d1O2 - d1d2 * (dot(newBeamDir, Vector(_ray.o)) - dot(newBeamDir, Vector(p1)))) / d1d2SqrMinus1;

    // Out of range on ray 1.
    if (w <= _ray.mint || w >= _ray.maxt)
      return false;

    v = (w + d1O1 - d1O2) / d1d2;

    // Out of range on ray 2.
    if (v <= 0.0 || std::isnan(v)) // || v >= getLength() infinite beam
      return false;

    const float sinThetaConst = std::sqrt(sinThetaSqr);

    u = std::abs(ad) / sinThetaConst;
    sinTheta = sinThetaConst;

    return true;
  }

  bool rayIntersectInternal1D(Float radius, const Ray &_ray, Float tminBeam, Float tmaxBeam,
                              Float &u, Float &v, Float &w,
                              Float &sinTheta) const {
    /*
     * Code taken from smallUPDT
     */

    // v: distance to intersection along this beam
    // w: distance to intersection along the other ray
    // u: distance between the two intersection points on this beam and the other ray

    // Compute the direction and the max distance
    // Of the photon beams

    const Vector d1d2c = cross(_ray.d, getDir());
    const float
        sinThetaSqr = dot(d1d2c, d1d2c); // Square of the sine between the two lines (||cross(d1, d2)|| = sinTheta).

    // Slower code to test if the lines are too far apart.
    // oDistance = absDot((O2 - O1), d1d2c) / d1d2c.size();
    // if(oDistance*oDistance >= maxDistSqr) return false;

    const float ad = dot((p1 - _ray.o), d1d2c);

    // Lines too far apart.
    if (ad * ad >= (radius * radius) * sinThetaSqr)
      return false;

    // Cosine between the two lines.
    const float d1d2 = dot(_ray.d, getDir());
    const float d1d2Sqr = d1d2 * d1d2;
    const float d1d2SqrMinus1 = d1d2Sqr - 1.0f;

    // Parallel lines?
    if (d1d2SqrMinus1 < 1e-5f && d1d2SqrMinus1 > -1e-5f)
      return false;

    const float d1O1 = dot(_ray.d, Vector(_ray.o));
    const float d1O2 = dot(_ray.d, Vector(p1));

    w = (d1O1 - d1O2 - d1d2 * (dot(getDir(), Vector(_ray.o)) - dot(getDir(), Vector(p1)))) / d1d2SqrMinus1;

    // Out of range on ray 1.
    if (w <= _ray.mint || w >= _ray.maxt)
      return false;

    v = (w + d1O1 - d1O2) / d1d2;

    // Out of range on ray 2.
    if (v <= 0.0 || v >= getLength() || std::isnan(v))
      return false;
    // Out of range due to beam section
    if (tminBeam >= v || tmaxBeam < v)
      return false;

    const float sinThetaConst = std::sqrt(sinThetaSqr);

    u = std::abs(ad) / sinThetaConst;
    sinTheta = sinThetaConst;

    return true;
  }

  /*
   * Sub beams creation
   */
  int nbSubBeams(const Float sizeSubbeams = 1.0) const {
    return (int) std::ceil(getLength() / sizeSubbeams);
  }

  std::string toString() const {
    std::ostringstream oss;
    oss << "Beam[" << endl
        << "  p1 = " << p1.toString() << "," << endl
        << "  p2 = " << p2.toString() << "," << endl
        << "  medium = " << (medium == 0 ? "NULL" : medium->toString()) << "," << endl
        << "  flux = " << flux.toString() << "," << endl
        << "  depth = " << depth << "," << endl
        << "]";
    return oss.str();
  }

  Float getRadius() const {
    return radius;
  }
};

MTS_NAMESPACE_END

#endif //MITSUBA_BEAMS_STRUCT_H
