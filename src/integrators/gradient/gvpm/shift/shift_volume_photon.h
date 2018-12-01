#pragma once

#include <mitsuba/mitsuba.h>

#include "shift_cameraPath.h"
#include "shift_utilities.h"

#include "../gvpm_accel.h"

#ifndef MITSUBA_SHIFT_VOLUME_PHOTON_H
#define MITSUBA_SHIFT_VOLUME_PHOTON_H

MTS_NAMESPACE_BEGIN

class AbstractVolumeGradientRecord {
public:
  // The information about the shift
  Spectrum mediumFlux;                    // volume flux
  Spectrum shiftedMediumFlux[4];
  Spectrum weightedMediumFlux[4];
protected:

  // The original informations
  Scene *scene;
  GatherPoint *baseGather;
  std::vector<ShiftGatherPoint> &shiftGPs;
  Ray baseRay;
  size_t currEdge;
  const Medium *medium;

  // Additional informations
  const GPMConfig &config;
  GPMThreadData &thdata;
  Sampler *sampler;

public:
  AbstractVolumeGradientRecord(Scene *scene, GatherPoint *gp,
                               const GPMConfig &config,
                               GPMThreadData &thdata,
                               std::vector<ShiftGatherPoint> &_shiftGPs,
                               size_t _currEdge, Sampler *_sampler) :
      scene(scene), baseGather(gp), shiftGPs(_shiftGPs),
      currEdge(_currEdge), config(config), thdata(thdata), sampler(_sampler) {
    medium = nullptr;
    clear();
  }
  virtual ~AbstractVolumeGradientRecord() = default;

  void clear() {
    for (int i = 0; i < 4; ++i) {
      shiftedMediumFlux[i] = Spectrum(0.f);
      weightedMediumFlux[i] = Spectrum(0.f);
    }
    mediumFlux = Spectrum(0.f);
  }
  void newRayBase(const Ray &_baseRay, const Medium* med) {
    baseRay = _baseRay;
    medium = med;
  }

protected:
  Point getShiftPos(const Ray &shiftRay,
                    Float radius,
                    Point basePhotonPos,
                    bool coherent = false);

  // The global function for the shifting
  bool shiftPhoton(const Point &offsetPos, const Path *lt, size_t currVertex, const ShiftGatherPoint &shiftGP,
                   const Ray &shiftRay, const MediumSamplingRecord& shiftMRec,
                   GradientSamplingResult &result, Float photonRadius,
                   Float pdfBaseRay = 1.f, Float pdfShiftRay = 1.f,
                   Float additionalJacobian = 1.f);

  /**
   * If wi points to a camera, the transport mode is ERadiance (like path tracing)
   * If wi points to a light source, the transport mode is EImportance (like photon tracing)
   */
  inline Spectrum getVolumePhotonContrib(const Spectrum &flux,
                                         const MediumSamplingRecord& shiftMRec,
                                         const Vector &wi,
                                         const Vector &wo,
                                         ETransportMode mode = ERadiance) {
    PhaseFunctionSamplingRecord pRec(shiftMRec, wi, wo, mode);
    return shiftMRec.sigmaS * flux * medium->getPhaseFunction()->eval(pRec);
  }

  bool shiftNull(const Spectrum &photonFlux, const Vector3 &photonWi, const ShiftGatherPoint &shiftGP,
                 const Ray &shiftRay, const MediumSamplingRecord& shiftMRec,
                 GradientSamplingResult &result,
                 Float pdfBaseRay = 1.f, Float pdfShiftRay = 1.f,
                 Float additionalJacobian = 1.f);

private:
  bool shiftPhotonDiffuse(const Point &offsetPos,
                          const Path *lt, size_t currVertex, const ShiftGatherPoint &shiftGP,
                          const Ray &shiftRay, const MediumSamplingRecord& shiftMRec,
                          GradientSamplingResult &result,
                          Float pdfBaseRay = 1.f, Float pdfShiftRay = 1.f,
                          Float additionalJacobian = 1.f);

  bool shiftPhotonMedium(const Point &offsetPos,
                         const Path *lt, size_t currVertex, const ShiftGatherPoint &shiftGP,
                         const Ray &shiftRay, const MediumSamplingRecord& shiftMRec,
                         GradientSamplingResult &result,
                         Float pdfBaseRay = 1.f, Float pdfShiftRay = 1.f,
                         Float additionalJacobian = 1.f);

  bool shiftPhotonManifold(int b, int c,
                           const Point &offsetPos, const Path *lt,
                           const ShiftGatherPoint &shiftGP,
                           const Ray &shiftRay, const MediumSamplingRecord& shiftMRec,
                           GradientSamplingResult &result, Float photonRadius,
                           Float pdfBaseRay = 1.f, Float pdfShiftRay = 1.f,
                           Float additionalJacobian = 1.f);

public:
#if HAVE_ADDITIONAL_STATS
  // To have some statistics about the shifts, per pixels
  ShiftStats shiftStats[4];
#endif

protected:
};

class VolumeGradientPositionQuery : public AbstractVolumeGradientRecord {
protected:
  MediumSamplingRecord shiftMRec[4];
  bool shiftMRecIntialized;
  bool validShiftDist[4];
  Float shiftDistCamera[4];

public:
  Float searchRadius;
  const MediumSamplingRecord *baseMRec;

public:
  VolumeGradientPositionQuery(Scene *scene, GatherPoint *gp, const GPMConfig &config, GPMThreadData &thdata,
                              std::vector<ShiftGatherPoint> &_shiftGPs, size_t _currEdge, Sampler *_sampler) :
      AbstractVolumeGradientRecord(scene, gp, config, thdata, _shiftGPs, _currEdge, _sampler),
      searchRadius(0), baseMRec(nullptr) {
    resetMediumRecCache();
  }

  void newRayBase(const Ray &_baseRay, const MediumSamplingRecord &_baseMRec, Float searchRadii) {
    AbstractVolumeGradientRecord::newRayBase(_baseRay, _baseMRec.medium);
    resetMediumRecCache();
    baseMRec = &_baseMRec;
    searchRadius = searchRadii;
  }

  void operator()(const GPhotonNodeKD &nodePhoton);

protected:
  virtual Float pdfBaseRay() const = 0;
  virtual Float pdfShiftRay(int id, const Medium *m, Float shiftDist) const = 0;

  void resetMediumRecCache() {
    for (int i = 0; i < 4; ++i) {
      validShiftDist[i] = false;
      shiftDistCamera[i] = -1.f;
    }
    shiftMRecIntialized = false;
  }
};

class VolumeGradientDistanceQuery : public VolumeGradientPositionQuery {
protected:
  Float baseDistPDF = 0.f;
  Float pdfSelSection = 1.f;
public:
  VolumeGradientDistanceQuery(Scene *scene, GatherPoint *gp, const GPMConfig &config, GPMThreadData &thdata,
                              std::vector<ShiftGatherPoint> &_shiftGPs, int _currEdge, Sampler *_sampler) :
      VolumeGradientPositionQuery(scene, gp, config, thdata, _shiftGPs, _currEdge, _sampler) {}

  void newRayBase(const Ray &_baseRay, const MediumSamplingRecord &_baseMRec, Float searchRadii, Float _baseDistPDF) {
    baseDistPDF = _baseDistPDF;
    VolumeGradientPositionQuery::newRayBase(_baseRay, _baseMRec, searchRadii);
  }

  // FIXME: Merge with new ray base as the information is redundant here
  void changeEdge(size_t _currEdge, Float _pdfSel) {
    currEdge = _currEdge;
    pdfSelSection = _pdfSel;
  }

protected:
  Float pdfBaseRay() const override {

    return baseDistPDF * pdfSelSection;
  }
  Float pdfShiftRay(int id, const Medium *m, Float shiftDist) const override {
    // We use the same PDF as the base path for the selection of the edge
    // because we want to keep the things coherent.
    return shiftMRec[id].pdfSuccess * pdfSelSection;
  }

};

class VolumeGradientBREQuery : public AbstractVolumeGradientRecord {
public:
  VolumeGradientBREQuery(Scene *scene, GatherPoint *gp,
                         const GPMConfig &config,
                         GPMThreadData &thdata,
                         std::vector<ShiftGatherPoint> &_shiftGPs, size_t _currEdge, Sampler *_sampler)
      : AbstractVolumeGradientRecord(scene, gp, config, thdata, _shiftGPs, _currEdge, _sampler) {
  }

  // The operator call when we do the gathering
  void operator()(const GPhotonNodeKD &nodePhoton, Float photonRadius, Float randValue);
};

MTS_NAMESPACE_END

#endif //MITSUBA_SHIFT_VOLUME_PHOTON_H
