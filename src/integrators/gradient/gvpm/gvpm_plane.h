//
// Created by beltegeuse on 8/30/17.
//

#ifndef MITSUBA_GVPM_PLANE_H
#define MITSUBA_GVPM_PLANE_H

#include <mitsuba/render/particleproc.h>
#include <mitsuba/bidir/mempool.h>
#include <mitsuba/bidir/path.h>

#include "../../photonmapper/plane_struct.h"

#include <list>
#include <tuple>

MTS_NAMESPACE_BEGIN

struct LTPhotonPlane : public PhotonPlane {
  int edgeID;
  const Path *path;
  int pathID;

  LTPhotonPlane() :
      PhotonPlane(),
      edgeID(-1), path(nullptr), pathID(-1) {}

  LTPhotonPlane(const Path *p, int i, int _pathID, const Vector &addDir, const Float addLength) :
      PhotonPlane(p->vertex(i)->getPosition(),
                  p->edge(i)->d,
                  p->edge(i)->length,
                  addDir,
                  addLength,
                  p->edge(i)->medium,
                  Spectrum(1.f), i),
      edgeID(i), path(p), pathID(_pathID) {
    // Define a beam from vertex i to vertex i + 1
    // Flux does not include the transmittance on the last edge yet.
    for (int k = 0; k < i; k++) {
      _flux *= p->vertex(k)->rrWeight *
          p->vertex(k)->weight[EImportance] *
          p->edge(k)->weight[EImportance];
    }
    _flux *= p->vertex(i)->weight[EImportance];
    _flux *= p->vertex(i)->rrWeight;
  }

  /**
   * Helper to transform photon beams into plane
   * @param beam
   * @return
   */
  static LTPhotonPlane transformBeam(const LTPhotonBeam &beam, Sampler *sampler) {
    // TODO: Assume homogenous media
    // So the direction and the location of sample distance
    // Doesn't really matter here.

    // 1) Sample a distance (without taking into account the visibility)
    // The transmittance value of the first edge
    Ray rayNew(beam.getOri(), beam.getDir(), 0.f);
    MediumSamplingRecord mRecNew;
    beam.medium->sampleDistance(rayNew, mRecNew, sampler); // This sample the distance without vis tests.

    // 2) Sample additional scattering direction
    PhaseFunctionSamplingRecord pRec(mRecNew, -beam.getDir(), EImportance);
    do {
      beam.medium->getPhaseFunction()->sample(pRec, sampler);
    } while (absDot(beam.getDir(), pRec.wo) == 1);

    // 3) Create the photon plane and return it
    return LTPhotonPlane(beam.path, beam.edgeID, beam.pathID,
                         pRec.wo, mRecNew.t);
  }

  Vector offset(Float t0) const {
    return _w0 * t0;
  }

};

MTS_NAMESPACE_END

#endif //MITSUBA_GVPM_PLANE_H
