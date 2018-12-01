#ifndef MITSUBA_GVPM_BEAMS_H
#define MITSUBA_GVPM_BEAMS_H

#include <mitsuba/render/particleproc.h>
#include <mitsuba/bidir/mempool.h>
#include <mitsuba/bidir/path.h>

#include "../../photonmapper/beams.h"

#include <list>
#include <tuple>

MTS_NAMESPACE_BEGIN

// Quasi the same structure as photon beams
// However, still manage the light path information
// to be able to get the entire path
struct LTPhotonBeam : public PhotonBeam {
  size_t edgeID;
  const Path *path;
  size_t pathID;

  LTPhotonBeam(const Path *p, size_t i, Float _radius, int _pathID) :
      PhotonBeam(p->vertex(i)->getPosition(),
                 p->edge(i)->medium, Spectrum(1.f), i),
      edgeID(i), path(p), pathID(_pathID) {
    // Define a beam from vertex i to vertex i + 1
    // Flux does not include the transmittance on the last edge yet.
    for (size_t k = 0; k < i; k++) {
      flux *= p->vertex(k)->rrWeight *
          p->vertex(k)->weight[EImportance] *
          p->edge(k)->weight[EImportance];
    }
    flux *= p->vertex(i)->weight[EImportance];
    flux *= p->vertex(i)->rrWeight;

    const PathVertex *lastVertex = p->vertex(i + 1);
    setEndPoint(lastVertex->getPosition());

    // Set the radius as it is now attached to the beam
    radius = _radius;
  }
};

class LTBeamMap : public BeamMap<LTPhotonBeam> {
public:
  LTBeamMap(const GPMConfig &config,
            const Point &sensorPos) :
      BeamMap<LTPhotonBeam>(config.volumePhotonCount),
      m_config(config), m_sensorPos(sensorPos), m_nbBeamSkip(0) {
    m_nbLightPathAdded = 0;
  }

  inline int tryAppendLT(int idWorker,
                         const Path *lt,
                         Float _beamRadius,
                         int minDepth) {
    if (outCapacity()) {
      return -1;
    }

    int nbAppendVol = 0;

    for (size_t i = std::max(minDepth, 1); i < lt->edgeCount(); i++) {
      if (lt->edge(i)->medium != nullptr) {
        // Create the beam information
        if (!cameraHit(lt->vertex(i)->getPosition(), lt->vertex(i + 1)->getPosition())) {
          LTPhotonBeam pbeam(lt, i, _beamRadius, m_nbLightPathAdded);
          if (!tryAppend(idWorker, pbeam)) {
            break; // Out capacity
          }
          nbAppendVol++;
        }
      }
    }

    if (nbAppendVol == 0) {
      return -1;
    } else {
      m_nbLightPathAdded += 1;
    }

    return nbAppendVol;
  }

  inline size_t nbSkip() const { return m_nbBeamSkip; }

private:
  inline bool cameraHit(const Point &prev, const Point &p) {
    if (m_config.cameraSphere != 0.f) {
      if (isIntersectedPoint(m_sensorPos, prev, p, m_config.cameraSphere)) {
        ++m_nbBeamSkip;
        return true;
      }
    }
    return false;
  }

protected:
  // Keep the light path information
  // For later memory release
  GPMConfig m_config;
  Point m_sensorPos;
  size_t m_nbBeamSkip;

  int m_nbLightPathAdded;
};

MTS_NAMESPACE_END

#endif //MITSUBA_GVPM_BEAMS_H
