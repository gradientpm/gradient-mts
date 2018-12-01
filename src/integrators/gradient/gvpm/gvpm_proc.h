#pragma once

#include <mitsuba/render/particleproc.h>
#include <mitsuba/bidir/mempool.h>

#include "gvpm_struct.h"
#include "gvpm_beams.h"
#include "gvpm_accel.h"

MTS_NAMESPACE_BEGIN

class GradientPhotonProcess : public ParticleProcess {
public:
  /**
   * Create a new process for parallel photon gathering
   * \param type
   *     Specifies the type of requested photons (surface/caustic/volume)
   * \param photonCount
   *     Specifies the number of requested photons
   * \param granularity
   *     Size of the internally used work units (in photons)
   * \param isLocal
   *     Should the parallel process only be executed locally? (sending
   *     photons over the network may be unnecessary and wasteful)
   * \param autoCancel
   *     Indicates if the gathering process should be canceled if there
   *     are not enough photons generated
   * \param progressReporterPayload
   *    Custom pointer payload to be delivered with progress messages
   */
  GradientPhotonProcess(const GPMConfig &config, const Point &sensorPos,
                        size_t granularity,
                        bool isLocal, bool autoCancel, const void *progressReporterPayload,
                        std::vector<ShootingThreadData> &shootingData,
                        Float beamRadius = 0.f);

  /**
   * Once the process has finished, this returns a reference
   * to the (still unbalanced) photon map
   */
  inline GPhotonMap *getPhotonMap() { return m_photonMap; }
  inline GPhotonMap *getPhotonVolumeMap() { return m_photonVolumeMap; }
  inline LTBeamMap *getBeamMap() { return m_beamMap; }
  /**
   * \brief Return the number of discarded photons
   *
   * Due to asynchronous processing, some excess photons
   * will generally be produced. This function returns the number
   * of excess photons that had to be discarded. If this is too
   * high, the granularity should be decreased.
   */
  inline size_t getExcessPhotons() const { return m_excess; }
  bool outVolumeCapacity() const {
    if (!m_storeBeams)
      return m_photonVolumeMap->outCapacity();
    else
      return m_beamMap->outCapacity();
  }
  size_t nbSkip() const {
    if (m_storeBeams)
      return m_beamMap->nbSkip();
    else
      return m_photonVolumeMap->nbSkip();
  }

  /**
   * \brief Lists the nuber of particles that had to be shot
   * in order to fill the photon map.
   */
  inline size_t getShotParticlesSurface() const { return m_numShotSurface; }
  inline size_t getShotParticlesVolume() const { return m_numShotVolume; }

  // ======================================================================
  /// @{ \name ParallelProcess implementation
  // ======================================================================

  bool isLocal() const override;
  ref<WorkProcessor> createWorkProcessor() const override;

  int pushSurfaceLT(int idWorker, const Path *lt);
  int pushVolumeLT(int idWorker, const Path *lt);

  void processResult(const WorkResult *wr, bool cancelled) override;
  EStatus generateWork(WorkUnit *unit, int worker) override;
  size_t size() const {
    size_t total = 0;
    if (m_photonMap != nullptr)
      total += m_photonMap->size();
    if (m_photonVolumeMap != nullptr)
      total += m_photonVolumeMap->size();
    if (m_beamMap != nullptr)
      total += m_beamMap->size();
    return total;
  }

  /// @}
  // ======================================================================

  MTS_DECLARE_CLASS()
protected:
  /// Virtual destructor
  ~GradientPhotonProcess() override = default;

  /**
   * \brief Checks if the configuration of needed, generated and shot
   * photons indicates an unsuccessful progress of the gathering. This
   * check is taken from PBRT.
   */
  inline bool unsuccessful(size_t needed, size_t gen, size_t shot) {
    return (gen < needed && (gen == 0 || gen < shot / 1024));
  }
protected:
  const GPMConfig &m_config;

  ref<GPhotonMap> m_photonMap;
  int m_nbSurfaceLT;
  ref<GPhotonMap> m_photonVolumeMap;
  int m_nbVolumeLT;

  bool m_isLocal;
  bool m_autoCancel;

  size_t m_excess;
  size_t m_numShotSurface, m_numShotVolume;

  // For the photon beams
  ref<LTBeamMap> m_beamMap;
  Float m_beamRadius;
  bool m_storeBeams;

  // Pools used to generate light paths
  std::vector<ShootingThreadData> &m_thdata;
  mutable int m_workerID;

  std::vector<bool> m_workGenerated;
};

MTS_NAMESPACE_END