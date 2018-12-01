#pragma once
#if !defined(__MITSUBA_RENDER_PHOTONBREPROC_H_)
#define __MITSUBA_RENDER_PHOTONBREPROC_H_

#include <mitsuba/render/particleproc.h>
#include <mitsuba/render/photonmap.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Process for parallel photon map construction
 *
 * Given a number and type (surface/caustic/volume) of photons, this
 * class distributes the work over an arbitrary number of machines.
 * \ingroup librender
 */
class GatherPhotonCustomProcess : public ParticleProcess {
public:
  enum EGatherType {
    /// Surface photons (indirect on diffuse surfaces, last bounce was not through a delta BSDF)
        ESurfacePhotons,
    /// Caustic photons (indirect on diffuse surfaces, last bounce was through a delta BSDF)
        ECausticPhotons,
    /// Surface photons (all of them, even direct illumination)
        EAllSurfacePhotons,
    /// Volumetric photons
        EVolumePhotons
  };

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
  GatherPhotonCustomProcess(EGatherType type, size_t photonCount,
                            size_t granularity, int maxDepth, int rrDepth, bool isLocal,
                            bool autoCancel, const void *progressReporterPayload,
                            bool adjointCompensation, Float m_cameraSphere, int minDepth);

  /**
   * Once the process has finished, this returns a reference
   * to the (still unbalanced) photon map
   */
  inline PhotonMap *getPhotonMap() { return m_photonMap; }

  /**
   * \brief Return the number of discarded photons
   *
   * Due to asynchronous processing, some excess photons
   * will generally be produced. This function returns the number
   * of excess photons that had to be discarded. If this is too
   * high, the granularity should be decreased.
   */
  inline size_t getExcessPhotons() const { return m_excess; }

  /**
   * \brief Lists the nuber of particles that had to be shot
   * in order to fill the photon map.
   */
  inline size_t getShotParticles() const { return m_numShot; }

  inline size_t nbSkip() const { return m_nbSkip; }

  // ======================================================================
  /// @{ \name ParallelProcess implementation
  // ======================================================================

  bool isLocal() const;
  ref<WorkProcessor> createWorkProcessor() const;
  void processResult(const WorkResult *wr, bool cancelled);
  EStatus generateWork(WorkUnit *unit, int worker);

  /// @}
  // ======================================================================

  MTS_DECLARE_CLASS()
protected:
  /// Virtual destructor
  virtual ~GatherPhotonCustomProcess() {}

  /**
   * \brief Checks if the configuration of needed, generated and shot
   * photons indicates an unsuccessful progress of the gathering. This
   * check is taken from PBRT.
   */
  inline bool unsuccessful(size_t needed, size_t gen, size_t shot) {
    return (gen < needed && (gen == 0 || gen < shot / 1024));
  }
protected:
  EGatherType m_type;
  ref<PhotonMap> m_photonMap;
  size_t m_photonCount;
  int m_maxDepth, m_minDepth;
  int m_rrDepth;
  bool m_isLocal;
  bool m_autoCancel;
  size_t m_excess, m_numShot;
  Float m_cameraSphere;
  size_t m_nbSkip;
};

MTS_NAMESPACE_END

#endif //MITSUBA_BRE_PROC_H
