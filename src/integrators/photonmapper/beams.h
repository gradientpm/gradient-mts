#if !defined(__BEAMS_H)
#define __BEAMS_H

#include <mitsuba/render/photonmap.h>
#include <mitsuba/render/gatherproc.h>

#include <fstream>

#include "beams_struct.h"
#include "beams_accel.h"
#include "../volume_utils.h"

#include "beams_3d_intersections.h"

MTS_NAMESPACE_BEGIN

template<typename T>
class BeamRadianceQuery {
public:
  BeamRadianceQuery(const Ray &ray_,
                    const Medium *medium_,
                    int maxDepth_,
                    int minDepth_,
                    ref<Sampler> sampler_,
                    EVolumeTechnique volTechnique_) :
      baseCameraRay(ray_), medium(medium_),
      maxDepth(maxDepth_), minDepth(minDepth_), Li(0.f), sampler(sampler_), volTechnique(volTechnique_) {}

  bool operator()(const T *beam, Float tmin = 0, Float tmax = std::numeric_limits<Float>::infinity()) {
    if (tmax > beam->getLength()) {
      tmax = beam->getLength();
    }

    if (maxDepth != -1 && beam->depth > maxDepth) {
      return false; // Skip this path
    }
    if (minDepth != 0 && beam->depth < minDepth) {
      return false; // Skip also this path
    }

    if (volTechnique == EBeamBeam1D) {
      Float u, v, w, sinTheta;
      if (beam->rayIntersect1D(baseCameraRay,
                               tmin, tmax,
                               u, v, w, sinTheta)) {
        // if beam differential is enable and beam has differential then we test intersection against it's varying radius
        // otherwise we use the global 'radius' constant that is also the same for all beams
        // if differential radius is smaller than shortest distance to beam, there is no intersection.
        if (beam->getRadius() <= u) {
          return false;
        }

        // Evaluate the medium
        Ray rayBeam(baseCameraRay);
        rayBeam.mint = Epsilon;
        rayBeam.maxt = w;
        MediumSamplingRecord mRec;
        medium->eval(rayBeam, mRec);

        Float radius = beam->getRadius();
        Float weightKernel = 0.5f / radius;
        Li += beam->getContrib(v, mRec, rayBeam.d) * weightKernel / sinTheta;
        return true;
      } else {
        return false;
      }
    } else if (volTechnique == EBeamBeam3D_Naive || volTechnique == EBeamBeam3D_EGSR
        || volTechnique == EBeamBeam3D_Optimized) {
      Float radius = beam->getRadius();

      // Optimized version of 3D beam-beam kernel.
      if (beam->isInvalid())
        return false;

      Float beamSegmentRand = -1.f, cameraSegmentRand = -1.f, invPDF = 0.f;

      if (volTechnique == EBeamBeam3D_Naive) {
        // When a camera ray intersects a sub-beam AABB, this query will be executed
        // Inside here we have to perform radius check


        // Sample a point on the subbeam [tb-, tb+] to put the kernel (uniform sampling)
        // This effectively slides the kernel along the subbeam (over many Monte Carlo estimations)
        beamSegmentRand = tmin + (tmax - tmin) * sampler->next1D();
        invPDF = tmax - tmin;

        // Put the kernel at the sampled point
        Point kernelCentroid = beam->getPos(beamSegmentRand);
//              Float kernelVol = ((4.0 / 3.0) * M_PI * std::pow(radius, 3));

        // Find the projection of the photon position to the camera ray
        Float distToProj = dot(kernelCentroid - baseCameraRay.o, baseCameraRay.d);
        Float distSqr = (baseCameraRay(distToProj) - kernelCentroid).lengthSquared();
        Float radSqr = radius * radius;
        if (distSqr >= radSqr)
          return false;    // ray does not intersect the kernel

        // Sample camera segment [tc-, tc+] (uniform sampling)
        Float deltaT = math::safe_sqrt(radSqr - distSqr); // equal to half of [tc-, tc+] range
        cameraSegmentRand = (distToProj - deltaT)
            + 2 * deltaT * sampler->next1D();     // the random point is always inside the kernel
        invPDF *= std::max(2.0 * deltaT, 0.0001);
      } else {
        // Pre-compute ray form for camera and beam
        // in order to call the intersection procedure
        Ray _cam = Ray(baseCameraRay(baseCameraRay.mint), baseCameraRay.d,
                       0, baseCameraRay.maxt - baseCameraRay.mint, 0.f);

        // In this sampling strategy, we consider the beam entirely.
        // However, as a beam can be represented by several subbeam
        // We discard the computation of non valid subbeams (to select only one subbeam at the end)
        Ray _beam = Ray(beam->getOri(), beam->getDir(),
                        0.f, beam->getLength(), 0.f); // TODO: MAke a function

        // 1) Compute tNear and tFar of the beam to camera "beam" (Cylinder intersection)
        double tNearBeam, tFarBeam;
        if (!cylinderIntersection(_cam, _beam, radius, tNearBeam, tFarBeam)) {
          return false; // No intersection
        }

        // 2) Check if it is the correct sub-beam (to avoid to integrate twice)
        if (tNearBeam < 0 && tmin <= Epsilon) {
          // Acceptable, it is the first subbeam
        } else if (tNearBeam > tmin && tNearBeam < tmax) {
          // The tNearBeam need to be inside the subbeam
        } else {
          return false; // Another subbeam will take care of it
        }

        // 3) Sample a point on the photon beam
        beamSegmentRand = tNearBeam + (tFarBeam - tNearBeam) * sampler->next1D();
        invPDF = std::max(tFarBeam - tNearBeam, 0.0001);
        if (beamSegmentRand < 0 || beamSegmentRand > beam->getLength()) {
          return false;
        }

        // 4) Determine the camera position
        if (volTechnique == EBeamBeam3D_EGSR) {
          // Compute tNear and tFar over the camera (Cylinder intersection)
          double tNearCam, tFarCam;
          if (!cylinderIntersection(_beam, _cam, radius, tNearCam, tFarCam)) {
            return false; // No intersection
          }

          // To double check
          SAssert(tFarBeam >= tNearBeam);
          SAssert(tFarCam >= tNearCam);

          cameraSegmentRand = tNearCam + (tFarCam - tNearCam) * sampler->next1D();
          invPDF *= std::max(tFarCam - tNearCam, 0.0001);
        } else {
          // Optimized: sample the point over the camera beam
          // on the overlapping region
          Point kernelCentroid = beam->getPos(beamSegmentRand);

          // Find the projection of the photon position to the camera ray
          Float distToProj = dot(kernelCentroid - baseCameraRay.o, baseCameraRay.d);
          Float distSqr = (baseCameraRay(distToProj) - kernelCentroid).lengthSquared();
          Float radSqr = radius * radius;
          if (distSqr >= radSqr) {
            return false;  // ray does not intersect the kernel
          }

          // Sample camera segment [tc-, tc+] (uniform sampling)
          Float deltaT = math::safe_sqrt(radSqr - distSqr);    // equal to half of [tc-, tc+] range
          cameraSegmentRand = distToProj - deltaT
              + 2 * deltaT * sampler->next1D(); // the random point is always inside the kernel
          invPDF *= std::max(2.0 * deltaT, 0.0001);
        }

        // Check if it is in a acceptable distance for the camera
        if (cameraSegmentRand < baseCameraRay.mint || cameraSegmentRand > baseCameraRay.maxt) {
          return false;
        }

        // Check that the photon is inside the kernel
        if (volTechnique == EBeamBeam3D_EGSR) {
          // "Optimized" generate always a photon inside the kernel
          // so we don't need to test it.
          Point kernelCentroid = beam->getPos(beamSegmentRand);
          Float distSqr = (baseCameraRay(cameraSegmentRand) - kernelCentroid).lengthSquared();
          Float radSqr = radius * radius;
          if (distSqr >= radSqr) {
            return false;  // ray does not intersect the kernel
          }
        }
      }

      // Compute final contribution
      // a) Transmittance value over the photon beam
      Ray rayTrans(beam->getOri(), beam->getDir(), 0.f); // TODO: Make a function
      rayTrans.mint = 0.f;
      rayTrans.maxt = beamSegmentRand;
      MediumSamplingRecord mRecBeam;
      medium->eval(rayTrans, mRecBeam);
      Spectrum segmentFlux = beam->flux * mRecBeam.transmittance;
      // b) Transmittance value over the camera beam
      Ray cameraRay(baseCameraRay);
      cameraRay.mint = Epsilon;
      cameraRay.maxt = cameraSegmentRand;
      MediumSamplingRecord mRecCamera;
      medium->eval(cameraRay, mRecCamera);
      // c) Phase function value
      PhaseFunctionSamplingRecord pRec(mRecCamera, -beam->getDir(), -cameraRay.d, EImportance);
      Float phaseTerm = mRecCamera.getPhaseFunction()->eval(pRec);

      // Total contribution
      Float kernelVol = ((4.0 / 3.0) * M_PI * std::pow(radius, 3));
      Spectrum beamContrib = segmentFlux * mRecCamera.sigmaS * mRecCamera.transmittance * phaseTerm
          * (invPDF / kernelVol);
      if (!beam->longBeams) {
        // need to incorporate the russian roulette
        // basically when the end point is not sampled by medium's sampleDistance function, use pdfFailure.
        beamContrib /= mRecBeam.pdfFailure;
      }
      // Add the contribution
      Li += beamContrib;
      return true;
    }

    SLog(EError, "Not supported kernel type");
    return false;
  }

  // The initial ray
  const Ray &baseCameraRay;
  const Medium *medium;
  int maxDepth;
  int minDepth;

  // Accumulate spectrum value
  Spectrum Li;

  ref<Sampler> sampler;
  EVolumeTechnique volTechnique;
};

template<typename T>
class BeamMap : public SerializableObject {
public:
  enum EBeamAccelStructure {
    ENoAccel = 1,
    EBVHAccel = 2
  };

  /**
   * \brief Create a Beam acceleration data structure
   */
  BeamMap(std::size_t beamCount) : m_beamBVH(0) {
    // Pre-allocate the memory for the beams
    m_beams.reserve(beamCount + 100);
  }

  /**
   * \brief Unserialize a BRE acceleration data structure from
   * a binary data stream
   */
  BeamMap(Stream *stream, InstanceManager *manager) {
    SLog(EError, "Beam map doesn't support serilisation yet");
  }

  /**
     * \brief Try to append a photon beam into the photon beam map
     *
     * \return \c false If the photon beam map is full
     */
  inline bool tryAppend(int idWorker, const T &beam) {
    if (m_beams.size() < m_beams.capacity()) {
      m_beams.push_back(std::pair<int, T>(idWorker, beam));
      return true;
    } else {
      return false;
    }
  }

  inline size_t size() const { return m_beams.size(); }

  inline bool outCapacity() const { return m_beams.size() >= m_beams.capacity(); }

  /// Serialize to a binary data stream
  void serialize(Stream *stream, InstanceManager *manager) const {
    SLog(EError, "Serialisation is not implemented");
  }

  /// Compute the beam radiance estimate for the given ray segment and medium
  template<typename Q>
  void query(Q &bRadQuery) const {
    if (m_accelType == ENoAccel) {
      // For now, just looping over all the photon beams
      for (size_t i = 0; i < m_beams.size(); i++) {
        // Use the ray min and max distance because we want to intersection along all the view ray
        bRadQuery(&m_beams[i].second);
      }
    } else if (m_accelType == EBVHAccel) {
      m_beamBVH->query(bRadQuery);
    } else {
      SLog(EError, "Unknown accel type building strategy");
    }
  }

  void build(EBeamAccelStructure accelType = ENoAccel, bool deterministic = false) {
    m_accelType = accelType;
    if (m_accelType == ENoAccel) {
      // Nothing to do
    } else if (m_accelType == EBVHAccel) {
      // In theory, BVH should not be affected by the order of the beams.
      // But in practice, to ensure deterministic behavior, beside ensuring threads taking up the same jobs,
      // we have to sort the beams based on worker ID.
      // This means there could have an issue in BVH construction or intersection test.
      if (deterministic) {
        typedef std::pair<int, T> BeamPair;
        std::stable_sort(m_beams.begin(), m_beams.end(),
                         [](const BeamPair &lhs, const BeamPair &rhs) {
                           return lhs.first < rhs.first;
                         });
      }

      m_beamBVH = new SubBeamBVH<T>(m_beams);
    } else {
      SLog(EError, "Unknown accel type building strategy");
    }
  }

  void convertToLongBeams(const Scene *scene) {
    // Normally we generate short beams, but by extend the beam (to the surface intersection)
    // We can change short beams to long one
    for (size_t i = 0; i < m_beams.size(); i++) {
      m_beams[i].second.convertToLong(scene);
    }
  }

  const std::vector<std::pair<int, T> > &getBeams() const {
    return m_beams;
  }

  /// Release all memory
  virtual ~BeamMap() {
    if (m_beamBVH != 0)
      delete m_beamBVH;
  }

protected:
  std::vector<std::pair<int, T> > m_beams;
  EBeamAccelStructure m_accelType;
  SubBeamBVH<T> *m_beamBVH;
};

class GatherPhotonBeamProcess : public ParticleProcess {
public:
  /**
   * Create a new process for parallel photon gathering
   * \param beamCount
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
  GatherPhotonBeamProcess(size_t beamCount, size_t granularity, int maxDepth, int rrDepth, Float beamRadius,
                          bool isLocal, bool autoCancel, const void *progressReporterPayload,
                          bool adjointCompensation, Float cameraSphere, int minDepth);

  inline BeamMap<PhotonBeam> *getBeamMap() { return m_beamMap; }

  /**
   * \brief Return the number of discarded photons
   *
   * Due to asynchronous processing, some excess photons
   * will generally be produced. This function returns the number
   * of excess photons that had to be discarded. If this is too
   * high, the granularity should be decreased.
   */
  inline size_t getExcessBeams() const { return m_excess; }

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
  virtual ~GatherPhotonBeamProcess() {}

  /**
   * \brief Checks if the configuration of needed, generated and shot
   * photons indicates an unsuccessful progress of the gathering. This
   * check is taken from PBRT.
   */
  inline bool unsuccessful(size_t needed, size_t gen, size_t shot) {
    return (gen < needed && (gen == 0 || gen < shot / 1024));
  }

protected:
  BeamMap<PhotonBeam> *m_beamMap;

  size_t m_beamCount;
  int m_maxDepth, m_minDepth;
  int m_rrDepth;
  bool m_isLocal;
  bool m_autoCancel;
  size_t m_excess, m_numShot;
  Float m_beamRadius;   // beam radius value is really used only if m_useBeamDifferential is False
  Float m_cameraSphere;
  size_t m_nbSkip;

  mutable int m_nextWorkerID;   // for keeping track of which beam is spawn by which worker if parallelization is used
};

// OR/Opt: See the line kd-tree (but double check the bias)

MTS_NAMESPACE_END

#endif /* __BEAMS_H */