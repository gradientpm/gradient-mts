#include "beams.h"

// Acceleration structure

MTS_NAMESPACE_BEGIN

class PhotonBeamVector : public WorkResult {
public:
  PhotonBeamVector(int workerID) : m_workerID(workerID) {
    m_nbBeamSkip = 0;
  }

  inline void nextParticle() {
    m_particleIndices.push_back((uint32_t) m_beams.size());
  }

  inline void put(const PhotonBeam &p) {
    m_beams.push_back(p);
  }

  inline size_t size() const {
    return m_beams.size();
  }

  inline size_t getParticleCount() const {
    return m_particleIndices.size() - 1;
  }

  inline size_t getParticleIndex(size_t idx) const {
    return m_particleIndices.at(idx);
  }

  inline void clear() {
    m_beams.clear();
    m_particleIndices.clear();
    m_nbBeamSkip = 0;
  }

  inline const PhotonBeam &operator[](size_t index) const {
    return m_beams[index];
  }

  inline void incSkip() {
    m_nbBeamSkip += 1;
  }
  inline const size_t nbSkip() const {
    return m_nbBeamSkip;
  }

  void load(Stream *stream) {
    clear();
    size_t count = (size_t) stream->readUInt();
    m_particleIndices.resize(count);
    stream->readUIntArray(&m_particleIndices[0], count);
    count = (size_t) stream->readUInt();
    m_beams.resize(count);
    /*for (size_t i=0; i<count; ++i)
      m_photons[i] = Photon(stream);*/
    SLog(EError, "Impossible to load serialized photon beams results.");
  }

  void save(Stream *stream) const {
    stream->writeUInt((uint32_t) m_particleIndices.size());
    stream->writeUIntArray(&m_particleIndices[0], m_particleIndices.size());
    stream->writeUInt((uint32_t) m_beams.size());
    /*for (size_t i=0; i<m_photons.size(); ++i)
      m_photons[i].serialize(stream);*/
    SLog(EError, "Impossible to serialize photon beams for now");
  }

  std::string toString() const {
    std::ostringstream oss;
    oss << "PhotonBeamVector[size=" << m_beams.size() << "]";
    return oss.str();
  }

  MTS_DECLARE_CLASS()
protected:
  // Virtual destructor
  virtual ~PhotonBeamVector() {}
private:
  std::vector<PhotonBeam> m_beams;
  std::vector<uint32_t> m_particleIndices;
  size_t m_nbBeamSkip;
public:
  int m_workerID;
};

inline RayDifferential updateToCylinder(RayDifferential &ray, const Point &newP, const Vector &newD) {
  Float dist = distance(newP, ray.o);
  Point xP = ray.rxOrigin + ray.rxDirection * dist;
  Point yP = ray.ryOrigin + ray.ryDirection * dist;

  Float xRadius = distance(newP, xP);
  Float yRadius = distance(newP, yP);

  // Update the ray differential
  FrameCoherent frame(ray.d);
  ray = RayDifferential(newP, newD, ray.time);
  ray.rxOrigin = newP + xRadius * frame.s;
  ray.ryOrigin = newP + yRadius * frame.t;
  ray.rxDirection = ray.d; // Cylinder
  ray.ryDirection = ray.d; // Cylinder

  return ray;
}

class BeamPhotonWorker : public ParticleTracer {
public:

  BeamPhotonWorker(int workerID,
                   size_t beamCount,
                   size_t granularity,
                   int maxDepth,
                   int rrDepth,
                   Float radius, Float cameraSphere,
                   int minDepth) :
      ParticleTracer(maxDepth, rrDepth, false),
      haveBeam(false),
      m_beamCount(beamCount),
      m_granularity(granularity), m_beamRadius(radius),
      m_cameraSphere(cameraSphere),
      m_minDepth(minDepth),
      m_workerID(workerID) {}

  BeamPhotonWorker(Stream *stream, InstanceManager *manager)
      : ParticleTracer(stream, manager) {
    m_granularity = stream->readSize();
    m_beamCount = stream->readSize();
    SLog(EError, "No serialization is implemented for photon mapping");
  }

  ref<WorkProcessor> clone() const {
    SLog(EError, "Cloning of beam worker not allowed.");

    return new BeamPhotonWorker(m_workerID,
                                m_beamCount,
                                m_granularity,
                                m_maxDepth,
                                m_rrDepth,
                                m_beamRadius,
                                m_cameraSphere,
                                m_minDepth);
  }

  void serialize(Stream *stream, InstanceManager *manager) const {
    ParticleTracer::serialize(stream, manager);
    stream->writeSize(m_granularity);
    stream->writeSize(m_beamCount);
    SLog(EError, "No serialization is implemented for photon mapping");
  }

  ref<WorkResult> createWorkResult() const {
    return new PhotonBeamVector(m_workerID);
  }

  void process(const WorkUnit *workUnit, WorkResult *workResult,
               const bool &stop) {
    // Get the work results and clear it
    m_workResult = static_cast<PhotonBeamVector *>(workResult);
    m_workResult->clear();
    haveBeam = false;

    // Get the sensor position
    m_sensorPos = m_scene->getSensor()->getWorldTransform()->operator()(0.f, Point(0.f, 0.f, 0.f));

    // Run the shooting
    ParticleTracer::process(workUnit, workResult, stop);

    // Finish and free the smart pointer
    m_workResult->nextParticle();
    m_workResult = NULL;
  }

  void handleCastNewParticule(const RayDifferential &ray,
                              const Medium *medium,
                              const Spectrum &power,
                              const Emitter *emitter) override {
    // Check if there is a medium
    m_workResult->nextParticle();
    haveBeam = false;

    if (m_minDepth > 0) {
      return; // Do not create the beam
    }

    if (medium != NULL) {
      // If the luminaire is inside a medium,
      // Start to create a photon beams
      currBeam = PhotonBeam(ray.o, medium, power, 1, m_beamRadius);
      currBeam.setDifferential(ray);
      haveBeam = true;
    }
  }

  void handleSurfaceInteraction(int depth, int nullInteractions, bool delta,
                                const Intersection &its, const Medium *medium,
                                const Spectrum &weight) {
    // We was creating a beam, just create it
    if (haveBeam && currBeam.isInvalid()) {
      // In this case, create the beam
      currBeam.setEndPoint(its.p);
      if (!hitCamera()) {
        m_workResult->put(currBeam);
      }
      haveBeam = false;
    }
  }

  void handleMediumInteraction(int depth, int nullInteractions, bool delta,
                               const MediumSamplingRecord &mRec, const Medium *medium,
                               const Vector &wi, const Spectrum &weight) {
    if (medium != 0) {
      // We was creating a beam, just create it
      if (haveBeam && currBeam.isInvalid()) {
        // In this case, create the beam
        currBeam.setEndPoint(mRec.p);
        if (!hitCamera()) {
          m_workResult->put(currBeam);
        }
        haveBeam = false;
      }
    }
  }

  void handledBounce(int depth_, int nullInteractions,
                     bool delta, const Point &p, const Medium *medium, bool wasNullInteraction,
                     const Spectrum &weight) {
    SAssert(haveBeam == false);

    int depth = depth_;
    if (depth < m_minDepth) {
      return; // Do not create the beam
    }

    // Do not create beam if there is no valid weight
    if (medium != 0 && !weight.isZero()) {
      // If there is a medium, create a beam anyway
      currBeam = PhotonBeam(p, medium, weight, depth, m_beamRadius);
      haveBeam = true;
    }
  }

  void handleSetRayDifferentialFromEmitter(RayDifferential &ray,
                                           const PositionSamplingRecord &pRec,
                                           const DirectionSamplingRecord &dRec) {
    // Only handle non discrete case
    // The only valid setup is area light for now
    ray.hasDifferentials = false;
    if (pRec.measure == EDiscrete || dRec.measure == EDiscrete)
      return;


    // dA = r^2*M_PI = 1 / (N * pdf(x)) -> bRadius = r
    Float bRadius = std::sqrt(1.f / (m_beamCount * pRec.pdf * M_PI));
    // Use the pdf as it mention in Sec 8.1 "A Comprehensive Theory of Volumetric Radiance ..."
    Float beamSolidAngle = 1.f / (m_beamCount * dRec.pdf);
    // Project angle in 2D
    Float angle2D = std::asin(std::sqrt(beamSolidAngle / M_PI));
    // Compute the apex (Fig. 9)
    Float dist = bRadius / std::tan(angle2D);
    Point apexVertex = ray(-dist);

    // Construct a frame to fill ray differentials
    FrameCoherent frame(ray.d);
    ray.rxOrigin = ray.o + bRadius * frame.s;
    ray.ryOrigin = ray.o + bRadius * frame.t;
    ray.rxDirection = normalize(ray.rxOrigin - apexVertex);
    ray.ryDirection = normalize(ray.ryOrigin - apexVertex);
    ray.hasDifferentials = true;

    // Scale up the ray differential if needed
    ray.scaleDifferential(m_beamRadius);
  }

  void handleMediumInteractionScattering(const MediumSamplingRecord &mRec, const PhaseFunctionSamplingRecord &pRec,
                                         RayDifferential &ray) override {
    ray = RayDifferential(mRec.p, pRec.wo, ray.time);
    ray.hasDifferentials = false;
    ray.mint = 0;
  }

  void handleSurfaceInteractionScattering(const BSDFSamplingRecord &bRec,
                                          RayDifferential &ray) override {
    ray.setOrigin(bRec.its.p);
    ray.setDirection(bRec.its.toWorld(bRec.wo));
    ray.hasDifferentials = false;
    ray.mint = Epsilon; // to avoid self intersection
  }

  void handleSetRayDifferential(const RayDifferential &ray) {
    currBeam.setDifferential(ray);
  }

  void handleFinishParticule() {
    if (haveBeam && currBeam.isInvalid()) {
      // Just ignore it,
      SLog(EError, "Impossible if there is a limit to the participating media: %s",
           currBeam.toString().c_str());
    }
  }

  bool hitCamera() {
    if (m_cameraSphere != 0.f) {
      if (haveBeam && currBeam.isInvalid()) {
        if (isIntersectedPoint(m_sensorPos, currBeam.getOri(), currBeam.getEndPoint(), m_cameraSphere)) {
          m_workResult->incSkip();
          return true;
        }
      }
    }
    return false;
  }

  MTS_DECLARE_CLASS()
protected:
  /// Virtual destructor
  virtual ~BeamPhotonWorker() {}
protected:
  bool haveBeam;
  PhotonBeam currBeam;
  size_t m_granularity = 0;
  size_t m_beamCount = 0;
  Float m_beamRadius;
  ref<PhotonBeamVector> m_workResult;
  Float m_cameraSphere;
  int m_minDepth;

  Point m_sensorPos;
  int m_workerID;
};

/*
* SubBeam kd tree
*/
PhotonSubBeam::PhotonSubBeam(const Point &pos, const PhotonBeam *beam, Float t1, Float t2) {
  position = pos;
  data.beam = beam;
  data.t1 = t1;
  data.t2 = t2;
  flags = 0;
}

/*
* Parallel process
*/
GatherPhotonBeamProcess::GatherPhotonBeamProcess(size_t beamCount, size_t granularity, int maxDepth, int rrDepth,
                                                 Float beamRadius,
                                                 bool isLocal, bool autoCancel,
                                                 const void *progressReporterPayload,
                                                 bool adjointCompensation,
                                                 Float cameraSphere, int minDepth)
    : ParticleProcess(ParticleProcess::EGather, beamCount, granularity, "Shooting Photon Beams",
                      progressReporterPayload, adjointCompensation),
      m_beamCount(beamCount), m_maxDepth(maxDepth), m_minDepth(minDepth),
      m_rrDepth(rrDepth), m_isLocal(isLocal), m_autoCancel(autoCancel),
      m_excess(0), m_numShot(0),
      m_beamRadius(beamRadius), m_cameraSphere(cameraSphere), m_nbSkip(0), m_nextWorkerID(0) {

  m_beamMap = new BeamMap<PhotonBeam>(beamCount);
}

bool GatherPhotonBeamProcess::isLocal() const {
  return m_isLocal;
}

ref<WorkProcessor> GatherPhotonBeamProcess::createWorkProcessor() const {
  // It is sort of a hack to use mutable keyword for m_nextWorkerID here since the function is const.
  std::cout << "Spawn SPPM beam work processor. ID " << m_nextWorkerID << std::endl;
  return new BeamPhotonWorker(m_nextWorkerID++, m_beamCount, m_granularity,
                              m_maxDepth, m_rrDepth,
                              m_beamRadius, m_cameraSphere, m_minDepth);
}

void GatherPhotonBeamProcess::processResult(const WorkResult *wr, bool cancelled) {
  if (cancelled)
    return;
  const PhotonBeamVector &vec = *static_cast<const PhotonBeamVector *>(wr);
  LockGuard lock(m_resultMutex);
  m_nbSkip += vec.nbSkip();

  size_t nParticles = 0;      // Each particle marks a light path (with several photons on it)
  for (size_t i = 0; i < vec.getParticleCount(); ++i) {
    size_t start = vec.getParticleIndex(i),
        end = vec.getParticleIndex(i + 1);
    ++nParticles;
    bool full = false;
    for (size_t j = start; j < end; ++j) {
      if (vec[j].isInvalid()) {
        SLog(EError, "Try to push an invalid beam");
      }
      if (!m_beamMap->tryAppend(vec.m_workerID, vec[j])) {
        m_excess += vec.size() - j;
        full = true;
        break;
      }
    }
    if (full)
      break;
  }
  m_numShot += nParticles;
  increaseResultCount(vec.size());
}

ParallelProcess::EStatus GatherPhotonBeamProcess::generateWork(WorkUnit *unit, int worker) {
  /* Use the same approach as PBRT for auto canceling */
  LockGuard lock(m_resultMutex);
  if (m_autoCancel && m_numShot > 100000
      && unsuccessful(m_beamCount, m_beamMap->size(), m_numShot)) {
    Log(EInfo, "Not enough photons could be collected, giving up");
    return EFailure;
  }

  return ParticleProcess::generateWork(unit, worker);
}

// Worker/Parallel structure for beams shooting
MTS_IMPLEMENT_CLASS(GatherPhotonBeamProcess, false, ParticleProcess)
MTS_IMPLEMENT_CLASS_S(BeamPhotonWorker, false, ParticleTracer)
MTS_IMPLEMENT_CLASS(PhotonBeamVector, false, WorkResult)

MTS_NAMESPACE_END
