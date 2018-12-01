#include "photon_proc.h"

// For the statistics
#include <mitsuba/core/statistics.h>

// For isIntersectedPoint
#include "../volume_utils.h"

MTS_NAMESPACE_BEGIN

/**
 * \brief This work result implementation stores a sequence of photons, which can be
 * sent over the wire as needed.
 *
 * It is used to implement parallel networked photon tracing passes.
 */
class PhotonVector : public WorkResult {
public:
  PhotonVector() {
    m_nbSkip = 0;
  }

  inline void nextParticle() {
    m_particleIndices.push_back((uint32_t) m_photons.size());
  }

  inline void put(const Photon &p) {
    m_photons.push_back(p);
  }

  inline void incSkip() {
    m_nbSkip += 1;
  }

  inline size_t nbSkip() const {
    return m_nbSkip;
  }

  inline size_t size() const {
    return m_photons.size();
  }

  inline size_t getParticleCount() const {
    return m_particleIndices.size() - 1;
  }

  inline size_t getParticleIndex(size_t idx) const {
    return m_particleIndices.at(idx);
  }

  inline void clear() {
    m_photons.clear();
    m_particleIndices.clear();
    m_nbSkip = 0;
  }

  inline const Photon &operator[](size_t index) const {
    return m_photons[index];
  }

  void load(Stream *stream) {
    clear();
    size_t count = (size_t) stream->readUInt();
    m_particleIndices.resize(count);
    stream->readUIntArray(&m_particleIndices[0], count);
    count = (size_t) stream->readUInt();
    m_photons.resize(count);
    for (size_t i = 0; i < count; ++i)
      m_photons[i] = Photon(stream);
  }

  void save(Stream *stream) const {
    stream->writeUInt((uint32_t) m_particleIndices.size());
    stream->writeUIntArray(&m_particleIndices[0], m_particleIndices.size());
    stream->writeUInt((uint32_t) m_photons.size());
    for (size_t i = 0; i < m_photons.size(); ++i)
      m_photons[i].serialize(stream);
  }

  std::string toString() const {
    std::ostringstream oss;
    oss << "PhotonVector[size=" << m_photons.size() << "]";
    return oss.str();
  }

  MTS_DECLARE_CLASS()
protected:
  // Virtual destructor
  virtual ~PhotonVector() {}
private:
  std::vector<Photon> m_photons;
  std::vector<uint32_t> m_particleIndices;
  size_t m_nbSkip;
};

/**
 * This class does the actual photon tracing work
 */
class GatherPhotonBREWorker : public ParticleTracer {
public:
  GatherPhotonBREWorker(GatherPhotonCustomProcess::EGatherType type, size_t granularity,
                        int maxDepth, int rrDepth, Float cameraSphere, int minDepth) :
      ParticleTracer(maxDepth, rrDepth, false),
      m_type(type), m_granularity(granularity), m_cameraSphere(cameraSphere), m_minDepth(minDepth) {}

  GatherPhotonBREWorker(Stream *stream, InstanceManager *manager)
      : ParticleTracer(stream, manager) {
    m_type = (GatherPhotonCustomProcess::EGatherType) stream->readInt();
    m_granularity = stream->readSize();
  }

  ref<WorkProcessor> clone() const {
    return new GatherPhotonBREWorker(m_type, m_granularity, m_maxDepth,
                                     m_rrDepth, m_cameraSphere, m_minDepth);
  }

  void serialize(Stream *stream, InstanceManager *manager) const {
    ParticleTracer::serialize(stream, manager);
    stream->writeInt(m_type);
    stream->writeSize(m_granularity);
  }

  ref<WorkResult> createWorkResult() const {
    return new PhotonVector();
  }

  void process(const WorkUnit *workUnit, WorkResult *workResult,
               const bool &stop) {
    m_workResult = static_cast<PhotonVector *>(workResult);
    m_workResult->clear();
    m_sensorPos = m_scene->getSensor()->getWorldTransform()->operator()(0.f, Point(0.f, 0.f, 0.f));
    ParticleTracer::process(workUnit, workResult, stop);
    m_workResult->nextParticle();
    m_workResult = NULL;
  }

  void handleNewParticle() {
    m_workResult->nextParticle();
  }

  void handleCastNewParticule(const RayDifferential &ray,
                              const Medium *medium,
                              const Spectrum &power,
                              const Emitter *emitter) override {
    m_previousPos = ray.o;
  }

  void handleSurfaceInteraction(int depth_, int nullInteractions, bool delta,
                                const Intersection &its, const Medium *medium,
                                const Spectrum &weight) {
    int bsdfType = its.getBSDF()->getType(), depth = depth_;
    if (!(bsdfType & BSDF::EDiffuseReflection) && !(bsdfType & BSDF::EGlossyReflection)) {
      m_previousPos = its.p;
      return;
    }

    if (depth < m_minDepth) {
      m_previousPos = its.p;
      return;
    }

    if ((m_type == GatherPhotonCustomProcess::ECausticPhotons && depth > 1 && delta)
        || (m_type == GatherPhotonCustomProcess::ESurfacePhotons && depth > 1 && !delta)
        || (m_type == GatherPhotonCustomProcess::EAllSurfacePhotons)) {

      if (!cameraHit(its.p)) {
        m_workResult->put(Photon(its.p, its.geoFrame.n, -its.toWorld(its.wi), weight, depth));
      }

    }
    m_previousPos = its.p;
  }

  void handleMediumInteraction(int depth_, int nullInteractions, bool delta,
                               const MediumSamplingRecord &mRec, const Medium *medium,
                               const Vector &wi, const Spectrum &weight) {
    int depth = depth_;
    if (depth < m_minDepth) {
      m_previousPos = mRec.p;
      return;
    }

    if (m_type == GatherPhotonCustomProcess::EVolumePhotons) {
      if (!cameraHit(mRec.p)) {
        m_workResult->put(Photon(mRec.p, Normal(0.0f, 0.0f, 0.0f),
                                 -wi, weight, depth));
      }
    }
    m_previousPos = mRec.p;
  }

  bool cameraHit(const Point &p) {
    if (m_cameraSphere != 0.f) {
      if (isIntersectedPoint(m_sensorPos, m_previousPos, p, m_cameraSphere)) {
        m_workResult->incSkip();
        return true;
      }
    }
    return false;
  }

  MTS_DECLARE_CLASS()
protected:
  /// Virtual destructor
  virtual ~GatherPhotonBREWorker() {}

protected:
  GatherPhotonCustomProcess::EGatherType m_type;
  size_t m_granularity;
  Float m_cameraSphere;
  ref<PhotonVector> m_workResult;
  int m_minDepth;

  // To dermine the creation or not of the photon
  Point m_previousPos;
  Point m_sensorPos;
};

GatherPhotonCustomProcess::GatherPhotonCustomProcess(EGatherType type, size_t photonCount,
                                                     size_t granularity, int maxDepth, int rrDepth, bool isLocal,
                                                     bool autoCancel,
                                                     const void *progressReporterPayload,
                                                     bool adjointCompensation, Float cameraSphere, int minDepth)
    : ParticleProcess(ParticleProcess::EGather, photonCount, granularity, "Gathering photons",
                      progressReporterPayload, adjointCompensation), m_type(type), m_photonCount(photonCount),
      m_maxDepth(maxDepth), m_minDepth(minDepth),
      m_rrDepth(rrDepth), m_isLocal(isLocal),
      m_autoCancel(autoCancel), m_excess(0), m_numShot(0),
      m_cameraSphere(cameraSphere), m_nbSkip(0) {
  m_photonMap = new PhotonMap(photonCount);
}

bool GatherPhotonCustomProcess::isLocal() const {
  return m_isLocal;
}

ref<WorkProcessor> GatherPhotonCustomProcess::createWorkProcessor() const {
  return new GatherPhotonBREWorker(m_type,
                                   m_granularity,
                                   m_maxDepth,
                                   m_rrDepth,
                                   m_cameraSphere,
                                   m_minDepth);
}

void GatherPhotonCustomProcess::processResult(const WorkResult *wr, bool cancelled) {
  if (cancelled)
    return;
  const PhotonVector &vec = *static_cast<const PhotonVector *>(wr);
  LockGuard lock(m_resultMutex);
  m_nbSkip += vec.nbSkip();

  size_t nParticles = 0;      // Each particle marks a light path (with several photons on it)
  for (size_t i = 0; i < vec.getParticleCount(); ++i) {
    size_t start = vec.getParticleIndex(i),
        end = vec.getParticleIndex(i + 1);
    ++nParticles;
    bool full = false;
    for (size_t j = start; j < end; ++j) {
      if (!m_photonMap->tryAppend(vec[j])) {
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

ParallelProcess::EStatus GatherPhotonCustomProcess::generateWork(WorkUnit *unit, int worker) {
  /* Use the same approach as PBRT for auto canceling */
  LockGuard lock(m_resultMutex);
  if (m_autoCancel && m_numShot > 100000
      && unsuccessful(m_photonCount, m_photonMap->size(), m_numShot)) {
    Log(EInfo, "Not enough photons could be collected, giving up");
    return EFailure;
  }

  return ParticleProcess::generateWork(unit, worker);
}

MTS_IMPLEMENT_CLASS(GatherPhotonCustomProcess, false, ParticleProcess)
MTS_IMPLEMENT_CLASS_S(GatherPhotonBREWorker, false, ParticleTracer)
MTS_IMPLEMENT_CLASS(PhotonVector, false, WorkResult)
MTS_NAMESPACE_END
