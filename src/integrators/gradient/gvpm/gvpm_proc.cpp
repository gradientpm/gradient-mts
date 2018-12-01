#include <mitsuba/bidir/path.h>
#include <mitsuba/render/range.h>

#include "gvpm_proc.h"

MTS_NAMESPACE_BEGIN

/**
 * \brief This work result implementation stores a sequence of photons, which can be
 * sent over the wire as needed.
 *
 * It is used to implement parallel networked photon tracing passes.
 */
class PhotonPathVector : public WorkResult {
public:
  PhotonPathVector() {
    idWorker = -1;
  }

  inline void addPath(Path *p) {
    m_lightPaths.push_back(p);
  }
  inline void addVoidPath() {
    m_lightPaths.push_back(nullptr);
  }

  inline void clear() {
    m_lightPaths.clear();
  }

  inline const size_t getLightPathCount() const {
    return m_lightPaths.size();
  }

  inline const Path *operator[](size_t index) const {
    return m_lightPaths[index];
  }

  void load(Stream *stream) override {
    SLog(EError, "Impossible to stream GPM");
  }

  void save(Stream *stream) const override {
    SLog(EError, "Impossible to stream GPM");
  }

  std::string toString() const override {
    std::ostringstream oss;
    oss << "PhotonVector[size=" << m_lightPaths.size() << "]";
    return oss.str();
  }

  int idWorker;

  MTS_DECLARE_CLASS()
protected:
  // Virtual destructor
  virtual ~PhotonPathVector() = default;
private:
  std::vector<Path *> m_lightPaths;
};

/**
 * This class does the actual photon tracing work
 */
class GradientPhotonWorker : public WorkProcessor {
public:
  GradientPhotonWorker(const GPMConfig &config,
                       int workerID, std::vector<ShootingThreadData> &thdata,
                       Float beamRadius) :
      m_config(config), m_workerID(workerID),
      m_myThdata(&thdata[workerID]), m_beamRadius(beamRadius) {}

  GradientPhotonWorker(Stream *stream, InstanceManager *manager)
      : WorkProcessor(stream, manager), m_config(GPMConfig()) {
    SLog(EError, "No serialization for GPM");
  }

  ref<WorkProcessor> clone() const override {
    // This is not possible to ensure the uniqueness
    // of the worker ID
    Log(EError, "No suport of worker cloning ... ");
    return nullptr;
  }

  void serialize(Stream *stream, InstanceManager *manager) const override {
    SLog(EError, "No serialization GPM");
  }

  ref<WorkResult> createWorkResult() const override {
    return new PhotonPathVector();
  }

  ref<WorkUnit> createWorkUnit() const override {
    return new RangeWorkUnit();
  }

  void prepare() override {
    Scene *scene = static_cast<Scene *>(getResource("scene"));
    m_scene = new Scene(scene);
    Sensor *newSensor = static_cast<Sensor *>(getResource("sensor"));
    m_scene->removeSensor(scene->getSensor());
    m_scene->addSensor(newSensor);
    m_scene->setSensor(newSensor);
    m_scene->initializeBidirectional();
  }

  bool haveVolumeInteraction(bool beamGeneration, const Path *lPath) const {
    if (beamGeneration) {
      for (size_t idEdge = std::max(m_config.minDepth, 1);
           idEdge < lPath->edgeCount(); idEdge++) {
        if (lPath->edge(idEdge)->medium != nullptr)
          return true;
      }
    } else {
      for (size_t idVertex = std::max(m_config.minDepth, 1) + 1;
           idVertex < lPath->vertexCount(); idVertex++) {
        if (lPath->vertex(idVertex)->isMediumInteraction())
          return true;
      }
    }
    return false;
  }

  bool generatePath(bool beamGeneration, bool &needCreateNew, Path **lPath) {
    // Generate light path
    if (needCreateNew) {
      (*lPath) = m_myThdata->getPath();
    }
    (*lPath)->initialize(m_scene, 0.f, EImportance, m_myThdata->pool);
    (*lPath)->randomWalk(m_scene, m_myThdata->sampler,
                         m_config.maxDepth,
                         m_config.rrDepth,
                         EImportance,
                         m_myThdata->pool,
                         m_config.adjointCompensation,
                         true);
    // Check the PDF
    for (size_t i = 1; i < (*lPath)->vertexCount() - 1; i++) {
      if ((*lPath)->vertex(i)->pdf[EImportance] == 0.f) {
        SLog(EWarn, "Importance pdf is 0 for the vertex id: %i", i);
        return false; // Reject this type of paths
      }
    }

    return haveVolumeInteraction(beamGeneration, *lPath);
  }

  void process(const WorkUnit *workUnit, WorkResult *workResult,
               const bool &stop) override {
    // Get the work task
    const RangeWorkUnit *range = static_cast<const RangeWorkUnit *>(workUnit);
    bool beamGeneration = EVolumeTechniqueHelper::needBeams(m_config.volTechnique);

    // Get the work results and clear it from previous
    m_workResult = static_cast<PhotonPathVector *>(workResult);
    m_workResult->clear();
    m_workResult->idWorker = m_workerID; // Attach the worker ID

    bool needCreateNew = true;
    Path *lPath = nullptr;
    if (m_config.deterministic) {
      // In deterministic, make sure to generate
      size_t index = range->getRangeStart();
      while (index <= range->getRangeEnd() && !stop) {
        if (!generatePath(beamGeneration, needCreateNew, &lPath)) {
          m_workResult->addVoidPath();
          lPath->release(m_myThdata->pool);
          needCreateNew = false;
          continue;
        }

        // Advance based on the number of particule generated
        if (beamGeneration) {
          index += numberOfBeams(lPath, m_config.minDepth);
        } else {
          index += numberOfVolumePhotons(lPath, m_config.minDepth);
        }

        if (index > range->getRangeEnd()) {
          // Don't count the last path (to not reach maximum of beams)
          // This will create a small but acceptable bias
          lPath->release(m_myThdata->pool);
        } else {
          // This is a valid path, we just add it to (later) push it inside the KD-tree.
          m_workResult->addPath(lPath);
          needCreateNew = true;
        }

      }
    } else {
      // Just do the loop to how many light path generated
      for (size_t index = range->getRangeStart(); index <= range->getRangeEnd() && !stop; ++index) {
        bool haveBeam = generatePath(beamGeneration, needCreateNew, &lPath);

        if (!haveBeam && !m_config.needSurfaceRendering()) {
          // No interactions, add a empty path.
          // We do that in order to avoid new path creation and reuse the previous one.
          m_workResult->addVoidPath();
          lPath->release(m_myThdata->pool);
          needCreateNew = false;
        } else {
          // This is a valid path, we just add it to (later) push it inside the KD-tree.
          m_workResult->addPath(lPath);
          needCreateNew = true;
        }
      }
    }

    m_workResult = nullptr;
  }

  MTS_DECLARE_CLASS()
protected:
  /// Virtual destructor
  virtual ~GradientPhotonWorker() = default;
protected:
  const GPMConfig &m_config;
  ref<Scene> m_scene;

  int m_workerID;
  ref<PhotonPathVector> m_workResult;
  ShootingThreadData *m_myThdata;
  Float m_beamRadius;
};

GradientPhotonProcess::GradientPhotonProcess(const GPMConfig &config, const Point &sensorPos,
                                             size_t granularity, bool isLocal, bool autoCancel,
                                             const void *progressReporterPayload,
                                             std::vector<ShootingThreadData> &shootingData, Float beamRadius)
    : ParticleProcess(config.deterministic ? ParticleProcess::ETrace : ParticleProcess::EGather,
                      config.volumePhotonCount + config.photonCount, granularity, "Gathering photons",
                      progressReporterPayload), m_config(config), m_photonMap(nullptr), m_nbSurfaceLT(0),
      m_photonVolumeMap(nullptr), m_nbVolumeLT(0),
      m_isLocal(isLocal), m_autoCancel(autoCancel), m_excess(0),
      m_numShotSurface(0), m_numShotVolume(0), m_beamMap(nullptr), m_beamRadius(beamRadius),
      m_thdata(shootingData) {
  if (m_config.needSurfaceRendering()) {
    m_photonMap = new GPhotonMap(m_config.photonCount,
                                 true,
                                 sensorPos,
                                 m_config.cameraSphere);
  }

  if (!EVolumeTechniqueHelper::needBeams(m_config.volTechnique)) {
    m_photonVolumeMap = new GPhotonMap(m_config.volumePhotonCount,
                                       false,
                                       sensorPos,
                                       m_config.cameraSphere);
    m_storeBeams = false;
  } else {
    if (m_beamRadius == 0) {
      SLog(EError, "Need to set the beam radius before tracing");
    }
    m_beamMap = new LTBeamMap(m_config, sensorPos);
    m_storeBeams = true;
  } // Force to 0 the volume photons

  m_workerID = 0;
  if (m_config.deterministic) {
    ref<Scheduler> sched = Scheduler::getInstance();
    size_t nWorkers = sched->getWorkerCount();
    m_workGenerated.resize(nWorkers);
    for (size_t i = 0; i < nWorkers; ++i)
      m_workGenerated[i] = false;
  }
}

bool GradientPhotonProcess::isLocal() const {
  return m_isLocal;
}

ref<WorkProcessor> GradientPhotonProcess::createWorkProcessor() const {
  //std::cout << "Spawn GVPM beam work processor. ID " << m_workerID << std::endl;
  // NOTE: m_workerID is not necessarily equal to Mitsuba's local worker ID. Not a big deal 
  return new GradientPhotonWorker(m_config, m_workerID++, m_thdata, m_beamRadius);
}

int GradientPhotonProcess::pushSurfaceLT(int idWorker, const Path *lt) {
  if (m_photonMap == nullptr) {
    return -1;
  }

  // Add regular photon if it is needed
  int nPPushed = m_photonMap->tryAppend(idWorker, lt, m_config.minDepth);
  if (nPPushed != -1) {
    m_numShotSurface += 1;
  } else {
    m_excess++;         // Count excessive light paths
    // In this case, the path was excessive for volume and surface
  }

  return nPPushed;
}

int GradientPhotonProcess::pushVolumeLT(int idWorker, const Path *lt) {
  if (m_photonVolumeMap == nullptr && m_beamMap == nullptr) {
    return -1;
  }

  // Check the capacity
  bool isOutCapacity = outVolumeCapacity();
  if (isOutCapacity || lt == nullptr) {
    if (!isOutCapacity) {
      // Need to count this path
      m_numShotVolume += 1;
      return 0;
    } else {
      m_excess++;
      return -1;
    }
  }

  int nPPushed = 0;
  if (!m_storeBeams) {
    nPPushed = m_photonVolumeMap->tryAppend(idWorker,
                                            lt,
                                            m_config.minDepth);
  } else {
    nPPushed = m_beamMap->tryAppendLT(idWorker,
                                      lt,
                                      m_beamRadius,
                                      m_config.minDepth);
  }

  // Add regular photon if it is needed
  if (nPPushed != -1) {
    m_numShotVolume += 1;
  } else {
    m_excess++;         // Count excessive light paths
    // In this case, the path was excessive for volume and surface
  }

  return nPPushed;
}

void GradientPhotonProcess::processResult(const WorkResult *wr, bool cancelled) {
  if (cancelled)
    return;
  const PhotonPathVector &vec = *static_cast<const PhotonPathVector *>(wr);
  LockGuard lock(m_resultMutex);

  size_t nPhotons = 0;
  for (size_t i = 0; i < vec.getLightPathCount(); ++i) {
    int nbSurfPushed = pushSurfaceLT(vec.idWorker, vec[i]);
    int nbVolPushed = pushVolumeLT(vec.idWorker, vec[i]);
    nPhotons += std::max(0, nbSurfPushed) + std::max(0, nbVolPushed);
  }

  increaseResultCount(nPhotons);
}

ParallelProcess::EStatus GradientPhotonProcess::generateWork(WorkUnit *unit, int worker) {
  /* Use the same approach as PBRT for auto canceling */
  LockGuard lock(m_resultMutex);
  size_t numShot = (m_numShotVolume + m_numShotSurface);
  size_t photonCount = m_config.volumePhotonCount + m_config.photonCount;
  if (m_autoCancel && numShot > 100000
      && unsuccessful(photonCount, size(), numShot)) {
    Log(EInfo, "Not enough photons could be collected, giving up");
    return EFailure;
  }

  if (m_config.deterministic) {
    /* In this mode we only allow one work unit per worker. No worker should take over other worker's work, 
     * which might prevent other worker from being spawn. When this happens determinism is unknown.  
     */
    if (ParticleProcess::hasMoreWork() && m_workGenerated[worker]) {
      unit->clear();
      SLog(EInfo, "Skip regenerating work unit on a previous worker.");
      return ParallelProcess::EStatus::ESuccess;
    }

    // Proceed to generate work for this worker
    m_workGenerated[worker] = true;
  }

  return ParticleProcess::generateWork(unit, worker);
}

MTS_IMPLEMENT_CLASS(GradientPhotonProcess, false, ParticleProcess)
MTS_IMPLEMENT_CLASS_S(GradientPhotonWorker, false, ParticleTracer)
MTS_IMPLEMENT_CLASS(PhotonPathVector, false, WorkResult)

MTS_NAMESPACE_END