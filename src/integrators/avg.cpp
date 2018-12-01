// STL includes
#include <fstream>

// MTS includes
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

class AvgIntegrator : public Integrator {
public:
  AvgIntegrator(const Properties &props) : Integrator(props) {
    m_maxPasses = props.getInteger("maxPasses");
    m_maxRenderingTime = props.getInteger("maxRenderingTime");
    m_dumpIteration = props.getInteger("dumpIteration", 5); // Every 5 frames
    m_finiteGrad = props.getBoolean("finiteGrad", false);
    m_twostepalgo = props.getBoolean("twostepalgo", false);
    m_iterationMult = props.getInteger("iterationMult", 1);
  }

  bool preprocess(const Scene *scene, RenderQueue *queue,
      const RenderJob *job, int sceneResID, int sensorResID,
      int samplerResID) {
    if(!m_twostepalgo)
      m_subIntegrator->preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID);
    m_stop = false;
    return true;
  }

  bool render(Scene *scene, RenderQueue *queue, const RenderJob *job,
      int sceneResID, int sensorResID, int samplerResID) {

    /// Get data
    ref<Scheduler> sched = Scheduler::getInstance();
    ref<Sensor> sensor = scene->getSensor();
    ref<Film> film = sensor->getFilm();
    size_t nCores = sched->getCoreCount();
    Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT " %s, " SSE_STR ") ..",
      film->getCropSize().x, film->getCropSize().y,
      nCores, nCores == 1 ? "core" : "cores");

    Vector2i cropSize = film->getCropSize();

    //FIXME: No support of cropping
//    Point2i cropOffset = film->getCropOffset();
//    int blockSize = scene->getBlockSize();

    /// Create bitmaps
    ref<Bitmap> bitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, cropSize);
    ref<Bitmap> bitmapAvg = bitmap->clone();
    bitmapAvg->clear();

    /// Create all struct
    std::string timeFilename = scene->getDestinationFile().string()
            + "_time.csv";
    std::ofstream timeFile(timeFilename.c_str());
    ref<Timer> renderingTimer = new Timer;
    Float cumulativeTime = 0.f;

    for(int it = 1; it < m_maxPasses && (!m_stop) && (cumulativeTime < m_maxRenderingTime); it++) {
      bitmap->clear(); //< Clear bitmap this pass
      film->clear(); //< and setup

      /// Recall preprocess if it is a two step algorithm
      if(m_twostepalgo)
        m_subIntegrator->preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID);

      /// Render
      m_subIntegrator->render(scene, queue, job, sceneResID, sensorResID, samplerResID);
      film->develop(Point2i(0,0), cropSize, Point2i(0,0) ,bitmap);

      /// Avg results
      if(!m_stop) {
        bitmapAvg->scale(it-1);
        bitmapAvg->accumulate(bitmap);
        bitmapAvg->scale(1.f/(it));

        /// Write image
        if((it*m_iterationMult) % m_dumpIteration == 0){
          /// Path computation
          std::stringstream ss;
          ss << scene->getDestinationFile().string() << "_pass_" << it*m_iterationMult;
          std::string path = ss.str();

          /// Develop image
          film->setBitmap(bitmapAvg);
          film->setDestinationFile(path,0);
          film->develop(scene,0.f);

          /// Develop GX Finite
          if(m_finiteGrad) {
            ref<Bitmap> gXBitmap = bitmapAvg->clone();
            ref<Bitmap> gYBitmap = bitmapAvg->clone();

            computeGradientFinite(cropSize, *bitmapAvg, gXBitmap.get(), gYBitmap.get());

            develop(scene, film, gXBitmap, it*m_iterationMult, "_dxAbs_");
            develop(scene, film, gYBitmap, it*m_iterationMult, "_dyAbs_");
          }

          /// Revert destination file
          film->setDestinationFile(scene->getDestinationFile(),0);
        }

        /// Time it
        unsigned int milliseconds = renderingTimer->getMilliseconds();
        for(int localTime = 0; localTime < m_iterationMult; localTime++) {
            timeFile << (milliseconds / (1000.f*m_iterationMult)) << ",\n";
        }
        timeFile.flush();
        Log(EInfo, "Rendering time: %i, %i", milliseconds / 1000,
          milliseconds % 1000);
        cumulativeTime += (milliseconds / 1000.f);

        renderingTimer->reset();

        // === Print the statistic at each step of the rendering
        // to see the algorithm behaviors.
        Statistics::getInstance()->printStats();
      }
    }

    return true;
  }

  void cancel() {
    m_subIntegrator->cancel();
    m_stop = true;
  }

  ///////////////
  // Config subintergrator
  ///////////////
  void configureSampler(const Scene *scene, Sampler *sampler) {
    m_subIntegrator->configureSampler(scene, sampler);
  }

  void bindUsedResources(ParallelProcess *proc) const {
    m_subIntegrator->bindUsedResources(proc);
  }

  void wakeup(ConfigurableObject *parent,
      std::map<std::string, SerializableObject *> &params) {
    m_subIntegrator->wakeup(this, params);
  }

  void addChild(const std::string &name, ConfigurableObject *child) {
    const Class *cClass = child->getClass();

    if (cClass->derivesFrom(MTS_CLASS(Integrator))) {
      m_subIntegrator = static_cast<Integrator *>(child);
      m_subIntegrator->setParent(this);
    } else {
      Integrator::addChild(name, child);
    }
  }

  void computeGradientFinite(const Vector2i &cropSize, const Bitmap& pixel, Bitmap *gXBitmap, Bitmap *gYBitmap, bool useAbs = true) {
    Spectrum *targetGX = (Spectrum *) gXBitmap->getUInt8Data();
    Spectrum *targetGY = (Spectrum *) gYBitmap->getUInt8Data();
    for (int y = 1; y < cropSize.y - 1; ++y) {
      for (int x = 1; x < cropSize.x - 1; ++x) {
        Spectrum curr = pixel.getPixel(Point2i(x,y));
        Spectrum right = pixel.getPixel(Point2i(x+1,y));
        Spectrum top = pixel.getPixel(Point2i(x,y+1));

        Spectrum gX(0.f);
        Spectrum gY(0.f);


        gX += (right- curr);
        gY += (top - curr);


        if (useAbs) {
          targetGX[y * pixel.getWidth() + x] = gX.abs();
          targetGY[y * pixel.getWidth() + x] = gY.abs();
        } else {
          targetGX[y * pixel.getWidth() + x] = gX;
          targetGY[y * pixel.getWidth() + x] = gY;
        }
      }
    }
  }

  void develop(Scene *scene, Film *film, Bitmap *bitmap,
               int currentIteration, const std::string &suffixName = "_") {
    std::stringstream ss;
    ss << scene->getDestinationFile().string() << suffixName
       << currentIteration;
    std::string path = ss.str();

    film->setBitmap(bitmap);
    film->setDestinationFile(path, 0);
    film->develop(scene, 0.f);

  }

  std::string toString() const {
    std::ostringstream oss;
    oss << "AvgIntegrator[" << endl
      << "  subIntegrator = " << indent(m_subIntegrator->toString()) << "," << endl
      << "]";
    return oss.str();
  }

  MTS_DECLARE_CLASS()

protected:
  int m_maxPasses;
  float m_maxRenderingTime;
  ref<Integrator> m_subIntegrator;
  bool m_stop;
  bool m_twostepalgo;
  int m_dumpIteration;
  bool m_finiteGrad;
  int m_iterationMult;
};

MTS_IMPLEMENT_CLASS(AvgIntegrator, false, Integrator)
MTS_EXPORT_PLUGIN(AvgIntegrator, "Avg pass integrator");
MTS_NAMESPACE_END
