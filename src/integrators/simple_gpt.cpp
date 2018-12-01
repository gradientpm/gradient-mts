#include <mitsuba/core/plugin.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/fstream.h>

#include <omp.h>

#include "export.h"

MTS_NAMESPACE_BEGIN

class ReplaySampler: public  Sampler {

public:
  ReplaySampler(ref<Sampler> otherSampler):
    Sampler(Properties{}),
    m_sampler(otherSampler), m_sampleIndex(0) {
    m_sampleCount = otherSampler->getSampleCount();
  }

  Float next1D() override {
    if(m_random_numbers.size() > m_sampleIndex) {
      Float rand = m_random_numbers[m_sampleIndex];
      m_sampleIndex+= 1;
      return rand;
    }

    Float rand = m_sampler->next1D();
    m_random_numbers.push_back(rand);
    m_sampleIndex += 1;
    return rand;
  }

  Point2 next2D() override {
    return Point2(next1D(), next1D());
  }

  void advance() override {
    m_sampler->advance();
  }

  void clear() {
    m_random_numbers.clear();
    m_sampleIndex = 0;
  }

  void reset() {
    m_sampleIndex = 0;
  }

protected:
  ref<Sampler> m_sampler;
  std::vector<Float> m_random_numbers;
  size_t m_sampleIndex;
};

class SimpleGPTIntegrator : public Integrator {
public:
  SimpleGPTIntegrator(const Properties &props) : Integrator(props) {
  }

  bool preprocess(const Scene *scene, RenderQueue *queue,
                  const RenderJob *job, int sceneResID, int sensorResID,
                  int samplerResID) {
    m_subIntegrator->preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID);
    m_stop = false;
    return true;
  }

  bool render(Scene *scene, RenderQueue *queue, const RenderJob *job,
              int sceneResID, int sensorResID, int samplerResID) {

    /// Get data from Mitsuba
    auto sensor = scene->getSensor();
    auto film = sensor->getFilm();
    Log(EInfo, "Starting render job (%ix%i)",
        film->getCropSize().x, film->getCropSize().y);
    Vector2i cropSize = film->getCropSize();
    size_t sampleCount = sensor->getSampler()->getSampleCount();
    auto pixel_index = [&](size_t x, size_t y) -> size_t {
      return y * cropSize.x + x;
    };

    /// Create bitmaps
    ref <Bitmap> bitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, cropSize);
    ref <Bitmap> bitmapGx = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, cropSize);
    ref <Bitmap> bitmapGy = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, cropSize);
    bitmap->clear();
    bitmapGx->clear();
    bitmapGx->clear();

    /// Create the random number
    /// generator that will capture the random sequence
    auto sampler = new ReplaySampler(sensor->getSampler());

    // Rendering loop.

    film->clear(); //< and setup
    Spectrum *primal = (Spectrum *) bitmap->getUInt8Data();
    Spectrum *Gx = (Spectrum *) bitmapGx->getUInt8Data();
    Spectrum *Gy = (Spectrum *) bitmapGy->getUInt8Data();

    // Rendering the image
    for (size_t y = 0; y < cropSize.y; y++) {
      for (size_t x = 0; x < cropSize.x; x++) {

        for (size_t j = 0; j < sampler->getSampleCount(); j++) {
          // Clear the capture random number sequence.
          sampler->clear();
          auto base_contrib = Li(Point2i(x,y), sampler, sensor, scene);
          primal[pixel_index(x,y)] += base_contrib;

          if(x > 0) {
            sampler->reset();
            Gx[pixel_index(x-1,y)] += 0.5*(base_contrib - Li(Point2i(x-1,y), sampler, sensor, scene));
          }
          if(y > 0) {
            sampler->reset();
            Gy[pixel_index(x,y-1)] += 0.5*(base_contrib - Li(Point2i(x,y-1), sampler, sensor, scene));
          }
          if(x < cropSize.x - 1) {
            sampler->reset();
            Gx[pixel_index(x,y)] += 0.5*(Li(Point2i(x+1,y), sampler, sensor, scene) - base_contrib);
          }
          if(x < cropSize.y - 1) {
            sampler->reset();
            Gy[pixel_index(x,y)] += 0.5*(Li(Point2i(x,y+1), sampler, sensor, scene) - base_contrib);
          }

          sampler->advance();
        }
      }
    }

    /// Dump images
    develop(scene, film, bitmap, "_primal");
    develop(scene, film, bitmapGx, "_Gx");
    develop(scene, film, bitmapGy, "_Gy");

    // Do the reconstruction
    ref <Bitmap> reconstructionBitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, cropSize);
    reconstructionBitmap->clear();

    // Initialize the recons
    ref<Bitmap> nextReconsBitmap = reconstructionBitmap->clone();
    auto next_recons = (Spectrum*)nextReconsBitmap->getData();
    auto recons = (Spectrum*)reconstructionBitmap->getData();
    for(size_t t = 0; t < cropSize.x * cropSize.y; t++) {
      recons[t] = primal[t];
    }

    for(size_t iter = 0; iter < 50; iter++) {
      for (int i = 0; i < cropSize.y; i++) {
        for (int j = 0; j < cropSize.x; ++j) {
          int cur = i * cropSize.x + j;

          // Current pixel
          int nbAccum = 1;
          Spectrum next = recons[cur];

          // Bottom
          if (i != 0) {
            int bottom = (i - 1) * cropSize.x + j;
            next += (Gy[bottom] + recons[bottom]);
            nbAccum += 1;
          }

          // Up
          if (i != cropSize.y - 1) {
            int top = (i + 1) * cropSize.x + j;
            next += (-Gy[cur] + recons[top]);
            nbAccum += 1;
          }

          // Left
          if (j != 0) {
            int left = i * cropSize.x + j - 1;
            next += (Gx[left] + recons[left]);
            nbAccum += 1;
          }

          // Right
          if (j != cropSize.x - 1) {
            int right = i * cropSize.x + j + 1;
            next += (-Gx[cur] + recons[right]);
            nbAccum += 1;
          }

          next_recons[cur] = next * (1.0 / (Float) nbAccum);
        }
      }

      // Recopy the informations
      for(size_t t = 0; t < cropSize.x * cropSize.y; t++) {
        recons[t] = next_recons[t];
      }
    }
    develop(scene, film, reconstructionBitmap, "_recons");

    return true;
  }

  Spectrum Li(const Point2i pix, Sampler* sampler, const Sensor* sensor, const Scene *scene) {
    // Get information about the sensor
    Float diffScaleFactor = 1.0f /
        std::sqrt((Float) sampler->getSampleCount());
    bool needsApertureSample = sensor->needsApertureSample();
    bool needsTimeSample = sensor->needsTimeSample();

    RadianceQueryRecord rRec(scene, sampler);
    Point2 apertureSample(0.5f);
    Float timeSample = 0.5f;
    RayDifferential sensorRay;
    uint32_t queryType = RadianceQueryRecord::ESensorRay;
    if (!sensor->getFilm()->hasAlpha()) /* Don't compute an alpha channel if we don't have to */
      queryType &= ~RadianceQueryRecord::EOpacity;

    // Sample the camera ray
    rRec.newQuery(queryType, sensor->getMedium());
    Point2 samplePos(Point2(pix) + Vector2(rRec.nextSample2D()));
    if (needsApertureSample)
      apertureSample = rRec.nextSample2D();
    if (needsTimeSample)
      timeSample = rRec.nextSample1D();
    // Inside spec there is the importance of the sensor
    Spectrum spec = sensor->sampleRayDifferential(
        sensorRay, samplePos, apertureSample, timeSample);
    sensorRay.scaleDifferential(diffScaleFactor);

    return spec * m_subIntegrator->Li(sensorRay, rRec) / (Float) sampler->getSampleCount();
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
      m_subIntegrator = dynamic_cast<MonteCarloIntegrator *>(child);
      if (m_subIntegrator == 0) {
        SLog(EError, "The type of the sub integrator need to be MonteCarlo type");
      }
      m_subIntegrator->setParent(this);
    } else {
      Integrator::addChild(name, child);
    }
  }

  void develop(Scene *scene, Film *film, Bitmap *bitmap, const std::string &suffixName = "_") {
    std::stringstream ss;
    ss << scene->getDestinationFile().string() << suffixName;
    std::string path = ss.str();

    saveExr(bitmap, path + ".exr");

  }

  std::string toString() const {
    std::ostringstream oss;
    oss << "SimpleGPTIntegrator[" << endl
        << "  subIntegrator = " << indent(m_subIntegrator->toString()) << "," << endl
        << "]";
    return oss.str();
  }

  MTS_DECLARE_CLASS()

protected:
  int m_maxPasses;
  ref <MonteCarloIntegrator> m_subIntegrator;
  bool m_stop;
};

MTS_IMPLEMENT_CLASS(SimpleGPTIntegrator,
false, Integrator)
MTS_EXPORT_PLUGIN(SimpleGPTIntegrator,
"Naive BCD integrator");
MTS_NAMESPACE_END
