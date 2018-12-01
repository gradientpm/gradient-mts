/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/render/gatherproc.h>
#include <mitsuba/render/renderqueue.h>

// BRE implementation
#include "bre.h"

// Photon Beams implementation
#include "beams.h"

// Photon plane implementation
#include "plane_struct.h"
#include "plane_accel.h"

// Other utils

MTS_NAMESPACE_BEGIN

struct Beam {
  Point p1;
  Point p2;
  const Medium *medium;
  Spectrum weight; // Initial weight of the beam
  Float t; // Initial distance for the scale reduction
  int depth;
  Spectrum transmittance; // The transmittance along the ray

  // Just a protection if p2 is not set...
  bool invalid;

  Beam(const Point &p, const Medium *m, const Spectrum &w, Float dist, int curDepth) :
      p1(p), medium(m), weight(w), t(dist), depth(curDepth), transmittance(0.0) {
    invalid = true;
  }

  bool sampleDistance(Float distance, const Vector& direction) {
    Ray rayTrans(p1, direction, 0.f);
    rayTrans.mint = Epsilon;
    rayTrans.maxt = distance;

    MediumSamplingRecord mRec;
    bool succeed = medium->sampleDistance(rayTrans, mRec, nullptr, EDistanceLong);

    transmittance = mRec.transmittance;
    p2 = rayTrans(mRec.t);
    invalid = false;

    return succeed;
  }

  bool isInvalid() const {
    return invalid;
  }
};

/// Represents one individual PPM gather point including relevant statistics
struct GatherPoint {
  // Surface informations
  Intersection its;
  Float radius;
  Spectrum weight;
  Spectrum flux;
  Spectrum emission;
  Float N;
  int depth;
  Point2i pos;
  Float scale;

  // BSDF informations
  int sampledComponent;
  float pdfComponent;

  // Volume informations
  std::vector<Beam> beams;
  Float NVol;
  Spectrum fluxVol;
  Float scaleVol;

  inline GatherPoint() : radius(0.0f), weight(0.0f), flux(0.0f), emission(0.0f), N(0.0f), depth(-1), scale(1.0f),
                         sampledComponent(-1), pdfComponent(1.0f),
                         NVol(0.f), fluxVol(0.f), scaleVol(1.f) {}

  void resetTemp() {
    // Nothing to do?
  }

  inline void rescaleRadii(Float v) {
    flux *= v;
  }

};

MTS_NAMESPACE_END

#include "utilities/initializeRadius.h"
#include "photon_proc.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{sppm}{Stochastic progressive photon mapping integrator}
* \order{8}
* \parameters{
*     \parameter{maxDepth}{\Integer}{Specifies the longest path depth
*         in the generated output image (where \code{-1} corresponds to $\infty$).
*	       A value of \code{1} will only render directly visible light sources.
*	       \code{2} will lead to single-bounce (direct-only) illumination,
*	       and so on. \default{\code{-1}}
*	   }
*     \parameter{photonCount}{\Integer}{Number of photons to be shot per iteration\default{250000}}
*     \parameter{initialRadius}{\Float}{Initial radius of gather points in world space units.
*         \default{0, i.e. decide automatically}}
*     \parameter{alpha}{\Float}{Radius reduction parameter \code{alpha} from the paper\default{0.7}}
*     \parameter{granularity}{\Integer}{
    Granularity of photon tracing work units for the purpose
    of parallelization (in \# of shot particles) \default{0, i.e. decide automatically}
*     }
*	   \parameter{rrDepth}{\Integer}{Specifies the minimum path depth, after
*	      which the implementation will start to use the ``russian roulette''
*	      path termination criterion. \default{\code{5}}
*	   }
*     \parameter{maxPasses}{\Integer}{Maximum number of passes to render (where \code{-1}
*        corresponds to rendering until stopped manually). \default{\code{-1}}}
* }
* This plugin implements stochastic progressive photon mapping by Hachisuka et al.
* \cite{Hachisuka2009Stochastic}. This algorithm is an extension of progressive photon
* mapping (\pluginref{ppm}) that improves convergence
* when rendering scenes involving depth-of-field, motion blur, and glossy reflections.
*
* Note that the implementation of \pluginref{sppm} in Mitsuba ignores the sampler
* configuration---hence, the usual steps of choosing a sample generator and a desired
* number of samples per pixel are not necessary. As with \pluginref{ppm}, once started,
* the rendering process continues indefinitely until it is manually stopped.
*
* \remarks{
*    \item Due to the data dependencies of this algorithm, the parallelization is
*    limited to the local machine (i.e. cluster-wide renderings are not implemented)
*    \item This integrator does not handle participating media
*    \item This integrator does not currently work with subsurface scattering
*    models.
* }
*/
class SPPMIntegrator : public Integrator {
public:
  SPPMIntegrator(const Properties &props) : Integrator(props) {

    m_initialScale = props.getFloat("initialScale", 1.f);
    m_initialScaleVolume = props.getFloat("initialScaleVolume", 1.f);
    /* Alpha parameter from the paper (influences the speed, at which the photon radius is reduced) */
    m_alpha = props.getFloat("alpha", .7);
    /* Number of photons to shoot in each iteration */
    m_photonSurfaceCount = props.getInteger("photonCount", 250000);
    m_photonVolumeCount = props.getInteger("volumePhotonCount", 250000);
    /* Granularity of the work units used in parallelizing the
       particle tracing task (default: choose automatically). */
    m_granularity = props.getInteger("granularity", 0);
    /* Longest visualized path length (<tt>-1</tt>=infinite). When a positive value is
       specified, it must be greater or equal to <tt>2</tt>, which corresponds to single-bounce
       (direct-only) illumination */
    m_maxDepth = props.getInteger("maxDepth", -1);
    m_minDepth = props.getInteger("minDepth", 0);
    /* Depth to start using russian roulette */
    m_rrDepth = props.getInteger("rrDepth", 3);
    /* Indicates if the gathering steps should be canceled if not enough photons are generated. */
    m_autoCancelGathering = props.getBoolean("autoCancelGathering", false);
    /* Maximum number of passes to render. -1 renders until the process is stopped. */
    m_maxPasses = props.getInteger("maxPasses", -1);
    m_mutex = new Mutex();
    if (m_maxDepth <= 1 && m_maxDepth != -1)
      Log(EError, "Maximum depth must be set to \"2\" or higher!");
    if (m_maxPasses <= 0 && m_maxPasses != -1)
      Log(EError, "Maximum number of Passes must either be set to \"-1\" or \"1\" or higher!");

    // Create the instance to generate the gather point in the similar way for all
    // PPM/SPPM codes
    m_gpManager = new RadiusInitializer(props);

    m_maxRenderingTime = props.getInteger("maxRenderingTime", INT_MAX);
    m_dumpIteration = props.getInteger("dumpIteration", 1);
    m_nbInternalIteration = props.getInteger("nbInternalIteration", 10);
    if (m_maxRenderingTime != INT_MAX && m_maxPasses != INT_MAX) {
      Log(EError, "Max pass and time is incompatible!");
    }

    // Debug options / Developement proposes
    m_computeGrad = props.getBoolean("computeGrad", true);
    m_directTracing = props.getBoolean("directTracing", true);

    m_surfaceRendering = props.getBoolean("surfaceRendering", true);
    m_volumeRendering = props.getBoolean("volumeRendering", true);

    std::string volRenderingTech = props.getString("volTechnique", "distance");
    m_volTechnique = parseVolumeTechnique(volRenderingTech);
    m_convertLong = props.getBoolean("convertLong", false);
    m_deterministic = props.getBoolean("deterministic", false);

    // Special props for distance sampling strategy
    m_nbCameraSamples = props.getInteger("nbCameraSamples", 40);

    // Error checking:
    if (m_surfaceRendering && m_photonSurfaceCount == 0)
      Log(EError,
          "No surface photons and need to render surfaces LT");
    if (m_volumeRendering && m_photonVolumeCount == 0)
      Log(EError,
          "No volume photons/beams and need to render volume LT");
    m_nanCheck = props.getBoolean("nanCheck", false);

    // Camera debug
    m_minCameraDepth = [&props]() -> size_t {
      auto v = props.getInteger("minCameraDepth", 0);
      if(v < 0) { SLog(EError, "minCamera depth need to be null or positive"); }
      return v;
    }();
    m_maxCameraDepth = props.getInteger("maxCameraDepth", -1);

    // Add variance but need to be able to reproduce
    // results of path tracing.
    m_adjointCompensation = true;

    // Based on the "optimization" for independence
    // between number of samples and alpha.
    // See "Progressive Photon Beams" section 5.2
    m_independentScale = false; // Put at false as the scaling this too strong
    m_cameraSphere = props.getFloat("cameraSphere", 1.f);
    m_forceAPA = props.getString("forceAPA", "");
  }

  SPPMIntegrator(Stream *stream, InstanceManager *manager)
      : Integrator(stream, manager) {}

  void serialize(Stream *stream, InstanceManager *manager) const {
    Integrator::serialize(stream, manager);
    Log(EError, "Network rendering is not supported!");
  }

  void cancel() {
    m_running = false;
  }

  void scaleVolumeAPA(int it) {
    // Scale the gathering size for the next iteration

    it -= 1; // Fix the bug as it == 1 at the first iteration.
    double ratioVolAPA = 1.0;

    if (!m_independentScale) {
      ratioVolAPA = (it + m_alpha) / (it + 1);
    } else {
      for (int i = 0; i < m_photonVolumeCount; i++) {
        size_t newIt = m_photonVolumeCount * it + i;
        ratioVolAPA *= double(newIt + m_alpha) / double(newIt + 1.f);
      }
    }

    if (m_forceAPA == "") {
      // If not force APA is provided use the classical reduction rates
      if (EVolumeTechniqueHelper::use3DKernel(m_volTechnique)) {
        globalScaleVolume = globalScaleVolume * std::cbrt(ratioVolAPA);
      } else if (EVolumeTechniqueHelper::use2DKernel(m_volTechnique)) {
        globalScaleVolume = globalScaleVolume * std::sqrt(ratioVolAPA);
      } else {
        globalScaleVolume = globalScaleVolume * ratioVolAPA;
      }
    } else {
      if (m_forceAPA == "1D") {
        globalScaleVolume = globalScaleVolume * ratioVolAPA;
      } else if (m_forceAPA == "2D") {
        globalScaleVolume = globalScaleVolume * std::sqrt(ratioVolAPA);
      } else if (m_forceAPA == "3D") {
        globalScaleVolume = globalScaleVolume * std::cbrt(ratioVolAPA);
      } else {
        SLog(EError, "No Force APA: %s", m_forceAPA.c_str());
      }
    }
  }

  bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
                  int sceneResID, int sensorResID, int samplerResID) {
    Integrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID);
    const_cast<Scene *>(scene)->weightEmitterFlux();

    if(scene->getMedia().size() == 0) {
        m_volumeRendering = false;
        m_photonVolumeCount = 0;
    } else {
        // The bounding box of the smoke
        m_smokeAABB = max_AABB_medium(scene);
    }

    // Check if there is one surface with a BSDF
    // if there is no surface like this, disable surface rendering
    bool haveValidBSDF = false;
    for (auto mesh: scene->getMeshes()) {
      if (!mesh->isMediumTransition() && !mesh->isEmitter())
        haveValidBSDF = true;
    }
    if (!haveValidBSDF)
      m_surfaceRendering = false;

    if (!m_surfaceRendering || m_photonSurfaceCount == 0) {
      Log(EInfo, "Only volume rendering, change probability to sample the medium");
      auto medium = scene->getMedia().begin();
      while (medium != scene->getMedia().end()) {
        const_cast<Medium *>(medium->get())->computeOnlyVolumeInteraction();
        medium++;
      }
    }

    // Update size of exclusion for the camera sphere
    m_cameraSphere = m_smokeAABB.getBSphere().radius * m_cameraSphere * POURCENTAGE_BS;

    return true;
  }

  bool render(Scene *scene, RenderQueue *queue,
              const RenderJob *job, int sceneResID, int sensorResID, int unused) {
    ref<Scheduler> sched = Scheduler::getInstance();
    ref<Sensor> sensor = scene->getSensor();
    ref<Film> film = sensor->getFilm();
    size_t nCores = sched->getCoreCount();
    Log(EInfo, "Starting render job (%ix%i, "
        SIZE_T_FMT
        " %s, "
        SSE_STR
        ") ..",
        film->getCropSize().x, film->getCropSize().y,
        nCores, nCores == 1 ? "core" : "cores");

    Vector2i cropSize = film->getCropSize();

    m_gatherBlocks.clear();
    m_running = true;
    m_totalEmitted = 0;
    m_totalPhotons = 0;

    m_totalEmittedVol = 0;
    m_totalPhotonsVol = 0;

    ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
        createObject(MTS_CLASS(Sampler), Properties("independent")));

    int blockSize = scene->getBlockSize();

    /* Allocate memory */
    m_bitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, cropSize);
    m_bitmap->clear();
    for (int yofs = 0; yofs < cropSize.y; yofs += blockSize) {
      for (int xofs = 0; xofs < cropSize.x; xofs += blockSize) {
        m_gatherBlocks.push_back(std::vector<GatherPoint>());
        m_offset.push_back(Point2i(xofs, yofs));
        std::vector<GatherPoint> &gatherPoints = m_gatherBlocks[m_gatherBlocks.size() - 1];
        int nPixels = std::min(blockSize, cropSize.y - yofs)
            * std::min(blockSize, cropSize.x - xofs);
        gatherPoints.resize(nPixels);
        for (int i = 0; i < nPixels; ++i) {
          gatherPoints[i].radius = 0;
          gatherPoints[i].scale = m_initialScale;
          gatherPoints[i].scaleVol = m_initialScaleVolume;
        }
      }
    }
    globalScaleVolume = m_initialScaleVolume;

    // Initialize the class responsible to the GP generation
    // and the radii initialization
    m_gpManager->init(scene, m_maxDepth, m_gatherBlocks, m_offset);

    /* Create a sampler instance for every core */
    std::vector<SerializableObject *> samplers(sched->getCoreCount());
    std::vector<Sampler *> samplersConcrete(sched->getCoreCount());
    for (size_t i = 0; i < sched->getCoreCount(); ++i) {
      ref<Sampler> clonedSampler = sampler->clone();
      clonedSampler->incRef();
      samplers[i] = clonedSampler.get();
      samplersConcrete[i] = clonedSampler.get();
    }

    int samplerResID = sched->registerMultiResource(samplers);

#ifdef MTS_DEBUG_FP
    enableFPExceptions();
#endif

    // Create the files to dump information about the rendering
    // Also create timer to track algorithm performance
    std::string timeFilename = scene->getDestinationFile().string()
        + "_time.csv";
    std::ofstream timeFile(timeFilename.c_str());
    ref<Timer> renderingTimer = new Timer;

    // Rendering loop
    int it = 1;
    int global_it = 1;
    while (m_running && (m_maxPasses == -1 || global_it < m_maxPasses)) {
      for(int local_it = 0; local_it < m_nbInternalIteration; local_it++) {
          Log(EInfo, "Regenerating gather points positions and radius!");
          m_gpManager->regeneratePositionAndRadius();
          Log(EInfo, "Done regenerating!");
          m_gpManager->rescaleFlux();

          if (it == 1) {
            m_imgGP.resize(cropSize.x);
            for (int x = 0; x < cropSize.x; x += 1) {
              m_imgGP[x].resize(cropSize.y);
            }
            for (size_t idBlock = 0; idBlock < m_gatherBlocks.size(); idBlock++) {
              std::vector<GatherPoint> &currBlock = m_gatherBlocks[idBlock];
              for (GatherPoint &gp: currBlock) {
                // The position of the gatherpoint doesn't have cropOffset
                int xcrop = gp.pos.x;
                int ycrop = gp.pos.y;
                m_imgGP[xcrop][ycrop] = &gp;
              }
            }
          }

          photonMapPass(it, queue, job, film, sceneResID,
                        sensorResID, samplerResID, scene,
                        samplersConcrete);

          if (m_nanCheck) {
              for (int blockIdx = 0; blockIdx < (int) m_gatherBlocks.size(); ++blockIdx) {
                  std::vector<GatherPoint> &gatherPoints = m_gatherBlocks[blockIdx];

                  Spectrum *target = (Spectrum *) m_bitmap->getUInt8Data();
                  for (size_t i = 0; i < gatherPoints.size(); ++i) {
                      GatherPoint &gp = gatherPoints[i];
                      Spectrum value = target[gp.pos.y * m_bitmap->getWidth() + gp.pos.x];
                      if (value.isNaN()) {
                          SLog(EWarn, "Pixel (%d, %d) is NaN", gp.pos.x, gp.pos.y);
                          target[gp.pos.y * m_bitmap->getWidth() + gp.pos.x] = Spectrum(0.f);
                      }
                  }
              }
          }

          ++it;
      }


      // Write down some
      if ((global_it % m_dumpIteration) == 0) {
        develop(scene, film, m_bitmap, global_it, "_pass_");
        if (m_computeGrad) {
          // Compute gradient using finite differences
          ref<Bitmap> gXBitmap = m_bitmap->clone();
          ref<Bitmap> gYBitmap = m_bitmap->clone();

          computeGradientFinite(cropSize, gXBitmap.get(), gYBitmap.get());

          develop(scene, film, gXBitmap, global_it, "_dxAbs_");
          develop(scene, film, gYBitmap, global_it, "_dyAbs_");
        }
      }

      // === Update the log time
      unsigned int milliseconds = renderingTimer->getMilliseconds();
      timeFile << (milliseconds / 1000.f) << ",\n";
      timeFile.flush();
      Log(EInfo, "Rendering time: %i, %i", milliseconds / 1000,
          milliseconds % 1000);
      m_maxRenderingTime -= (milliseconds / 1000);
      if (m_maxRenderingTime < 0) {
        m_running = false;
        Log(EInfo, "Max time reaching !");
      }

      // === update variables && free memory
      renderingTimer->reset();

      ++global_it;
    }

#ifdef MTS_DEBUG_FP
    disableFPExceptions();
#endif

    for (size_t i = 0; i < samplers.size(); ++i)
      samplers[i]->decRef();

    sched->unregisterResource(samplerResID);
    return true;
  }

  void photonMapPass(int it, RenderQueue *queue, const RenderJob *job,
                     Film *film, int sceneResID, int sensorResID, int samplerResID, Scene *scene,
                     std::vector<Sampler *> &samplers) {
    Log(EInfo, "Performing a photon mapping pass %i ("
        SIZE_T_FMT
        " photons so far)",
        it, m_totalPhotons);
    ref<Scheduler> sched = Scheduler::getInstance();

    // Clear the bitmap to accumulate
    // contribution from the surface and the volume
    m_bitmap->clear();

    ////////////////// Surface Only
    if (m_surfaceRendering) {
      /* Generate the global photon map */
      ref<GatherPhotonProcess> proc = new GatherPhotonProcess(
          GatherPhotonProcess::EAllSurfacePhotons, m_photonSurfaceCount,
          m_granularity, m_maxDepth == -1 ? -1 : m_maxDepth - 1, m_rrDepth, true,
          m_autoCancelGathering, job, m_adjointCompensation, m_minDepth);

      proc->bindResource("scene", sceneResID);
      proc->bindResource("sensor", sensorResID);
      proc->bindResource("sampler", samplerResID);

      sched->schedule(proc);
      sched->wait(proc);

      ref<PhotonMap> photonMap = proc->getPhotonMap();
      photonMap->build();
      Log(EDebug, "Photon map full. Shot "
          SIZE_T_FMT
          " particles, excess photons due to parallelism: "
          SIZE_T_FMT, proc->getShotParticles(), proc->getExcessPhotons());

      Log(EInfo, "Gathering Surface...");
      m_totalEmitted += proc->getShotParticles();
      m_totalPhotons += photonMap->size();
      film->clear();

      size_t nCores = sched->getCoreCount();
      BlockScheduler blockSched(m_gatherBlocks.size(), nCores);
      blockSched.run([&](int blockIdx, int tid) {

        std::vector<GatherPoint> &gatherPoints = m_gatherBlocks[blockIdx];

        Spectrum *target = (Spectrum *) m_bitmap->getUInt8Data();
        for (size_t i = 0; i < gatherPoints.size(); ++i) {
          GatherPoint &gp = gatherPoints[i];
          Spectrum contrib(0.0f);

          // Surface estimation
          if (gp.depth != -1) {

            Spectrum flux(0.0f);

            Float M = (Float) photonMap->estimateRadianceGP(
                gp.its, gp.radius, flux, m_maxDepth == -1 ? INT_MAX : m_maxDepth - gp.depth,
                gp.sampledComponent);

            // Update GP
            Float N = gp.N;
            if (N + M != 0) {
              Float ratio = (N + m_alpha * M) / (N + M);
              gp.scale = gp.scale * std::sqrt(ratio);
              gp.flux = (gp.flux + gp.weight * flux);
              gp.flux /= gp.radius * gp.radius * M_PI; // become radiance
              gp.N = N + m_alpha * M;
            }
          }

          // When depth == -1, gp.flux is not rescaled and is still radiance value. No need to divide radius.
          contrib = gp.flux / ((Float) m_totalEmitted);
          if (m_directTracing) {
            contrib += gp.emission / it; // Monte carlo estimator
          }

          target[gp.pos.y * m_bitmap->getWidth() + gp.pos.x] += contrib;
        }
      });
    } // End surface rendering


    ////////////// Volume only
    if (m_volumeRendering) {
      switch (m_volTechnique) {
      case EVolBRE2D:
      case EVolBRE3D:volumePhotonPassBRE(it, queue, job, film, sceneResID, sensorResID, samplerResID, scene, samplers);
        break;
      case EBeamBeam1D:
      case EBeamBeam3D_Naive:
      case EBeamBeam3D_EGSR:
      case EBeamBeam3D_Optimized:
        volumePhotonBeamPass(it,
                             queue,
                             job,
                             film,
                             sceneResID,
                             sensorResID,
                             samplerResID,
                             scene,
                             samplers);
        break;
      case EVolPlane0D:
        volumePhotonPlanePass(it,
                              queue,
                              job,
                              film,
                              sceneResID,
                              sensorResID,
                              samplerResID,
                              scene,
                              samplers);
        break;
      case EDistance:
        volumePhotonPassDistance(it,
                                 queue,
                                 job,
                                 film,
                                 sceneResID,
                                 sensorResID,
                                 samplerResID,
                                 scene,
                                 samplers);
        break;
      default:SLog(EError, "Unknow volume rendering technique");
      }

    }

    film->setBitmap(m_bitmap);
    queue->signalRefresh(job);
  }

  void volumePhotonPlanePass(int it, RenderQueue *queue, const RenderJob *job,
                             Film *film, int sceneResID, int sensorResID, int samplerResID, Scene *scene,
                             std::vector<Sampler *> &samplers) {
    ref<Scheduler> sched = Scheduler::getInstance();

    const Float beamInitSize = m_smokeAABB.getBSphere().radius * globalScaleVolume * POURCENTAGE_BS;
    Log(EInfo, "Initial size of the beams: %f", beamInitSize);

    // FIXME: Warning, minDepth and maxDepth are not compatible at all
    // with photon beams.

    // Reduce the number of beam to a small amount
    ref<GatherPhotonBeamProcess> proc = new GatherPhotonBeamProcess(
        m_photonVolumeCount,
        m_granularity,
        m_maxDepth == -1 ? -1 : m_maxDepth - 2,
        m_rrDepth,
        beamInitSize,   //beam initial radius
        true,
        m_autoCancelGathering, job,
        m_adjointCompensation,
        m_cameraSphere,
        std::min(m_minDepth, m_minDepth - 2));
    // Change the minDepth to render the same
    // minDepth of the other rendering techniques

    // Launch the computation
    proc->bindResource("scene", sceneResID);
    proc->bindResource("sensor", sensorResID);
    proc->bindResource("sampler", samplerResID);
    sched->schedule(proc);
    sched->wait(proc);

    Log(EInfo, "Number of volume photon/beam skip: %i", proc->nbSkip());

    Log(EDebug, "Beam map full. Shot "
        SIZE_T_FMT
        " particles, excess photons due to parallelism: "
        SIZE_T_FMT, proc->getShotParticles(), proc->getExcessBeams());

    // Just contains all plane and beams ...
    BeamMap<PhotonBeam> *beamMap = proc->getBeamMap();
    auto beams = beamMap->getBeams();
    std::vector<PhotonPlane> planes(beams.size());
    for (size_t i = 0; i < beams.size(); i++) {
      planes[i] = PhotonPlane::transformBeam(beams[i].second, samplers[0]);
    }

    // Construct the photon plane BVH
    ref<PhotonPlaneBVH<PhotonPlane>> planeBVH = new PhotonPlaneBVH<PhotonPlane>(planes);
    planeBVH->topNodeAABBDump();

    m_totalEmittedVol += proc->getShotParticles();
    m_totalPhotonsVol += beamMap->size();

    Log(EInfo, "Gathering...");

    size_t nCores = sched->getCoreCount();
    BlockScheduler blockSched(m_gatherBlocks.size(), nCores);
    blockSched.run([&](int blockIdx, int tid) {

      std::vector<GatherPoint> &gatherPoints = m_gatherBlocks[blockIdx];

      // Get a sampler to sample a random dist for the step marching
      Sampler *sampler = samplers[tid];

      Spectrum *target = (Spectrum *) m_bitmap->getUInt8Data();
      for (size_t i = 0; i < gatherPoints.size(); ++i) {
        GatherPoint &gp = gatherPoints[i];

        // The volume radiance gathered during this iteration
        Spectrum fluxVolIter(0.f);

        // Volume estimation
        if (!gp.beams.empty()) {
          // For each beams, need to do
          for (Beam &beam: gp.beams) {
            if (beam.isInvalid()) {
              SLog(EError, "Invalid beam");
            }

            Vector d = beam.p2 - beam.p1;
            Float distTotal = d.length();
            d /= distTotal;

            // Query the BRE by casting a ray inside it
            Ray ray(beam.p1, d, Epsilon, distTotal - Epsilon, 0.f);
            PhotonPlaneQuery qPlane(ray, beam.medium,
                                    m_maxDepth == -1 ? -1 : m_maxDepth - beam.depth,
                                    std::max(0, m_minDepth - beam.depth),
                                    sampler,
                                    m_volTechnique);
#if 0
            for(int idPlane = 0; idPlane < planes.size(); idPlane++) {
                qPlane(&planes[idPlane]);
            }
#else
            planeBVH->query<PhotonPlaneQuery>(qPlane);
#endif
            fluxVolIter += qPlane.Li * beam.weight;
          }
        } else {
          // For rendering participating media, we always a bounding sphere for the scene.
          // If there is no such a sphere, or the ray does not hit the sphere, or ray intersection tmax is too small,
          // this could result in no beams and thus black patches in the image.
          SLog(EInfo, "No primary hit point. Check bounding sphere, ray tmax, or far plane.");
        }

        // Normalization of this pass
        fluxVolIter /= proc->getShotParticles();

        // In APA, we always average the rendering
        // even if there is no photon collected
        gp.fluxVol = (gp.fluxVol * (it - 1) + fluxVolIter) / it; // APA estimator
        //gp.fluxVol = fluxVolIter;

        // No normalization about how many photon generated
        // It already carried by the scaling factor inside the photon map
        target[gp.pos.y * m_bitmap->getWidth() + gp.pos.x] += gp.fluxVol;
      }
    });

    // Scale the gathering size for the next iteration
    scaleVolumeAPA(it);

    // Free the memory
    delete beamMap;
  }

  void volumePhotonBeamPass(int it, RenderQueue *queue, const RenderJob *job,
                            Film *film, int sceneResID, int sensorResID, int samplerResID, Scene *scene,
                            std::vector<Sampler *> &samplers) {
    ref<Scheduler> sched = Scheduler::getInstance();

    const Float beamInitSize = m_smokeAABB.getBSphere().radius * globalScaleVolume * POURCENTAGE_BS;
    Log(EInfo, "Initial size of the beams: %f", beamInitSize);

    // Reduce the number of beam to a small amount
    ref<GatherPhotonBeamProcess> proc = new GatherPhotonBeamProcess(
        m_photonVolumeCount,
        m_granularity,
        m_maxDepth == -1 ? -1 : m_maxDepth - 1,
        m_rrDepth,
        beamInitSize,   //beam initial radius
        true,
        m_autoCancelGathering, job,
        m_adjointCompensation,
        m_cameraSphere,
        m_minDepth);

    // Launch the computation
    proc->bindResource("scene", sceneResID);
    proc->bindResource("sensor", sensorResID);
    proc->bindResource("sampler", samplerResID);
    sched->schedule(proc);
    sched->wait(proc);

    Log(EInfo, "Number of volume photon/beam skip: %i", proc->nbSkip());

    BeamMap<PhotonBeam> *beamMap = proc->getBeamMap();
    // If it is requested, generate long beams
    if (m_convertLong) {
      Log(EInfo, "Converting short beams into long one");
      beamMap->convertToLongBeams(scene);
    }

    // Build the acceleration structure
    beamMap->build(BeamMap<PhotonBeam>::EBVHAccel, m_deterministic);

    Log(EDebug, "Beam map full. Shot "
        SIZE_T_FMT
        " particles, excess photons due to parallelism: "
        SIZE_T_FMT, proc->getShotParticles(), proc->getExcessBeams());

    Log(EDebug, "Number of beam stored in the Beam map "
        SIZE_T_FMT,
        beamMap->size());

    m_totalEmittedVol += proc->getShotParticles();
    m_totalPhotonsVol += beamMap->size();

    size_t nCores = sched->getCoreCount();
    BlockScheduler blockSched(m_gatherBlocks.size(), nCores);
    blockSched.run([&](int blockIdx, int tid) {

      std::vector<GatherPoint> &gatherPoints = m_gatherBlocks[blockIdx];

      // Get a sampler to sample a random dist for the step marching
      Sampler *sampler = samplers[tid];

      Spectrum *target = (Spectrum *) m_bitmap->getUInt8Data();
      for (size_t i = 0; i < gatherPoints.size(); ++i) {
        GatherPoint &gp = gatherPoints[i];

        // The volume radiance gathered during this iteration
        Spectrum fluxVolIter(0.f);

        // Volume estimation
        if (!gp.beams.empty()) {
          // For each beams, need to do
          for (Beam &beam: gp.beams) {
            if(m_minCameraDepth > beam.depth)
              continue;
            if(m_maxCameraDepth != -1 && m_maxCameraDepth < beam.depth)
              continue;

            if (beam.isInvalid()) {
              SLog(EError, "Invalid beam");
            }

            Vector d = beam.p2 - beam.p1;
            Float distTotal = d.length();
            d /= distTotal;

            // Query the BRE by casting a ray inside it
            Ray ray(beam.p1, d, Epsilon, distTotal - Epsilon, 0.f);
            BeamRadianceQuery<PhotonBeam> bRadQuery(ray, beam.medium,
                                                    m_maxDepth == -1 ? -1 : m_maxDepth - beam.depth,
                                                    std::max(0, m_minDepth - beam.depth),
                                                    sampler,
                                                    m_volTechnique);
            beamMap->query(bRadQuery);
            fluxVolIter += bRadQuery.Li * beam.weight;
          }
        }

        // Normalization of this pass
        fluxVolIter /= proc->getShotParticles();

        // In APA, we always average the rendering
        // even if there is no photon collected
        gp.fluxVol = (gp.fluxVol * (it - 1) + fluxVolIter) / it; // APA estimator

        // No normalization about how many photon generated
        // It already carried by the scaling factor inside the photon map
        target[gp.pos.y * m_bitmap->getWidth() + gp.pos.x] += gp.fluxVol;
      }
    });

    // Scale the gathering size for the next iteration
    scaleVolumeAPA(it);

    // Free the memory
    delete beamMap;
  }

  void volumePhotonPassBRE(int it, RenderQueue *queue, const RenderJob *job,
                           Film *film, int sceneResID, int sensorResID, int samplerResID, Scene *scene,
                           std::vector<Sampler *> &samplers) {
    {

    }

    ref<Timer> photonShootingTime = new Timer;
    ref<Scheduler> sched = Scheduler::getInstance();

    /* Generate the global photon map */
    ref<GatherPhotonCustomProcess> proc = new GatherPhotonCustomProcess(
        GatherPhotonCustomProcess::EVolumePhotons, m_photonVolumeCount,
        m_granularity, m_maxDepth == -1 ? -1 : m_maxDepth - 1, m_rrDepth, true,
        m_autoCancelGathering, job, m_adjointCompensation, m_cameraSphere, m_minDepth);

    proc->bindResource("scene", sceneResID);
    proc->bindResource("sensor", sensorResID);
    proc->bindResource("sampler", samplerResID);

    sched->schedule(proc);
    sched->wait(proc);
    Log(EInfo, "Volume Shooting time: %i ms", photonShootingTime->getMilliseconds());

    Log(EInfo, "Number of volume photon/beam skip: %i", proc->nbSkip());

    // Build the first acceleration structure
    ref<PhotonMap> photonMap = proc->getPhotonMap();
    photonMap->build();
    Log(EDebug, "Photon map full. Shot "
        SIZE_T_FMT
        " particles, excess photons due to parallelism: "
        SIZE_T_FMT, proc->getShotParticles(), proc->getExcessPhotons());

    m_totalEmittedVol += proc->getShotParticles();
    m_totalPhotonsVol += photonMap->size();

    // Build the second acceleration structure
    Log(EInfo, "Build BRE ...");
    photonMap->setScaleFactor(1 / (Float) proc->getShotParticles());

    // In BRE, the radius managed is directly the reduced size
    const Float breInitSize = m_smokeAABB.getBSphere().radius * globalScaleVolume * POURCENTAGE_BS;
    Log(EInfo, "BRE photon size: %f", breInitSize);
    ref<BeamRadianceEstimator> bre = new BeamRadianceEstimator(photonMap, 120, breInitSize, true);

    // Do the gathering
    ref<Timer> GatheringTime = new Timer;

    size_t nCores = sched->getCoreCount();
    BlockScheduler blockSched(m_gatherBlocks.size(), nCores);
    blockSched.run([&](int blockIdx, int tid) {

      std::vector<GatherPoint> &gatherPoints = m_gatherBlocks[blockIdx];

      // Get a sampler to sample a random dist for the step marching
      Sampler *sampler = samplers[tid];

      Spectrum *target = (Spectrum *) m_bitmap->getUInt8Data();
      for (size_t i = 0; i < gatherPoints.size(); ++i) {
        GatherPoint &gp = gatherPoints[i];
        Spectrum contrib(0.0f);

        // The volume radiance gathered during this iteration
        Spectrum fluxVolIter(0.f);

        // Volume estimation
        if (!gp.beams.empty()) {
          // For each beams, need to do
          for (size_t idBeam = 0; idBeam < gp.beams.size();
               idBeam++) {
            Beam &beam = gp.beams[idBeam];
            if (m_minCameraDepth > beam.depth) {
              continue;
            }
            if (m_maxCameraDepth != -1 && beam.depth > m_maxCameraDepth) {
              break; // Skip this path
            }

            if (beam.isInvalid()) {
              SLog(EError, "Invalid beam");
            }

            if (beam.depth < m_minCameraDepth) {

            }

            // Compute the direction and the total distance
            Vector d = beam.p2 - beam.p1;
            Float distTotal = d.length();
            d /= distTotal;

            // Query the BRE by casting a ray inside it
            Ray ray(beam.p1, d, Epsilon, distTotal - Epsilon, 0.f);
            fluxVolIter += bre->query(ray,
                                      beam.medium,
                                      (m_maxDepth == -1 ? -1 : m_maxDepth - beam.depth),
                                      EVolumeTechniqueHelper::use3DKernel(m_volTechnique),
                                      sampler) * beam.weight;
          }
        }

        // In APA, we always average the rendering
        // even if there is no photon collected
        gp.fluxVol = (gp.fluxVol * (it - 1) + fluxVolIter) / it; // APA estimator

        // No normalization about how many photon generated
        // It already carried by the scaling factor inside the photon map
        contrib = gp.fluxVol;
        target[gp.pos.y * m_bitmap->getWidth() + gp.pos.x] += contrib;
      }
    });

    // In APA, we just have a global scaling factor
    // Update it after each iteration
    scaleVolumeAPA(it);

    Log(EInfo, "Volume Gathering time: %i ms", GatheringTime->getMilliseconds());

  }

  void volumePhotonPassDistance(int it, RenderQueue *queue, const RenderJob *job,
                                Film *film, int sceneResID, int sensorResID, int samplerResID, Scene *scene,
                                std::vector<Sampler *> &samplers) {
    ref<Timer> photonShootingTime = new Timer;
    const Float BBPourcentageCONST = m_smokeAABB.getBSphere().radius * POURCENTAGE_BS;

    ref<Scheduler> sched = Scheduler::getInstance();

    /* Generate the global photon map */
    ref<GatherPhotonCustomProcess> proc = new GatherPhotonCustomProcess(
        GatherPhotonCustomProcess::EVolumePhotons, m_photonVolumeCount,
        m_granularity, m_maxDepth == -1 ? -1 : m_maxDepth - 1, m_rrDepth, true,
        m_autoCancelGathering, job, m_adjointCompensation, m_cameraSphere, m_minDepth);

    proc->bindResource("scene", sceneResID);
    proc->bindResource("sensor", sensorResID);
    proc->bindResource("sampler", samplerResID);

    sched->schedule(proc);
    sched->wait(proc);
    Log(EInfo, "Volume Shooting time: %i ms", photonShootingTime->getMilliseconds());

    Log(EInfo, "Number of volume photon/beam skip: %i", proc->nbSkip());

    // Build the photon map
    ref<PhotonMap> photonMap = proc->getPhotonMap();
    photonMap->build();
    Log(EDebug, "Photon map full. Shot "
        SIZE_T_FMT
        " particles, excess photons due to parallelism: "
        SIZE_T_FMT, proc->getShotParticles(), proc->getExcessPhotons());

    // Do the gathering
    ref<Timer> GatheringTime = new Timer;
    Log(EInfo, "Gathering Volume ...");
    m_totalEmittedVol += proc->getShotParticles();
    m_totalPhotonsVol += photonMap->size();
    film->clear();

    size_t nCores = sched->getCoreCount();
    BlockScheduler blockSched(m_gatherBlocks.size(), nCores);
    blockSched.run([&](int blockIdx, int tid) {

      std::vector<GatherPoint> &gatherPoints = m_gatherBlocks[blockIdx];

      // Get a sampler to sample a random dist for the step martching
      Sampler *sampler = samplers[tid];

      Spectrum *target = (Spectrum *) m_bitmap->getUInt8Data();
      for (size_t i = 0; i < gatherPoints.size(); ++i) {
        GatherPoint &gp = gatherPoints[i];
        Spectrum contrib(0.0f);

        // Volume estimation
        if (!gp.beams.empty()) {
          // Generate the CDF
          Float totalBeamDist = 0.0;
          int nbBeams = gp.beams.size();
          DiscreteDistribution selBeam(nbBeams);
          for (int i = 0; i < nbBeams; i++) {
            if (gp.beams[i].isInvalid()) {
              SLog(EError, "Invalid beam");
            }
            selBeam.append(gp.beams[i].weight.max());

            Vector d = gp.beams[i].p2 - gp.beams[i].p1;
            totalBeamDist += d.length();

          }

          // Make the same size as the ray marching (twice bigger by default)
          const Float querySize = BBPourcentageCONST * gp.scaleVol;
          const Float kernelVol = (4.0 / 3.0) * M_PI * std::pow(querySize, 3);

          // No participating media section
          // go to the next GP
          if (totalBeamDist == 0.f) {
            continue;
          }

          selBeam.normalize();

          gp.fluxVol *= kernelVol;
          const Float MCNorm = 1.f / m_nbCameraSamples;
          Float MVol = 0.f;
          for (int k = 0; k < m_nbCameraSamples; k++) {
            Float randSample = sampler->next1D();
            int beamIndex = selBeam.sampleReuse(randSample);
            const Beam &beam = gp.beams[beamIndex];
            if(m_minCameraDepth > beam.depth)
              continue;
            if(m_maxCameraDepth != -1 && m_maxCameraDepth < beam.depth)
              continue;

            Vector d = beam.p2 - beam.p1;
            Float distTotal = d.length();
            d /= distTotal;
            // Sample an distance inside the beam
            MediumSamplingRecord mRec;
            Ray ray(beam.p1, d, Epsilon, distTotal, 0.f);
            if (beam.medium->sampleDistance(ray, mRec, sampler, EDistanceAlwaysValid, randSample)) {
              Spectrum resVolRad(0.f);
              MVol += photonMap->estimateVolumeRadiance(mRec, ray.d, querySize, resVolRad,
                                                        m_maxDepth == -1 ? INT_MAX : m_maxDepth
                                                            - beam.depth);
              // Add the contribution
              gp.fluxVol += MCNorm * resVolRad * beam.weight * mRec.transmittance
                  / (mRec.pdfSuccess * selBeam[beamIndex]);
            }
          }
          gp.fluxVol /= kernelVol;

          // Do only this at the final step
          if (MVol + gp.NVol > 0) {
            Float ratioVol = (gp.NVol + m_alpha * MVol) / (gp.NVol + MVol);
            gp.scaleVol = gp.scaleVol * std::cbrt(ratioVol);
            gp.NVol = gp.NVol + m_alpha * MVol;
          }
        }

        // When depth == -1, gp.flux is not rescaled and is still radiance value. No need to divide radius.
        contrib = gp.fluxVol / ((Float) m_totalEmittedVol);
        target[gp.pos.y * m_bitmap->getWidth() + gp.pos.x] += contrib;
      }
    });
    Log(EInfo, "Volume Gathering time: %i ms", GatheringTime->getMilliseconds());
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

  void computeGradientFinite(const Vector2i &cropSize, Bitmap *gXBitmap, Bitmap *gYBitmap, bool useAbs = true) {
    Spectrum *targetGX = (Spectrum *) gXBitmap->getUInt8Data();
    Spectrum *targetGY = (Spectrum *) gYBitmap->getUInt8Data();
    for (int y = 1; y < cropSize.y - 1; ++y) { //FIXME: Check the border conditions !!!!
      for (int x = 1; x < cropSize.x - 1; ++x) {
        GatherPoint *curr = m_imgGP[x][y];
        GatherPoint *right = m_imgGP[x + 1][y];
        GatherPoint *top = m_imgGP[x][y + 1];

        Spectrum gX(0.f);
        Spectrum gY(0.f);

        if (m_surfaceRendering) {
          gX += (right->flux - curr->flux) / ((Float) m_totalEmitted);
          gY += (top->flux - curr->flux) / ((Float) m_totalEmitted);
        }
        if (m_volumeRendering) {
          if (!EVolumeTechniqueHelper::useAPA(m_volTechnique)) {
            gX += (right->fluxVol - curr->fluxVol) / ((Float) m_totalEmittedVol);
            gY += (top->fluxVol - curr->fluxVol) / ((Float) m_totalEmittedVol);
          } else {
            gX += (right->fluxVol - curr->fluxVol);
            gY += (top->fluxVol - curr->fluxVol);
          }
        }

        if (useAbs) {
          targetGX[y * m_bitmap->getWidth() + x] = gX.abs();
          targetGY[y * m_bitmap->getWidth() + x] = gY.abs();
        } else {
          targetGX[y * m_bitmap->getWidth() + x] = gX;
          targetGY[y * m_bitmap->getWidth() + x] = gY;
        }
      }
    }
  }

  std::string toString() const {
    std::ostringstream oss;
    oss << "SPPMIntegrator[" << endl
        << "  maxDepth = " << m_maxDepth << "," << endl
        << "  rrDepth = " << m_rrDepth << "," << endl
        << "  alpha = " << m_alpha << "," << endl
        << "  photonCount (surface) = " << m_photonSurfaceCount << "," << endl
        << "  photonCount (volume) = " << m_photonVolumeCount << "," << endl
        << "  granularity = " << m_granularity << "," << endl
        << "  maxPasses = " << m_maxPasses << endl
        << "]";
    return oss.str();
  }

  MTS_DECLARE_CLASS()
private:
  std::vector<std::vector<GatherPoint> > m_gatherBlocks;
  RadiusInitializer *m_gpManager;

  std::vector<Point2i> m_offset;
  ref<Mutex> m_mutex;
  ref<Bitmap> m_bitmap;
  Float m_initialScale, m_initialScaleVolume;
  Float m_alpha;
  int m_photonSurfaceCount, m_photonVolumeCount, m_granularity;
  int m_maxDepth, m_rrDepth, m_minDepth;
  size_t m_totalEmitted, m_totalPhotons;
  size_t m_totalEmittedVol, m_totalPhotonsVol;
  bool m_running;
  bool m_autoCancelGathering;
  int m_maxPasses;

  Float globalScaleVolume;
  bool m_independentScale;

  // A structure to easy gather neighbors pixels
  std::vector<std::vector<GatherPoint *>> m_imgGP;
  int m_dumpIteration;
  int m_nbInternalIteration;
  bool m_computeGrad;
  int m_maxRenderingTime;
  bool m_directTracing;

  bool m_surfaceRendering;
  bool m_volumeRendering;

  EVolumeTechnique m_volTechnique;
  bool m_deterministic;

  bool m_convertLong; /// Convert to long beams
  bool m_nanCheck;

  size_t m_minCameraDepth;
  int m_maxCameraDepth;

  bool m_adjointCompensation;
  Float m_cameraSphere;

  // Distance sampling
  int m_nbCameraSamples;

  // AABB of the smoke
  AABB m_smokeAABB;
  std::string m_forceAPA = "";
};

MTS_IMPLEMENT_CLASS_S(SPPMIntegrator, false, Integrator)

MTS_EXPORT_PLUGIN(SPPMIntegrator, "Stochastic progressive photon mapper");
MTS_NAMESPACE_END
