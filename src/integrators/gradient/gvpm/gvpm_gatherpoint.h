#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/scene.h>

#include "../../photonmapper/utilities/block_sched.h"
#include "gvpm_geoOps.h"

MTS_NAMESPACE_BEGIN

/**
 * Attempt to trace a path that ends on the first diffuse vertex.
 *
 * When maximum path length is reached, the last vertex might or might not be diffused.
 *
 * As this is used for generating gather point, any participating medium encountered is skipped.
 * They will be handled during photon gathering.
 *
 * \param nSteps    Maximum vertex ID allowed (= expected depth + 1).
 * \return          Vertex ID of the last vertex (= traced depth + 1).
 */
inline int randomWalkFromPixelToFirstDiffuse(Path &path, const Scene *scene, Sampler *sampler,
                                             int nSteps, const Point2i &pixelPosition, MemoryPool &pool,
                                             Float &pdfComponent) {

  PathVertex *v1 = pool.allocVertex(), *v2 = pool.allocVertex();
  PathEdge *e0 = pool.allocEdge(), *e1 = pool.allocEdge();

  /* Use a special sampling routine for the first two sensor vertices so that
  the resulting subpath passes through the specified pixel position */
  int number_intersection = path.vertex(0)->sampleSensor(scene,
                                                         sampler, pixelPosition, e0, v1, e1, v2, true,
                                                         Point2(-1), true); // Use long beams

  // number_intersection is latest vertex ID
  if (number_intersection < 1) {
    pool.release(e0);
    pool.release(v1);
    pool.release(e1);
    pool.release(v2);
    return 0;
  }
  path.append(e0, v1);

  if (number_intersection < 2) {
    pool.release(e1);
    pool.release(v2);
    return 1;
  }
  path.append(e1, v2);

  PathVertex *predVertex = v1, *curVertex = v2;
  PathEdge *predEdge = e1;

  // Compute the throughput
  Spectrum throughput(1.0f);
  throughput *= v1->rrWeight * v1->weight[ERadiance] * e1->weight[ERadiance];

  for (; number_intersection < nSteps || nSteps == -1; ++number_intersection) {

    int componentSelected = -1;
    if (curVertex->isSurfaceInteraction()) {

      const BSDF *bsdf = curVertex->getIntersection().getBSDF();
      pdfComponent = 1.f;

      // Hit emitter with null BSDF, break
      if (curVertex->getIntersection().isEmitter()) {
        return -2; // Invalid GP
      }

      if (bsdf->getComponentCount() == 0) {
        // Case that we hit a medium interface
        // So we continue until we hit nothing...
        curVertex->sampledComponentIndex = 0;
      } else {
        // Select one of the BSDF component randomly
        Point2 randomBSDFSample = sampler->next2D();
        BSDFSamplingRecord bRec(curVertex->getIntersection(), sampler);
        componentSelected = bsdf->sampleComponent(bRec, pdfComponent,
                                                  randomBSDFSample,
                                                  VertexClassifier::roughnessThreshold);

        // Check its roughness
        int componentQuery = (componentSelected == -1 ? 0 : componentSelected);
        bool roughFunction = bsdf->getRoughness(curVertex->getIntersection(),
                                                componentQuery) > VertexClassifier::roughnessThreshold;
        //roughFunction &= !(curVertex->componentType & BSDF::EDelta);
        if (roughFunction) {
          // Take one component, if it is enough smooth,
          // Stop here and store the sampledComponentIndex
          // We will use this index later to decide the vertex classification
          curVertex->sampledComponentIndex = componentSelected;
          return number_intersection;
        }
      }
    } else if (curVertex->isMediumInteraction()) {
      return -1;
    } else {
      SLog(EError, "Unknow type");
      return -1;
    }

    PathVertex *succVertex = pool.allocVertex();
    PathEdge *succEdge = pool.allocEdge();

    // Use long beam sampling, all photons stay at surfaces.
    // The edge weight is equal to the transmittance.
    PathVertex::SamplingOption option;
    option.russianRoulette = number_intersection > 10; // Only do RR in this case
    option.componentSel = componentSelected;
    option.roughness = VertexClassifier::roughnessThreshold;
    option.longBeam = true;

    if (!curVertex->sampleNext(scene,
                               sampler,
                               predVertex,
                               predEdge,
                               succEdge,
                               succVertex,
                               ERadiance,
                               option,
                               &throughput)) {
      pool.release(succVertex);
      pool.release(succEdge);
      return -1; // Invalid GP
    }

    path.append(succEdge, succVertex);

    predVertex = curVertex;
    curVertex = succVertex;
    predEdge = succEdge;
  }

  // If we arrive here, it's mean that we have reach the max steps
  if (curVertex->isSurfaceInteraction()) {
    const BSDF *bsdf = curVertex->getIntersection().getBSDF();
    pdfComponent = 1.f;

    // Hit emitter with null BSDF, break
    if (curVertex->getIntersection().isEmitter()) {
      return -2;
    }

    // Case that we hit a medium interface
    if (bsdf->getComponentCount() == 0) {
      curVertex->sampledComponentIndex = 0;
      return number_intersection;
    }

    Point2 randomBSDFSample = sampler->next2D();
    BSDFSamplingRecord bRec(curVertex->getIntersection(), sampler);
    int componentSelected = bsdf->sampleComponent(bRec, pdfComponent,
                                                  randomBSDFSample,
                                                  VertexClassifier::roughnessThreshold);
    // Check its roughness
    int componentQuery = (componentSelected == -1 ? 0 : componentSelected);
    bool roughFunction = bsdf->getRoughness(curVertex->getIntersection(),
                                            componentQuery) > VertexClassifier::roughnessThreshold;
    //roughFunction &= !(curVertex->componentType & BSDF::EDelta);
    if (roughFunction) {
      // Take one component, if it is enough smooth,
      // Stop here and store the sampledComponentIndex
      // We will use this index later to decide the vertex classification
      curVertex->sampledComponentIndex = componentSelected;
      return number_intersection;
    } else {
      return -1;
    }
  } else {
    return -1;
  }
}

class GVPMRadiusInitializer {
public:
  GVPMRadiusInitializer()
      : scene(0),
        m_maxDepth(-1) 
  {    
  }

  virtual ~GVPMRadiusInitializer() = default;

  virtual void init(const GPMConfig &config, Scene *scene, 
                    std::vector<std::vector<GatherPoint>> &gatherPoints,
                    const std::vector<Point2i> &offset) {
    this->m_config = config;
    SLog(EInfo, "Initialize GP generation object");

    // Copy the object
    // These will be used to generate
    // the gather points through the scene
    this->scene = scene;
    this->m_maxDepth = config.maxDepth;
    this->gatherBlocks = &gatherPoints;
    this->offset = &offset;

    // And then initialise the samplers
    // without shift in the random
    // sequence
    generateSamplers(0);

    m_pools.resize(offset.size());
  }

  void generateSamplers(int shiftRandomNumber) {
    // Initialize samplers for generate GP
    // Each pixel blocks will have is own
    // sampler to generate the same set
    // of gather points

    Properties props("independent");
    if (m_config.referenceMod) {
      props.setBoolean("randInit", true);
    }

    // --- Sampler to generate seed of other
    ref<Sampler> samplerIndependent =
        static_cast<Sampler *>(PluginManager::getInstance()->createObject(
            MTS_CLASS(Sampler), props));

    if (shiftRandomNumber != 0) {
      SLog(EInfo, "Make an shift of random number generator: %i (%i)",
           shiftRandomNumber, offset->size());
      // Make the shift by advancing in the
      // sequence by calling the
      for (size_t i = 0; i < (offset->size() * shiftRandomNumber); ++i) {
        samplerIndependent->next2D();
      }
    }

    // --- Create all samplers
    samplers.resize(offset->size());  //< Create sampler
    for (size_t i = 0; i < offset->size(); ++i) {
      ref<Sampler> clonedIndepSampler = samplerIndependent->clone();
      samplers[i] = clonedIndepSampler;
    }
  }

  virtual Float getRadiusRayDifferential(RayDifferential &ray, Float totalDist) {

    if (ray.hasDifferentials) {  // nbComponentExtra == 0 &&
      ray.scaleDifferential(1.f);  // Scale to one (TODO: Is it necessary ??? )
      Point posProj = ray.o + ray.d * totalDist;
      Point rX = ray.rxOrigin + ray.rxDirection * totalDist;
      Point rY = ray.ryOrigin + ray.ryDirection * totalDist;
      Float dX = (rX - posProj).length();
      Float dY = (rY - posProj).length();

      Float r = std::max(dX, dY);
      if (r > 100) {
        SLog(EError, "Infinite radius %f", r);
      }
      return r;
    } else {
      SLog(EError, "No ray differential");
      return 0.f;
    }
  }

  void regeneratePositionAndRadius() {
    ref<Timer> timer = new Timer();
    ref<Sensor> sensor = scene->getSensor();

    bool needsApertureSample = sensor->needsApertureSample();
    bool needsTimeSample = sensor->needsTimeSample();
    ref<Film> film = sensor->getFilm();
    Vector2i cropSize = film->getCropSize();
    int blockSize = scene->getBlockSize();

    ref<Scheduler> sched = Scheduler::getInstance();
    size_t nCores = sched->getCoreCount();
    BlockScheduler blockSched(gatherBlocks->size(), nCores);
    blockSched.run([&](int blockIdx, int tid) {
      // Get the sampler associated
      // to the block
      ref<Sampler> sampler = samplers[blockIdx];
      MemoryPool &pool = m_pools[blockIdx];

      // For all the gather points int the block
      // be carefull of the offset of the image plane
      std::vector<GatherPoint> &gb = (*gatherBlocks)[blockIdx];
      int xofs = (*offset)[blockIdx].x, yofs = (*offset)[blockIdx].y;
      int index = 0;
      for (int yofsInt = 0; yofsInt < blockSize; ++yofsInt) {
        if (yofsInt + yofs >= cropSize.y)
          continue;
        for (int xofsInt = 0; xofsInt < blockSize; ++xofsInt) {
          if (xofsInt + xofs >= cropSize.x)
            continue;

          // Get the gather point and clear it
          // (temp data erased) + prepare data
          // to sample the gp position
          GatherPoint &gatherPoint = (gb[index++]);

          // Initialize the GP plane position
          gatherPoint.pixel = Point2i(xofs + xofsInt, yofs + yofsInt);
          sampler->generate(gatherPoint.pixel);

          // Invalid - for now
          gatherPoint.clear();

          // Initialize the camera path
          gatherPoint.path.initialize(scene, 0.f, ERadiance, pool);
          Float pdfComponent = 1.f;

          // Generate the camera path
          int nSteps = (m_maxDepth - 1) + 1;   // nSteps parameter
          int lastVertexID = randomWalkFromPixelToFirstDiffuse(gatherPoint.path,
                                                               scene,
                                                               sampler,
                                                               nSteps,
                                                               gatherPoint.pixel,
                                                               pool,
                                                               pdfComponent);

          for (size_t k = 2; k < gatherPoint.path.vertexCount() - 1; ++k) {
            // Check that there is only surface interactions
            // and no volume interactions
            if (gatherPoint.path.vertex(k)->getType() == PathVertex::EMediumInteraction) {
              SLog(EError, "There is a medium interaction catched during the camera"
                           "path generation");
            }

            // Relax the pureSpecular check by allowing high glossy path.
            // Note that this is a double-edged sword. It can leave some highlight regions as is (no reconstruction)
            // if bounceRoughness is set too high.
            bool diffuse = !gatherPoint.path.vertex(k)->degenerate;
            if (m_config.nearSpecularDirectImage) {
              diffuse = (VertexClassifier::type(*gatherPoint.path.vertex(k),
                                                gatherPoint.path.vertex(k)->sampledComponentIndex)
                  == VERTEX_TYPE_DIFFUSE);
            }
            if (diffuse) {
              gatherPoint.pureSpecular = false;
              break;
            }
          }

          // Compute PDF, weights for each vertex of the camera path
          gatherPoint.generateVertexInfo();

          // Emitter hit: we assume that emitter always has a NullBSDF.
          // The tracing will stop when a NullBSDF is encountered.
          if (lastVertexID == -2) {
            gatherPoint.depth = -1;

            const Path &pt = gatherPoint.path;
            if (m_config.directTracing && pt.vertexCount() > 2) {
              // Complete the contribution with light emission
              const PathVertex *last = pt.vertex((int) pt.vertexCount() - 1);
              if (last->getType() == PathVertex::ESurfaceInteraction) {
                const Intersection &itsLocal = last->getIntersection();
                if (itsLocal.isEmitter()) {
                  int depth = pt.vertexCount() - 2;
                  if (depth >= m_config.minDepth) {
                    Spectrum emission = gatherPoint.surfInfo().weight
                        * itsLocal.Le(pt.edge(pt.vertexCount() - 2)->d);
                    if (!gatherPoint.pureSpecular) {
                      gatherPoint.currEmission += emission;
                      gatherPoint.emission += emission;
                    } else {
                      // Separate emission due to pure specular path
                      // Usually, these causes bright edges and hence artifacts in L2 reconstruction
                      // We treat them as direct image and add it after reconstruction.
                      gatherPoint.directEmission += emission;
                    }
                  }
                }
              }
            }

            sampler->advance();
            continue;
          }

          // Sanity check
          PathVertex *last = gatherPoint.path.vertex((int) gatherPoint.path.vertexCount() - 1);

          // Strictly no glossy gather point
          if (last->getType() != PathVertex::ESurfaceInteraction ||
              VertexClassifier::type(*last, last->sampledComponentIndex) != VERTEX_TYPE_DIFFUSE
              ) {

            if (m_config.directTracing && last->getType() != PathVertex::ESurfaceInteraction) {
              // If no intersection, check the envmap
              //gatherPoint.emission += scene->evalEnvironment(ray); // No weight => bounce over scene
            }
            gatherPoint.depth = -1;

          } else {
            // Compute PDF values and weights for all the camera path
            gatherPoint.depth = lastVertexID - 1;

            // Sample the end point one more time to generate the BSDF component type
            // and component index at the gather point
            int componentSelected = last->sampledComponentIndex;
            gatherPoint.sampledComponent = componentSelected;
            gatherPoint.pdfComponent = pdfComponent;
            // pdf of the component will be accounted during gathering

            // Check the PDF for sample an component
            // There is a problem when it's equal to 0
            // Arrives when ray leaking in RoughtPlastic
            if (pdfComponent == 0) {
              if (componentSelected != -1) {
                SLog(EError, "Not the good component is returned for component selection");
              }
              break;
            }

            // Sanity check
            // -1 for component selected = same behavior for all the bsdf component
            // in this case, the probability to pick all need to be 1
            if (componentSelected == -1 && pdfComponent != 1.f) {
              SLog(EError, "All component is selected but the pdf is not 1");
            }


            // Compute the radius using the ray differential if this is the first gather point
            // Note that this radius is rescale by the scale factor of the gather point
            Point2 apertureSample;
            Float timeSample = 0.0f;
            if (needsApertureSample)
              apertureSample = Point2(0.5f);
            if (needsTimeSample)
              timeSample = 0.5f;

            if (gatherPoint.path.vertexCount() < 2)
              SLog(EError, "Problem of initialize the path");
            // Now estimate initial radius
            const Point2 pixel = gatherPoint.path.vertex(1)->getSamplePosition();
            RayDifferential rayCamera;
            sensor->sampleRayDifferential(rayCamera, pixel, apertureSample, timeSample);

            if (gatherPoint.scale == 0.f) {
              SLog(EError, "Zero scale on valid gather point");
            }

            Float traveledDistance = 0.f;
            for (size_t i = 0; i < gatherPoint.path.vertexCount() - 1; ++i) {
              traveledDistance += gatherPoint.path.edge(i)->length;
            }

            // Should not be zero
            Float estimatedRadius = m_config.initialRadius;
            if(estimatedRadius == 0.0) {
              estimatedRadius = getRadiusRayDifferential(rayCamera, traveledDistance);
            }
            estimatedRadius = std::max<Float>(estimatedRadius, m_config.epsilonRadius);
            // Radius of the current pass is the initial radius times the reduction ratio
            gatherPoint.radius = estimatedRadius * gatherPoint.scale;

            // Finish, we have generated a gather point
          }

          sampler->advance();
        }
      }
    });  // End of the sampling

    unsigned int milliseconds = timer->getMilliseconds();
    SLog(EInfo, "Radius lookup time: %i, %i", milliseconds / 1000, milliseconds % 1000);

    //////////////////////////////////
    // Display extra informations
    // And in case not initialized GP
    // gets the maximum radii
    //////////////////////////////////
    // === Find max Size
    Float radiusMax = 0;
    Float radiusMin = std::numeric_limits<Float>::max();
    Float radiusAvg = 0;
    int k = 0;
    for (size_t j = 0; j < gatherBlocks->size(); j++) {
      std::vector<GatherPoint> &gb = (*gatherBlocks)[j];
      for (size_t i = 0; i < gb.size(); ++i) {
        GatherPoint &gp = gb[i];
        radiusMax = std::max(radiusMax, gp.radius);
        radiusMin = std::min(radiusMin, gp.radius);
        radiusAvg += gp.radius;
        k++;
      }
    }
    radiusAvg /= k;
    SLog(EInfo, "Max radius: %f", radiusMax);
    SLog(EInfo, "Min radius: %f", radiusMin);
    SLog(EInfo, "Avg radius: %f", radiusAvg);
  }

  void rescaleFlux() {
    for (size_t j = 0; j < gatherBlocks->size(); j++) {
      std::vector<GatherPoint> &gb = (*gatherBlocks)[j];
      for (size_t i = 0; i < gb.size(); ++i) {
        GatherPoint &gp = gb[i];
        if (gp.depth != -1 && gp.lastIts().isValid()) {
          if (gp.radius == 0) {
            // No radius => Error because we will lose the flux
            SLog(EError, "Valid GP with null radius");
          } else {
            gp.rescaleRadii(gp.radius * gp.radius * M_PI);
          }

        } else {
          // No rescale, gp.flux is still radiance.
        }
      }
    }
  }

  ref<Sampler> getSamplerBlock(int i) {
    return samplers[i];
  }

protected:

  void resetInitialRadius(Float initialRadius) {
    SLog(EInfo, "Reset Initial radius to: %f", initialRadius);
    for (size_t j = 0; j < gatherBlocks->size(); j++) {
      std::vector<GatherPoint> &gb = (*gatherBlocks)[j];

      // === Create & Initialize all gather points in the block
      for (size_t i = 0; i < gb.size(); ++i) {
        GatherPoint &gp = gb[i];
        gp.radius = initialRadius;
      }
    }
  }

protected:
  // === Config attributs
  Scene *scene;
  int m_maxDepth;

  // === Gather blocks attributes
  std::vector<std::vector<GatherPoint>> *gatherBlocks;
  const std::vector<Point2i> *offset;

  // === Sampler attributes
  ref_vector<Sampler> samplers;
  std::vector<MemoryPool> m_pools;

  GPMConfig m_config;
};

MTS_NAMESPACE_END
