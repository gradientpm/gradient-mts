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

#include <mitsuba/render/scene.h>
#include <mitsuba/core/statistics.h>

#include "../volume_utils.h"

MTS_NAMESPACE_BEGIN

// For testing, make possible to disable
// each type of transport (without compromise MIS)
#define MEDIUM_DIRECT_TRANS 1
#define MEDIUM_INTER_TRANS 1
#define SURFACE_DIRECT_TRANS 1
#define SURFACE_INTER_TRANS 1

static StatsCounter avgPathLength("Volumetric path tracer", "Average path length", EAverage);

enum EVertexType {
  VERTEX_TYPE_GLOSSY,      ///< "Specular" vertex that requires the half-vector duplication shift.
  VERTEX_TYPE_DIFFUSE,     ///< "Non-specular" vertex that is rough enough for the reconnection shift.
  VERTEX_TYPE_INVALID
};

struct VertexClassifier {
  static Float roughnessThreshold;
  /**
  * \brief Classify a pathvertex into diffuse or specular.
  *
  * \param bsdfComponentType         The component type (ESmooth or EDelta) being considered.
  *                                  If a BSDF has a smooth and a delta component,
  *                                  the delta will be ignored in the classification
  *                                  if the smooth component is being considered.
  */

  static EVertexType type(const BSDF *bsdf, const Intersection &its, int component) {
    component = component == -1 ? 0 : component;
    if (component >= bsdf->getComponentCount()) {
      return VERTEX_TYPE_INVALID; // This can be a quite bit drastic
    } else {
      return getVertexTypeByRoughness(bsdf->getRoughness(its, component));
    }
  }

  static bool enoughRough(const Intersection &its, const BSDF *bsdf) {
    if (bsdf->getComponentCount() == 0) {
      // No BSDF, just said no
      return false;
    }

    for (int i = 0; i < bsdf->getComponentCount(); i++) {
      if (getVertexTypeByRoughness(bsdf->getRoughness(its, i)) == VERTEX_TYPE_DIFFUSE)
        return true;
    }
    // All the component are "glossy", return false
    return false;
  }

  /// Returns the vertex type of a vertex by its roughness value.
  static EVertexType getVertexTypeByRoughness(Float roughness) {
    if (roughness <= roughnessThreshold) {
      return VERTEX_TYPE_GLOSSY;
    } else {
      return VERTEX_TYPE_DIFFUSE;
    }
  }
};
Float VertexClassifier::roughnessThreshold;

/*!\plugin{volpath}{Extended volumetric path tracer}
 * \order{4}
 * \parameters{
 *     \parameter{maxDepth}{\Integer}{Specifies the longest path depth
 *         in the generated output image (where \code{-1} corresponds to $\infty$).
 *	       A value of \code{1} will only render directly visible light sources.
 *	       \code{2} will lead to single-bounce (direct-only) illumination,
 *	       and so on. \default{\code{-1}}
 *	   }
 *	   \parameter{rrDepth}{\Integer}{Specifies the minimum path depth, after
 *	      which the implementation will start to use the ``russian roulette''
 *	      path termination criterion. \default{\code{5}}
 *	   }
 *     \parameter{strictNormals}{\Boolean}{Be strict about potential
 *        inconsistencies involving shading normals? See
 *        page~\pageref{sec:strictnormals} for details.
 *        \default{no, i.e. \code{false}}
 *     }
 *     \parameter{hideEmitters}{\Boolean}{Hide directly visible emitters?
 *        See page~\pageref{sec:hideemitters} for details.
 *        \default{no, i.e. \code{false}}
 *     }
 * }
 *
 * This plugin provides a volumetric path tracer that can be used to
 * compute approximate solutions of the radiative transfer equation.
 * Its implementation makes use of multiple importance sampling to
 * combine BSDF and phase function sampling with direct illumination
 * sampling strategies. On surfaces, it behaves exactly
 * like the standard path tracer.
 *
 * This integrator has special support for \emph{index-matched} transmission
 * events (i.e. surface scattering events that do not change the direction
 * of light). As a consequence, participating media enclosed by a stencil shape (see
 * \secref{shapes} for details) are rendered considerably more efficiently when this
 * shape has \emph{no}\footnote{this is what signals to Mitsuba that the boundary is
 * index-matched and does not interact with light in any way. Alternatively,
 * the \pluginref{mask} and \pluginref{thindielectric} BSDF can be used to specify
 * index-matched boundaries that involve some amount of interaction.} BSDF assigned
 * to it (as compared to, say, a \pluginref{dielectric} or \pluginref{roughdielectric} BSDF).
 *
 * \remarks{
 *    \item This integrator will generally perform poorly when rendering
 *      participating media that have a different index of refraction compared
 *      to the surrounding medium.
 *    \item This integrator has poor convergence properties when rendering
 *      caustics and similar effects. In this case, \pluginref{bdpt} or
 *      one of the photon mappers may be preferable.
 * }
 */
class VolumetricPathTracer : public MonteCarloIntegrator {
public:
  VolumetricPathTracer(const Properties &props) : MonteCarloIntegrator(props) {
    m_directLightVisibility = props.getBoolean("directLight", true);
    m_minDepth = props.getInteger("minDepth", 0);

    /* This option is specific for GVPM project:
     * The only propose is to compute a very specific light transport when computing volume interaction
     * This option does not have (really) meaning for path tracing but it does for photon mapping-based
     * approaches. The idea is to restrict the light transport only if the path have bounce over
     * classified glossy surface.
     */
    m_minCameraDepth = props.getInteger("minCameraDepth", 0);

    /*
    In photon mapping,
    - For all2surface, gathering is done at the first diffuse vertex on camera path.
    - For all2media, gathering is done at the vertex generated during ray marching, which is always a medium vertex.
    It is important to make sure that the interaction mode check here is the same
    as those used in photon mapping.
    For example,
    - ESML path is a media2media path (assume eye in medium). Basically S can be ignored.
    - EDML or ES*DML path is a media2surface path.
    */
    std::string lightingInteractionMode = props.getString("lightingInteractionMode", "");
    if (lightingInteractionMode == "")
      lightingInteractionMode = props.getString("interactionMode", "all2all");
    m_lightingInteractionMode = parseMediaInteractionMode(lightingInteractionMode);

    /* Setup the roughtness classifier */
    VertexClassifier::roughnessThreshold = props.getFloat("bounceRoughness", 0.05);

    if(m_minCameraDepth > 0 && needRenderSurface(m_lightingInteractionMode)) {
      SLog(EError, "minCameraDepth is only meaningful when only computing volume interactions.");
    }
  }

  /// Unserialize from a binary data stream
  VolumetricPathTracer(Stream *stream, InstanceManager *manager)
      : MonteCarloIntegrator(stream, manager) {}

  bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
                  int sceneResID, int sensorResID, int samplerResID) {
    MonteCarloIntegrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID);

    for (auto media: scene->getMedia()) {
      const Medium *m = media.get();
      bool foundAttachedMesh = false;
      // Check the meshes
      for (auto mesh: scene->getMeshes()) {
        if (mesh->getExteriorMedium() != NULL) {
          foundAttachedMesh = foundAttachedMesh || (mesh->getExteriorMedium() == m);
        }
        if (mesh->getInteriorMedium() != NULL) {
          foundAttachedMesh = foundAttachedMesh || (mesh->getInteriorMedium() == m);
        }
        if (foundAttachedMesh) {
          Log(EInfo, "Medium: %s attached to %s", m->toString().c_str(), mesh->toString().c_str());
          break; // Stop the loop
        }
      }
      if (foundAttachedMesh) {
        continue;// Do not need to check the shapes.
      }
      // Check the shapes
      for (auto mesh: scene->getShapes()) {
        if (mesh->getExteriorMedium() != NULL) {
          foundAttachedMesh = foundAttachedMesh || (mesh->getExteriorMedium() == m);
        }
        if (mesh->getInteriorMedium() != NULL) {
          foundAttachedMesh = foundAttachedMesh || (mesh->getInteriorMedium() == m);
        }
        if (foundAttachedMesh) {
          Log(EInfo, "Medium: %s attached to %s", m->toString().c_str(), mesh->toString().c_str());
          break; // Stop the loop
        }
      }

      if (!foundAttachedMesh) {
        Log(EError, "No shape define interior/exterior limit for the media: %s", m->toString().c_str());
      }
    }

    if (!(m_lightingInteractionMode & EAll2Surf)) {
      Log(EInfo, "Only volume rendering, change probability to sample the medium");
      auto medium = scene->getMedia().begin();
      while (medium != scene->getMedia().end()) {
        const_cast<Medium *>(medium->get())->computeOnlyVolumeInteraction();
        medium++;
      }
    }

    return true;
  }

  Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
    /* Some aliases and local variables */
    const Scene *scene = rRec.scene;
    Intersection &its = rRec.its;
    MediumSamplingRecord mRec;
    RayDifferential ray(r);
    Spectrum Li(0.0f);
    Float eta = 1.0f;

    /* Perform the first ray intersection (or ignore if the
       intersection has already been provided). */
    rRec.rayIntersect(ray);

    // This boolean controls if we need to scatter inside the media first
    // in order to add the path contribution.
    bool mediumScattered = needRenderSurface(m_lightingInteractionMode);

    Spectrum throughput(1.0f);
    bool scattered = false;
    int nbSurfaceInteractions = 0;

    while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
      /* ==================================================================== */
      /*                 Radiative Transfer Equation sampling                 */
      /* ==================================================================== */
      if (rRec.medium && rRec.medium->sampleDistance(Ray(ray, 0, its.t), mRec, rRec.sampler)) {
        /* Sample the integral
         \int_x^y tau(x, x') [ \sigma_s \int_{S^2} \rho(\omega,\omega') L(x,\omega') d\omega' ] dx'
      */
        if(nbSurfaceInteractions < m_minCameraDepth) {
          return Li; // Early kill the path
        }
        mediumScattered = true;
        const PhaseFunction *phase = mRec.getPhaseFunction();

        if (rRec.depth >= m_maxDepth && m_maxDepth != -1) // No more scattering events allowed
          break;

        throughput *= mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;

        /* ==================================================================== */
        /*                          Luminaire sampling                          */
        /* ==================================================================== */

        /* Estimate the single scattering component if this is requested */
        DirectSamplingRecord dRec(mRec.p, mRec.time);
        if (mediumScattered) {
          if (rRec.type & RadianceQueryRecord::EDirectMediumRadiance) {
            int interactions = m_maxDepth - rRec.depth - 1;

            Spectrum value = scene->sampleAttenuatedEmitterDirect(
                dRec, rRec.medium, interactions,
                rRec.nextSample2D(), rRec.sampler);

            if (!value.isZero()) {
              const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

              /* Evaluate the phase function */
              PhaseFunctionSamplingRecord pRec(mRec, -ray.d, dRec.d);
              Float phaseVal = phase->eval(pRec);

              if (phaseVal != 0) {
                /* Calculate prob. of having sampled that direction using
                   phase function sampling */
                Float phasePdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle)
                                 ? phase->pdf(pRec) : (Float) 0.0f;

                /* Weight using the power heuristic */
                const Float weight = miWeight(dRec.pdf, phasePdf);
                if (rRec.depth + 1 >= m_minDepth && nbSurfaceInteractions >= m_minCameraDepth) {
#if MEDIUM_DIRECT_TRANS
                  Li += throughput * value * phaseVal * weight;
#endif
                }
              }
            }
          }
        }

        /* ==================================================================== */
        /*                         Phase function sampling                      */
        /* ==================================================================== */

        Float phasePdf;
        PhaseFunctionSamplingRecord pRec(mRec, -ray.d);
        Float phaseVal = phase->sample(pRec, phasePdf, rRec.sampler);
        if (phaseVal == 0)
          break;
        throughput *= phaseVal;

        /* Trace a ray in this direction */
        ray = Ray(mRec.p, pRec.wo, ray.time);
        ray.mint = 0;

        Spectrum value(0.0f);
        rayIntersectAndLookForEmitter(scene, rRec.sampler, rRec.medium,
                                      m_maxDepth - rRec.depth - 1, ray, its, dRec, value);

        /* If a luminaire was hit, estimate the local illumination and
           weight using the power heuristic */
        if (mediumScattered) {
          if (!value.isZero() && (rRec.type & RadianceQueryRecord::EDirectMediumRadiance)) {
            const Float emitterPdf = scene->pdfEmitterDirect(dRec);
            if (rRec.depth + 1 >= m_minDepth && nbSurfaceInteractions >= m_minCameraDepth) {
#if MEDIUM_INTER_TRANS
              Li += throughput * value * miWeight(phasePdf, emitterPdf);
#endif
            }
          }
        }

        /* ==================================================================== */
        /*                         Multiple scattering                          */
        /* ==================================================================== */

        /* Stop if multiple scattering was not requested */
        if (!(rRec.type & RadianceQueryRecord::EIndirectMediumRadiance))
          break;
        rRec.type = RadianceQueryRecord::ERadianceNoEmission;
      } else {
        /* Sample
            tau(x, y) (Surface integral). This happens with probability mRec.pdfFailure
            Account for this and multiply by the proper per-color-channel transmittance.
        */
        if (rRec.medium)
          throughput *= mRec.transmittance / mRec.pdfFailure;

        if(!mediumScattered) {
          nbSurfaceInteractions += 1;
        }

        if (!its.isValid()) {
          if (!m_directLightVisibility && rRec.depth == 1 || !mediumScattered) { // Direct visibility on the light
            break;
          }
          /* If no intersection could be found, possibly return
             attenuated radiance from a background luminaire */
          if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
              && (!m_hideEmitters || scattered)) {
            Spectrum value = throughput * scene->evalEnvironment(ray);
            if (rRec.medium)
              value *= rRec.medium->evalTransmittance(ray, rRec.sampler);
            if (rRec.depth + 1 >= m_minDepth && nbSurfaceInteractions >= m_minCameraDepth) {
              Li += value;
            }
          }

          break;
        }

        /* Possibly include emitted radiance if requested */
        if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
            && ((!m_hideEmitters || scattered) || mediumScattered)) {
          if (!m_directLightVisibility && rRec.depth == 1) {
            // We do not add the light contribution
          } else {
            if (rRec.depth + 1 >= m_minDepth && nbSurfaceInteractions >= m_minCameraDepth) {
//              SLog(EError, "Impossible !");
              Li += throughput * its.Le(-ray.d);
            }
          }
        }


        /* Include radiance from a subsurface integrator if requested */
        if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
          Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

        if (rRec.depth >= m_maxDepth && m_maxDepth != -1)
          break;

        /* Prevent light leaks due to the use of shading normals */
        Float wiDotGeoN = -dot(its.geoFrame.n, ray.d),
            wiDotShN = Frame::cosTheta(its.wi);
        if (wiDotGeoN * wiDotShN < 0 && m_strictNormals)
          break;

        /* ==================================================================== */
        /*                          Luminaire sampling                          */
        /* ==================================================================== */

        const BSDF *bsdf = its.getBSDF(ray);

        // Now proceed for luminaire sampling
        DirectSamplingRecord dRec(its);
        if (mediumScattered) {
          /* Estimate the direct illumination if this is requested */
          if ((rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance) &&
              (bsdf->getType() & BSDF::ESmooth)) {
            int interactions = m_maxDepth - rRec.depth - 1;

            Spectrum value = scene->sampleAttenuatedEmitterDirect(
                dRec, its, rRec.medium, interactions,
                rRec.nextSample2D(), rRec.sampler);

            if (!value.isZero()) {
              const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

              /* Evaluate BSDF * cos(theta) */
              BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));
              const Spectrum bsdfVal = bsdf->eval(bRec);

              Float woDotGeoN = dot(its.geoFrame.n, dRec.d);

              /* Prevent light leaks due to the use of shading normals */
              if (!bsdfVal.isZero() && (!m_strictNormals ||
                  woDotGeoN * Frame::cosTheta(bRec.wo) > 0)) {
                /* Calculate prob. of having generated that direction
                   using BSDF sampling */
                Float bsdfPdf = (emitter->isOnSurface()
                    && dRec.measure == ESolidAngle)
                                ? bsdf->pdf(bRec) : (Float) 0.0f;

                /* Weight using the power heuristic */
                const Float weight = miWeight(dRec.pdf, bsdfPdf);
                if (rRec.depth + 1 >= m_minDepth && nbSurfaceInteractions >= m_minCameraDepth) {
#if SURFACE_DIRECT_TRANS
                  Li += throughput * value * bsdfVal * weight;
#endif
                }
              }
            }
          }
        }


        /* ==================================================================== */
        /*                            BSDF sampling                             */
        /* ==================================================================== */

        // We do BSDF sampling before luminaire sampling so we know which BSDF component
        // to use for vertex classification
        /* Sample BSDF * cos(theta) */
        BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
        Float bsdfPdf;
        Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());

        // This part of the code might be responsible to the difference between G-PT
        // However, the difference might be rather minor and does not appears inside the scenes
        // used inside the paper.
        if (VertexClassifier::type(bsdf, its, bRec.sampledComponent) == VERTEX_TYPE_DIFFUSE) {
          if (!mediumScattered) {
            // In this case the path is imediately stop
            break;
          }
        }

        // BSDF sampled before. Now just proceed.
        if (bsdfWeight.isZero())
          break;

        /* Prevent light leaks due to the use of shading normals */
        const Vector wo = its.toWorld(bRec.wo);
        Float woDotGeoN = dot(its.geoFrame.n, wo);
        if (woDotGeoN * Frame::cosTheta(bRec.wo) <= 0 && m_strictNormals)
          break;

        /* Trace a ray in this direction */
        ray = Ray(its.p, wo, ray.time);

        /* Keep track of the throughput, medium, and relative
           refractive index along the path */
        throughput *= bsdfWeight;
        eta *= bRec.eta;
        if (its.isMediumTransition())
          rRec.medium = its.getTargetMedium(ray.d);

        /* Handle index-matched medium transitions specially */
        if (bRec.sampledType == BSDF::ENull) {
          if (!(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
            break;
          rRec.type = scattered ? RadianceQueryRecord::ERadianceNoEmission
                                : RadianceQueryRecord::ERadiance;
          scene->rayIntersect(ray, its);
          rRec.depth++;
          continue;
        }

        Spectrum value(0.0f);
        rayIntersectAndLookForEmitter(scene, rRec.sampler, rRec.medium,
                                      m_maxDepth - rRec.depth - 1, ray, its, dRec, value);

        /* If a luminaire was hit, estimate the local illumination and
           weight using the power heuristic */
        if (mediumScattered) {
          if (!value.isZero() && (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
            const Float emitterPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
                                     scene->pdfEmitterDirect(dRec) : 0;
            if (rRec.depth + 1 >= m_minDepth && nbSurfaceInteractions >= m_minCameraDepth) {
#if SURFACE_INTER_TRANS
              Li += throughput * value * miWeight(bsdfPdf, emitterPdf);
#endif
            }
          }
        }
        /* ==================================================================== */
        /*                         Indirect illumination                        */
        /* ==================================================================== */

        /* Stop if indirect illumination was not requested */
        if (!(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
          break;

        rRec.type = RadianceQueryRecord::ERadianceNoEmission;
      }

      if (rRec.depth++ >= m_rrDepth) {
        /* Russian roulette: try to keep path weights equal to one,
           while accounting for the solid angle compression at refractive
           index boundaries. Stop with at least some probability to avoid
           getting stuck (e.g. due to total internal reflection) */

        Float q = std::min(throughput.max() * eta * eta, (Float) 0.95f);
        if (rRec.nextSample1D() >= q)
          break;
        throughput /= q;
      }

      scattered = true;
    }
    avgPathLength.incrementBase();
    avgPathLength += rRec.depth;
    return Li;
  }

  /**
   * This function is called by the recursive ray tracing above after
   * having sampled a direction from a BSDF/phase function. Due to the
   * way in which this integrator deals with index-matched boundaries,
   * it is necessarily a bit complicated (though the improved performance
   * easily pays for the extra effort).
   *
   * This function
   *
   * 1. Intersects 'ray' against the scene geometry and returns the
   *    *first* intersection via the '_its' argument.
   *
   * 2. It checks whether the intersected shape was an emitter, or if
   *    the ray intersects nothing and there is an environment emitter.
   *    In this case, it returns the attenuated emittance, as well as
   *    a DirectSamplingRecord that can be used to query the hypothetical
   *    sampling density at the emitter.
   *
   * 3. If current shape is an index-matched medium transition, the
   *    integrator keeps on looking on whether a light source eventually
   *    follows after a potential chain of index-matched medium transitions,
   *    while respecting the specified 'maxDepth' limits. It then returns
   *    the attenuated emittance of this light source, while accounting for
   *    all attenuation that occurs on the wya.
   */
  void rayIntersectAndLookForEmitter(const Scene *scene, Sampler *sampler,
                                     const Medium *medium, int maxInteractions, Ray ray, Intersection &_its,
                                     DirectSamplingRecord &dRec, Spectrum &value) const {
    Intersection its2, *its = &_its;
    Spectrum transmittance(1.0f);
    bool surface = false;
    int interactions = 0;

    while (true) {
      surface = scene->rayIntersect(ray, *its);

      if (medium)
        transmittance *= medium->evalTransmittance(Ray(ray, 0, its->t), sampler);

      if (surface && (interactions == maxInteractions ||
          !(its->getBSDF()->getType() & BSDF::ENull) ||
          its->isEmitter())) {
        /* Encountered an occluder / light source */
        break;
      }

      if (!surface)
        break;

      if (transmittance.isZero())
        return;

      if (its->isMediumTransition())
        medium = its->getTargetMedium(ray.d);

      Vector wo = its->shFrame.toLocal(ray.d);
      BSDFSamplingRecord bRec(*its, -wo, wo, ERadiance);
      bRec.typeMask = BSDF::ENull;
      transmittance *= its->getBSDF()->eval(bRec, EDiscrete);

      ray.o = ray(its->t);
      ray.mint = Epsilon;
      its = &its2;

      if (++interactions > 100) { /// Just a precaution..
        Log(EWarn, "rayIntersectAndLookForEmitter(): round-off error issues?");
        return;
      }
    }

    if (surface) {
      /* Intersected something - check if it was a luminaire */
      if (its->isEmitter()) {
        dRec.setQuery(ray, *its);
        value = transmittance * its->Le(-ray.d);
      }
    } else {
      /* Intersected nothing -- perhaps there is an environment map? */
      const Emitter *env = scene->getEnvironmentEmitter();

      if (env && env->fillDirectSamplingRecord(dRec, ray))
        value = transmittance * env->evalEnvironment(RayDifferential(ray));
    }
  }

  inline Float miWeight(Float pdfA, Float pdfB) const {
    pdfA *= pdfA;
    pdfB *= pdfB;
    return pdfA / (pdfA + pdfB);
  }

  void serialize(Stream *stream, InstanceManager *manager) const {
    MonteCarloIntegrator::serialize(stream, manager);
  }

  std::string toString() const {
    std::ostringstream oss;
    oss << "VolumetricPathTracer[" << endl
        << "  maxDepth = " << m_maxDepth << "," << endl
        << "  rrDepth = " << m_rrDepth << "," << endl
        << "  strictNormals = " << m_strictNormals << endl
        << "]";
    return oss.str();
  }

  MTS_DECLARE_CLASS()

private:
  int m_minDepth;
  int m_minCameraDepth;
  ELightingEffects m_lightingInteractionMode;
  bool m_directLightVisibility;
};

MTS_IMPLEMENT_CLASS_S(VolumetricPathTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(VolumetricPathTracer, "Volumetric path tracer");
MTS_NAMESPACE_END
