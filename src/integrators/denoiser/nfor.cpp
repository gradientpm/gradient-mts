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
#include <mitsuba/render/renderproc.h>
#include <mitsuba/core/plugin.h>

#include "../denoiser/aux_buffer.h"
#include "../denoiser/nfor/nfor.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{multichannel}{Multi-channel integrator}
 * \order{16}
 * \parameters{
 *     \parameter{\Unnamed}{\Integrator}{One or more sub-integrators whose output
 *     should be rendered into a combined multi-channel image}
 * }
 *
 * The multi-channel integrator groups several sub-integrators together
 * and invokes them at the same time for each pixel; the result from each
 * integrator is written into a separate channel of the output image.
 * This could include things like surface normals or the distance
 * from the camera (via the \pluginref{field} plugin) or ambient occlusion
 * (via the \pluginref{ao} plugin).
 * In this way, this integrator can be a powerful tool for unusual applications
 * of Mitsuba, e.g. to create reference data for computer vision algorithms. Currently, it only
 * works with a subset of the other plugins---see the red box for details.
 *
 * Thee \code{multichannel} plugin also disables certain checks for negative or infinite
 * radiance values during rendering that normally cause warnings to be emitted.
 * This is simply to process extracted fields for which it is fine
 * to take on such values.
 *
 * The following example contains a typical setup for rendering an 7 channel EXR image:
 * 3 for a path traced image (RGB), 3 for surface normals
 * (encoded as RGB), and 1 channel for the ray distance measured from the camera.
 *
 * \vspace{2mm}
 * \begin{xml}
 * <scene>
 *     <integrator type="multichannel">
 *         <integrator type="path"/>
 *         <integrator type="field">
 *             <string name="field" value="shNormal"/>
 *         </integrator>
 *         <integrator type="field">
 *             <string name="field" value="distance"/>
 *         </integrator>
 *     </integrator>
 *
 *     <sensor type="perspective">
 *         <sampler type="halton">
 *             <integer name="sampleCount" value="32"/>
 *         </sampler>
 *         <film type="hdrfilm">
 *             <string name="pixelFormat" value="rgb, rgb, luminance"/>
 *             <string name="channelNames" value="color, normal, distance"/>
 *         </film>
 *     </sensor>
 *     <!-- **** scene contents **** -->
 * </scene>
 * \end{xml}
 *
 * \remarks{
 * \item Requires the \pluginref{hdrfilm} or \pluginref{tiledhdrfilm}.
 * \item All nested integrators must
 * conform to Mitsuba's basic \emph{SamplingIntegrator} interface.
 * Currently, only a few of them do this, including:
 * \pluginref{field}, \pluginref{ao}, \pluginref{direct}, \pluginref{path},
 * \pluginref{volpath}, \pluginref[volpathsimple]{volpath\_simple},
 * and \pluginref{irrcache}.
 * }
 */

#define MEDIUM_DIRECT_TRANS 1
#define MEDIUM_INTER_TRANS 1
#define SURFACE_DIRECT_TRANS 1
#define SURFACE_INTER_TRANS 1


class NFORIntegrator : public SamplingIntegrator {
public:
    NFORIntegrator(const Properties &props) : SamplingIntegrator(props) {
        m_auxBuffer = 0;
        m_throughputBuffer = 0;

        /* Depth to begin using russian roulette */
        m_rrDepth = props.getInteger("rrDepth", 5);

        /* Longest visualized path depth (\c -1 = infinite).
           A value of \c 1 will visualize only directly visible light sources.
           \c 2 will lead to single-bounce (direct-only) illumination, and so on. */
        m_maxDepth = props.getInteger("maxDepth", -1);

        /**
         * This parameter specifies the action to be taken when the geometric
         * and shading normals of a surface don't agree on whether a ray is on
         * the front or back-side of a surface.
         *
         * When \c strictNormals is set to \c false, the shading normal has
         * precedence, and rendering proceeds normally at the risk of
         * introducing small light leaks (this is the default).
         *
         * When \c strictNormals is set to \c true, the random walk is
         * terminated when encountering such a situation. This may
         * lead to silhouette darkening on badly tesselated meshes.
         */
        m_strictNormals = props.getBoolean("strictNormals", false);

        /**
         * When this flag is set to true, contributions from directly
         * visible emitters will not be included in the rendered image
         */
        m_hideEmitters = props.getBoolean("hideEmitters", false);
        m_maxPasses = 2000;

        if (m_rrDepth <= 0)
            Log(EError, "'rrDepth' must be set to a value greater than zero!");

        if (m_maxDepth <= 0 && m_maxDepth != -1)
            Log(EError, "'maxDepth' must be set to -1 (infinite) or a value greater than zero!");
    }

    NFORIntegrator(Stream *stream, InstanceManager *manager)
            : SamplingIntegrator(stream, manager) {
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        SamplingIntegrator::serialize(stream, manager);
    }

    bool preprocess(const Scene *scene, RenderQueue *queue,
                    const RenderJob *job, int sceneResID, int sensorResID,
                    int samplerResID) {
        m_stop = false;
        if (!SamplingIntegrator::preprocess(scene, queue, job, sceneResID,
                                            sensorResID, samplerResID))
            return false;
        return true;
    }

    bool render(Scene *scene,
                RenderQueue *queue, const RenderJob *job,
                int sceneResID, int sensorResID, int samplerResID) {
        ref<Scheduler> sched = Scheduler::getInstance();
        ref<Sensor> sensor = static_cast<Sensor *>(sched->getResource(sensorResID));
        ref<Film> film = sensor->getFilm();
        auto cropSize = film->getCropSize();

        scene->getSensor()->getFilm(); // We can replace it??


        m_auxBuffer = new AuxiliaryBuffer(film->getSize(), true, true);
        m_auxBuffer->clear();
        m_throughputBuffer = new AccumBuffer(1, film->getSize(), true, true);

        size_t nCores = sched->getCoreCount();
        const Sampler *sampler = static_cast<const Sampler *>(sched->getResource(samplerResID, 0));
        size_t sampleCount = sampler->getSampleCount();

        // Time for rendering
        std::string timeFilename = scene->getDestinationFile().string()
                                   + "_time.csv";
        std::ofstream timeFile(timeFilename.c_str());
        ref<Timer> renderingTimer = new Timer;
        Float cumulativeTime = 0.f;

        for (int it = 1; it < m_maxPasses && (!m_stop); it++) {
            Log(EInfo, "Starting render job (%ix%i, "
                    SIZE_T_FMT
                    " %s, "
                    SIZE_T_FMT
                    " %s, "
                    SSE_STR
                    ") ..", film->getCropSize().x, film->getCropSize().y,
                sampleCount, sampleCount == 1 ? "sample" : "samples", nCores,
                nCores == 1 ? "core" : "cores");

            /* This is a sampling-based integrator - parallelize */
            ref<BlockedRenderProcess> proc = new BlockedRenderProcess(job,
                                                                      queue, scene->getBlockSize());
            int integratorResID = sched->registerResource(this);
            proc->bindResource("integrator", integratorResID);
            proc->bindResource("scene", sceneResID);
            proc->bindResource("sensor", sensorResID);
            proc->bindResource("sampler", samplerResID);
            scene->bindUsedResources(proc);
            bindUsedResources(proc);
            sched->schedule(proc);

            m_process = proc;
            sched->wait(proc);
            m_process = NULL;
            sched->unregisterResource(integratorResID);

            // Output information buffers
            Vector2i size = scene->getFilm()->getSize();
            Properties pHDRFilm("hdrfilm");
            pHDRFilm.setInteger("width", size.x);
            pHDRFilm.setInteger("height", size.y);
            pHDRFilm.setString("fileFormat", "rgbe");
            // Create new film based on HDR film
            ref<Film> hdrFilm = dynamic_cast<Film *> (PluginManager::getInstance()->
                    createObject(MTS_CLASS(Film), pHDRFilm));

            ref<Bitmap> throughputBitmap = m_throughputBuffer->getBuffer(0);
            ref<Bitmap> normalBitmap = m_auxBuffer->getBuffer(EAuxBuffer::NormalBuffer);
            ref<Bitmap> positionBitmap = m_auxBuffer->getBuffer(EAuxBuffer::PositionBuffer);
            ref<Bitmap> albedoBitmap = m_auxBuffer->getBuffer(EAuxBuffer::AlbedoBuffer);
            ref<Bitmap> nforBitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, cropSize);

            /// Time it
            unsigned int milliseconds = renderingTimer->getMilliseconds();
            timeFile.flush();
            Log(EInfo, "Rendering time: %i, %i", milliseconds / 1000,
                milliseconds % 1000);
            timeFile << (milliseconds / (1000.f)) << ",\n";
            cumulativeTime += (milliseconds / 1000.f);
            renderingTimer->reset();

            // Exclude from the rendering time the reconstruction
            auto nfor = new Nfor();
            nfor->denoise(m_throughputBuffer.get(), nforBitmap, m_auxBuffer.get());
            unsigned int millisecondsRecons = renderingTimer->getMilliseconds();
            Log(EInfo, "Reconstruction time: %i, %i", millisecondsRecons / 1000,
                millisecondsRecons % 1000);
            renderingTimer->reset();

            develop(scene, hdrFilm, nforBitmap, it, "_recons_");
            develop(scene, hdrFilm, throughputBitmap, it, "_th_");
            develop(scene, hdrFilm, positionBitmap, it, "_position_");
            develop(scene, hdrFilm, albedoBitmap, it, "_albedo_");
            develop(scene, hdrFilm, normalBitmap, it, "_normal_");
            if (proc->getReturnStatus() != ParallelProcess::ESuccess) {
                return false;
            }
        }
        return true;
    }

    void develop(Scene *scene, ref<Film> film, ref<Bitmap> bitmap,
                 int currentIteration, const std::string &suffixName = "_") {
        std::stringstream ss;
        ss << scene->getDestinationFile().string() << suffixName
           << currentIteration;
        std::string path = ss.str();

        // bool exr = false
        film->setBitmap(bitmap);
        film->setDestinationFile(path, 0);
        film->develop(scene, 0.f);
    }

    void renderBlock(const Scene *scene,
                     const Sensor *sensor, Sampler *sampler, ImageBlock *block,
                     const bool &stop, const std::vector<TPoint2<uint8_t> > &points) const {
        Float diffScaleFactor = 1.0f /
                                std::sqrt((Float) sampler->getSampleCount());

        bool needsApertureSample = sensor->needsApertureSample();
        bool needsTimeSample = sensor->needsTimeSample();

        RadianceQueryRecord rRec(scene, sampler);
        Point2 apertureSample(0.5f);
        Float timeSample = 0.5f;
        RayDifferential sensorRay;

        auto throughputBuffer = const_cast<AccumBuffer *>(m_throughputBuffer.get());
        block->clear();

        uint32_t queryType = RadianceQueryRecord::ESensorRay;

        if (!sensor->getFilm()->hasAlpha()) /* Don't compute an alpha channel if we don't have to */
            queryType &= ~RadianceQueryRecord::EOpacity;

        for (size_t i = 0; i < points.size(); ++i) {
            Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
            if (stop)
                break;

            sampler->generate(offset);

            for (size_t j = 0; j < sampler->getSampleCount(); j++) {
                rRec.newQuery(queryType, sensor->getMedium());
                Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

                if (needsApertureSample)
                    apertureSample = rRec.nextSample2D();
                if (needsTimeSample)
                    timeSample = rRec.nextSample1D();

                Spectrum spec = sensor->sampleRayDifferential(
                        sensorRay, samplePos, apertureSample, timeSample);
                sensorRay.scaleDifferential(diffScaleFactor);
                if (m_noRayDiff) {
                    sensorRay.hasDifferentials = false;
                }

                spec *= LiAux(sensorRay, rRec, const_cast<AuxiliaryBuffer *>(m_auxBuffer.get()), samplePos);
                block->put(samplePos, spec, rRec.alpha);
                throughputBuffer->addSample(samplePos.x, samplePos.y, spec);
                sampler->advance();
            }
        }
    }


    bool hasSmoothComponent(const BSDF *bsdf,
                            Intersection &its) const {

        bool found_smooth = false;
        bool found_dirac = false;
        for (int i = 0, component_count = bsdf->getComponentCount(); i < component_count; ++i) {
            Float component_roughness = bsdf->getRoughness(its, i);

            if (component_roughness == Float(0)) {
                found_dirac = true;
            } else {
                found_smooth = true;
            }
        }
        return found_smooth;
    }

    Spectrum LiAux(const RayDifferential &r, RadianceQueryRecord &rRec,
                   AuxiliaryBuffer *auxBuffer, Point2 samplePosition) const {
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
        bool mediumScattered = true;

        Spectrum throughput(1.0f);
        bool scattered = false;
        int nbSurfaceInteractions = 0;

        bool isAlbedoStored = false;
        {
            const int x = (int) samplePosition.x;
            const int y = (int) samplePosition.y;

            const Transform &camera2world = scene->getSensor()->getWorldTransform()->eval(0.f);
            const Transform world2camera = camera2world.inverse();
            const Point cameraCoord = world2camera(its.p);

            // add position (depth)
            auxBuffer->addSample(x, y, Vector(cameraCoord.z, cameraCoord.z, cameraCoord.z), EAuxBuffer::PositionBuffer);

            // add normal
            const auto nc = world2camera(its.shFrame.n);
            auxBuffer->addSample(x, y, Vector(nc.x, nc.y, nc.z), EAuxBuffer::NormalBuffer);

            //revert incoming dir (-wi)
            auxBuffer->addSample(x, y, ray.d, EAuxBuffer::RayW);
        }

        while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
            /* ==================================================================== */
            /*                 Radiative Transfer Equation sampling                 */
            /* ==================================================================== */
            if (rRec.medium && rRec.medium->sampleDistance(Ray(ray, 0, its.t), mRec, rRec.sampler)) {
                isAlbedoStored = true;
                /* Sample the integral
                 \int_x^y tau(x, x') [ \sigma_s \int_{S^2} \rho(\omega,\omega') L(x,\omega') d\omega' ] dx'
              */
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
#if MEDIUM_DIRECT_TRANS
                                Li += throughput * value * phaseVal * weight;
#endif

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
#if MEDIUM_INTER_TRANS
                        Li += throughput * value * miWeight(phasePdf, emitterPdf);
#endif
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

                if (!mediumScattered) {
                    nbSurfaceInteractions += 1;
                }

                if (!its.isValid()) {
                    /* If no intersection could be found, possibly return
                       attenuated radiance from a background luminaire */
                    if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
                        && (!m_hideEmitters || scattered)) {
                        Spectrum value = throughput * scene->evalEnvironment(ray);
                        if (rRec.medium)
                            value *= rRec.medium->evalTransmittance(ray, rRec.sampler);
                        Li += value;
                    }

                    break;
                }

                /* Possibly include emitted radiance if requested */
                if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
                    && ((!m_hideEmitters || scattered) || mediumScattered)) {

                    Li += throughput * its.Le(-ray.d);

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
#if SURFACE_DIRECT_TRANS
                                Li += throughput * value * bsdfVal * weight;
#endif
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

                // BSDF sampled before. Now just proceed.
                if (bsdfWeight.isZero())
                    break;

                if (!isAlbedoStored) {
                    bool hasDiffuse = hasSmoothComponent(bsdf, its);
                    // We want to be consistent by collecting the albedo and shading normal at the same shading point
                    if (hasDiffuse) {
                        const int x = (int) samplePosition.x;
                        const int y = (int) samplePosition.y;

                        const Spectrum reflectance = bsdf->getDiffuseReflectance(its);
                        auxBuffer->addSample(x, y, Vector(reflectance[0], reflectance[1], reflectance[2]),
                                             EAuxBuffer::AlbedoBuffer);
                        isAlbedoStored = true;

                        // Have to call getFrame() just in case we have bumpmap
                        Frame shFrame = bsdf->getFrame(its);
                        auxBuffer->addSample(x, y, Vector(shFrame.n.x, shFrame.n.y, shFrame.n.z),
                                             EAuxBuffer::ShadingNormalW);
                    }
                }

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
#if SURFACE_INTER_TRANS
                        Li += throughput * value * miWeight(bsdfPdf, emitterPdf);
#endif
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
        return Li;
    }

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

    Spectrum LiAuxSurface(const RayDifferential &r, RadianceQueryRecord &rRec,
                          AuxiliaryBuffer *auxBuffer, Point2 samplePosition) const {
        /* Some aliases and local variables */
        const Scene *scene = rRec.scene;
        Intersection &its = rRec.its;
        RayDifferential ray(r);
        bool scattered = false;

        /* Perform the first ray intersection (or ignore if the
           intersection has already been provided). */
        rRec.rayIntersect(ray);
        ray.mint = Epsilon;

        Spectrum Li(0.f);
        Spectrum throughput(1.0f);
        Float eta = 1.0f;
        Float depth = rRec.its.t;

        if (!its.isValid()) {
            if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
                && (!m_hideEmitters || scattered))
                Li += throughput * scene->evalEnvironment(ray);
        }

        bool isAlbedoStored = false;
        {
            const int x = (int) samplePosition.x;
            const int y = (int) samplePosition.y;

            const Transform &camera2world = scene->getSensor()->getWorldTransform()->eval(0.f);
            const Transform world2camera = camera2world.inverse();
            const Point cameraCoord = world2camera(its.p);

            // add position (depth)
            auxBuffer->addSample(x, y, Vector(cameraCoord.z, cameraCoord.z, cameraCoord.z), EAuxBuffer::PositionBuffer);

            // add normal
            const auto nc = world2camera(its.shFrame.n);
            auxBuffer->addSample(x, y, Vector(nc.x, nc.y, nc.z), EAuxBuffer::NormalBuffer);

            //revert incoming dir (-wi)
            auxBuffer->addSample(x, y, ray.d, EAuxBuffer::RayW);
        }

        while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
            if (!its.isValid()) {
                /* If no intersection could be found, potentially return
                   radiance from a environment luminaire if it exists */
                if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
                    && (!m_hideEmitters || scattered))
                    Li += throughput * scene->evalEnvironment(ray);
                break;
            }

            const BSDF *bsdf = its.getBSDF(ray);

            /* Possibly include emitted radiance if requested */
            if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
                && (!m_hideEmitters || scattered))
                Li += throughput * its.Le(-ray.d);

            /* Include radiance from a subsurface scattering model if requested */
            if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
                Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

            if ((rRec.depth >= m_maxDepth && m_maxDepth > 0)
                || (m_strictNormals && dot(ray.d, its.geoFrame.n)
                                       * Frame::cosTheta(its.wi) >= 0)) {

                /* Only continue if:
                   1. The current path length is below the specifed maximum
                   2. If 'strictNormals'=true, when the geometric and shading
                      normals classify the incident direction to the same side */
                break;
            }

            /* ==================================================================== */
            /*                     Direct illumination sampling                     */
            /* ==================================================================== */

            /* Estimate the direct illumination if this is requested */
            DirectSamplingRecord dRec(its);

            if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
                (bsdf->getType() & BSDF::ESmooth)) {
                Spectrum value = scene->sampleEmitterDirect(dRec, rRec.nextSample2D());
                if (!value.isZero()) {
                    const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

                    /* Allocate a record for querying the BSDF */
                    BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);

                    /* Evaluate BSDF * cos(theta) */
                    const Spectrum bsdfVal = bsdf->eval(bRec);

                    /* Prevent light leaks due to the use of shading normals */
                    if (!bsdfVal.isZero() && (!m_strictNormals
                                              || dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0)) {

                        /* Calculate prob. of having generated that direction
                           using BSDF sampling */
                        Float bsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle)
                                        ? bsdf->pdf(bRec) : 0;

                        /* Weight using the power heuristic */
                        Float weight = miWeight(dRec.pdf, bsdfPdf);
                        Li += throughput * value * bsdfVal * weight;
                    }
                }
            }

            /* ==================================================================== */
            /*                            BSDF sampling                             */
            /* ==================================================================== */

            /* Sample BSDF * cos(theta) */
            Float bsdfPdf;
            BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
            Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
            if (bsdfWeight.isZero())
                break;

            scattered |= bRec.sampledType != BSDF::ENull;
            if (!isAlbedoStored) {
                bool hasDiffuse = hasSmoothComponent(bsdf, its);
                // We want to be consistent by collecting the albedo and shading normal at the same shading point
                if (hasDiffuse) {
                    const int x = (int) samplePosition.x;
                    const int y = (int) samplePosition.y;

                    const Spectrum reflectance = bsdf->getDiffuseReflectance(its);
                    auxBuffer->addSample(x, y, Vector(reflectance[0], reflectance[1], reflectance[2]),
                                         EAuxBuffer::AlbedoBuffer);
                    isAlbedoStored = true;

                    // Have to call getFrame() just in case we have bumpmap
                    Frame shFrame = bsdf->getFrame(its);
                    auxBuffer->addSample(x, y, Vector(shFrame.n.x, shFrame.n.y, shFrame.n.z),
                                         EAuxBuffer::ShadingNormalW);
                }
            }

            /* Prevent light leaks due to the use of shading normals */
            const Vector wo = its.toWorld(bRec.wo);
            Float woDotGeoN = dot(its.geoFrame.n, wo);
            if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
                break;

            bool hitEmitter = false;
            Spectrum value;

            /* Trace a ray in this direction */
            ray = Ray(its.p, wo, ray.time);
            if (scene->rayIntersect(ray, its)) {
                depth += its.t;
                /* Intersected something - check if it was a luminaire */
                if (its.isEmitter()) {
                    value = its.Le(-ray.d);
                    dRec.setQuery(ray, its);
                    hitEmitter = true;
                }
            } else {
                /* Intersected nothing -- perhaps there is an environment map? */
                const Emitter *env = scene->getEnvironmentEmitter();

                if (env) {
                    if (m_hideEmitters && !scattered)
                        break;

                    value = env->evalEnvironment(ray);
                    if (!env->fillDirectSamplingRecord(dRec, ray))
                        break;
                    hitEmitter = true;
                } else {
                    break;
                }
            }

            /* Keep track of the throughput and relative
               refractive index along the path */
            throughput *= bsdfWeight;
            eta *= bRec.eta;

            /* If a luminaire was hit, estimate the local illumination and
               weight using the power heuristic */
            if (hitEmitter &&
                (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
                /* Compute the prob. of generating that direction using the
                   implemented direct illumination sampling technique */
                const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
                                     scene->pdfEmitterDirect(dRec) : 0;
                Li += throughput * value * miWeight(bsdfPdf, lumPdf);
            }

            /* ==================================================================== */
            /*                         Indirect illumination                        */
            /* ==================================================================== */

            /* Set the recursive query type. Stop if no surface was hit by the
               BSDF sample or if indirect illumination was not requested */
            if (!its.isValid() || !(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
                break;
            rRec.type = RadianceQueryRecord::ERadianceNoEmission;

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
        }

        return Li;
    }

    inline Float miWeight(Float pdfA, Float pdfB) const {
        pdfA *= pdfA;
        pdfB *= pdfB;
        return pdfA / (pdfA + pdfB);
    }


    void bindUsedResources(ParallelProcess *proc) const {
        SamplingIntegrator::bindUsedResources(proc);
    }

    void wakeup(ConfigurableObject *parent, std::map<std::string, SerializableObject *> &params) {
        SamplingIntegrator::wakeup(parent, params);
    }


    void configureSampler(const Scene *scene, Sampler *sampler) {
        SamplingIntegrator::configureSampler(scene, sampler);
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(SamplingIntegrator))) {
            Log(EError, "Do not support nested integrators");
        } else {
            SamplingIntegrator::addChild(name, child);
        }
    }

    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        NotImplementedError("Li");
    }

    void cancel() {
        m_stop = true;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "MultiChannelIntegrator[" << endl
            << "  integrators = {" << endl;
        oss << "  }" << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()

private:
    // Path tracer configuration
    int m_maxDepth;
    int m_rrDepth;
    bool m_strictNormals;
    bool m_hideEmitters;

    ref<AuxiliaryBuffer> m_auxBuffer;
    ref<AccumBuffer> m_throughputBuffer;
    size_t m_maxPasses;
    bool m_stop;
};

MTS_IMPLEMENT_CLASS_S(NFORIntegrator, false, SamplingIntegrator)

MTS_EXPORT_PLUGIN(NFORIntegrator, "NFOR with Path tracing integrator");
MTS_NAMESPACE_END
