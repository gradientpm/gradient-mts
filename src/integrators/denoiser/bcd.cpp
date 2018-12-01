#include <mitsuba/core/plugin.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/render/renderqueue.h>

#include <omp.h>

#include "bcd/SamplesAccumulator.h"
#include "bcd/Utils.h"
#include "bcd/SpikeRemovalFilter.h"
#include "bcd/Denoiser.h"
#include "bcd/MultiscaleDenoiser.h"
#include "bcd/IDenoiser.h"

// For testing, make possible to disable
// each type of transport (without compromise MIS)
#define MEDIUM_DIRECT_TRANS 1
#define MEDIUM_INTER_TRANS 1
#define SURFACE_DIRECT_TRANS 1
#define SURFACE_INTER_TRANS 1

MTS_NAMESPACE_BEGIN

class BCDIntegrator : public Integrator {
public:
    BCDIntegrator(const Properties &props) : Integrator(props) {
        m_maxPasses = props.getInteger("maxPasses", 200000);
        m_maxRenderingTime = props.getInteger("maxRenderingTime", 200000000);
        m_dumpIteration = props.getInteger("dumpIteration", 1);

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
        m_noRayDiff = true;
    }

    bool preprocess(const Scene *scene, RenderQueue *queue,
                    const RenderJob *job, int sceneResID, int sensorResID,
                    int samplerResID) {
        m_stop = false;
        return true;
    }

    bool render(Scene *scene, RenderQueue *queue, const RenderJob *job,
                int sceneResID, int sensorResID, int samplerResID) {

        /// Get data from Mitsuba
        ref<Scheduler> sched = Scheduler::getInstance();
        ref<Sensor> sensor = scene->getSensor();
        ref<Film> film = sensor->getFilm();
        auto sampler = sensor->getSampler();
        size_t nCores = sched->getCoreCount();
        Log(EInfo, "Starting render job (%ix%i, "
                SIZE_T_FMT
                " %s, "
                SSE_STR
                ") ..",
            film->getCropSize().x, film->getCropSize().y,
            nCores, nCores == 1 ? "core" : "cores");
        Vector2i cropSize = film->getCropSize();
        size_t sampleCount = sensor->getSampler()->getSampleCount();

        /// Create bitmaps
        ref<Bitmap> bitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, cropSize);
        ref<Bitmap> bitmapAvg = bitmap->clone();
        ref<Bitmap> bitmapVeryDirect = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, cropSize);
        ref<Bitmap> bitmapVeryDirectAvg = bitmap->clone();
        bitmapAvg->clear();
        bitmapVeryDirectAvg->clear();

        /// Create all struct
        std::string timeFilename = scene->getDestinationFile().string()
                                   + "_time.csv";
        std::ofstream timeFile(timeFilename.c_str());
        ref<Timer> renderingTimer = new Timer;
        Float cumulativeTime = 0.f;

        bcd::HistogramParameters histoParams;
        histoParams.m_nbOfBins = 20;
        histoParams.m_gamma = 2.2f;
        histoParams.m_maxValue = 2.5f;
        bcd::SamplesAccumulator samplesAccumulator(film->getCropSize().x, film->getCropSize().y, histoParams);

        // Generates all the buffer to store all the samples.
        // note that the current implementation is very memory demanding.
        // streaming directly the samples by using BCD as a library might be
        // a better idea
        m_values.reserve(cropSize.x * cropSize.y);
        for (size_t t = 0; t < cropSize.x * cropSize.y; t++) {
            std::vector<Spectrum> pixel_values;
            pixel_values.reserve(sampleCount);
            m_values.emplace_back(pixel_values);
        }
        // We also need to generate all the samplers.
        // indeed, the rendering loop in BCD is a little bit different
        // as we want to keep all the samples. It is possible to integrate
        // BCD deeper inside Mistuba. However, for early testing, this approach is good enough.
        std::vector<Sampler *> samplers(sched->getCoreCount());
        for (size_t i = 0; i < sched->getCoreCount(); ++i) {
            ref<Sampler> clonedSampler = sampler->clone();
            clonedSampler->incRef();
            samplers[i] = clonedSampler.get();
        }

        // Rendering loop.
        // We split the samples over time to be able to generate convergence curves for this type
        // of technique.
        for (int it = 1; it < m_maxPasses && (!m_stop) && (cumulativeTime < m_maxRenderingTime); it++) {
            bitmap->clear(); //< Clear bitmap this pass
            bitmapVeryDirect->clear();
            film->clear(); //< and setup
            Spectrum *target = (Spectrum *) bitmap->getUInt8Data();
            Spectrum *targetVeryDirect = (Spectrum *) bitmapVeryDirect->getUInt8Data();

            // Only store the samples at this given iteration
            // The different iterations will be concatenated with a python script
            for (auto &v: m_values) {
                v.clear();
            }

            // Rendering the image
            SLog(EInfo, "Rendering...");
#pragma omp parallel for num_threads(sched->getCoreCount())
            for (size_t y = 0; y < cropSize.y; y++) {
                int tid = omp_get_thread_num();
                auto local_sampler = samplers[tid];
                for (size_t x = 0; x < cropSize.x; x++) {
                    auto results = renderPixel(Point2i(x, y),
                                               scene, sensor,
                                               local_sampler,
                                               m_stop, m_values[y * cropSize.x + x]);
                    target[y * cropSize.x + x] = results.Li;
                    targetVeryDirect[y * cropSize.x + x] = results.LiVery;
                }
            }

            /// Write the raw information for BCD in a raw format
            /// This part of the process might be long depending to the disk performance
            /// So it might be more fair to remove this part when we time the technique.
            if (!m_stop) {
                SLog(EInfo, "Accumulating...");
                for (size_t y = 0; y < cropSize.y; y++) {
                    for (size_t x = 0; x < cropSize.x; x++) {
                        const std::vector<Spectrum> &values = m_values[y * cropSize.x + x];
                        if (values.size() != sampleCount) {
                            SLog(EError, "Wrong number of samples stored. Abord (%i != %i)", values.size(),
                                 sampleCount * it);
                        }
                        for (size_t n = 0; n < values.size(); n++) {
                            Float r, g, b;
                            values[n].toLinearRGB(r, g, b);
                            samplesAccumulator.addSample(y, x, r, g, b);
                        }
                    }
                }
            }

            /// Average the image results (for the visualization)
            if (!m_stop) {
                bitmapAvg->scale(it - 1);
                bitmapAvg->accumulate(bitmap);
                bitmapAvg->scale(1.f / (it));

                bitmapVeryDirectAvg->scale(it - 1);
                bitmapVeryDirectAvg->accumulate(bitmapVeryDirect);
                bitmapVeryDirectAvg->scale(1.f / (it));

                /// Write image
                if (it % m_dumpIteration == 0) {
                    develop(scene, film, bitmapAvg, it, "_pass_");
                    develop(scene, film, bitmapVeryDirectAvg, it, "_very_direct_");
                }

                /// Time it
                unsigned int milliseconds = renderingTimer->getMilliseconds();
                timeFile << (milliseconds / (1000.f)) << ",\n";
                timeFile.flush();
                Log(EInfo, "Rendering time: %i, %i", milliseconds / 1000,
                    milliseconds % 1000);
                cumulativeTime += (milliseconds / 1000.f);

                // Note: Exclude filtering in the rendering time
                /// Extract information
                SLog(EInfo, "Filtering...");
                bcd::SamplesStatisticsImages samplesStats = samplesAccumulator.getSamplesStatistics();
                if (true) {
                    bcd::SpikeRemovalFilter::filter(
                            samplesStats.m_meanImage,
                            samplesStats.m_nbOfSamplesImage,
                            samplesStats.m_histoImage,
                            samplesStats.m_covarImage,
                            2.f);
                }

                bcd::DenoiserInputs inputs;
                bcd::DenoiserOutputs outputs;
                bcd::DenoiserParameters parameters;

                inputs.m_pColors = &(samplesStats.m_meanImage);
                inputs.m_pNbOfSamples = &(samplesStats.m_nbOfSamplesImage);
                inputs.m_pHistograms = &(samplesStats.m_histoImage);
                inputs.m_pSampleCovariances = &(samplesStats.m_covarImage);

                bcd::Deepimf outputDenoisedColorImage(samplesStats.m_meanImage);
                outputs.m_pDenoisedColors = &outputDenoisedColorImage;

                parameters.m_histogramDistanceThreshold = 1.f;
                parameters.m_patchRadius = 1;
                parameters.m_searchWindowRadius = 6;
                parameters.m_minEigenValue = 1.e-8f;
                parameters.m_useRandomPixelOrder = true;
                parameters.m_markedPixelsSkippingProbability = 1.f;
                parameters.m_nbOfCores = nCores;
                parameters.m_useCuda = false;

                SLog(EInfo, "Reconstructing...");
                std::unique_ptr<bcd::IDenoiser> uDenoiser = nullptr;
                uDenoiser.reset(new bcd::MultiscaleDenoiser(3)); // 3 scales
                uDenoiser->setInputs(inputs);
                uDenoiser->setOutputs(outputs);
                uDenoiser->setParameters(parameters);
                uDenoiser->denoise();

                unsigned int millisecondsRecons = renderingTimer->getMilliseconds();
                Log(EInfo, "Reconstruction time: %i, %i", millisecondsRecons / 1000,
                    millisecondsRecons % 1000);
                renderingTimer->reset();

                ref<Bitmap> bitmapRecons = bitmap->clone();
                Spectrum *targetRecons = (Spectrum *) bitmapRecons->getUInt8Data();
                Spectrum *targetVeryAvg = (Spectrum *) bitmapVeryDirectAvg->getUInt8Data();
                for (size_t y = 0; y < cropSize.y; y++) {
                    for (size_t x = 0; x < cropSize.x; x++) {
                        Float value[3] = {
                                std::max(outputs.m_pDenoisedColors->get(bcd::PixelPosition(y, x), 0), 0.f),
                                std::max(outputs.m_pDenoisedColors->get(bcd::PixelPosition(y, x), 1), 0.f),
                                std::max(outputs.m_pDenoisedColors->get(bcd::PixelPosition(y, x), 2), 0.f)
                        };
                        targetRecons[y * bitmapRecons->getWidth() + x] = Spectrum(value)
                                                                         + targetVeryAvg[y * bitmapRecons->getWidth() +
                                                                                         x];
                    }
                }
                develop(scene, film, bitmapRecons, it, "_recons_");


                // === Print the statistic at each step of the rendering
                // to see the algorithm behaviors.
                Statistics::getInstance()->printStats();

                film->setBitmap(bitmapRecons);
                queue->signalRefresh(job);
            }
        }

        return true;
    }

    struct Results {
        Spectrum Li = Spectrum(0.f);
        Spectrum LiVery = Spectrum(0.f);
    };

    Results renderPixel(const Point2i offset,
                        const Scene *scene,
                        const Sensor *sensor, Sampler *sampler,
                        const bool &stop,
                        std::vector<Spectrum> &values) const {

        Results result;

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

            auto local_result = Li(sensorRay, rRec);
            local_result.Li *= spec;
            local_result.LiVery *= spec;
            values.emplace_back(local_result.Li);
            result.Li += local_result.Li;
            result.LiVery += local_result.LiVery;

            sampler->advance();
        }


        result.Li /= (Float) sampler->getSampleCount();
        result.LiVery /= (Float) sampler->getSampleCount();


        return result;
    }

    bool hasSmoothComponent(const BSDF *bsdf,
                            Intersection &its) const {

        bool found_smooth = false;
        for (int i = 0, component_count = bsdf->getComponentCount(); i < component_count; ++i) {
            Float component_roughness = bsdf->getRoughness(its, i);
            if (component_roughness != Float(0)) {
                found_smooth = true;
            }
        }
        return found_smooth;
    }

    void cancel() {
        m_stop = true;
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        const Class *cClass = child->getClass();

        if (cClass->derivesFrom(MTS_CLASS(Integrator))) {
            SLog(EError, "Do not support neasted integrators");
        } else {
            Integrator::addChild(name, child);
        }
    }

    Results Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        /* Some aliases and local variables */
        const Scene *scene = rRec.scene;
        Intersection &its = rRec.its;
        MediumSamplingRecord mRec;
        RayDifferential ray(r);
        Results result;
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

        while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
            /* ==================================================================== */
            /*                 Radiative Transfer Equation sampling                 */
            /* ==================================================================== */
            if (rRec.medium && rRec.medium->sampleDistance(Ray(ray, 0, its.t), mRec, rRec.sampler)) {
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
                                result.Li += throughput * value * phaseVal * weight;
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
                        result.Li += throughput * value * miWeight(phasePdf, emitterPdf);
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

                if(!mediumScattered) {
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
                        result.Li += value;
                    }

                    break;
                }

                /* Possibly include emitted radiance if requested */
                if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
                    && ((!m_hideEmitters || scattered) || mediumScattered)) {
                    result.Li += throughput * its.Le(-ray.d);
                }


                /* Include radiance from a subsurface integrator if requested */
                if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
                    result.Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

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
                                    result.Li += throughput * value * bsdfVal * weight;
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
                        result.Li += throughput * value * miWeight(bsdfPdf, emitterPdf);
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
        return result;
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

    inline Float miWeight(Float pdfA, Float pdfB) const {
        pdfA *= pdfA;
        pdfB *= pdfB;
        return pdfA / (pdfA + pdfB);
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
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()

protected:
    int m_maxPasses;
    float m_maxRenderingTime;
    bool m_stop;
    int m_dumpIteration;

    std::vector<std::vector<Spectrum>> m_values;

    // Path tracer configuration
    int m_maxDepth;
    int m_rrDepth;
    bool m_strictNormals;
    bool m_hideEmitters;
    bool m_noRayDiff;
};

MTS_IMPLEMENT_CLASS(BCDIntegrator, false, Integrator)

MTS_EXPORT_PLUGIN(BCDIntegrator, "Naive BCD integrator");
MTS_NAMESPACE_END