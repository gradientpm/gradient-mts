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

#include <mitsuba/bidir/vertex.h>
#include <mitsuba/bidir/edge.h>
#include "gdvcm_proc.h"
#include "mitsuba/bidir/util.h"
#include "../reconstruction.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{gbdpt}{Gradient-Domain Bidirectional path tracer}
 * \order{5}
 * \parameters{
 *     \parameter{maxDepth}{\Integer}{Specifies the longest path depth
 *         in the generated output image (where \code{-1} corresponds to $\infty$).
 *	       A value of \code{1} will only render directly visible light sources.
 *	       \code{2} will lead to single-bounce (direct-only) illumination,
 *	       and so on. \default{\code{-1}}
 *	   }
 *	   \parameter{lightImage}{\Boolean}{Include sampling strategies that connect
 *	      paths traced from emitters directly to the camera? (i.e. what \pluginref{ptracer} does)
 *	      This improves the effectiveness of bidirectional path tracing
 *	      but severely increases the local and remote communication
 *	      overhead, since large \emph{light images} must be transferred between threads
 *	      or over the network. See the text below for a more detailed explanation.
 *	      \default{include these strategies, i.e. \code{true}}
 *     }
 *	   \parameter{rrDepth}{\Integer}{Specifies the minimum path depth, after
 *	      which the implementation will start to use the ``russian roulette''
 *	      path termination criterion. \default{\code{5}}
 *	   }
 *	   \parameter{reconstructL1}{\Boolean}{If set, the rendering method reconstructs the final image using a reconstruction method
 *           that efficiently kills many image artifacts. The reconstruction is slightly biased, but the bias will go away by increasing sample count. \default{\code{true}}
 *     }
 *	   \parameter{reconstructL2}{\Boolean}{If set, the rendering method reconstructs the final image using a reconstruction method that is unbiased,
 *			but sometimes introduces severe dipole artifacts. \default{\code{false}}
 *     }
 *	   \parameter{shiftThreshold}{\Float}{Specifies the roughness threshold for classifying materials as 'diffuse', in contrast to 'specular',
 *			for the purposes of constructing paths pairs for estimating pixel differences. This value should usually be somewhere between 0.0005 and 0.01.
 *			If the result image has noise similar to standard path tracing, increasing or decreasing this value may sometimes help. This implementation assumes that this value is small.\default{\code{0.001}}
 *	   }
 *	   \parameter{reconstructAlpha}{\Float}{
 *			Higher value makes the reconstruction trust the noisy color image more, giving less weight to the usually lower-noise gradients.
 *			The optimal value tends to be around 0.2, but for scenes with much geometric detail at sub-pixel level a slightly higher value such as 0.3 or 0.4 may be tried.\default{\code{0.2}}
 }
 * }
 *
 *
 * This plugin implements a gradient-domain bidirectional path tracer (short: G-BDPT)
 * as described in the paper "Gradient-Domain Bidirectional Path Tracing" by Manzi et al.. 
 * It supports classical materials like diffuse, specular and glossy materials, depth-of-field, 
 * motion blur and low discrepancy samplers. However it is still an experimental 
 * implementation that hasn't been tested with all of Mitsubas features.
 * Notably there is no support yet for any kind of participating media or pixel filters 
 * beside the box filter. Also smart sampling of the direct illumination is not 
 * implemented (i.e. no sampleDirect option as in BDPT). G-BDPT only works with \pluginref{multifilm}
 * since multiple output files are generated. Note that in the GUI MultiFilm is automatically 
 * chosen if G-BDPT is selected.
 *
 */
class GDVCMIntegrator : public VCMIntegratorBase {
public:

    GDVCMIntegrator(const Properties &props) : VCMIntegratorBase(props) {
        /* Load the parameters / defaults */
        m_config.maxDepth = props.getInteger("maxDepth", -1);
        m_config.rrDepth = props.getInteger("rrDepth", 5);
        m_config.initialRadius = props.getFloat("initialRadius", 0.f);
        m_config.radiusReductionAlpha = props.getFloat("radiusReductionAlpha", 0.9f);
        m_config.lightImage = props.getBoolean("lightImage", true);
        m_config.phExponent = props.getFloat("phExponent", 1.0);
        m_config.m_shiftThreshold = props.getFloat("shiftThreshold", Float(0.001));
        m_config.directTracing = props.getBoolean("directTracing", true);

        m_config.m_reconstructL1 = props.getBoolean("reconstructL1", false); //no matter what is chosen, both reconstructions are written to disc. this only changes which one is shown in the GUI
        m_config.m_reconstructL2 = props.getBoolean("reconstructL2", true);
		m_config.m_reconstructUni = props.getBoolean("reconstructUni", false);
        m_config.m_reconstructAlpha = (Float) props.getFloat("reconstructAlpha", Float(0.2));
        m_config.mergeOnly = props.getBoolean("mergeOnly", false);
        m_config.forceBlackPixels = props.getBoolean("forceBlackPixels", true);
		m_config.maxManifoldIterations = props.getInteger("maxManifoldIterations", -1);;

        if (m_config.mergeOnly) m_config.lightImage = false;

        if (m_config.m_reconstructAlpha <= 0.0f)
            Log(EError, "'reconstructAlpha' must be set to a value greater than zero!");

        m_config.nNeighbours = 4; //fixed at the moment

        if (m_config.rrDepth <= 0)
            Log(EError, "'rrDepth' must be set to a value greater than zero!");

        if (m_config.maxDepth <= 0 && m_config.maxDepth != -1)
            Log(EError, "'maxDepth' must be set to -1 (infinite) or a value greater than zero!");

        if (m_config.maxDepth == -1) {
            Log(EWarn, "maxDepth is unlimited, set to 12!");
            m_config.maxDepth = 12;
        }
    }

    /// Unserialize from a binary data stream

    GDVCMIntegrator(Stream *stream, InstanceManager *manager)
    : VCMIntegratorBase(stream, manager) {
        m_config = GDVCMConfiguration(stream);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Integrator::serialize(stream, manager);
        m_config.serialize(stream);
    }

    bool preprocess(const Scene *scene, RenderQueue *queue,
            const RenderJob *job, int sceneResID, int sensorResID,
            int samplerResID) {
        Integrator::preprocess(scene, queue, job, sceneResID,
                sensorResID, samplerResID);

        if (scene->getSubsurfaceIntegrators().size() > 0)
            Log(EError, "Subsurface integrators are not supported "
                "by the bidirectional path tracer!");
        if (m_config.initialRadius == 0) {
            /* Guess an initial radius if not provided
              (scene width / horizontal or vertical pixel count) * 5 */
            Float rad = scene->getBSphere().radius;
            Vector2i filmSize = scene->getSensor()->getFilm()->getSize();

            m_config.initialRadius = std::min(rad / filmSize.x, rad / filmSize.y);
        }

        return true;
    }

    void cancel() {
        m_cancelled = true;
        Scheduler::getInstance()->cancel(m_process);
    }

    void configureSampler(const Scene *scene, Sampler *sampler) {
        /* Prepare the sampler for tile-based rendering */
        sampler->setFilmResolution(scene->getFilm()->getCropSize(), true);
    }

    bool render(Scene *scene, RenderQueue *queue, const RenderJob *job,
            int sceneResID, int sensorResID, int samplerResID) {

        ref<Scheduler> scheduler = Scheduler::getInstance();
        ref<Sensor> sensor = scene->getSensor();
        ref<Film> film = sensor->getFilm();
        size_t sampleCount = scene->getSampler()->getSampleCount();
        size_t nCores = scheduler->getCoreCount();

		/* Create CSV file to dump all the rendering timings */
		// Also create an timer
		std::string timeFilename = scene->getDestinationFile().string()
								   + "_time.csv";
		std::ofstream timeFile(timeFilename.c_str());
		ref<Timer> renderingTimer = new Timer;

        Log(EDebug, "Size of data structures: PathVertex=%i bytes, PathEdge=%i bytes",
                (int) sizeof (PathVertex), (int) sizeof (PathEdge));

        Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT " samples, " SIZE_T_FMT
                " %s, " SSE_STR ") ..", film->getCropSize().x, film->getCropSize().y,
                sampleCount, nCores, nCores == 1 ? "core" : "cores");

        m_config.blockSize = scene->getBlockSize();
        m_config.cropSize = film->getCropSize();
        m_config.extraBorder = 0;
        m_config.dump();

        /* setup MultiFilm */
        std::vector<std::string> outNames;

        outNames.push_back(m_config.m_reconstructL1 ? "-L1" : "-L2");
        outNames.push_back("-gradientNegY");
        outNames.push_back("-gradientNegX");
        outNames.push_back("-gradientPosX");
        outNames.push_back("-gradientPosY");
#if defined(SEPARATE_DIRECT)
        outNames.push_back("-direct");
#else
        outNames.push_back(m_config.m_reconstructL1 ? "-L2" : "-L1");
#endif
        outNames.push_back("-primal");
        
        if (!film->setBuffers(outNames)) {
            Log(EError, "Cannot render image! G-VCM has been called without MultiFilm.");
            return false;
        }

        ref<Scene> bidir_scene = new Scene(scene);
        bidir_scene->initializeBidirectional();
        int bidirSceneResID = scheduler->registerResource(bidir_scene);

        /* run job */
        ref<GDVCMProcess> process = new GDVCMProcess(job, queue, m_config);
        m_process = process;
        process->bindResource("scene", bidirSceneResID);
        process->bindResource("sensor", sensorResID);
        process->bindResource("sampler", samplerResID);

#if defined(MTS_OPENMP)
        Thread::initializeOpenMP(nCores);
#endif
        m_cancelled = false;
		int currentIteration = 1;
		while(true) {
			// Change the number of sample count because we want to be able
			// normalize correctly the light image.
			process->sampleCount = currentIteration*sampleCount;
			size_t offset_iteration = (currentIteration-1)*sampleCount;
			for (size_t i = 0; i < sampleCount; i++) {
				if (!iterateVCM(process, sensorResID, i + offset_iteration)) {
					return false;
				}
			}
			process->develop();

			/* data buffers for solver */
			int len = 3 * film->getCropSize().x * film->getCropSize().y;
			std::vector<float> imgf(len, 0.f);
			std::vector<float> dxf(len, 0.f);
			std::vector<float> dyf(len, 0.f);
			std::vector<float> direct(len, 0.f);
			std::vector<Float *> grad(m_config.nNeighbours);

			/* allocate temporary bitmaps for solvers */
			ref<Bitmap> imgBaseBuff, recBuff, errBuff;
			std::vector<ref<Bitmap> > gradBuff;
			imgBaseBuff = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getCropSize());
			recBuff = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getCropSize());
			errBuff = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getCropSize());
			for (int j = 0; j < m_config.nNeighbours; j++) {
				gradBuff.push_back(new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getCropSize()));
			}

			/* develop primal and gradient data into bitmaps */
			film->developMulti(Point2i(0, 0), film->getCropSize(), Point2i(0, 0), imgBaseBuff, 0);
			for (int j = 1; j <= m_config.nNeighbours; j++)
				film->developMulti(Point2i(0, 0), film->getCropSize(), Point2i(0, 0), gradBuff[j - 1], j);

			/* read data from bitmaps */
			for (int i = 0; i < m_config.nNeighbours; i++)
				grad[i] = gradBuff[i]->getFloatData();

			/* prepare data for solvers */
			prepareDataForSolver(1.f, &imgf[0], imgBaseBuff->getFloatData(), len, NULL, 0);
			prepareDataForSolver(1.f, &dyf[0], grad[3], len, grad[0], film->getCropSize().x);
			prepareDataForSolver(1.f, &dxf[0], grad[2], len, grad[1], 1);
#if defined(SEPARATE_DIRECT)
			if(m_config.directTracing) {
				ref<Bitmap> imgDirectBuff = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getCropSize());
				film->developMulti(Point2i(0, 0), film->getCropSize(), Point2i(0, 0), imgDirectBuff, 5);
				prepareDataForSolver(1.f, &direct[0], imgDirectBuff->getFloatData(), len, NULL, 0);
			}
#endif
			// Allocated bitmaps
			ref<Bitmap> dxBitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getCropSize());
			ref<Bitmap> dyBitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getCropSize());
			ref<Bitmap> imgBitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getCropSize());

			// Copy values
			setBitmapFromArray(dxBitmap, &dxf[0], true); // use abs
			setBitmapFromArray(dyBitmap, &dyf[0], true); // use abs
			setBitmapFromArray(imgBitmap, &imgf[0]);

			// Dump extra information
			Properties pHDRFilm("hdrfilm");
			pHDRFilm.setInteger("width", film->getCropSize().x);
			pHDRFilm.setInteger("height", film->getCropSize().y);
			ref<Film> hdrFilm = static_cast<Film *> (PluginManager::getInstance()->
					createObject(MTS_CLASS(Film), pHDRFilm));
			develop(scene, hdrFilm, dxBitmap, currentIteration, "_dxAbs_");
			develop(scene, hdrFilm, dyBitmap, currentIteration, "_dyAbs_");
			develop(scene, hdrFilm, imgBitmap, currentIteration, "_throughput_");

			/* need to put primal img back into film such that it can be written to disc */
			film->setBitmapMulti(imgBaseBuff, 1, m_config.nNeighbours + 2);
			unsigned int milliseconds = renderingTimer->getMilliseconds();
			timeFile << (milliseconds / (1000.f)) << ",\n";
			timeFile.flush();
			Log(EInfo, "Rendering time: %i, %i", milliseconds / 1000,
				milliseconds % 1000);
			renderingTimer->reset();

			////////////////////////////////////////////////////
			////////////////////////////////////////////////////
			// Perform the reconstruction
			Reconstruction rec {};
			rec.reconstructL1 = m_config.m_reconstructL1;
			rec.reconstructL2 = m_config.m_reconstructL2;
			rec.reconstructUni = m_config.m_reconstructUni;
			rec.alpha = (float)m_config.m_reconstructAlpha;
			{
				auto rec_results = rec.reconstruct(film->getCropSize(),
												   imgf, dxf, dyf, direct,
                                                   Reconstruction::Variance {},
												   PostProcessOption{
														   forceBlackPixels: m_config.forceBlackPixels,
														   clampingValues: true
												   });
				renderingTimer->reset();
				if (rec_results.size() == 1) {
					develop(scene, hdrFilm, rec_results[0].img, currentIteration, "_recons_");
				} else {
					for (auto &result: rec_results) {
						develop(scene, hdrFilm, result.img, currentIteration, "_"+result.name+"_");
					}
				}
			}

			// In case something wrong happens
			if (!process->getReturnStatus() == ParallelProcess::ESuccess)
				return false;

			++currentIteration;
		}

		scheduler->unregisterResource(bidirSceneResID);

        return true;
    }

    void prepareDataForSolver(float w, float* out, Float * data, int len, Float *data2, int offset) {

        for (int i = 0; i < len; i++)
            out[i] = w * float(data[i]);

        //used to merge inverse directions into one buffer
        if (data2 != NULL) {
            int io;
            for (int i = 0; i < len; i++) {
                io = i + 3 * offset;
                if (io >= 0 && io < len) {
                    out[i] *= 0.5;
                    out[i] -= 0.5 * w * float(data2[io]);
                }
            }
        }
    }

	void setBitmapFromArray(ref<Bitmap> &bitmap, float *img, bool useAbs = false, bool useClamp = false)
	{
		int x, y;
		Float tmp[3];
		for (int i = 0; i < bitmap->getSize().x*bitmap->getSize().y; i++){
			y = i / bitmap->getSize().x;
			x = i - y*bitmap->getSize().x;
			if(useAbs) {
				tmp[0] = std::abs(Float(img[3 * i]));
				tmp[1] = std::abs(Float(img[3 * i + 1]));
				tmp[2] = std::abs(Float(img[3 * i + 2]));
			} else if(useClamp) {
				tmp[0] = std::max(Float(img[3 * i]), Float(0.f));
				tmp[1] = std::max(Float(img[3 * i + 1]), Float(0.f));
				tmp[2] = std::max(Float(img[3 * i + 2]), Float(0.f));
			} else {
				tmp[0] = Float(img[3 * i]);
				tmp[1] = Float(img[3 * i + 1]);
				tmp[2] = Float(img[3 * i + 2]);
			}
			bitmap->setPixel(Point2i(x, y), Spectrum(tmp));
		}
	}
	void develop(Scene* scene, Film* film, Bitmap* bitmap,
				 int currentIteration, const std::string& suffixName = "_") {
		std::stringstream ss;
		ss << scene->getDestinationFile().string() << suffixName
		   << currentIteration;
		std::string path = ss.str();

		film->setBitmap(bitmap);
		film->setDestinationFile(path, 0);
		film->develop(scene, 0.f);

	}

    MTS_DECLARE_CLASS()
private:
    ref<ParallelProcess> m_process;
    GDVCMConfiguration m_config;
};

MTS_IMPLEMENT_CLASS_S(GDVCMIntegrator, false, Integrator)
MTS_EXPORT_PLUGIN(GDVCMIntegrator, "Gradient Domain VCM");
MTS_NAMESPACE_END
