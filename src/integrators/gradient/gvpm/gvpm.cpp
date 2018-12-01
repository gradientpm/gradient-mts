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
#include <mitsuba/core/statistics.h>
#include <mitsuba/render/gatherproc.h>
#include <mitsuba/render/renderqueue.h>
#include <mitsuba/bidir/manifold.h>

#include "gvpm_proc.h"
#include "gvpm_gatherpoint.h"

#include "gvpm_plane.h"
#include "../../photonmapper/plane_accel.h"

#include "shift/shift_surface.h"
#include "shift/shift_volume_photon.h"
#include "shift/shift_volume_beams.h"
#include "shift/shift_volume_planes.h"

// Solver
#include "../reconstruction.h"
#include "../../denoiser/nfor/nfor.h"

// For get the same ray diff heuristic
#include "../../photonmapper/utilities/tools.h"

// Include tools shared between SPPM and GVPM

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
Float VertexClassifier::roughnessThreshold;

class GPMIntegrator : public Integrator {
public:
	explicit GPMIntegrator(const Properties &props) : Integrator(props) {
		m_config.load(props);

		if (m_config.volTechnique == EBeamBeam1D) {
			m_config.newShiftBeam = true; // Use the correct fix here
		}

		// Generate the gather point in the similar way for all PPM/SPPM codes
		m_mutex = new Mutex();		

		// Based on the "optimization" for independence
		// between number of samples and alpha.
		// See "Progressive Photon Beams" section 5.2
		m_independentScale = false; // Put at false as the scaling this too strong
	}

	GPMIntegrator(Stream *stream, InstanceManager *manager)
			: Integrator(stream, manager) {}

	~GPMIntegrator() override {
		
	}

	void serialize(Stream *stream, InstanceManager *manager) const override {
		Integrator::serialize(stream, manager);
		Log(EError, "Network rendering is not supported!");
	}

	void cancel() override {
		m_running = false;
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
					int sceneResID, int sensorResID, int samplerResID) override {
		Integrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID);

		if (m_config.fluxLightDensity)
			const_cast<Scene *>(scene)->weightEmitterFlux();

		if (scene->getMedia().size() == 0) {
			m_config.volumePhotonCount = 0;
			m_config.lightingInteractionMode = ESurf2Surf;
		} else {
			m_smokeAABB = max_AABB_medium(scene);
		}

		if (!m_config.needSurfaceRendering()) {
			Log(EInfo, "Only volume rendering, change probability to sample the medium");
			auto medium = scene->getMedia().begin();
			while (medium != scene->getMedia().end()) {
				const_cast<Medium *>(medium->get())->computeOnlyVolumeInteraction();
				medium++;
			}
		}

#if 1
		///////////////////////////
		// Check scene requirement

		// 1) Check if all the lights are Area ones
		if (m_config.minDepth < 2) {
			auto emitters = scene->getEmitters();
			for (auto e: emitters) {
				std::string emitterType = e->getClass()->getName();
				SLog(EInfo, "Emitter type: %s", emitterType.c_str());
				if (emitterType != "AreaLight") {
					SLog(EError, "Unsupported light type: %s", emitterType.c_str());
				}
			}
		}
#endif

		// Scale update the cameraSphere
		m_config.cameraSphere = m_smokeAABB.getBSphere().radius * m_config.cameraSphere * POURCENTAGE_BS;

		if (m_config.isPlaneEstimator()) {
			if (m_config.minDepth < 2) {
				SLog(EError, "Impossible to use plane with minDepth smaller than 2");
			}

			// Remove the dimentionality for the planes
			// TODO: This not working if the camera path can bounce
			m_config.minDepth -= 1;
			if (m_config.maxDepth) {
				m_config.maxDepth -= 1;
			}

		}

		return true;
	}

	void scaleVolumeAPA(int it) {
		// Scale the gathering size for the next iteration

		it -= 1; // Fix the bug as it == 1 at the first iteration.
		double ratioVolAPA = 1.0;

		if (!m_independentScale) {
			ratioVolAPA = (it + m_config.alpha) / (it + 1);
		} else {
			for (int i = 0; i < m_config.volumePhotonCount; i++) {
				size_t newIt = m_config.volumePhotonCount * it + i;
				ratioVolAPA *= (newIt + m_config.alpha) / (newIt + 1);
			}
		}

		if (m_config.forceAPA.empty()) {
			if (EVolumeTechniqueHelper::use3DKernel(m_config.volTechnique) || m_config.use3DKernelReduction) {
				globalScaleVolume[m_bufferIndex] = globalScaleVolume[m_bufferIndex] * std::cbrt(ratioVolAPA);
			} else if (EVolumeTechniqueHelper::use2DKernel(m_config.volTechnique)) {
				globalScaleVolume[m_bufferIndex] = globalScaleVolume[m_bufferIndex] * std::sqrt(ratioVolAPA);
			} else {
				globalScaleVolume[m_bufferIndex] = globalScaleVolume[m_bufferIndex] * ratioVolAPA;
			}
		} else {
			if (m_config.forceAPA == "1D") {
				globalScaleVolume[m_bufferIndex] = globalScaleVolume[m_bufferIndex] * ratioVolAPA;
			} else if (m_config.forceAPA == "2D") {
				globalScaleVolume[m_bufferIndex] = globalScaleVolume[m_bufferIndex] * std::sqrt(ratioVolAPA);
			} else if (m_config.forceAPA == "3D") {
				globalScaleVolume[m_bufferIndex] = globalScaleVolume[m_bufferIndex] * std::cbrt(ratioVolAPA);
			} else {
				SLog(EError, "No Force APA: %s", m_config.forceAPA.c_str());
			}
		}
	}

	void allocOffset(Scene *scene, std::vector<Point2i> &offset) {
		offset.clear();
		int blockSize = scene->getBlockSize();
		Vector2i size = scene->getFilm()->getSize();		
		for (int yofs = 0; yofs < size.y; yofs += blockSize) {
			for (int xofs = 0; xofs < size.x; xofs += blockSize) {				
				m_offset.emplace_back(Point2i(xofs, yofs));
			}
		}
	}

	void allocGatherBlocks(Scene *scene, GatherBlocks &gatherBlocks) {
		int blockSize = scene->getBlockSize();
		Vector2i size = scene->getFilm()->getSize();
		gatherBlocks.clear();
		for (int yofs = 0; yofs < size.y; yofs += blockSize) {
			for (int xofs = 0; xofs < size.x; xofs += blockSize) {
				gatherBlocks.emplace_back(std::vector<GatherPoint>());				
				std::vector<GatherPoint> &gatherPoints = gatherBlocks[gatherBlocks.size() - 1];
				int nPixels = std::min(blockSize, size.y - yofs)
							* std::min(blockSize, size.x - xofs);
				gatherPoints.reserve(nPixels);

				// Generate the gather points
				for (int i = 0; i < nPixels; ++i) {
					GatherPoint newGP;
					newGP.radius = 0;
					newGP.scale = m_config.initialScale;
					newGP.scaleVol = m_config.initialScaleVolume;
					gatherPoints.emplace_back(std::move(newGP));
				}
			}
		}
	}

	void allocImageBlocks(Scene *scene, GatherBlocks &gatherBlocks, std::vector<std::vector<GatherPoint *>> &imgGP) {
		// 2D maps that allows iteration over all pixels.
		// We need this because gatherBlocks are stored as blocks, not row-major pixels in memory.
		// We can only build this map once after the gather blocks are generated. 				
		Vector2i size = scene->getFilm()->getSize();
		imgGP.resize(size.x);
		for (int x = 0; x < size.x; ++x) {
			imgGP[x].resize(size.y);
		}
		for (size_t idBlock = 0; idBlock < gatherBlocks.size(); idBlock++) {
			std::vector<GatherPoint> &currBlock = gatherBlocks[idBlock];
			for (GatherPoint &gp: currBlock) {
				imgGP[gp.pixel.x][gp.pixel.y] = &gp;
			}
		}
	}

	void averageAndVariance(Bitmap *bitmapA, Bitmap *bitmapB, Bitmap *meanBitmap, Bitmap *varianceBitmap) {
		Vector2i size = bitmapA->getSize();
		Spectrum *mean = (Spectrum *)meanBitmap->getFloatData();
		Spectrum *variance = (Spectrum *)varianceBitmap->getFloatData();
		Spectrum *A = (Spectrum *)bitmapA->getFloatData();
		Spectrum *B = (Spectrum *)bitmapB->getFloatData();
		for (int y = 0; y < size.y; ++y) {
			for (int x = 0; x < size.x; ++x) {
				int index = y * size.x + x;
				mean[index] = (A[index] + B[index]) * 0.5;
				variance[index] = ((A[index] - B[index]) * (A[index] - B[index])) * 0.25;
			}
		}
	}

	void auxiliaryPass(Scene *scene, int N) {
		ref<Scheduler> sched = Scheduler::getInstance();
		size_t nCores = sched->getCoreCount();
		BlockScheduler blockSched(m_gatherBlocks->size(), nCores);

		blockSched.run([this, scene, N](int blockIdx, int tid) {
			GPMThreadData *threadData = m_threadData[tid];
			std::vector<GatherPoint> &gatherPoints = (*m_gatherBlocks)[blockIdx];
			for (size_t i = 0; i < gatherPoints.size(); ++i) {
				GatherPoint &gp = gatherPoints[i];
				int lastVertexID = gp.path.vertexCount() - 1;
				if (lastVertexID > 1) {
					PathVertex *vertex = gp.path.vertex(lastVertexID);
					Intersection &its = vertex->getIntersection();
					if (!its.isEmitter()) {
						const Spectrum reflectance = its.getBSDF()->getDiffuseReflectance(its);
						gp.albedo = (gp.albedo * (N - 1) + reflectance) / N;
					}

					Intersection &primaryIts = gp.path.vertex(2)->getIntersection();
					const Transform &camera2world = scene->getSensor()->getWorldTransform()->eval(0.f);
					const Transform world2camera = camera2world.inverse();
					const Point cameraCoord = world2camera(primaryIts.p);
					gp.position = (gp.position * (N - 1) + Vector(cameraCoord.z, cameraCoord.z, cameraCoord.z)) / N; // depth only

					if (!primaryIts.isEmitter()) {
						Frame shFrame = primaryIts.getBSDF()->getFrame(primaryIts);
						gp.shadingNormal = (gp.shadingNormal * (N - 1) + Vector(shFrame.n.x, shFrame.n.y, shFrame.n.z)) / N;
					}
				}
			}
		});
	}

	bool render(Scene *scene, RenderQueue *queue,
				const RenderJob *job, int sceneResID, int sensorResID, int unused) override {
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
		
		m_running = true;
		m_totalEmittedSurface[0] = m_totalEmittedSurface[1] = 0;
		m_totalEmittedVolume[0] = m_totalEmittedVolume[1] = 0;
		m_totalPhotons[0] = m_totalPhotons[1] = 0;

		ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(Sampler), Properties("independent")));
		
		/* Create CSV file to dump all the rendering timings */
		// Also create an timer
		std::string timeFilename = scene->getDestinationFile().string()
								   + "_time.csv";
		SLog(EInfo, "Time file %s", timeFilename.c_str());
		std::ofstream timeFile(timeFilename.c_str());
		ref<Timer> renderingTimer = new Timer;

		// Accumulation buffers with variance tracking		
		bufferAB = (m_config.reconstructNfor || m_config.reconstructWeighted);
				
		// Allocate the special accumulation buffer
		// This buffer will contains light transport that we do not have added
		// in the gradient computation
		m_extraRendering = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, cropSize);
		m_extraRendering->clear();

		// Allocate gather points for each image block
		allocOffset(scene, m_offset);				
		allocGatherBlocks(scene, m_gatherBlocksAB[0]);
		allocGatherBlocks(scene, m_gatherBlocksAB[1]);
		
		// Initialize the class responsible to the GP generation
		// and the radii initialization		
		m_gpManagerAB[0].init(m_config, scene, m_gatherBlocksAB[0], m_offset);
		m_gpManagerAB[1].init(m_config, scene, m_gatherBlocksAB[1], m_offset);

		globalScaleVolume[0] = globalScaleVolume[1] = m_config.initialScaleVolume;

		// Initialize memory for thread objects
		m_threadData.resize(nCores);
		for (size_t i = 0; i < nCores; i++) {
			m_threadData[i] = new GPMThreadData(scene, m_config);
		}

		/* Create a sampler instance for every core */
		/* and a memory pool */
		const size_t preinitLightPath = (m_config.volumePhotonCount / nCores);
		std::vector<ShootingThreadData> shootingData(nCores);
		for (size_t i = 0; i < nCores; ++i) {
			shootingData[i].sampler = sampler->clone();
			shootingData[i].paths.reserve(preinitLightPath);
			for (size_t k = 0; k < preinitLightPath; k++) {
				shootingData[i].paths.push_back(new Path);
			}
		}

		int samplesPerPixel = sampler->getSampleCount();

#ifdef MTS_DEBUG_FP
		enableFPExceptions();
#endif
		int itGlobal = 1;
		int it = 1;
		while (m_running && (m_config.maxPasses == -1 || it < m_config.maxPasses)) {			
			for (int local_it = 0; local_it < samplesPerPixel; local_it++) {								
				m_bufferIndex = bufferAB ? (itGlobal - 1) % 2 : 0;
				m_gatherBlocks = &m_gatherBlocksAB[m_bufferIndex];
				m_gpManager = &m_gpManagerAB[m_bufferIndex];

				Log(EInfo, "Pass %d / %d", it, m_config.maxPasses);
				Log(EInfo, "Regenerating gather points positions and radius!");
				m_gpManager->regeneratePositionAndRadius();
				Log(EInfo, "Done regenerating!");
				m_gpManager->rescaleFlux();

				if (itGlobal == 1) {
					// FIXME: Change this weird behavior. Indeed gp.pixel is set during the gp tracing...
					/* Create a copy of the gather point using image space,
					 * to be able to retrieve them easily
					 */
					allocImageBlocks(scene, m_gatherBlocksAB[0], m_imgGPAB[0]);
				} else if (itGlobal == 2) {
					allocImageBlocks(scene, m_gatherBlocksAB[1], m_imgGPAB[1]);
				}

				m_imgGP = &m_imgGPAB[m_bufferIndex];

				// GPM
				photonMapPass(scene, bufferAB ? (itGlobal + 1) / 2 : itGlobal, queue, job, film, sceneResID,
							  sensorResID, shootingData);

				// Collect auxiliary buffer 
				if (m_config.estimateGradients == false && m_config.reconstructNfor)
					auxiliaryPass(scene, bufferAB ? (itGlobal + 1) / 2 : itGlobal);
				
				//film->setBitmap(m_bitmap);
				queue->signalRefresh(job);

				++itGlobal;
			}

			if (it % m_config.dumpIteration == 0) {
				// Update the log time
				unsigned int milliseconds = renderingTimer->getMilliseconds();
				timeFile << (milliseconds / 1000.f) << ",\n";
				timeFile.flush();
				Log(EInfo, "Rendering time: %i, %i", milliseconds / 1000,
					milliseconds % 1000);				
				renderingTimer->reset();

				int currentIteration = it;

				ref<Film> hdrFilm = newHDRFilm(scene);

				if (m_config.estimateGradients == false) {
					// Primal only mode
					
					ref<Bitmap> throughputBitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, cropSize);					
					computeThroughput(m_bufferIndex, itGlobal-1, cropSize, throughputBitmap.get());
					develop(scene, film, throughputBitmap, it);

					if (m_config.reconstructNfor) {
						// Perform NFOR with aux buffer
						Vector2i size = film->getSize();
						ref<Bitmap> throughputBitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, size);					
						ref<Bitmap> throughputBitmapVariance = throughputBitmap->clone();	
						ref<Bitmap> throughputBitmapA = throughputBitmap->clone();
						ref<Bitmap> throughputBitmapB = throughputBitmap->clone();
						computeThroughput(0, itGlobal / 2, size, throughputBitmapA);
						computeThroughput(1, itGlobal / 2, size, throughputBitmapB);
						averageAndVariance(throughputBitmapA, throughputBitmapB, throughputBitmap, throughputBitmapVariance);
						ref<Bitmap> directBitmap = throughputBitmap->clone();
						computeVeryDirect(m_bufferIndex, itGlobal / 2, size, directBitmap);	
						
						ref<AccumBuffer> throughputBuffer = new AccumBuffer(
										{throughputBitmap},
										{throughputBitmapVariance},
										{throughputBitmapA},
										{throughputBitmapB},
										currentIteration
						);

						ref<Bitmap> albedoBitmapA = throughputBitmap->clone();
						ref<Bitmap> albedoBitmapB = throughputBitmap->clone();
						ref<Bitmap> shadingNormalBitmapA = throughputBitmap->clone();
						ref<Bitmap> shadingNormalBitmapB = throughputBitmap->clone();
						ref<Bitmap> positionBitmapA = throughputBitmap->clone();
						ref<Bitmap> positionBitmapB = throughputBitmap->clone();
						computeAuxBitmaps(0, albedoBitmapA, positionBitmapA, shadingNormalBitmapA);
						computeAuxBitmaps(1, albedoBitmapB, positionBitmapB, shadingNormalBitmapB);
						ref<Bitmap> albedoBitmap = throughputBitmap->clone();
						ref<Bitmap> albedoVarianceBitmap = throughputBitmap->clone();
						ref<Bitmap> shadingNormalBitmap = throughputBitmap->clone();
						ref<Bitmap> shadingNormalVarianceBitmap = throughputBitmap->clone();
						ref<Bitmap> positionBitmap = throughputBitmap->clone();
						ref<Bitmap> positionVarianceBitmap = throughputBitmap->clone();
						averageAndVariance(albedoBitmapA, albedoBitmapB, albedoBitmap, albedoVarianceBitmap);
						averageAndVariance(positionBitmapA, positionBitmapB, positionBitmap, positionVarianceBitmap);
						averageAndVariance(shadingNormalBitmapA, shadingNormalBitmapB, shadingNormalBitmap, shadingNormalVarianceBitmap);
						
						ref<AccumBuffer> auxBuffer = new AccumBuffer(
										{albedoBitmap, positionBitmap, shadingNormalBitmap},
										{albedoVarianceBitmap, positionVarianceBitmap, shadingNormalVarianceBitmap},
										{albedoBitmapA, positionBitmapA, shadingNormalBitmapA},
										{albedoBitmapB, positionBitmapB, shadingNormalBitmapB},
										currentIteration
						);

						// Use all buffers
						std::vector<int> bufferIndices(auxBuffer->getBufferCount());
						for (int i = 0; i < bufferIndices.size(); ++i) bufferIndices[i] = i;

						{
							ref<Nfor> nfor = new Nfor();
							ref<Bitmap> nforBitmap = throughputBitmap->clone();
							nfor->denoise(throughputBuffer, nforBitmap, auxBuffer, bufferIndices);							
							if (m_config.directTracing) {
								nforBitmap->accumulate(directBitmap);
							}
							develop(scene, hdrFilm, nforBitmap, currentIteration, "_nfor_");
							develop(scene, hdrFilm, albedoBitmap, it, "_albedo_");
							develop(scene, hdrFilm, positionBitmap, it, "_position_");
							develop(scene, hdrFilm, shadingNormalBitmap, it, "_shadingNormal_");							
						}

					}

				} else {
								
					// Traditional reconstruction
					if (m_config.reconstructL1 || m_config.reconstructL2 ||
						m_config.reconstructUni)
					{
						Reconstruction rec {};
						rec.reconstructL1 = m_config.reconstructL1;
						rec.reconstructL2 = m_config.reconstructL2;
						rec.reconstructUni = m_config.reconstructUni;
						rec.reconstructWeighted = false;
						rec.alpha = (float)m_config.reconstructAlpha;
						
						/* Create bitmap to display gradient images */
						ref<Bitmap> throughputBitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, cropSize);					
						computeThroughput(m_bufferIndex, itGlobal-1, cropSize, throughputBitmap.get());
						ref<Bitmap> gXBitmap = throughputBitmap->clone();
						ref<Bitmap> gYBitmap = throughputBitmap->clone();
						computeGradient(m_bufferIndex, itGlobal-1, cropSize, gXBitmap.get(), gYBitmap.get());
						ref<Bitmap> directBitmap(new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getSize()));
						computeVeryDirect(m_bufferIndex, itGlobal-1, cropSize, directBitmap);					
						auto directVector = bitmap2vec(directBitmap);
						if (!m_config.directTracing) {
							auto subPixelCount = size_t(3 * directBitmap->getSize().x * directBitmap->getSize().y);
							directVector = std::vector<float>(subPixelCount, 0.f);
						}
						
						// Develop all films
						develop(scene, film, throughputBitmap, it);     // Save throughput + direct image
						develop(scene, hdrFilm, directBitmap, currentIteration, "_direct_");
						ref<Bitmap> gXBitmapAbs = gXBitmap->clone();
						ref<Bitmap> gYBitmapAbs = gYBitmap->clone();					
						computeAbs(cropSize, gXBitmapAbs);
						computeAbs(cropSize, gYBitmapAbs);					
						develop(scene, film, gXBitmapAbs, it, "_dxAbs_");
						develop(scene, film, gYBitmapAbs, it, "_dyAbs_");					

						if (m_config.nanCheck) {
							/*
							for (size_t i = 0; i < subPixelCount; ++i) {
								if (!std::isfinite(dxVector[i])) {
									int y = (i / 3) / film->getCropSize().x;
									int x = (i / 3) % film->getCropSize().x;
									int c = i % 3;
									SLog(EWarn, "Gradient x at pixel (%d, %d, %d) is NaN", x, y, c);
									dxVector[i] = 0.f;
								}

								if (!std::isfinite(dyVector[i])) {
									int y = (i / 3) / film->getCropSize().x;
									int x = (i / 3) % film->getCropSize().x;
									int c = i % 3;
									SLog(EWarn, "Gradient y at pixel (%d, %d, %d) is NaN", x, y, c);
									dyVector[i] = 0.f;
								}

								if (!std::isfinite(throughputVector[i])) {
									int y = (i / 3) / film->getCropSize().x;
									int x = (i / 3) % film->getCropSize().x;
									int c = i % 3;
									SLog(EWarn, "Throughput at pixel (%d, %d, %d) is NaN", x, y, c);
									throughputVector[i] = 0.f;
								}
							}*/
						}

						/* Reconstruct. */
						{
							auto throughputVector = bitmap2vec(throughputBitmap);
							auto dxVector = bitmap2vec(gXBitmap);
							auto dyVector = bitmap2vec(gYBitmap);
							auto rec_results = rec.reconstruct(film->getCropSize(),
															throughputVector, 
															dxVector, dyVector,
															directVector,
															Reconstruction::Variance {},
															PostProcessOption{
																	forceBlackPixels: m_config.forceBlackPixels,
																	clampingValues: true
															});
							
							Bitmap *bitmapRecons = NULL;
							if (rec_results.size() == 1) {
								develop(scene, film, rec_results[0].img, it, "_recons_");
								bitmapRecons = rec_results[0].img;
							} else {
								for (auto &result: rec_results) {
									develop(scene, film, result.img, it, "_" + result.name + "_");
								}
								bitmapRecons = rec_results[0].img;
							}

							// Show the results if it is needed
							film->setBitmap(bitmapRecons);
							queue->signalRefresh(job);
						}
					}

					if ((m_config.reconstructNfor || m_config.reconstructWeighted) && (itGlobal - 1) % 2 == 0) 
					{	
						Reconstruction rec {};
						rec.reconstructL1 = m_config.reconstructL1;
						rec.reconstructL2 = m_config.reconstructL2;
						rec.reconstructUni = m_config.reconstructUni;
						rec.reconstructWeighted = m_config.reconstructWeighted;
						rec.alpha = (float)m_config.reconstructAlpha;

						// Estimate variance maps using AB buffer
						Vector2i size = film->getSize();
						ref<Bitmap> throughputBitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, size);					
						ref<Bitmap> throughputBitmapVariance = throughputBitmap->clone();	

						ref<Bitmap> throughputBitmapA = throughputBitmap->clone();
						ref<Bitmap> throughputBitmapB = throughputBitmap->clone();
						computeThroughput(0, itGlobal / 2, size, throughputBitmapA);
						computeThroughput(1, itGlobal / 2, size, throughputBitmapB);
						averageAndVariance(throughputBitmapA, throughputBitmapB, throughputBitmap, throughputBitmapVariance);

						ref<Bitmap> dxBitmap = throughputBitmap->clone();
						ref<Bitmap> dyBitmap = throughputBitmap->clone();
						ref<Bitmap> dxBitmapA = throughputBitmap->clone();
						ref<Bitmap> dyBitmapA = throughputBitmap->clone();
						ref<Bitmap> dxBitmapB = throughputBitmap->clone();
						ref<Bitmap> dyBitmapB = throughputBitmap->clone();
						ref<Bitmap> dxBitmapVariance = throughputBitmap->clone();
						ref<Bitmap> dyBitmapVariance = throughputBitmap->clone();
						computeGradient(0, itGlobal / 2, size, dxBitmapA, dyBitmapA);
						computeGradient(1, itGlobal / 2, size, dxBitmapB, dyBitmapB);
						averageAndVariance(dxBitmapA, dxBitmapB, dxBitmap, dxBitmapVariance);
						averageAndVariance(dyBitmapA, dyBitmapB, dyBitmap, dyBitmapVariance);

						ref<Bitmap> directBitmap = throughputBitmap->clone();
						computeVeryDirect(m_bufferIndex, itGlobal / 2, size, directBitmap);		// No need average since very direct is less noisy			
						auto directVector = bitmap2vec(directBitmap);
						if (!m_config.directTracing) {
							auto subPixelCount = size_t(3 * directBitmap->getSize().x * directBitmap->getSize().y);
							directVector = std::vector<float>(subPixelCount, 0.f);
						}

						if (m_config.reconstructNfor && !m_config.reconstructL1 && !m_config.reconstructL2 && !m_config.reconstructUni) {
							// When no specific technique for reconstruction is requested, use weighted control variates.
							rec.reconstructWeighted = true;
						}

						// Feature buffer
						auto throughputVector = bitmap2vec(throughputBitmap);
						auto gradientXVector = bitmap2vec(dxBitmap);
						auto gradientYVector = bitmap2vec(dyBitmap);
						auto results = rec.reconstruct(film->getSize(),
													throughputVector,
													gradientXVector,
													gradientYVector,
													directVector,
													Reconstruction::Variance{// Variance
															primal: bitmap2vec(throughputBitmapVariance),
															dx: bitmap2vec(dxBitmapVariance),
															dy: bitmap2vec(dyBitmapVariance)
													},
													PostProcessOption{
															forceBlackPixels: m_config.forceBlackPixels,
															clampingValues: true
													});
						
						develop(scene, hdrFilm, throughputBitmap, currentIteration, "_throughput_");
						develop(scene, hdrFilm, dxBitmap, currentIteration, "_dx_");
						develop(scene, hdrFilm, dyBitmap, currentIteration, "_dy_");
						if (results.size() == 1) {
							// Only save weighted because the others are saved before
							if (m_config.reconstructWeighted)
								develop(scene, hdrFilm, results[0].img, currentIteration, "_" + results[0].name + "_");
						} else {
							for (int i = 0; i < results.size(); ++i) {							
								if (results[i].name == "Weighted")
									develop(scene, hdrFilm, results[i].img, currentIteration, "_" + results[i].name + "_");
							}
						}

						if (m_config.reconstructNfor) 
						{					
							auto throughputVectorA = bitmap2vec(throughputBitmapA);
							auto gradientXVectorA = bitmap2vec(dxBitmapA);
							auto gradientYVectorA = bitmap2vec(dyBitmapA);
							auto resultsA = rec.reconstruct(film->getSize(),
															throughputVectorA,
															gradientXVectorA,
															gradientYVectorA,
															directVector,
															Reconstruction::Variance{// reuse variance of the main buffer
																	primal: bitmap2vec(throughputBitmapVariance),
																	dx: bitmap2vec(dxBitmapVariance),
																	dy: bitmap2vec(dyBitmapVariance)
															},
															PostProcessOption{
																	forceBlackPixels: m_config.forceBlackPixels,
																	clampingValues: true
															});

							auto throughputVectorB = bitmap2vec(throughputBitmapB);
							auto gradientXVectorB = bitmap2vec(dxBitmapB);
							auto gradientYVectorB = bitmap2vec(dyBitmapB);
							auto resultsB = rec.reconstruct(film->getSize(),
															throughputVectorB,
															gradientXVectorB,
															gradientYVectorB,
															directVector,
															Reconstruction::Variance{// reuse variance of the main buffer
																	primal: bitmap2vec(throughputBitmapVariance),
																	dx: bitmap2vec(dxBitmapVariance),
																	dy: bitmap2vec(dyBitmapVariance)
															},
															PostProcessOption{
																	forceBlackPixels: m_config.forceBlackPixels,
																	clampingValues: true
															});

							// Do NFOR for all reconstructed images
							for (int i = 0; i < results.size(); ++i) {
								ref<Bitmap> reconstructionBitmap = results[i].img;
								ref<Bitmap> reconstructionBitmapA = resultsA[i].img;
								ref<Bitmap> reconstructionBitmapB = resultsB[i].img;
								ref<Bitmap> reconstructionVarianceBitmap = throughputBitmap->clone();
								// Approximate variance of reconstruction by the reconstruction of bufferA and bufferB
								{
									Spectrum *recon = (Spectrum *) reconstructionBitmap->getFloatData();
									Spectrum *bufferA = (Spectrum *) reconstructionBitmapA->getFloatData();
									Spectrum *bufferB = (Spectrum *) reconstructionBitmapB->getFloatData();
									Spectrum *variance = (Spectrum *) reconstructionVarianceBitmap->getFloatData();
									Vector2i size = throughputBitmap->getSize();
									for (int i = 0; i < size.x * size.y; ++i) {
										// This recon is usually more accurate than averaging two recons from bufferA and B.
										// So we will get better variance approximation.
										Spectrum mean = recon[i];
										variance[i] = 0.5 * ((bufferA[i] - mean) * (bufferA[i] - mean) +
															(bufferB[i] - mean) * (bufferB[i] - mean));
									}
								}

								ref<AccumBuffer> throughputBuffer = new AccumBuffer(
										{throughputBitmap},
										{throughputBitmapVariance},
										{throughputBitmapA},
										{throughputBitmapB},
										currentIteration
								);

								ref<AccumBuffer> reconstructionBuffer = new AccumBuffer(
										{reconstructionBitmap},
										{reconstructionVarianceBitmap},
										{reconstructionBitmapA},
										{reconstructionBitmapB},
										currentIteration
								);
								
								// Custom NFOR
								{
									ref<Nfor> nfor = new Nfor();
									ref<Bitmap> nforBitmap = reconstructionBitmap->clone();
									nfor->denoise(throughputBuffer, nforBitmap, reconstructionBuffer, {0}, true);							
									if (m_config.directTracing) {
										nforBitmap->accumulate(directBitmap);
									}

									//if (m_config.reconstructL1 || m_config.reconstructL2 || m_config.reconstructUni || m_config.reconstructWeighted)
										develop(scene, hdrFilm, nforBitmap, currentIteration, "_" + results[i].name + "_nfor_custom_");
									//else {
									//	develop(scene, hdrFilm, nforBitmap, currentIteration, "_" + results[i].name + "_");
									//}							
								}
							}
						}
					}
				} // end of primal/gradient mode

				// Print statistics
				Statistics::getInstance()->printStats();

				renderingTimer->reset();
			}
			++it;
		}

		// Before stopping, display the results
		//film->setBitmap(m_bitmapRecons);
		queue->signalRefresh(job);

#ifdef MTS_DEBUG_FP
		disableFPExceptions();
#endif

		timeFile.close();

		// Free memory
		for (GPMThreadData *thData : m_threadData) {
			delete thData;
		}

		return true;
	}

	void photonMapPass(Scene *scene, int it,
					   RenderQueue *queue, const RenderJob *job,
					   Film *film, int sceneResID, int sensorResID,
					   std::vector<ShootingThreadData> &shootingData) {
		ref<Timer> photonShootingTime = new Timer;

		Log(EInfo, "Performing a photon mapping pass %i ("
				SIZE_T_FMT
				" photons so far for buffer %i)",
			it, m_totalPhotons[m_bufferIndex], m_bufferIndex);
		ref<Scheduler> sched = Scheduler::getInstance();

		/* Re-get the image size useful for testing if
		 * a gather point is at the image border */
		Vector2i cropSize = film->getCropSize();

		size_t granularity = 0;
		if (m_config.deterministic) {
			granularity = m_config.volumePhotonCount / sched->getWorkerCount();
			// To make sure that the we will have one task per proc
			m_config.volumePhotonCount = granularity * sched->getWorkerCount();
			SLog(EInfo,
				 "Deterministic mode. Granularity %d. Volume photon count adjusted to %d. Worker count %d",
				 granularity,
				 m_config.volumePhotonCount,
				 sched->getWorkerCount());
		}

		/* Generate the global photon map */
		const Float beamInitSize = m_smokeAABB.getBSphere().radius * globalScaleVolume[m_bufferIndex] * POURCENTAGE_BS;
		ref<GradientPhotonProcess> proc = new GradientPhotonProcess(m_config, sensorPosition(scene),
																	granularity,
																	true,
																	true, // Auto cancel the job is needed
																	job, shootingData,
																	beamInitSize);

		proc->bindResource("scene", sceneResID);
		proc->bindResource("sensor", sensorResID);

		sched->schedule(proc);
		sched->wait(proc);
		// Display the shooting time without considering building the acceleration structure
		Log(EInfo, "Shooting time: %i ms", photonShootingTime->getMilliseconds());
		Log(EDebug, "Photon map full. Shot "
				SIZE_T_FMT
				" paths, excess paths due to parallelism: "
				SIZE_T_FMT, std::max(proc->getShotParticlesSurface(), proc->getShotParticlesVolume()),
			proc->getExcessPhotons());
		Log(EInfo, "Number of volume photon/beam skip: %i", proc->nbSkip());
		m_totalEmittedSurface[m_bufferIndex] += proc->getShotParticlesSurface();
		m_totalEmittedVolume[m_bufferIndex] += proc->getShotParticlesVolume();

		// Gathering on surface
		ref<Timer> gatheringTime = new Timer;
		Log(EInfo, "Gathering ..");
		film->clear();
		if (m_config.needSurfaceRendering()) {
			ref<GPhotonMap> gradientPhotonMap = proc->getPhotonMap();
			gradientPhotonMap->build();
			computeSurfaceGradient(it, scene, gradientPhotonMap);
			m_totalPhotons[m_bufferIndex] += gradientPhotonMap->size();
		}
		Log(EInfo, "Surface gathering time: %i ms", gatheringTime->getMilliseconds());
		gatheringTime->reset();
		Log(EInfo, "Surface gathering done.");

		if (m_config.needVolumeRendering()) {
			ref<GPhotonMap> gradientPhotonMap = proc->getPhotonVolumeMap();
			ref<LTBeamMap> beamMap = proc->getBeamMap();
			if (!EVolumeTechniqueHelper::needBeams(m_config.volTechnique)) {
				gradientPhotonMap->build();
			}
			switch (m_config.volTechnique) {
				case EDistance:
					computeVolumeGradientPhoton(scene, gradientPhotonMap.get());
					break;
				case EVolBRE2D:
				case EVolBRE3D:
					computeVolumeGradientPhotonBRE(it,
												   scene,
												   gradientPhotonMap.get(),
												   proc->getShotParticlesVolume());
					break;
				case EBeamBeam1D:
				case EBeamBeam3D_Naive:
				case EBeamBeam3D_EGSR:
				case EBeamBeam3D_Optimized:
					computeVolumeGradientBeams(it, scene, beamMap, proc->getShotParticlesVolume());
					break;
				case EVolPlane0D:
					computeVolumeGradientPlanes(it, scene, beamMap, proc->getShotParticlesVolume());
					break;
				default:
					SLog(EError, "Unknow volume rendering technique");
			}
			Log(EInfo, "Volume gathering time: %i ms", gatheringTime->getMilliseconds());
			Log(EInfo, "Volume gathering done.");
		}

#if HAVE_ADDITIONAL_STATS
		// Acceptance
		ref<Bitmap> accBit = m_bitmap->clone();
		generateStatsBitmapAcc(cropSize, accBit);
		develop(scene, film, accBit, it, "_acc_");

		// The typeof shift
		ref<Bitmap> typeShift = m_bitmap->clone();
		generateStatsBitmapType(cropSize, typeShift);
		develop(scene, film, typeShift, it, "_type_");
#endif


		// Clean memory
		for (ShootingThreadData &d : shootingData) {
			d.clear();
		}
		for (GPMThreadData *thData : m_threadData) {
			if (!thData->pool.unused()) {
				SLog(EError, thData->pool.toString().c_str());
				SLog(EError, "Memory leak detected (th data)");
			}
		}
	}

	void computeSurfaceGradient(int it, Scene *scene, GPhotonMap *surfaceMap) {
		ref<Scheduler> sched = Scheduler::getInstance();
		size_t nCores = sched->getCoreCount();
		BlockScheduler blockSched(m_gatherBlocks->size(), nCores);

		blockSched.run([this, scene, surfaceMap](int blockIdx, int tid) {
			GPMThreadData *threadData = m_threadData[tid];
			std::vector<GatherPoint> &gatherPoints = (*m_gatherBlocks)[blockIdx];
			// Calculate flux arriving at the gather points
			for (size_t i = 0; i < gatherPoints.size(); ++i) {
				GatherPoint &gp = gatherPoints[i];
				if (m_config.directTracing) {
					if (!gp.currEmission.isZero()) {
						std::vector<ShiftGatherPoint> shiftGPs(4);
						if (m_config.estimateGradients == false) shiftGPs.clear();

						SurfaceGradientRecord gRec(scene, &gp, m_config, *threadData, shiftGPs);
						gRec.estimateGradEmitter(); // directly add the radiance inside the gather point
						// Free memory
						for (auto shiftGP: shiftGPs) {
							shiftGP.path.release(threadData->pool);
						}
					}
				}

				if (gp.depth != -1 && m_config.needSurfaceRendering()) {
					std::vector<ShiftGatherPoint> shiftGPs(4);
					if (m_config.estimateGradients == false) shiftGPs.clear();

					// Compute the surface radiance and taking into account gradients
					SurfaceGradientRecord gRec(scene, &gp, m_config, *threadData, shiftGPs);
					auto M = (Float) surfaceMap->evaluate(gRec, gp.lastIts().p, gp.radius);

					// Add the collected flux to surface flux data
					if (gp.N + M != 0) {
						Float ratio = (gp.N + m_config.alpha * M) / (gp.N + M);
						gp.scale = gp.scale * std::sqrt(ratio);

						// Update flux
						gp.baseFlux += gRec.baseFlux;
						gp.baseFlux /= gp.radius * gp.radius * M_PI; // Normalize it

						// Update flux used for gradient calculation
						for (int k = 0; k < shiftGPs.size(); ++k) {
							gp.shiftedFlux[k] += gRec.shiftedFlux[k];
							gp.shiftedFlux[k] /= gp.radius * gp.radius * M_PI;

							gp.weightedBaseFlux[k] += gRec.weightedBaseFlux[k];
							gp.weightedBaseFlux[k] /= gp.radius * gp.radius * M_PI;
						}
						gp.N = gp.N + m_config.alpha * M;
					}

					// Free memory
					for (auto shiftGP: shiftGPs) {
						shiftGP.path.release(threadData->pool);
					}
				} // End case for the surfaces
			}
		});
	}

	void computeVolumeGradientPlanes(int it, Scene *scene, LTBeamMap *beamMap, size_t nbPathBeams) {
		{
			auto sensor = scene->getSensor();
			if (sensor->getMedium() == nullptr) {
				SLog(EError, "Planes does not support camera outside the medium");
			}
		}

		auto beams = beamMap->getBeams();

		// Convert beams to planes
		Sampler *beamSampler = m_gpManager->getSamplerBlock(0);
		std::vector<LTPhotonPlane> planes(beams.size());
		for (size_t i = 0; i < beams.size(); i++) {
			planes[i] = LTPhotonPlane::transformBeam(beams[i].second, beamSampler);
		}

		// Construct the photon plane BVH
		ref<PhotonPlaneBVH<LTPhotonPlane>> planeBVH = new PhotonPlaneBVH<LTPhotonPlane>(planes);

		ref<Scheduler> sched = Scheduler::getInstance();
		size_t nCores = sched->getCoreCount();
		BlockScheduler blockSched(m_gatherBlocks->size(), nCores);
		blockSched.run([this, it, scene, beamMap, nbPathBeams, planeBVH](int blockIdx, int tid) {

			GPMThreadData *threadData = m_threadData[tid];
			std::vector<GatherPoint> &gatherPoints = (*m_gatherBlocks)[blockIdx];

			for (auto &gp : gatherPoints) {
				// The volume radiance gathered during this iteration
				Spectrum fluxVolIter(0.f);
				Spectrum shiftedMediumFluxIter[4];
				Spectrum weightedMediumFluxIter[4];
				for (int k = 0; k < 4; k++) {
					shiftedMediumFluxIter[k] = weightedMediumFluxIter[k] = Spectrum(0.f);
				}

				// Shift paths ...
				std::vector<ShiftGatherPoint> shiftGPs(4);

				// Count the number of edges that interact with the participating media
				for (int idEdge = 1; idEdge < int(gp.path.edgeCount()); idEdge++) {
					if (m_config.maxCameraDepth != -1 && idEdge > m_config.maxCameraDepth + 1) {
						break; // Skip this path
					}
					if (gp.path.edge(idEdge)->medium != nullptr) {
						// DEBUGGING
						gp.haveSmoke = true;
						if (int(gp.path.vertexCount()) < idEdge + 2)
						SLog(EError, "Problem of path");
						const PathEdge *currE = gp.path.edge(idEdge);
						Vector d = -currE->d;     // Camera ray towards light
						Float distTotal = currE->length;

						// Query the BRE by casting a ray inside it
						Ray ray(gp.path.vertex(idEdge)->getPosition(), d, Epsilon, distTotal - Epsilon, 0.f);
						PlaneGradRadianceQuery radQuery(scene, &gp, shiftGPs, ray, currE->medium,
														m_config, *threadData,
														idEdge);
						planeBVH->query(radQuery);
						fluxVolIter += radQuery.mediumFlux;

						for (int k = 0; k < 4; ++k) {
							shiftedMediumFluxIter[k] += radQuery.shiftedMediumFlux[k];
							weightedMediumFluxIter[k] += radQuery.weightedMediumFlux[k];
						}
					}
				}

				// Normalization of this pass
				fluxVolIter /= nbPathBeams;
				for (int k = 0; k < 4; k++) {
					shiftedMediumFluxIter[k] /= nbPathBeams;
					weightedMediumFluxIter[k] /= nbPathBeams;
				}

				// In APA, we always average the rendering
				// even if there is no photon collected
				gp.mediumFlux = (gp.mediumFlux * (it - 1) + fluxVolIter) / it; // APA estimator
				for (int k = 0; k < 4; k++) {
					gp.shiftedMediumFlux[k] = (gp.shiftedMediumFlux[k] * (it - 1)
											   + shiftedMediumFluxIter[k]) / it;
					gp.weightedMediumFlux[k] = (gp.weightedMediumFlux[k] * (it - 1)
												+ weightedMediumFluxIter[k]) / it;
				}

				// Free memory
				for (auto shiftGP: shiftGPs) {
					shiftGP.path.release(threadData->pool);
				}
			}
		});

		scaleVolumeAPA(it);

		SLog(EInfo, "Gathering done");
	}

	void computeVolumeGradientBeams(int it, Scene *scene, LTBeamMap *beamMap, size_t nbPathBeams) {
		const Float beamInitSize = m_smokeAABB.getBSphere().radius * globalScaleVolume[m_bufferIndex] * POURCENTAGE_BS;
		Log(EInfo, "Initial size of the beams: %f", beamInitSize);

		if (m_config.convertLong) {
			Log(EInfo, "Convert to long beams ...");
			beamMap->convertToLongBeams(scene);
		}

		const Vector2i &cropSize = scene->getFilm()->getCropSize();

		SLog(EInfo, "Build the photon beam map");
		beamMap->build(LTBeamMap::EBVHAccel, m_config.deterministic);

		SLog(EInfo, "Gathering...");

		ref<Scheduler> sched = Scheduler::getInstance();
		size_t nCores = sched->getCoreCount();
		BlockScheduler blockSched(m_gatherBlocks->size(), nCores);
		blockSched.run([&](int blockIdx, int tid) {
			GPMThreadData *threadData = m_threadData[tid];
			Sampler *sampler = m_gpManager->getSamplerBlock(blockIdx);

			std::vector<GatherPoint> &gatherPoints = (*m_gatherBlocks)[blockIdx];

			for (auto &gp : gatherPoints) {
				size_t idPix = gp.pixel.y * cropSize.y + gp.pixel.x;

				// The volume radiance gathered during this iteration
				Spectrum fluxVolIter(0.f);
				Spectrum fluxExtraVolIter(0.f);
				Spectrum shiftedMediumFluxIter[4];
				Spectrum weightedMediumFluxIter[4];
				for (int k = 0; k < 4; k++) {
					shiftedMediumFluxIter[k] = weightedMediumFluxIter[k] = Spectrum(0.f);
				}

				// Shift paths ...
				std::vector<ShiftGatherPoint> shiftGPs(4);

				// Count the number of edges that interact with the participating media
				for (int idEdge = 1; idEdge < int(gp.path.edgeCount()); idEdge++) {
					if (m_config.minDepth > idEdge) {
						continue;
					}
					if (m_config.maxCameraDepth != -1 && idEdge > m_config.maxCameraDepth + 1) {
						break; // Skip this path
					}

					if (gp.path.edge(idEdge)->medium != nullptr) {
						const PathEdge *currE = gp.path.edge(idEdge);
						Vector d = -currE->d;     // Camera ray towards light
						Float distTotal = currE->length;
						gp.haveSmoke = true;

						// Query the BRE by casting a ray inside it
						Ray ray(gp.path.vertex(idEdge)->getPosition(), d, Epsilon, distTotal - Epsilon, 0.f);
						BeamGradRadianceQuery radQuery(scene, &gp, shiftGPs, ray, currE->medium,
													   m_config, *threadData,
													   idEdge, sampler);
						beamMap->query(radQuery);
						fluxVolIter += radQuery.mediumFlux;
						fluxExtraVolIter += radQuery.extraFlux;

						for (int k = 0; k < 4; ++k) {
							shiftedMediumFluxIter[k] += radQuery.shiftedMediumFlux[k];
							weightedMediumFluxIter[k] += radQuery.weightedMediumFlux[k];
						}

#if HAVE_ADDITIONAL_STATS
						for (int k = 0; k < 4; ++k) {
						  gp.stats.add(radQuery.shiftStats[k]);
						}
#endif
					}
				}


				// Normalization of this pass
				fluxVolIter /= nbPathBeams;
				fluxExtraVolIter /= nbPathBeams;
				for (int k = 0; k < 4; k++) {
					shiftedMediumFluxIter[k] /= nbPathBeams;
					weightedMediumFluxIter[k] /= nbPathBeams;
				}

				// In APA, we always average the rendering
				// even if there is no photon collected
				gp.mediumFlux = (gp.mediumFlux * (it - 1) + fluxVolIter) / it; // APA estimator
				for (int k = 0; k < 4; k++) {
					gp.shiftedMediumFlux[k] = (gp.shiftedMediumFlux[k] * (it - 1)
											   + shiftedMediumFluxIter[k]) / it;
					gp.weightedMediumFlux[k] = (gp.weightedMediumFlux[k] * (it - 1)
												+ weightedMediumFluxIter[k]) / it;
				}

				// Free memory
				for (auto shiftGP: shiftGPs) {
					shiftGP.path.release(threadData->pool);
				}
			}
		});

		scaleVolumeAPA(it);

		SLog(EInfo, "Gathering done");
	}

	void computeVolumeGradientPhotonBRE(int it, Scene *scene, GPhotonMap *gradientPhotonMap, size_t nbPathVolume) {
		const Float breInitSize = m_smokeAABB.getBSphere().radius * globalScaleVolume[m_bufferIndex] * POURCENTAGE_BS;
		Log(EInfo, "BRE photon size: %f", breInitSize);

		// Construct the photon size by looking to the camera distance.
		// For now we use APA estimator
		ref<GradientBeamRadianceEstimator> bre = new GradientBeamRadianceEstimator(gradientPhotonMap,
																				   breInitSize);

		ref<Scheduler> sched = Scheduler::getInstance();
		size_t nCores = sched->getCoreCount();
		BlockScheduler blockSched(m_gatherBlocks->size(), nCores);
		blockSched.run([this, it, scene, gradientPhotonMap, nbPathVolume, bre](int blockIdx, int tid) {

			GPMThreadData *threadData = m_threadData[tid];
			Sampler *sampler = m_gpManager->getSamplerBlock(blockIdx);

			// The base image
			std::vector<GatherPoint> &gatherPoints = (*m_gatherBlocks)[blockIdx];
			// Calculate flux arriving at the gather points
			for (auto &gp : gatherPoints) {
				std::vector<ShiftGatherPoint> shiftGPs(4);
				if (m_config.estimateGradients == false) shiftGPs.clear();

				Spectrum fluxVolIter(0.f);
				Spectrum shiftedMediumFluxIter[4];
				Spectrum weightedMediumFluxIter[4];
				for (int k = 0; k < shiftGPs.size(); k++) {
					shiftedMediumFluxIter[k] = weightedMediumFluxIter[k] = Spectrum(0.f);
				}

				// Count the number of edges that interact with the participating media
				for (int idEdge = 1; idEdge < int(gp.path.edgeCount()); idEdge++) {
					if (m_config.minCameraDepth > idEdge) {
						continue; // Skip this path
					}
					if (m_config.maxCameraDepth != -1 && idEdge > m_config.maxCameraDepth + 1) {
						break; // Skip this path
					}

					if (gp.path.edge(idEdge)->medium != nullptr) {
						// If we are here, its mean that we have a volume segment
						const Medium *currentMedium = gp.path.edge(idEdge)->medium;
						gp.haveSmoke = true;

						// Recreate the beam inside the medium
						Point oBeam = gp.path.vertex(idEdge)->getPosition();
						// TODO: Warning: This is not correct we handle sky beams
						Vector dBeam = gp.path.vertex(idEdge + 1)->getPosition() - oBeam;
						Float beamDist = dBeam.length();
						dBeam /= beamDist;

						Ray ray(oBeam, dBeam, Epsilon, beamDist - Epsilon, 0.f);
						VolumeGradientBREQuery gRec(scene, &gp, m_config, *threadData, shiftGPs, idEdge, sampler);
						gRec.newRayBase(ray, gp.path.edge(idEdge)->medium);
						gRec.clear();
						bre->query(ray, currentMedium, gRec, sampler->next1D());

						// Accumulate the values
						// Note that beam weight is already include
						fluxVolIter += gRec.mediumFlux;
						for (int k = 0; k < shiftGPs.size(); k++) {
							shiftedMediumFluxIter[k] += gRec.shiftedMediumFlux[k];
							weightedMediumFluxIter[k] += gRec.weightedMediumFlux[k];
						}
					}
				}

				// Normalization of this pass
				fluxVolIter /= nbPathVolume;
				for (int k = 0; k < shiftGPs.size(); k++) {
					shiftedMediumFluxIter[k] /= nbPathVolume;
					weightedMediumFluxIter[k] /= nbPathVolume;
				}

				// In APA, we always average the rendering
				// even if there is no photon collected
				gp.mediumFlux = (gp.mediumFlux * (it - 1) + fluxVolIter) / it; // APA estimator
				for (int k = 0; k < shiftGPs.size(); k++) {
					gp.shiftedMediumFlux[k] = (gp.shiftedMediumFlux[k] * (it - 1)
											   + shiftedMediumFluxIter[k]) / it;
					gp.weightedMediumFlux[k] = (gp.weightedMediumFlux[k] * (it - 1)
												+ weightedMediumFluxIter[k]) / it;
				}

				// Free memory
				for (auto &shiftGP : shiftGPs) {
					shiftGP.path.release(threadData->pool);
				}
			}
		});

		scaleVolumeAPA(it);
	}

	void computeVolumeGradientPhoton(Scene *scene, GPhotonMap *gradientPhotonMap) {
		const Float BBPourcentageCONST = m_smokeAABB.getBSphere().radius * POURCENTAGE_BS;

		ref<Scheduler> sched = Scheduler::getInstance();
		size_t nCores = sched->getCoreCount();
		BlockScheduler blockSched(m_gatherBlocks->size(), nCores);
		blockSched.run([this, scene, gradientPhotonMap, BBPourcentageCONST](int blockIdx, int tid) {

			GPMThreadData *threadData = m_threadData[tid];
			Sampler *sampler = m_gpManager->getSamplerBlock(blockIdx);

			// The base image
			std::vector<GatherPoint> &gatherPoints = (*m_gatherBlocks)[blockIdx];
			// Calculate flux arriving at the gather points
			for (auto &gp : gatherPoints) {
				std::vector<ShiftGatherPoint> shiftGPs(4);
				Float MVol = 0;

				// Count how many participating media edges
				int nbEdgePM = 0;
				for (size_t idEdge = 1; idEdge < gp.path.edgeCount(); idEdge++) {
					if (gp.path.edge(idEdge)->medium != nullptr) {
						nbEdgePM += 1;
					}
				}

				// Skip this gp, no participating media.
				if (nbEdgePM == 0)
					continue;

				// Compute CDF and all the distance
				Float totalBeamDist = 0.0;
				DiscreteDistribution selBeam(nbEdgePM);
				std::vector<size_t> beamIndex;
				beamIndex.reserve(nbEdgePM);

				Spectrum weightBeam(1.f);
				for (size_t idEdge = 1; idEdge < gp.path.edgeCount(); idEdge++) {
					if (gp.path.edge(idEdge)->medium != nullptr) {
						selBeam.append(weightBeam.max());
						totalBeamDist += gp.path.edge(idEdge)->length;
						beamIndex.push_back(idEdge);
					}
					// Account to the changes
					// FIXME: RR????
					weightBeam *= gp.path.vertex(idEdge + 1)->weight[EImportance];
					weightBeam *= gp.path.edge(idEdge)->weight[EImportance];
				}
				selBeam.normalize();

				// Make the same size as the ray marching (twice bigger by default)
				const Float querySize = BBPourcentageCONST * gp.scaleVol;
				Float normalization = 1.f / m_config.nbCameraSamples;

				// Create a global volume query
				VolumeGradientDistanceQuery gRec(scene, &gp, m_config, *threadData, shiftGPs, 0, sampler);

				for (int idSamples = 0; idSamples < m_config.nbCameraSamples; idSamples++) {
					Float randSample = sampler->next1D();
					if (m_config.stratified) {
						randSample = idSamples * normalization + randSample * normalization;
					}
					size_t sampleIndex = selBeam.sampleReuse(randSample);
					size_t idEdge = beamIndex[sampleIndex];

					// Skipping the max/min camera depth
					if (m_config.maxCameraDepth != -1 && idEdge > m_config.maxCameraDepth + 1) {
						continue;
					}
					if (m_config.minCameraDepth != 0 && idEdge < m_config.minCameraDepth + 1) {
						continue;
					}
					gp.haveSmoke = true;
					gRec.changeEdge(idEdge, selBeam[sampleIndex]);

					// Recreate the beam inside the medium
					Point oBeam = gp.path.vertex(idEdge)->getPosition();
					Vector dBeam = gp.path.vertex(idEdge + 1)->getPosition() - oBeam; // FIXME: Possible to use edge
					Float beamDist = dBeam.length(); // FIXME: Possible to use edge
					dBeam /= beamDist; // FIXME: Possible to use edge

					// If we are here, its mean that we have a volume segment
					const Medium *currentMedium = gp.path.edge(idEdge)->medium;

					// Sample a distance inside the participating media
					Ray ray(oBeam, dBeam, Epsilon, beamDist, 0.f);
					MediumSamplingRecord mRec;
					if (currentMedium->sampleDistance(ray, mRec, sampler, EDistanceAlwaysValid, randSample)) {
						// Update the ray max t because it is used in grad. volume computation
						ray.maxt = mRec.t;

						gRec.newRayBase(ray, mRec, querySize, mRec.pdfSuccess);
						gRec.clear();
						MVol += gradientPhotonMap->evaluate(gRec, ray.o + mRec.t * ray.d, querySize);

						gp.mediumFlux += (gRec.mediumFlux * normalization);

						for (int k = 0; k < 4; ++k) {
							gp.shiftedMediumFlux[k] += (gRec.shiftedMediumFlux[k] * normalization);
							gp.weightedMediumFlux[k] += (gRec.weightedMediumFlux[k] * normalization);
						}
					}
				} // All the samples

#if HAVE_ADDITIONAL_STATS
				for (int k = 0; k < 4; ++k) {
				  gp.stats.add(gRec.shiftStats[k]);
				}
#endif

				if (MVol + gp.NVol != 0) {
					Float ratioVol = (gp.NVol + m_config.alpha * MVol) / (gp.NVol + MVol);
					gp.scaleVol = gp.scaleVol * std::cbrt(ratioVol);
					gp.NVol = gp.NVol + m_config.alpha * MVol;
				}

				// Free memory
				for (auto shiftGP: shiftGPs) {
					shiftGP.path.release(threadData->pool);
				}
			}
		});
	}

	void computeThroughput(int bufferIndex, int N, const Vector2i &cropSize, Bitmap *bitmap) {
		GatherBlocks &blocks = m_gatherBlocksAB[bufferIndex];		
		int totalEmittedSurface = m_totalEmittedSurface[bufferIndex];
		int totalEmittedVolume = m_totalEmittedVolume[bufferIndex];
		
		if (! m_config.reusePrimal) 
		{
			for (auto &gatherPoints : blocks) {
				auto *target = (Spectrum *) bitmap->getUInt8Data();

				// Calculate flux arriving at the gather points
				for (size_t i = 0; i < gatherPoints.size(); ++i) {
					GatherPoint &gp = gatherPoints[i];
					// Fill throughput
					Spectrum fluxMedia(0.f);
					if (totalEmittedVolume != 0) {
						if (m_config.isAPAVolumeEstimator())
							fluxMedia = gp.mediumFlux;
						else
							fluxMedia = gp.mediumFlux / ((Float) totalEmittedVolume);
					}
					Spectrum fluxSurface(0.f);
					if (totalEmittedSurface != 0) {
						fluxSurface = gp.baseFlux / ((Float) totalEmittedSurface);
					}
					target[gp.pixel.y * bitmap->getWidth() + gp.pixel.x] = fluxMedia + fluxSurface + (gp.emission / N);
				}
			}
		} 
		else
		{
			if (m_config.estimateGradients == false) {
				SLog(EError, "reusePrimal not possible when no gradient is estimated.");
			}
			std::vector<std::vector<GatherPoint *>> &imgGP = m_imgGPAB[bufferIndex];
			
			auto *targetPrimal = (Spectrum *) bitmap->getUInt8Data();
			for (int y = 0; y < cropSize.y; ++y) {
				for (int x = 0; x < cropSize.x; ++x) {
					Spectrum TVol(0.f);
					Spectrum TSurf(0.f);
					Spectrum TEmitter(0.f);
					if (x != cropSize.x - 1) {
						TVol += imgGP[x + 1][y]->shiftedMediumFlux[ELeft];
						TSurf += imgGP[x + 1][y]->shiftedFlux[ELeft];
						if(m_config.directTracing) {
							TEmitter += imgGP[x + 1][y]->shiftedEmitterFlux[ELeft];
						}
					}
					if (x != 0) {
						TVol += imgGP[x - 1][y]->shiftedMediumFlux[ERight];
						TSurf += imgGP[x - 1][y]->shiftedFlux[ERight];
						if(m_config.directTracing) {
							TEmitter += imgGP[x - 1][y]->shiftedEmitterFlux[ERight];
						}
					}
					if (y != cropSize.y - 1) {
						TVol += imgGP[x][y + 1]->shiftedMediumFlux[EBottom];
						TSurf += imgGP[x][y + 1]->shiftedFlux[EBottom];
						if(m_config.directTracing) {
							TEmitter += imgGP[x][y + 1]->shiftedEmitterFlux[EBottom];
						}
					}
					if (y != 0) {
						TVol += imgGP[x][y - 1]->shiftedMediumFlux[ETop];
						TSurf += imgGP[x][y - 1]->shiftedFlux[ETop] ;
						if(m_config.directTracing) {
							TEmitter += imgGP[x][y - 1]->shiftedEmitterFlux[ETop];
						}
					}
					TVol += imgGP[x][y]->weightedMediumFlux[EBottom] + imgGP[x][y]->weightedMediumFlux[ETop] +
						 imgGP[x][y]->weightedMediumFlux[ERight] + imgGP[x][y]->weightedMediumFlux[ELeft];
					TSurf += imgGP[x][y]->weightedBaseFlux[EBottom] + imgGP[x][y]->weightedBaseFlux[ETop] +
							 imgGP[x][y]->weightedBaseFlux[ERight] + imgGP[x][y]->weightedBaseFlux[ELeft];
					if(m_config.directTracing) {
						TEmitter += imgGP[x][y]->weightedEmitterFlux[EBottom] + imgGP[x][y]->weightedEmitterFlux[ETop] +
									imgGP[x][y]->weightedEmitterFlux[ERight] + imgGP[x][y]->weightedEmitterFlux[ELeft];
					}
					targetPrimal[y * bitmap->getWidth() + x] = Spectrum(0.f);
					if(m_config.directTracing) {
						targetPrimal[y * bitmap->getWidth() + x] += TEmitter / ( 4.0 * N);
					}
					if(m_config.needSurfaceRendering()) {
						targetPrimal[y * bitmap->getWidth() + x] += (TSurf / 4.0) / ((Float) totalEmittedSurface);
					}
					if(m_config.needVolumeRendering()) {
						if (m_config.isAPAVolumeEstimator()) {
							targetPrimal[y * bitmap->getWidth() + x] += (TVol / 4.0);
						} else {
							targetPrimal[y * bitmap->getWidth() + x] += (TVol / 4.0) / ((Float) totalEmittedVolume);
						}
					}
				} // end pixel loop
			}
		}
	}

	void computeVeryDirect(int bufferIndex, int N, const Vector2i &cropSize, Bitmap *bitmap) {
		GatherBlocks &blocks = m_gatherBlocksAB[bufferIndex];		
		int totalEmittedSurface = m_totalEmittedSurface[bufferIndex];
		int totalEmittedVolume = m_totalEmittedVolume[bufferIndex];

		for (auto &gatherPoints : blocks)
			for (const auto &gp : gatherPoints) {
				if (!gp.directEmission.isZero()) {
					bitmap->setPixel(gp.pixel, gp.directEmission / N);
				} else {
					bitmap->setPixel(gp.pixel, Spectrum(0.f));
				}
			}
	}

	void computeAuxBitmaps(int bufferIndex, Bitmap *albedoBitmap, Bitmap *positionBitmap, Bitmap *shadingNormalBitmap) {
		GatherBlocks &blocks = m_gatherBlocksAB[bufferIndex];		
		
		for (auto &gatherPoints : blocks)
			for (const auto &gp : gatherPoints) {
				albedoBitmap->setPixel(gp.pixel, gp.albedo);
				Float pos[] = {gp.position.x, gp.position.y, gp.position.z};
				positionBitmap->setPixel(gp.pixel, Spectrum(pos));
				Float normal[]= {gp.shadingNormal.x, gp.shadingNormal.y, gp.shadingNormal.z};
				shadingNormalBitmap->setPixel(gp.pixel, Spectrum(normal));
			}
	}

	void computeGradient(int bufferIndex, int N, const Vector2i &cropSize, Bitmap *gXBitmap, Bitmap *gYBitmap) {
		std::vector<std::vector<GatherPoint *>> &imgGP = m_imgGPAB[bufferIndex];
		int totalEmittedSurface = m_totalEmittedSurface[bufferIndex];
		int totalEmittedVolume = m_totalEmittedVolume[bufferIndex];
		
		auto *targetGX = (Spectrum *) gXBitmap->getUInt8Data();
		auto *targetGY = (Spectrum *) gYBitmap->getUInt8Data();

		const bool DEBUG_GRAD = false;

		for (int y = 0; y < cropSize.y; ++y) {
			for (int x = 0; x < cropSize.x; ++x) {
				GatherPoint *curr = imgGP[x][y];

				targetGX[y * gXBitmap->getWidth() + x] = [&]() -> Spectrum {
					Spectrum gXsurface(0.0f), gXvolume(0.0f);
					if (totalEmittedSurface != 0) {
						// If the surface rendering is enable
						if (x == cropSize.x - 1) {
							gXsurface = (curr->shiftedFlux[ERight] - curr->weightedBaseFlux[ERight]);
						} else {
							GatherPoint *right = imgGP[x + 1][y];
							gXsurface = (curr->shiftedFlux[ERight] - curr->weightedBaseFlux[ERight]) +
										(right->weightedBaseFlux[ELeft] - right->shiftedFlux[ELeft]);
						}
						gXsurface /= ((Float) totalEmittedSurface);
					}

					if (m_config.directTracing) {
						// Gradient caused by non-density estimation path when gather point hits the light
						if (x == cropSize.x - 1) {
							gXsurface += (curr->shiftedEmitterFlux[ERight] - curr->weightedEmitterFlux[ERight]) / N;
						} else {
							GatherPoint *right = imgGP[x + 1][y];
							gXsurface += ((curr->shiftedEmitterFlux[ERight] - curr->weightedEmitterFlux[ERight]) +
										  (right->weightedEmitterFlux[ELeft] - right->shiftedEmitterFlux[ELeft])) / N;
						}
					}

					if (totalEmittedVolume != 0) {
						// If the volume rendering is enable
						if (x == cropSize.x - 1) {
							gXvolume = (curr->shiftedMediumFlux[ERight] - curr->weightedMediumFlux[ERight]);
							if (DEBUG_GRAD)
								gXvolume = 2.0 * (curr->shiftedMediumFlux[ERight]);
						} else {
							GatherPoint *right = imgGP[x + 1][y];
							gXvolume = (curr->shiftedMediumFlux[ERight] - curr->weightedMediumFlux[ERight]) +
									   (right->weightedMediumFlux[ELeft] - right->shiftedMediumFlux[ELeft]);
							if (DEBUG_GRAD)
								gXvolume = (curr->shiftedMediumFlux[ERight]) + right->shiftedMediumFlux[ELeft];
						}
						if (!m_config.isAPAVolumeEstimator())
							gXvolume /= ((Float) totalEmittedVolume);
					}
					return gXsurface + gXvolume;
				}();

				targetGY[y * gYBitmap->getWidth() + x] = [&]() -> Spectrum {
					Spectrum gYsurface(0.0f), gYvolume(0.0f);
					if (totalEmittedSurface != 0) {
						if (y == cropSize.y - 1) {
							gYsurface = (curr->shiftedFlux[ETop] - curr->weightedBaseFlux[ETop]);
						} else {
							GatherPoint *top = imgGP[x][y + 1];
							gYsurface = (curr->shiftedFlux[ETop] - curr->weightedBaseFlux[ETop]) +
										(top->weightedBaseFlux[EBottom] - top->shiftedFlux[EBottom]);
						}
						gYsurface /= ((Float) totalEmittedSurface);
					}

					if (m_config.directTracing) {
						if (y == cropSize.y - 1) {
							gYsurface += (curr->shiftedEmitterFlux[ETop] - curr->weightedEmitterFlux[ETop]) / N;
						} else {
							GatherPoint *top = imgGP[x][y + 1];
							gYsurface += ((curr->shiftedEmitterFlux[ETop] - curr->weightedEmitterFlux[ETop]) +
										  (top->weightedEmitterFlux[EBottom] - top->shiftedEmitterFlux[EBottom])) / N;
						}
					}

					if (totalEmittedVolume != 0) {
						if (y == cropSize.y - 1) {
							gYvolume = (curr->shiftedMediumFlux[ETop] - curr->weightedMediumFlux[ETop]);
						} else {
							GatherPoint *top = imgGP[x][y + 1];
							gYvolume = (curr->shiftedMediumFlux[ETop] - curr->weightedMediumFlux[ETop]) +
									   (top->weightedMediumFlux[EBottom] - top->shiftedMediumFlux[EBottom]);
						}
						if (!m_config.isAPAVolumeEstimator())
							gYvolume /= ((Float) totalEmittedVolume);
					}
					return gYsurface + gYvolume;
				}();
			}
		}
	}

	void computeAbs(const Vector2i &cropSize, Bitmap *bitmap) {
		auto *target = (Spectrum *) bitmap->getUInt8Data();
		for (int y = 0; y < cropSize.y; ++y) {
			for (int x = 0; x < cropSize.x; ++x) {
				target[y * bitmap->getWidth() + x] = target[y * bitmap->getWidth() + x].abs();
			}
		}
	}

	void computeShiftImage(int bufferIndex, int N,
						   const Vector2i &cropSize,
						   Bitmap *imgXPositive,
						   Bitmap *imgXNegative,
						   Bitmap *imgYPositive,
						   Bitmap *imgYNegative,
						   bool sameIntensityLevel = false) {
		
		std::vector<std::vector<GatherPoint *>> &imgGP = m_imgGPAB[bufferIndex];
		int totalEmittedSurface = m_totalEmittedSurface[bufferIndex];
		int totalEmittedVolume = m_totalEmittedVolume[bufferIndex];

		Spectrum *xPositive = (Spectrum *) imgXPositive->getUInt8Data();
		Spectrum *xNegative = (Spectrum *) imgXNegative->getUInt8Data();
		Spectrum *yPositive = (Spectrum *) imgYPositive->getUInt8Data();
		Spectrum *yNegative = (Spectrum *) imgYNegative->getUInt8Data();

		imgXPositive->clear();
		imgXNegative->clear();
		imgYPositive->clear();
		imgYNegative->clear();

		Float scale = 1.0;
		if (sameIntensityLevel) {
			if (m_config.useMIS == false)
				scale = 2.0;
			else {
				// Not yet supported. Have to add both shift directions together
			}
		}

		for (int y = 0; y < cropSize.y; ++y) {
			for (int x = 0; x < cropSize.x; ++x) {
				GatherPoint *curr = imgGP[x][y];

				Spectrum surX0(0.0f), surX1(0.0f), surY0(0.0f), surY1(0.0f);
				Spectrum volX0(0.0f), volX1(0.0f), volY0(0.0f), volY1(0.0f);
				Spectrum emitX0(0.0f), emitX1(0.0f), emitY0(0.0f), emitY1(0.0f);

				if (m_totalEmittedSurface != 0) {
					// If the surface rendering is enable
					surX0 = curr->shiftedFlux[ERight] * scale;
					surX1 = curr->shiftedFlux[ELeft] * scale;
					surY0 = curr->shiftedFlux[ETop] * scale;
					surY1 = curr->shiftedFlux[EBottom] * scale;

					surX0 /= ((Float) totalEmittedSurface);
					surX1 /= ((Float) totalEmittedSurface);
					surY0 /= ((Float) totalEmittedSurface);
					surY1 /= ((Float) totalEmittedSurface);
				}

				if (m_config.directTracing) {
					// Gradient caused by non-density estimation path when gather point hits the light
					emitX0 = curr->shiftedEmitterFlux[ERight] * scale;
					emitX1 = curr->shiftedEmitterFlux[ELeft] * scale;
					emitY0 = curr->shiftedEmitterFlux[ETop] * scale;
					emitY1 = curr->shiftedEmitterFlux[EBottom] * scale;

					emitX0 /= N;
					emitX1 /= N;
					emitY0 /= N;
					emitY1 /= N;
				}

				if (m_totalEmittedVolume != 0) {
					// If the volume rendering is enable
					volX0 = curr->shiftedMediumFlux[ERight] * scale;
					volX1 = curr->shiftedMediumFlux[ELeft] * scale;
					volY0 = curr->shiftedMediumFlux[ETop] * scale;
					volY1 = curr->shiftedMediumFlux[EBottom] * scale;

					if (!m_config.isAPAVolumeEstimator()) {
						volX0 /= ((Float) totalEmittedVolume);
						volX1 /= ((Float) totalEmittedVolume);
						volY0 /= ((Float) totalEmittedVolume);
						volY1 /= ((Float) totalEmittedVolume);
					}
				}

				Spectrum X0 = surX0 + emitX0 + volX0;
				Spectrum X1 = surX1 + emitX1 + volX1;
				Spectrum Y0 = surY0 + emitY0 + volY0;
				Spectrum Y1 = surY1 + emitY1 + volY1;

				xPositive[y * imgXPositive->getWidth() + x] = X0;
				xNegative[y * imgXNegative->getWidth() + x] = X1;
				yPositive[y * imgYPositive->getWidth() + x] = Y0;
				yNegative[y * imgYNegative->getWidth() + x] = Y1;
			}
		}
	}

#if HAVE_ADDITIONAL_STATS
	void generateStatsBitmapAcc(const Vector2i& cropSize, Bitmap* bitmap) {
	  Spectrum *target = (Spectrum *) bitmap->getUInt8Data();
	  for(int y = 0; y < cropSize.y; ++y) {
		for (int x = 0; x < cropSize.x; ++x) {
		  GatherPoint *curr = m_imgGP[x][y];

		  // Put the pourcentage of success
		  target[y * bitmap->getWidth() + x] = Spectrum(curr->stats.nbSuccessLightShifts /
			  ((Float)curr->stats.nbLightShifts));
		}
	  }
	}

	void generateStatsBitmapType(const Vector2i& cropSize, Bitmap* bitmap) {
	  Spectrum *target = (Spectrum *) bitmap->getUInt8Data();
	  for(int y = 0; y < cropSize.y; ++y) {
		for (int x = 0; x < cropSize.x; ++x) {
		  GatherPoint *curr = m_imgGP[x][y];

		  Float MEShift = curr->stats.MEShifts / ((Float)curr->stats.nbLightShifts);
		  Float DiffShift = curr->stats.DiffuseShifts / ((Float)curr->stats.nbLightShifts);
		  Float NullShift = curr->stats.nullShifts / ((Float)curr->stats.nbLightShifts);

		  Float rgb[3] = {MEShift, DiffShift, NullShift};

		  // Put the pourcentage of success
		  target[y * bitmap->getWidth() + x] = Spectrum(rgb);
		}
	  }
	}
#endif

	void clampNegative(const Vector2i &cropSize, Bitmap *bitmap) {
		auto *target = (Spectrum *) bitmap->getUInt8Data();
		for (int y = 0; y < cropSize.y; ++y) {
			for (int x = 0; x < cropSize.x; ++x) {
				target[y * bitmap->getWidth() + x].clampNegative();
			}
		}
	}

	const Point sensorPosition(const Scene *scene) {
		return scene->getSensor()->getWorldTransform()->operator()(0.f, Point(0.f, 0.f, 0.f));
	}

	void develop(Scene *scene, Film *film, const Bitmap *bitmap,
				 int currentIteration, const std::string &suffixName = "_") {
		std::stringstream ss;
		ss << scene->getDestinationFile().string() << suffixName
		   << currentIteration;
		std::string path = ss.str();
		film->setBitmap(bitmap);
		film->setDestinationFile(path, 0);
		film->develop(scene, 0.f);
		// in the last rendering pass, post process is called on scene which can call film develop once again
		// to prevent the curent file being overwritten with other bitmap (being displayed), we set the
		// file name to empty here
		film->setDestinationFile("", 0);
	}

	std::string toString() const override {
		std::ostringstream oss;
		oss << "SPPMIntegrator[" << endl
			<< "  maxDepth = " << m_config.maxDepth << "," << endl
			<< "  rrDepth = " << m_config.rrDepth << "," << endl
			<< "  alpha = " << m_config.alpha << "," << endl
			<< "  photonCount = " << m_config.photonCount << "," << endl
			<< "  maxPasses = " << m_config.maxPasses << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:	
	std::vector<GPMThreadData *> m_threadData;
		
	GatherBlocks m_gatherBlocksAB[2];	// two array of gather blocks for variance estimation
	GVPMRadiusInitializer m_gpManagerAB[2];
	GatherBlocks *m_gatherBlocks;		// point to one of the two arrays of the gather blocks
	GVPMRadiusInitializer *m_gpManager;
	std::vector<std::vector<GatherPoint *>> m_imgGPAB[2];
	std::vector<std::vector<GatherPoint *>> *m_imgGP;
	size_t m_totalEmittedSurface[2], m_totalEmittedVolume[2], m_totalPhotons[2];
	Float globalScaleVolume[2];
	int m_bufferIndex;
	bool bufferAB;
	
	std::vector<Point2i> m_offset;
	ref<Mutex> m_mutex;		
	ref<Bitmap> m_extraRendering;
	
	bool m_running;

	GPMConfig m_config;

	bool m_independentScale;
	AABB m_smokeAABB;
};

MTS_IMPLEMENT_CLASS_S(GPMIntegrator, false, Integrator)

MTS_EXPORT_PLUGIN(GPMIntegrator, "Gradient domain PM");
MTS_NAMESPACE_END
