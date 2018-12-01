/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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

#include <mitsuba/bidir/util.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/render/renderproc.h>
#include "mitsuba/core/plugin.h"

#include "gpt_proc.h"
#include "../reconstruction.h"
#include "shift_mapping/original.h"
#include "shift_mapping/explicit.h"
#include "shift_mapping/original_volume.h"
#include "../../denoiser/nfor/nfor.h"
#include "shift_mapping/random_replay.h"
#include "shift_mapping/random_replay_reconnect.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{gpt}{Gradient-domain path tracer}
* \order{5}
* \parameters{
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
* This plugin implements a gradient-domain path tracer (short: G-PT) as described in the paper "Gradient-Domain Path Tracing" by Kettunen et al. 
* It samples difference images in addition to the standard color image, and reconstructs the final image based on these.
* It supports classical materials like diffuse, specular and glossy materials, and area and point lights, depth-of-field, and low discrepancy samplers. 
* There is also experimental support for sub-surface scattering and motion blur. Note that this is still an experimental implementation of Gradient-Domain Path Tracing 
* that has not been tested with all of Mitsuba's features. Notably there is no support yet for any kind of participating media or directional lights. 
* Environment maps are supported, though. Does not support the 'hide emitters' option even though it is displayed.
*
*/

// NOTE: turn on CENTRAL_RADIANCE flag in shiftmapping.h. Do not turn it on here. 

StatsCounter shiftSuccessHVSurface("GPT",
                                   "Percentage of succesfull HVCopy (Surface) : ", EPercentage);
StatsCounter shiftSuccessHVVolume("GPT",
                                  "Percentage of succesfull HVCopy (Volume) : ", EPercentage);
StatsCounter shiftSuccessReconnectionVolume("GPT",
                                            "Percentage of succesfull reconnection (Volume) : ", EPercentage);
StatsCounter shiftSuccessReconnectionSurface("GPT",
                                             "Percentage of succesfull reconnection (Surface) : ", EPercentage);

StatsCounter shiftFailedVolume("GPT",
                               "Percentage of failure due to volume inconsistency: ", EPercentage);
StatsCounter shiftDied("GPT", "Percentage of shift died during the shift operation", EPercentage);

StatsCounter ReconnectionFailedVolume("GPT", "ReconnectionFailedVolume", ENumberValue);
StatsCounter ReconnectionFailedSurface("GPT", "ReconnectionFailedSurface", ENumberValue);
StatsCounter HVFailedVolume("GPT", "HVFailedVolume", ENumberValue);
StatsCounter HVFailedSurface("GPT", "HVFailedSurface", ENumberValue);
StatsCounter ShiftTransZero("GPT", "ShiftTransZero", ENumberValue);
StatsCounter TooShort("GPT", "TooShort", ENumberValue);
StatsCounter InconsistentBSDFShift("GPT", "InconsistentBSDFShift", ENumberValue);
StatsCounter HVSurfaceDisagreement("GPT", "HVSurfaceDisagreement", ENumberValue);
StatsCounter HVStepsSurface("GPT", "Number of HV over surfaces", ENumberValue);
StatsCounter ShiftZeroDensity("GPT", "Shift path kill 0 density", ENumberValue);

static StatsCounter avgPathLength("Gradient Path Tracer", "Average path length", EAverage);

// Output buffer names.
static const size_t BUFFER_FINAL = 0;       ///< Buffer index for the final image. Also used for preview.
static const size_t BUFFER_THROUGHPUT = 1;  ///< Buffer index for the noisy color image.
static const size_t BUFFER_DX = 2;          ///< Buffer index for the X gradients.
static const size_t BUFFER_DY = 3;          ///< Buffer index for the Y gradients.
static const size_t BUFFER_VERY_DIRECT = 4; ///< Buffer index for very direct light.

struct TheBuffers {
    ref<AccumBuffer> throughputBuffer;
    ref<AccumBuffer> throughputTmpBuffer;

    ref<AccumBuffer> gradientXBuffer;
    ref<AccumBuffer> gradientYBuffer;
    ref<AccumBuffer> gradientXTmpBuffer;
    ref<AccumBuffer> gradientYTmpBuffer;
};

enum EThroughputBuffer {
    Central, Left, Right, Top, Bottom, Count
};

/// The actual Gradient Path Tracer implementation.
class GradientPathTracer {
public:
    GradientPathTracer(const Scene *scene,
                       const Sensor *sensor,
                       Sampler *sampler,
                       TheBuffers *buffers,
                       const GradientPathTracerConfig *config)
            : m_scene(scene),
              m_sensor(sensor),
              m_sampler(sampler),
              m_config(config),
              m_buffers(buffers),
              m_shiftmapping(nullptr) {
        switch (m_config->shiftmapping) {
            case EShiftOriginal: {
                m_shiftmapping.reset(new ShiftMappingOriginal(m_config));
                break;
            }
            case EShiftOriginalVolume: {
                m_shiftmapping.reset(new ShiftMappingVolumeOriginal(m_config));
                break;
            }
            case EShiftExplicit: {
                m_shiftmapping.reset(new ShiftMappingExplicit(m_config));
                break;
            }
            case EShiftRandom: {
                m_shiftmapping.reset(new ShiftMappingRandom(m_config));
                break;
            }
            case EShiftRandomReconnect: {
                m_shiftmapping.reset(new ShiftMappingRandomReconnection(m_config));
                break;
            }
            default: {
                SLog(EError, "Not covered shift mapping");
                break;
            }
        }
    }

    /// Evaluates a sample at the given position.
    ///
    /// Outputs direct radiance to be added on top of the final image, the throughput to the central pixel, gradients to all neighbors,
    /// and throughput contribution to the neighboring pixels.
    void evaluatePoint(RadianceQueryRecord &rRec,
                       const Point2 &samplePosition,
                       const Point2 &apertureSample,
                       Float timeSample,
                       Float differentialScaleFactor,
                       Spectrum &out_very_direct,
                       Spectrum &out_throughput,
                       Spectrum *out_gradients,
                       Spectrum *out_neighborThroughputs) {
        // Initialize the base path.
        RayState mainRay;
        mainRay.throughput = m_sensor->sampleRayDifferential(mainRay.ray, samplePosition, apertureSample, timeSample);
        mainRay.ray.scaleDifferential(differentialScaleFactor);
        mainRay.ray.hasDifferentials = false;
        mainRay.rRec = rRec;
        mainRay.rRec.its = rRec.its;
        mainRay.rRec.medium = m_sensor->getMedium();

        // Initialize the offset paths.
        RayState shiftedRays[4];
        static const Vector2 pixelShifts[4] = {
                Vector2(1.0f, 0.0f),
                Vector2(0.0f, 1.0f),
                Vector2(-1.0f, 0.0f),
                Vector2(0.0f, -1.0f)
        };

        for (int i = 0; i < 4; ++i) {
            shiftedRays[i].throughput = m_sensor->sampleRayDifferential(shiftedRays[i].ray,
                                                                        samplePosition + pixelShifts[i],
                                                                        apertureSample,
                                                                        timeSample);
            shiftedRays[i].ray.scaleDifferential(differentialScaleFactor);
            shiftedRays[i].ray.hasDifferentials = false;
            shiftedRays[i].rRec = rRec;
            shiftedRays[i].rRec.its = rRec.its;
            shiftedRays[i].rRec.medium = m_sensor->getMedium();
        }

        // Evaluate the gradients. The actual algorithm happens here.
        Spectrum very_direct = Spectrum(0.0f);
        m_shiftmapping->evaluate(mainRay, shiftedRays, 4, very_direct);

        // Output results.
        out_very_direct = very_direct;
        out_throughput = mainRay.radiance;

        for (int i = 0; i < 4; i++) {
            out_gradients[i] = shiftedRays[i].gradient;
            out_neighborThroughputs[i] = shiftedRays[i].radiance;
        }
    }

    void renderBlockPathReuse(GPTWorkResult *block,
                              const bool &stop,
                              const Vector2i &rectSize) {
        AccumBuffer *throughputBuffer = m_buffers->throughputBuffer.get();
        AccumBuffer *gradientXBuffer = m_buffers->gradientXBuffer.get();
        AccumBuffer *gradientYBuffer = m_buffers->gradientYBuffer.get();

        bool needsApertureSample = m_sensor->needsApertureSample();
        bool needsTimeSample = m_sensor->needsTimeSample();

        // Original code from SamplingIntegrator.
        Float diffScaleFactor = 1.0f / std::sqrt((Float) m_sampler->getSampleCount());

        // Get ready for sampling.
        RadianceQueryRecord rRec(m_scene, m_sampler);

        Point2 apertureSample(0.5f);
        Float timeSample = 0.5f;
        RayDifferential sensorRay;
        block->clear();

        struct BlockInfo {
            const int number_pixels;
            const std::vector<Vector2> pixelShifts;
            std::vector<ShiftMapping::GradientInfo> gradients;
            std::vector<Spectrum> shiftedThroughputs;
            std::vector<Spectrum> centralVeryDirect;

            void clear() {
                for (auto &g: gradients) {
                    g.value = Spectrum(0.f);
                }
                for (auto &s: shiftedThroughputs) {
                    s = Spectrum(0.f);
                }
                for (auto &s: centralVeryDirect) {
                    s = Spectrum(0.f);
                }
            }
        };

        const bool avoid_overlap = true;
        const int tile_size = 4; // should be power of 2 here
        auto block_info = [&]() -> BlockInfo {
            const int number_pixels = tile_size * tile_size + 2 * tile_size;
            std::vector<ShiftMapping::GradientInfo> gradients;
            std::vector<Vector2> pixelShifts;
            pixelShifts.reserve(number_pixels);
            for (int iy = 0, ii = 0; iy < tile_size + 1; iy++) {
                for (int ix = 0; ix < tile_size + 1; ix++, ii++) {
                    // Exclude the last one
                    if (ix == tile_size && iy == tile_size)
                        continue;
                    pixelShifts.emplace_back(Vector2(ix, iy));
                    if (avoid_overlap && (ix == tile_size || iy == tile_size)) {
                        continue; // Do not add gradient between
                    }

                    if (ix != tile_size && (ii + 1) + 1 != (tile_size + 1) * (tile_size + 1)) {
                        gradients.emplace_back(
                                ShiftMapping::GradientInfo{
                                        Spectrum(0.f),
                                        true,
                                        Vector2(ix, iy),
                                        ii + 1, // Next id
                                        ii,
                                });
                    }
                    if (iy != tile_size && (ii + 1) + (tile_size + 1) != (tile_size + 1) * (tile_size + 1)) {
                        gradients.emplace_back(
                                ShiftMapping::GradientInfo{
                                        Spectrum(0.f),
                                        false,
                                        Vector2(ix, iy),
                                        ii + (tile_size + 1),
                                        ii,
                                });
                    }
                }
            }
            return BlockInfo{
                    number_pixels,
                    pixelShifts,
                    gradients,
                    std::vector<Spectrum>(number_pixels, Spectrum(0.f)),
                    std::vector<Spectrum>(number_pixels, Spectrum(0.f)),
            };
        }();

        // Sample at the given positions.
        for (int iy = 0; iy < rectSize.y; iy += tile_size) {
            for (int ix = 0; ix < rectSize.x; ix += tile_size) {
                if (stop) {
                    break;
                }

                Point2i offset = Point2i(ix, iy) + Vector2i(block->getOffset());
                m_sampler->generate(offset);
                for (size_t j = 0; j < m_sampler->getSampleCount(); ++j) {
                    if (stop) {
                        break;
                    }

                    block_info.clear();
                    for (int i = 0; i < block_info.number_pixels; ++i) {
                        // Get the initial ray to sample.
                        rRec.newQuery(RadianceQueryRecord::ESensorRay, m_sensor->getMedium());
                        Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));
                        if (needsApertureSample) {
                            apertureSample = rRec.nextSample2D();
                        }
                        if (needsTimeSample) {
                            timeSample = rRec.nextSample1D();
                        }

                        std::vector<RayState> shiftedRays(block_info.number_pixels, RayState{});
                        for (int k = 0; k < block_info.number_pixels; k++) {
                            shiftedRays[k].throughput = m_sensor->sampleRayDifferential(shiftedRays[k].ray,
                                                                                        samplePos +
                                                                                        block_info.pixelShifts[k],
                                                                                        apertureSample,
                                                                                        timeSample);
                            shiftedRays[k].ray.scaleDifferential(diffScaleFactor);
                            shiftedRays[k].ray.hasDifferentials = false;
                            shiftedRays[k].rRec = rRec;
                            shiftedRays[k].rRec.its = rRec.its;
                        }
                        m_shiftmapping->evaluateReuse(&shiftedRays[0], block_info.number_pixels,
                                                      block_info.centralVeryDirect[i], i,
                                                      block_info.gradients);
                        // Accumulate
                        for (int k = 0; k < block_info.number_pixels; k++) {
                            block_info.shiftedThroughputs[k] += shiftedRays[k].radiance;
                        }
                    }

                    // FIXME: Box filter is Ok
                    Point2 samplePos(Point2(offset) + Vector2(0.5));
                    {
                        for (int i = 0; i < block_info.number_pixels; i++) {
                            if (avoid_overlap
                                && (block_info.pixelShifts[i].x == tile_size ||
                                    block_info.pixelShifts[i].y == tile_size)) {
                                continue;
                            }

                            block->put(samplePos + block_info.pixelShifts[i],
                                       block_info.centralVeryDirect[i] + block_info.shiftedThroughputs[i],
                                       1.0f,
                                       1.0f,
                                       BUFFER_FINAL);
                        }
                    }

                    // Actual throughputs, with MIS between central and neighbor pixels for all neighbors.
                    // This can be replaced with a standard throughput sample without much loss of quality in most cases.
                    {
                        for (int i = 0; i < block_info.number_pixels; i++) {
                            if (avoid_overlap
                                && (block_info.pixelShifts[i].x == tile_size ||
                                    block_info.pixelShifts[i].y == tile_size)) {
                                continue;
                            }
                            block->put(samplePos + block_info.pixelShifts[i],
                                       block_info.shiftedThroughputs[i],
                                       1.0f,
                                       1.0f,
                                       BUFFER_THROUGHPUT);
                            block->put(samplePos + block_info.pixelShifts[i],
                                       block_info.centralVeryDirect[i],
                                       1.0f, 1.0f, BUFFER_VERY_DIRECT);

                            throughputBuffer->addSample(samplePos.x + block_info.pixelShifts[i].x,
                                                        samplePos.y + block_info.pixelShifts[i].y,
                                                        block_info.shiftedThroughputs[i]);
                        }
                    }

                    // Gradients.
                    {
                        for (auto &g: block_info.gradients) {
                            block->put(samplePos + g.offset, g.value, 1.0f, 1.0f, g.x_axis ? BUFFER_DX : BUFFER_DY);
                            if (g.x_axis) {
                                gradientXBuffer->addSample((samplePos + g.offset).x, (samplePos + g.offset).y, g.value);
                            } else {
                                gradientYBuffer->addSample((samplePos + g.offset).x, (samplePos + g.offset).y, g.value);
                            }
                        }
                    }

                }
            }
        }
    }

    void renderBlock(GPTWorkResult *block,
                     const bool &stop,
                     const std::vector<TPoint2<uint8_t>> &points) {

        AccumBuffer *throughputBuffer = m_buffers->throughputBuffer.get();
        AccumBuffer *throughputTmpBuffer = m_buffers->throughputTmpBuffer.get();
        AccumBuffer *gradientXTmpBuffer = m_buffers->gradientXTmpBuffer.get();
        AccumBuffer *gradientYTmpBuffer = m_buffers->gradientYTmpBuffer.get();

        bool needsApertureSample = m_sensor->needsApertureSample();
        bool needsTimeSample = m_sensor->needsTimeSample();

        // Original code from SamplingIntegrator.
        Float diffScaleFactor = 1.0f / std::sqrt((Float) m_sampler->getSampleCount());

        // Get ready for sampling.
        RadianceQueryRecord rRec(m_scene, m_sampler);

        Point2 apertureSample(0.5f);
        Float timeSample = 0.5f;
        RayDifferential sensorRay;

        block->clear();

        // Sample at the given positions.
        Spectrum gradients[4];
        Spectrum shiftedThroughputs[4];

        for (size_t i = 0; i < points.size(); ++i) {
            if (stop) {
                break;
            }

            Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
            m_sampler->generate(offset);

            for (size_t j = 0; j < m_sampler->getSampleCount(); ++j) {
                if (stop) {
                    break;
                }

                // Get the initial ray to sample.
                rRec.newQuery(RadianceQueryRecord::ESensorRay, m_sensor->getMedium());
                Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));
                if (needsApertureSample) {
                    apertureSample = rRec.nextSample2D();
                }
                if (needsTimeSample) {
                    timeSample = rRec.nextSample1D();
                }

                // Do the actual sampling.
                Spectrum centralVeryDirect = Spectrum(0.0f);
                Spectrum centralThroughput = Spectrum(0.0f);

                evaluatePoint(rRec,
                              samplePos,
                              apertureSample,
                              timeSample,
                              diffScaleFactor,
                              centralVeryDirect,
                              centralThroughput,
                              gradients,
                              shiftedThroughputs
                );

                // Accumulate results.
                const Point2 right_pixel = samplePos + Vector2(1.0f, 0.0f);
                const Point2 bottom_pixel = samplePos + Vector2(0.0f, 1.0f);
                const Point2 left_pixel = samplePos - Vector2(1.0f, 0.0f);
                const Point2 top_pixel = samplePos - Vector2(0.0f, 1.0f);
                const Point2 center_pixel = samplePos;

                static const int RIGHT = 0;
                static const int BOTTOM = 1;
                static const int LEFT = 2;
                static const int TOP = 3;


                // Note: Sampling differences and throughputs to multiple directions is essentially
                //       multiple importance sampling (MIS) between the pixels.
                //
                //       For a sample from a strategy participating in the MIS to be unbiased, we need to
                //       divide its result by the selection probability of that strategy.
                //
                //       As the selection probability is 0.5 for both directions (no adaptive sampling),
                //       we need to multiply the results by two.

                // Note: The central pixel is estimated as
                //               1/4 * (throughput estimate sampled from MIS(center, top)
                //                      + throughput estimate sampled from MIS(center, right)
                //                      + throughput estimate sampled from MIS(center, bottom)
                //                      + throughput estimate sampled from MIS(center, left)).
                //
                //       Variable centralThroughput is the sum of four throughput estimates sampled
                //       from each of these distributions, from the central pixel, so it's actually four samples,
                //       and thus its weight is 4.
                //
                //       The other samples from the MIS'd distributions will be sampled from the neighboring pixels,
                //       and their weight is 1.
                //
                //       If this feels too complicated, it should be OK to output a standard throughput sample from
                //       the path tracer.

                // Add the throughput image as a preview. Note: Preview and final buffers are shared.
                if (!m_config->reusePrimal) {
                    block->put(samplePos,
                               centralVeryDirect + centralThroughput,
                               1.0f,
                               1.0f,
                               BUFFER_FINAL); // Standard throughput estimate with direct.
                } else {
                    block->put(samplePos, (8 * centralVeryDirect) + (2 * centralThroughput), 4.0f, 4.0f,
                               0); // Adds very direct on top of the throughput image.
                    block->put(left_pixel, (2 * shiftedThroughputs[LEFT]), 1.0f, 1.0f,
                               BUFFER_FINAL);     // Negative x throughput.
                    block->put(right_pixel, (2 * shiftedThroughputs[RIGHT]), 1.0f, 1.0f,
                               BUFFER_FINAL);   // Positive x throughput.
                    block->put(top_pixel, (2 * shiftedThroughputs[TOP]), 1.0f, 1.0f,
                               BUFFER_FINAL);       // Negative y throughput.
                    block->put(bottom_pixel, (2 * shiftedThroughputs[BOTTOM]), 1.0f, 1.0f,
                               BUFFER_FINAL); // Positive y throughput.
                }

                // Actual throughputs, with MIS between central and neighbor pixels for all neighbors.
                // This can be replaced with a standard throughput sample without much loss of quality in most cases.
                if (!m_config->reusePrimal) {
                    block->put(samplePos, centralThroughput, 1.0f, 1.0f,
                               BUFFER_THROUGHPUT); // Standard throughput estimate.

                    throughputBuffer->addSample(samplePos.x, samplePos.y,
                                                centralThroughput); // Also track separately in our accumulation buffer for variance estimation.
                } else {
                    block->put(samplePos, (2 * centralThroughput), 4.0f, 4.0f,
                               BUFFER_THROUGHPUT); // Central throughput.
                    block->put(left_pixel, (2 * shiftedThroughputs[LEFT]), 1.0f, 1.0f,
                               BUFFER_THROUGHPUT);     // Negative x throughput.
                    block->put(right_pixel, (2 * shiftedThroughputs[RIGHT]), 1.0f, 1.0f,
                               BUFFER_THROUGHPUT);   // Positive x throughput.
                    block->put(top_pixel, (2 * shiftedThroughputs[TOP]), 1.0f, 1.0f,
                               BUFFER_THROUGHPUT);       // Negative y throughput.
                    block->put(bottom_pixel, (2 * shiftedThroughputs[BOTTOM]), 1.0f, 1.0f,
                               BUFFER_THROUGHPUT); // Positive y throughput.

                    // We keep track of all throughputs and combine them when the entire image is available 
                    throughputTmpBuffer->addSample(samplePos.x, samplePos.y, centralThroughput);
                    throughputTmpBuffer->addSample(samplePos.x, samplePos.y, shiftedThroughputs[LEFT],
                                                   EThroughputBuffer::Left);
                    throughputTmpBuffer->addSample(samplePos.x, samplePos.y, shiftedThroughputs[RIGHT],
                                                   EThroughputBuffer::Right);
                    throughputTmpBuffer->addSample(samplePos.x, samplePos.y, shiftedThroughputs[TOP],
                                                   EThroughputBuffer::Top);
                    throughputTmpBuffer->addSample(samplePos.x, samplePos.y, shiftedThroughputs[BOTTOM],
                                                   EThroughputBuffer::Bottom);
                }

                // Gradients.
                {
                    block->put(left_pixel,
                               -(2 * gradients[LEFT]), 1.0f, 1.0f, BUFFER_DX);    // Negative x gradient.
                    block->put(center_pixel, (
                            2 * gradients[RIGHT]), 1.0f, 1.0f, BUFFER_DX);  // Positive x gradient.
                    block->put(top_pixel,
                               -(2 * gradients[TOP]), 1.0f, 1.0f, BUFFER_DY);      // Negative y gradient.
                    block->put(center_pixel, (
                            2 * gradients[BOTTOM]), 1.0f, 1.0f, BUFFER_DY); // Positive y gradient.

                    // Store gradient into accumulation buffer and track variance for all computations related to this pixel. 
                    // To avoid confusion, do not write to other pixels. 
                    gradientXTmpBuffer->addSample(samplePos.x, samplePos.y, gradients[RIGHT], 0);
                    gradientXTmpBuffer->addSample(samplePos.x, samplePos.y, gradients[LEFT], 1);
                    gradientYTmpBuffer->addSample(samplePos.x, samplePos.y, gradients[BOTTOM], 0);
                    gradientYTmpBuffer->addSample(samplePos.x, samplePos.y, gradients[TOP], 1);
                }

                // Very direct.
                block->put(center_pixel, centralVeryDirect,
                           1.0f, 1.0f, BUFFER_VERY_DIRECT);
            }
        }
    }


private:
    const Scene *m_scene;
    const Sensor *m_sensor;
    Sampler *m_sampler;
    std::unique_ptr<ShiftMapping> m_shiftmapping;
    TheBuffers *m_buffers;
    const GradientPathTracerConfig *m_config;
};

GradientPathIntegrator::GradientPathIntegrator(const Properties &props)
        : MonteCarloIntegrator(props) {

    m_config.m_shiftThreshold = props.getFloat("shiftThreshold", Float(0.001));
    m_config.m_reconstructL1 = props.getBoolean("reconstructL1", true);
    m_config.m_reconstructL2 = props.getBoolean("reconstructL2", false);
    m_config.m_reconstructUni = props.getBoolean("reconstructUni", false);
    m_config.m_reconstructWeighted = props.getBoolean("reconstructWeighted", false);
    m_config.m_reconstructL2Weighted = props.getBoolean("reconstructL2Weighted", false);
    m_config.m_reconstructNfor = props.getBoolean("reconstructNfor", false);
    m_config.m_reconstructAlpha = (Float) props.getFloat("reconstructAlpha", Float(0.2));
    m_config.forceBlackPixels = props.getBoolean("forceBlackPixels", false);
    m_config.reusePrimal = props.getBoolean("reusePrimal", true);
    m_config.directTracing = props.getBoolean("directTracing", true);
    // For now, we compute variance only in this case.
    m_config.computeVariance =
            m_config.m_reconstructWeighted || m_config.m_reconstructL2Weighted || m_config.m_reconstructNfor;
    // Force to false as we do not support this mode yet
    m_strictNormals = false;

    m_config.shiftmapping = [&]() -> ShiftStrategy {
        auto name = props.getString("shiftmapping", "original");
        if (name == "original")
            return EShiftOriginal;
        else if (name == "original_volume")
            return EShiftOriginalVolume;
        else if (name == "explicit")
            return EShiftExplicit;
        else if (name == "random")
            return EShiftRandom;
        else if (name == "random_reconnect")
            return EShiftRandomReconnect;
        else
            SLog(EError, "incompatible shift mapping: %s", name.c_str());
    }();
    m_config.pathReuse = props.getBoolean("pathReuse", false);
    if ((!m_config.reusePrimal) && m_config.pathReuse) {
        SLog(EError, "Impossible to not reusePrimal and do path reuse");
    }

    //if (m_config.m_reconstructAlpha <= 0.0f)
    //    Log(EError, "'reconstructAlpha' must be set to a value greater than zero!");

    // Debug proposes
    m_showAbsGradients = true;
    m_config.m_minDepth = props.getInteger("minDepth", 0);
    m_config.rrShift = props.getFloat("rrShift", 0.1);
    if (m_config.rrShift <= 0.0 || m_config.rrShift > 1.0) {
        SLog(EError, "Wrong value for rrShift: %f", m_config.rrShift);
    }

    /* This option is specific for GVPM project:
       * The only propose is to compute a very specific light transport when computing volume interaction
       * This option does not have (really) meaning for path tracing but it does for photon mapping-based
       * approaches. The idea is to restrict the light transport only if the path have bounce over
       * classified glossy surface.
       */
    m_config.m_minCameraDepth = props.getInteger("minCameraDepth", 0);

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
    m_config.m_lightingInteractionMode = parseMediaInteractionMode(lightingInteractionMode);

    // Set to positive to dump a hdr files
    m_dumpIteration = props.getInteger("dumpIteration", 1);

    if (m_config.m_minCameraDepth > 0 && needRenderSurface(m_config.m_lightingInteractionMode)) {
        SLog(EError, "minCameraDepth is only meaningful when only computing volume interactions.");
    }
}

GradientPathIntegrator::GradientPathIntegrator(Stream *stream, InstanceManager *manager)
        : MonteCarloIntegrator(stream, manager) {
    m_config.m_shiftThreshold = stream->readFloat();
    m_config.m_reconstructL1 = stream->readBool();
    m_config.m_reconstructL2 = stream->readBool();
    m_config.m_reconstructAlpha = stream->readFloat();
    m_dumpIteration = 0;
}

void GradientPathIntegrator::renderBlock(const Scene *scene,
                                         const Sensor *sensor,
                                         Sampler *sampler,
                                         GPTWorkResult *block,
                                         const bool &stop,
                                         const std::vector<TPoint2<uint8_t>> &points,
                                         const Vector2i &rectSize) const {
    TheBuffers buffers;
    buffers.throughputBuffer = throughputBuffer;
    buffers.throughputTmpBuffer = throughputTmpBuffer;
    buffers.gradientXBuffer = gradientXBuffer;
    buffers.gradientYBuffer = gradientYBuffer;
    buffers.gradientXTmpBuffer = gradientXTmpBuffer;
    buffers.gradientYTmpBuffer = gradientYTmpBuffer;
    GradientPathTracer tracer(scene, sensor, sampler, &buffers, &m_config);
    if (m_config.pathReuse) {
        tracer.renderBlockPathReuse(block, stop, rectSize);
    } else {
        tracer.renderBlock(block, stop, points);
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

/// Custom render function that samples a number of paths for evaluating differences between pixels.
bool GradientPathIntegrator::render(Scene *scene,
                                    RenderQueue *queue, const RenderJob *job,
                                    int sceneResID, int sensorResID, int samplerResID) {

    m_hideEmitters = false; // Avoid the option
    if (m_hideEmitters) {
        /* Not supported! */
        Log(EError, "Option 'hideEmitters' not implemented for Gradient-Domain Path Tracing!");
    }

    /* Get config from the parent class. */
    m_config.m_maxDepth = m_maxDepth;
    m_config.m_rrDepth = m_rrDepth;
    m_config.m_strictNormals = m_strictNormals;

    /* Create CSV file to dump all the rendering timings */
    // Also create an timer
    std::string timeFilename = scene->getDestinationFile().string()
                               + "_time.csv";
    std::ofstream timeFile(timeFilename.c_str());
    ref<Timer> renderingTimer = new Timer;

    /* Code duplicated from SamplingIntegrator::Render. */
    ref<Scheduler> sched = Scheduler::getInstance();
    ref<Sensor> sensor = dynamic_cast<Sensor *>(sched->getResource(sensorResID));

    /* Set up MultiFilm. */
    ref<Film> film = sensor->getFilm();
    Vector2i cropSize = film->getCropSize();

    auto outNames = [&]() -> std::vector<std::string> {
        return {"-final", "-throughput", "-dx", "-dy", "-direct"};
    }();
    if (!film->setBuffers(outNames)) {
        Log(EError, "Cannot render image! G-PT has been called without MultiFilm.");
        return false;
    }


    size_t nCores = sched->getCoreCount();
    const auto *sampler = dynamic_cast<const Sampler *>(sched->getResource(samplerResID, 0));
    size_t sampleCount = sampler->getSampleCount();

    Log(EInfo, "Starting render job (GPT Custom::render) (%ix%i, "
            SIZE_T_FMT
            " %s, "
            SIZE_T_FMT
            " %s, "
            SSE_STR
            ") ..", film->getCropSize().x, film->getCropSize().y,
        sampleCount, sampleCount == 1 ? "sample" : "samples", nCores,
        nCores == 1 ? "core" : "cores");

    // Accumulative buffers
    bool bufferAB = (m_config.m_reconstructNfor == true);
    bool trackVariance = m_config.computeVariance;

    throughputBuffer = new AccumBuffer(1, film->getSize(), trackVariance, bufferAB);
    throughputTmpBuffer = new AccumBuffer(EThroughputBuffer::Count, film->getSize(), trackVariance, bufferAB);

    // Due to the calculation of gradient involves accessing current pixel and its neighbors, 
    // it is generally not convenient to track variance and buffer A, B in the accumulation buffer
    // The workaround is to keep track of all shift data at each pixel in the temporary buffers 
    // and consolidate them after each tracing pass
    gradientXBuffer = new AccumBuffer(1, film->getSize(), trackVariance, bufferAB);
    gradientYBuffer = new AccumBuffer(1, film->getSize(), trackVariance, bufferAB);
    // Temp buffers to record all forward and inverse shift
    gradientXTmpBuffer = new AccumBuffer(2, film->getSize(), trackVariance, bufferAB);
    gradientYTmpBuffer = new AccumBuffer(2, film->getSize(), trackVariance, bufferAB);

    //ref<AuxiliaryBuffer> auxBuffer = new AuxiliaryBuffer(film->getSize(), trackVariance, bufferAB);
    //auxBuffer->clear();

    int currentIteration = 1;

    while (true) {
        // Throughput tmp buffer caches only the current iteration 
        throughputTmpBuffer->clear();

        /* This is a sampling-based integrator - parallelize. */
        ref<BlockedRenderProcess> proc = new GPTRenderProcess(job, queue, scene->getBlockSize(), m_config);

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

        sched->unregisterResource(integratorResID);
        m_process = nullptr;

        if (proc->getReturnStatus() != ParallelProcess::ESuccess) {
            Log(EWarn, "There is a problem during the rendering process. Abord");
            break;
        }

        /* Manage throughput buffer in reuse primal case */
        if (m_config.reusePrimal && (!m_config.pathReuse)) {
            // In this case, the variance is only available from the 2nd iteration. 

            // TODO: parallelize this code 
            
            ref<Bitmap> throughputCentralBitmap = throughputTmpBuffer->getBuffer(EThroughputBuffer::Central);
            ref<Bitmap> throughputLeftBitmap = throughputTmpBuffer->getBuffer(EThroughputBuffer::Left);
            ref<Bitmap> throughputRightBitmap = throughputTmpBuffer->getBuffer(EThroughputBuffer::Right);
            ref<Bitmap> throughputTopBitmap = throughputTmpBuffer->getBuffer(EThroughputBuffer::Top);
            ref<Bitmap> throughputBottomBitmap = throughputTmpBuffer->getBuffer(EThroughputBuffer::Bottom);

            Vector2i size = scene->getFilm()->getSize();
            int width = size.x;

            Spectrum *throughputCentral = (Spectrum *) throughputCentralBitmap->getFloatData();
            Spectrum *throughputLeft = (Spectrum *) throughputLeftBitmap->getFloatData();
            Spectrum *throughputRight = (Spectrum *) throughputRightBitmap->getFloatData();
            Spectrum *throughputTop = (Spectrum *) throughputTopBitmap->getFloatData();
            Spectrum *throughputBottom = (Spectrum *) throughputBottomBitmap->getFloatData();

            for (int y = 0; y < size.y; ++y) {
                for (int x = 0; x < size.x; ++x) {
                    int centerIndex = y * width + x;
                    int leftIndex = y * width + x - 1;
                    int topIndex = (y - 1) * width + x;
                    int rightIndex = y * width + x + 1;
                    int bottomIndex = (y + 1) * width + x;

                    // To estimate pixel value, we have to add central with each shifted value as 
                    // they are both weighted by MIS for forward and backward mapping. 
                    // Since there are 4 directions, it results in 1/4 weight. 
                    // And note that the central value itself is already 4x the intensity so we only add it once. 
                    Spectrum val = 0.25 * throughputCentral[centerIndex];
                    if (x + 1 < size.x)
                        val += 0.25 * throughputLeft[rightIndex];
                    else
                        // when out of bound, we just replace the estimation with central 
                        // as central is 4x the intensity, we have to divide by 4 
                        val += 0.0625 * throughputCentral[centerIndex];

                    if (x > 0)
                        val += 0.25 * throughputRight[leftIndex];
                    else
                        val += 0.0625 * throughputCentral[centerIndex];

                    if (y + 1 < size.y)
                        val += 0.25 * throughputTop[bottomIndex];
                    else
                        val += 0.0625 * throughputCentral[centerIndex];

                    if (y > 0)
                        val += 0.25 * throughputBottom[topIndex];
                    else
                        val += 0.0625 * throughputCentral[centerIndex];

                    throughputBuffer->addSample(x, y, val);
                }
            }
        }

        auto fusion_gradient = [](Bitmap* destX, Bitmap* forwardX, Bitmap* inverseX,
                Bitmap* destY, Bitmap* forwardY, Bitmap* inverseY, int inverse = -1) -> void {
            destX->clear();
            destY->clear();

            auto dx = (Spectrum *) destX->getUInt8Data();
            auto dy = (Spectrum *) destY->getUInt8Data();
            auto dxForward = (Spectrum *) forwardX->getFloatData();
            auto dxInverse = (Spectrum *) inverseX->getFloatData();
            auto dyForward = (Spectrum *) forwardY->getFloatData();
            auto dyInverse = (Spectrum *) inverseY->getFloatData();

            Vector2i size = destX->getSize();
            int width = size.x;
            for (int y = 0; y < size.y; ++y) {
                for (int x = 0; x < size.x; ++x) {

                    int centerIndex = y * width + x;
                    int leftIndex = y * width + x - 1;
                    int topIndex = (y - 1) * width + x;

                    dx[centerIndex] += dxForward[centerIndex];
                    if (x > 0)
                        dx[leftIndex] += inverse*dxInverse[centerIndex];

                    dy[centerIndex] += dyForward[centerIndex];
                    if (y > 0)
                        dy[topIndex] += inverse*dyInverse[centerIndex];
                }
            }
        };

        /* Manage gradient buffers for variance and buffer A, B if needed */
        if (!m_config.pathReuse) {
            {
                ref<Bitmap> dxBitmap = gradientXBuffer->getBuffer(0);
                ref<Bitmap> dyBitmap = gradientYBuffer->getBuffer(0);
                ref<Bitmap> dxForwardBitmap = gradientXTmpBuffer->getBuffer(0);
                ref<Bitmap> dxInverseBitmap = gradientXTmpBuffer->getBuffer(1);
                ref<Bitmap> dyForwardBitmap = gradientYTmpBuffer->getBuffer(0);
                ref<Bitmap> dyInverseBitmap = gradientYTmpBuffer->getBuffer(1);
                fusion_gradient(dxBitmap.get(), dxForwardBitmap.get(), dxInverseBitmap.get(),
                        dyBitmap.get(), dyForwardBitmap.get(), dyInverseBitmap.get());
            }
            if (bufferAB) {
                ref<Bitmap> dxBitmap = gradientXBuffer->getBufferA(0);
                ref<Bitmap> dyBitmap = gradientYBuffer->getBufferA(0);
                ref<Bitmap> dxForwardBitmap = gradientXTmpBuffer->getBufferA(0);
                ref<Bitmap> dxInverseBitmap = gradientXTmpBuffer->getBufferA(1);
                ref<Bitmap> dyForwardBitmap = gradientYTmpBuffer->getBufferA(0);
                ref<Bitmap> dyInverseBitmap = gradientYTmpBuffer->getBufferA(1);
                fusion_gradient(dxBitmap.get(), dxForwardBitmap.get(), dxInverseBitmap.get(),
                                dyBitmap.get(), dyForwardBitmap.get(), dyInverseBitmap.get());
            }

            if (bufferAB) {
                ref<Bitmap> dxBitmap = gradientXBuffer->getBufferB(0);
                ref<Bitmap> dyBitmap = gradientYBuffer->getBufferB(0);
                ref<Bitmap> dxForwardBitmap = gradientXTmpBuffer->getBufferB(0);
                ref<Bitmap> dxInverseBitmap = gradientXTmpBuffer->getBufferB(1);
                ref<Bitmap> dyForwardBitmap = gradientYTmpBuffer->getBufferB(0);
                ref<Bitmap> dyInverseBitmap = gradientYTmpBuffer->getBufferB(1);
                fusion_gradient(dxBitmap.get(), dxForwardBitmap.get(), dxInverseBitmap.get(),
                                dyBitmap.get(), dyForwardBitmap.get(), dyInverseBitmap.get());
            }

            if (m_config.computeVariance) {
                // Estimate gradient variance: var(A - B) = var(A) + var(B) since forward and inverse shift are independent
                ref<Bitmap> dxVarianceBitmap = gradientXBuffer->getSampleVariance(0);
                ref<Bitmap> dyVarianceBitmap = gradientYBuffer->getSampleVariance(0);
                ref<Bitmap> dxForwardVarianceBitmap = gradientXTmpBuffer->getSampleVariance(0);
                ref<Bitmap> dxInverseVarianceBitmap = gradientXTmpBuffer->getSampleVariance(1);
                ref<Bitmap> dyForwardVarianceBitmap = gradientYTmpBuffer->getSampleVariance(0);
                ref<Bitmap> dyInverseVarianceBitmap = gradientYTmpBuffer->getSampleVariance(1);
                fusion_gradient(dxVarianceBitmap.get(), dxForwardVarianceBitmap.get(), dxInverseVarianceBitmap.get(),
                                dyVarianceBitmap.get(), dyForwardVarianceBitmap.get(), dyInverseVarianceBitmap.get(),
                                1);
            }
        }

        /* Reconstruct. */
        if (m_dumpIteration > 0 && currentIteration % m_dumpIteration == 0) {
            // === Update the log time
            unsigned int milliseconds = renderingTimer->getMilliseconds();
            timeFile << (milliseconds / 1000.f) << ",\n";
            timeFile.flush();
            Log(EInfo, "Rendering time: %i, %i", milliseconds / 1000,
                milliseconds % 1000);
            renderingTimer->reset();
            
            ref<Film> hdrFilm = newHDRFilm(scene);
            // Do not include reconstruction time
            // Perform the reconstruction
            Reconstruction rec{};
            rec.reconstructL1 = m_config.m_reconstructL1;
            rec.reconstructL2 = m_config.m_reconstructL2;
            rec.reconstructUni = m_config.m_reconstructUni;
            rec.reconstructWeighted = m_config.m_reconstructWeighted;
            rec.reconstructL2Weighted = m_config.m_reconstructL2Weighted;

            rec.alpha = (float) m_config.m_reconstructAlpha;

            ref<Bitmap> directBitmap(new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getCropSize()));
            film->developMulti(Point2i(0, 0), film->getCropSize(), Point2i(0, 0), directBitmap, BUFFER_VERY_DIRECT);
            develop(scene, hdrFilm, directBitmap, currentIteration, "_direct_");
            auto directVector = bitmap2vec(directBitmap);
            if (!m_config.directTracing) {
                auto subPixelCount = size_t(3 * directBitmap->getSize().x * directBitmap->getSize().y);
                directVector = std::vector<float>(subPixelCount, 0.f);
            }

            bool reconstructGradient = m_config.m_reconstructL1 || m_config.m_reconstructL2 || m_config.m_reconstructL2Weighted || 
                                       m_config.m_reconstructUni || m_config.m_reconstructWeighted;
            if (reconstructGradient)
            {
                /* Develop primal and gradient data into bitmaps. */
                ref<Bitmap> throughputBitmap(new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getCropSize()));
                ref<Bitmap> dxBitmap(new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getCropSize()));
                ref<Bitmap> dyBitmap(new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getCropSize()));

                film->developMulti(Point2i(0, 0), film->getCropSize(), Point2i(0, 0), throughputBitmap,
                                   BUFFER_THROUGHPUT);
                film->developMulti(Point2i(0, 0), film->getCropSize(), Point2i(0, 0), dxBitmap, BUFFER_DX);
                film->developMulti(Point2i(0, 0), film->getCropSize(), Point2i(0, 0), dyBitmap, BUFFER_DY);
                develop(scene, hdrFilm, throughputBitmap, currentIteration, "_throughput_");

                

                {
                    // For checking 
                    /*ref<Bitmap> bufferThroughputBitmap = throughputBuffer->getBuffer(0);
                    develop(scene, hdrFilm, bufferThroughputBitmap, currentIteration, "_buffer_throughput_");
                    ref<Bitmap> bufferThroughputVarianceBitmap = throughputBuffer->getSampleVariance(0);
                    develop(scene, hdrFilm, bufferThroughputVarianceBitmap, currentIteration,
                            "_buffer_throughput_variance_");

                    ref<Bitmap> bufferGradientXBitmap = gradientXBuffer->getBuffer(0);
                    develop(scene, hdrFilm, bufferGradientXBitmap, currentIteration, "_buffer_GX_");
                    ref<Bitmap> bufferGradientXVarianceBitmap = gradientXBuffer->getSampleVariance(0);
                    develop(scene, hdrFilm, bufferGradientXVarianceBitmap, currentIteration, "_buffer_GX_variance_");*/


                    auto throughputVector = bitmap2vec(throughputBitmap);
                    auto dxVector = bitmap2vec(dxBitmap);
                    auto dyVector = bitmap2vec(dyBitmap);

                    // And variance in our custom accumulation buffer
                    auto variance = [&]() -> Reconstruction::Variance {
                        if (m_config.computeVariance) {
                            return {
                                    primal: bitmap2vec(throughputBuffer->getSampleVariance(0)),
                                    dx: bitmap2vec(gradientXBuffer->getSampleVariance(0)),
                                    dy: bitmap2vec(gradientYBuffer->getSampleVariance(0))
                            };
                        } else {
                            return {};
                        }
                    }();

                    auto rec_results = rec.reconstruct(film->getCropSize(),
                                                       throughputVector, dxVector, dyVector, directVector,
                                                       variance,
                                                       PostProcessOption{
                                                               forceBlackPixels: m_config.forceBlackPixels,
                                                               clampingValues: true
                                                       });
                    if (rec_results.size() == 1) {
                        develop(scene, hdrFilm, rec_results[0].img, currentIteration, "_recons_");
                    } else {
                        for (auto &result: rec_results) {
                            develop(scene, hdrFilm, result.img, currentIteration, "_" + result.name + "_");
                        }
                    }
                }

                if (m_showAbsGradients) {
                    // Create new bitmap to store abs gradient values
                    ref<Bitmap> dxBitmapAbs = dxBitmap->clone();
                    ref<Bitmap> dyBitmapAbs = dyBitmap->clone();

                    // copy abs gradient values
                    auto *tDX = (Spectrum *) dxBitmap->getUInt8Data();
                    auto *tDXAbs = (Spectrum *) dxBitmapAbs->getUInt8Data();
                    auto *tDY = (Spectrum *) dyBitmap->getUInt8Data();
                    auto *tDYAbs = (Spectrum *) dyBitmapAbs->getUInt8Data();
                    for (int y = 0; y < cropSize.y; ++y) {
                        for (int x = 0; x < cropSize.x; ++x) {
                            int pixelID = y * cropSize.x + x;
                            tDXAbs[pixelID] = tDX[pixelID].abs();
                            tDYAbs[pixelID] = tDY[pixelID].abs();
                        }
                    }

                    develop(scene, hdrFilm, dxBitmapAbs, currentIteration, "_dxAbs_");
                    develop(scene, hdrFilm, dyBitmapAbs, currentIteration, "_dyAbs_");
                }
            }

            if (m_config.m_reconstructNfor && (!m_config.reusePrimal || (m_config.reusePrimal && currentIteration % 2 == 0)))
            {                
                // We only do NFOR at even iterations due to how we track bufferA and B in reusePrimal case
                ref<Bitmap> throughputBitmap = throughputBuffer->getBuffer(0);
                ref<Bitmap> dxBitmap = gradientXBuffer->getBuffer(0);
                ref<Bitmap> dyBitmap = gradientYBuffer->getBuffer(0);

                if (! reconstructGradient) {
                    // When no specific technique for reconstruction is requested, use weighted control variates.
                    rec.reconstructWeighted = true;
                }

                // Feature buffer
                auto throughputVector = bitmap2vec(throughputBuffer->getBuffer(0));
                auto gradientXVector = bitmap2vec(gradientXBuffer->getBuffer(0));
                auto gradientYVector = bitmap2vec(gradientYBuffer->getBuffer(0));
                auto results = rec.reconstruct(film->getSize(),
                                               throughputVector,
                                               gradientXVector,
                                               gradientYVector,
                                               directVector,
                                               Reconstruction::Variance{// Variance
                                                       primal: bitmap2vec(throughputBuffer->getSampleVariance(0)),
                                                       dx: bitmap2vec(gradientXBuffer->getSampleVariance(0)),
                                                       dy: bitmap2vec(gradientYBuffer->getSampleVariance(0))
                                               },
                                               PostProcessOption{
                                                       forceBlackPixels: m_config.forceBlackPixels,
                                                       clampingValues: true
                                               });

                auto throughputVectorA = bitmap2vec(throughputBuffer->getBufferA(0));
                auto gradientXVectorA = bitmap2vec(gradientXBuffer->getBufferA(0));
                auto gradientYVectorA = bitmap2vec(gradientYBuffer->getBufferA(0));
                auto resultsA = rec.reconstruct(film->getSize(),
                                                throughputVectorA,
                                                gradientXVectorA,
                                                gradientYVectorA,
                                                directVector,
                                                Reconstruction::Variance{// reuse variance of the main buffer
                                                        primal: bitmap2vec(throughputBuffer->getSampleVariance(0)),
                                                        dx: bitmap2vec(gradientXBuffer->getSampleVariance(0)),
                                                        dy: bitmap2vec(gradientYBuffer->getSampleVariance(0))
                                                },
                                                PostProcessOption{
                                                        forceBlackPixels: m_config.forceBlackPixels,
                                                        clampingValues: true
                                                });

                auto throughputVectorB = bitmap2vec(throughputBuffer->getBufferB(0));
                auto gradientXVectorB = bitmap2vec(gradientXBuffer->getBufferB(0));
                auto gradientYVectorB = bitmap2vec(gradientYBuffer->getBufferB(0));
                auto resultsB = rec.reconstruct(film->getSize(),
                                                throughputVectorB,
                                                gradientXVectorB,
                                                gradientYVectorB,
                                                directVector,
                                                Reconstruction::Variance{// reuse variance of the main buffer
                                                        primal: bitmap2vec(throughputBuffer->getSampleVariance(0)),
                                                        dx: bitmap2vec(gradientXBuffer->getSampleVariance(0)),
                                                        dy: bitmap2vec(gradientYBuffer->getSampleVariance(0))
                                                },
                                                PostProcessOption{
                                                        forceBlackPixels: m_config.forceBlackPixels,
                                                        clampingValues: true
                                                });

                // Do NFOR for all reconstructed images
                for (auto i = 0; i < results.size(); ++i) {
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
                            Spectrum mean = recon[i]; // TODO: This might be not optimal as the reconstruction are not
                            // TODO: independent
                            variance[i] = 0.5 * ((bufferA[i] - mean) * (bufferA[i] - mean) +
                                                 (bufferB[i] - mean) * (bufferB[i] - mean));
                        }
                    }

                    ref<AccumBuffer> reconstructionBuffer = new AccumBuffer(
                            {reconstructionBitmap},
                            {reconstructionVarianceBitmap},
                            {reconstructionBitmapA},
                            {reconstructionBitmapB},
                            currentIteration
                    );

                    // Filter the primal with reconstruction as a feature map. Note that for feature, only bufferA, B, and variance is used.
                    // The final feature is the average of filtered feature bufferA and B, as default by NFOR implementation.
                    // Obviously for gradient-domain reconstruction taking such average is sub-optimal. We should pass the reconstruction
                    // of full samples and filter this instead, but for now we just use the default implementation of NFOR.

                    // FIXME: throughput bitmap has no very direct, and reconstruction has very direct. Might affect the filtering

                    // Traditional NFOR
                    /*
                    {
                        ref<Nfor> nfor = new Nfor();
                        ref<Bitmap> nforBitmap = reconstructionBitmap->clone();
                        nfor->denoise(throughputBuffer, nforBitmap, reconstructionBuffer, {0});
                        if (m_config.directTracing) {
                            nforBitmap->accumulate(directBitmap);
                        }
                        develop(scene, hdrFilm, nforBitmap, currentIteration, "_" + results[i].name + "_nfor_");
                    }*/

                    // Custom NFOR
                    {
                        ref<Nfor> nfor = new Nfor();
                        ref<Bitmap> nforBitmap = reconstructionBitmap->clone();
                        nfor->denoise(throughputBuffer, nforBitmap, reconstructionBuffer, {0}, true);
                        if (m_config.directTracing) {
                            nforBitmap->accumulate(directBitmap);
                        }

                        if (reconstructGradient)
                            develop(scene, hdrFilm, nforBitmap, currentIteration, "_" + results[i].name + "_nfor_custom_");
                        else {
                            develop(scene, hdrFilm, nforBitmap, currentIteration, "_" + results[i].name + "_");
                        }
                    }
                }
            }

            // Reset timer because we exclude reconstruction timing (treat it as postprocess)
            renderingTimer->reset();
        }

       
        currentIteration += 1;
    }

    timeFile.close();

    return true; // Always true in this case
}

static Float miWeight(Float pdfA, Float pdfB) {
    pdfA *= pdfA;
    pdfB *= pdfB;
    return pdfA / (pdfA + pdfB);
}

Spectrum GradientPathIntegrator::Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
    // Duplicate of MIPathTracer::Li to support sub-surface scattering initialization.

    /* Some aliases and local variables */
    const Scene *scene = rRec.scene;
    Intersection &its = rRec.its;
    RayDifferential ray(r);
    Spectrum Li(0.0f);
    bool scattered = false;

    /* Perform the first ray intersection (or ignore if the
        intersection has already been provided). */
    rRec.rayIntersect(ray);
    ray.mint = Epsilon;

    Spectrum throughput(1.0f);
    Float eta = 1.0f;

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
                const auto *emitter = dynamic_cast<const Emitter *>(dRec.object);

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

void GradientPathIntegrator::serialize(Stream *stream, InstanceManager *manager) const {
    MonteCarloIntegrator::serialize(stream, manager);
    stream->writeFloat(m_config.m_shiftThreshold);
    stream->writeBool(m_config.m_reconstructL1);
    stream->writeBool(m_config.m_reconstructL2);
    stream->writeFloat(m_config.m_reconstructAlpha);
}

std::string GradientPathIntegrator::toString() const {
    std::ostringstream oss;
    oss << "GradientPathTracer[" << endl
        << "  maxDepth = " << m_maxDepth << "," << endl
        << "  rrDepth = " << m_rrDepth << "," << endl
        << "  shiftThreshold = " << m_config.m_shiftThreshold << endl
        << "  reconstructL1 = " << m_config.m_reconstructL1 << endl
        << "  reconstructL2 = " << m_config.m_reconstructL2 << endl
        << "  reconstructAlpha = " << m_config.m_reconstructAlpha << endl
        << "]";
    return oss.str();
}

MTS_IMPLEMENT_CLASS_S(GradientPathIntegrator, false, MonteCarloIntegrator)

MTS_EXPORT_PLUGIN(GradientPathIntegrator, "Gradient Path Integrator");
MTS_NAMESPACE_END
