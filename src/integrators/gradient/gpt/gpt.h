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

#if !defined(__GPT_H)
#define __GPT_H

#include <mitsuba/mitsuba.h>
#include "gpt_wr.h"

#include "../../volume_utils.h"
#include "../../denoiser/accum_buffer.h"

MTS_NAMESPACE_BEGIN


/* ==================================================================== */
/*                         Configuration storage                        */
/* ==================================================================== */

/// Configuration for the gradient path tracer.
enum ShiftStrategy {
    EShiftOriginal,
    EShiftOriginalVolume,
    EShiftExplicit,
    EShiftRandom,
    EShiftRandomReconnect,
};
struct GradientPathTracerConfig {
    int m_maxDepth;
    int m_minDepth;
    int m_minCameraDepth;
    int m_rrDepth;
    bool m_strictNormals;
    Float m_shiftThreshold;
    bool m_reconstructL1;
    bool m_reconstructL2;
    bool m_reconstructUni;
    bool m_reconstructWeighted;
    bool m_reconstructL2Weighted;
    bool m_reconstructNfor;
    Float m_reconstructAlpha;
    bool forceBlackPixels;
    ELightingEffects m_lightingInteractionMode;

    // If we want to compute the variance
    bool computeVariance;
    ShiftStrategy shiftmapping;
    bool pathReuse;
    Float rrShift;
    bool reusePrimal;
    bool directTracing;
};



/* ==================================================================== */
/*                         Integrator                         */
/* ==================================================================== */
class GradientPathIntegrator : public MonteCarloIntegrator {
public:
    explicit GradientPathIntegrator(const Properties &props);

    /// Unserialize from a binary data stream
    GradientPathIntegrator(Stream *stream, InstanceManager *manager);


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

        if (!(m_config.m_lightingInteractionMode & EAll2Surf)) {
            Log(EInfo, "Only volume rendering, change probability to sample the medium");
            auto medium = scene->getMedia().begin();
            while (medium != scene->getMedia().end()) {
                const_cast<Medium *>(medium->get())->computeOnlyVolumeInteraction();
                medium++;
            }
        }

        return true;
    }

    /// Starts the rendering process.
    bool render(Scene *scene,
                RenderQueue *queue, const RenderJob *job,
                int sceneResID, int sensorResID, int samplerResID);


    /// Renders a block in the image.
    void renderBlock(const Scene *scene, const Sensor *sensor, Sampler *sampler, GPTWorkResult *block,
                     const bool &stop, const std::vector<TPoint2<uint8_t> > &points, const Vector2i& rectSize) const;

    void serialize(Stream *stream, InstanceManager *manager) const;

    std::string toString() const;


    /// Used by Mitsuba for initializing sub-surface scattering.
    Spectrum Li(const RayDifferential &ray, RadianceQueryRecord &rRec) const;


    MTS_DECLARE_CLASS()

protected:


private:
    GradientPathTracerConfig m_config;

    // Debug options
    bool m_showAbsGradients;
    int m_dumpIteration;

    ref<AccumBuffer> throughputBuffer;
    ref<AccumBuffer> throughputTmpBuffer;
    ref<AccumBuffer> gradientXBuffer;
    ref<AccumBuffer> gradientYBuffer;
    ref<AccumBuffer> gradientXTmpBuffer;
    ref<AccumBuffer> gradientYTmpBuffer;    
};


MTS_NAMESPACE_END

#endif /* __GBDPT_H */
