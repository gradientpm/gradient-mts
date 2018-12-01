#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/bidir/path.h>

#ifndef MITSUBA_VOLUME_UTILS_H
#define MITSUBA_VOLUME_UTILS_H

MTS_NAMESPACE_BEGIN

enum EVolumeTechnique {
    EVolBRE2D,
    EVolBRE3D,
    EDistance,
    EBeamBeam1D,
    EBeamBeam3D_Naive,
    EBeamBeam3D_EGSR,
    EBeamBeam3D_Optimized,
    EVolPlane0D
};

struct EVolumeTechniqueHelper {
    static bool useAPA(EVolumeTechnique volTechnique) {
        return volTechnique == EVolPlane0D ||
               volTechnique == EBeamBeam1D ||
               volTechnique == EBeamBeam3D_Naive ||
               volTechnique == EBeamBeam3D_EGSR ||
               volTechnique == EBeamBeam3D_Optimized ||
               volTechnique == EVolBRE2D ||
               volTechnique == EVolBRE3D;
    }

    static bool use3DKernel(EVolumeTechnique volTechnique) {
        return volTechnique == EDistance ||
               volTechnique == EVolBRE3D ||
               volTechnique == EBeamBeam3D_Naive ||
               volTechnique == EBeamBeam3D_EGSR ||
               volTechnique == EBeamBeam3D_Optimized;
    }

    static bool use2DKernel(EVolumeTechnique volTechnique) {
        return volTechnique == EVolBRE2D;
    }

    static bool needBeams(EVolumeTechnique volTechnique) {
        return volTechnique == EBeamBeam1D ||
               volTechnique == EBeamBeam3D_Naive ||
               volTechnique == EBeamBeam3D_EGSR ||
               volTechnique == EBeamBeam3D_Optimized ||
               volTechnique == EVolPlane0D;
    }
};

inline EVolumeTechnique parseVolumeTechnique(std::string volRenderingTech) {
    EVolumeTechnique volTechnique;
    if (volRenderingTech == "distance") {
        volTechnique = EDistance;
    }
    else if (volRenderingTech == "bre") {
        volTechnique = EVolBRE3D;                   // default BRE is 3D
    }
    else if (volRenderingTech == "bre2d") {
        volTechnique = EVolBRE2D;
    }
    else if (volRenderingTech == "bre3d") {
        volTechnique = EVolBRE3D;
    }
    else if (volRenderingTech == "beam") {
        volTechnique = EBeamBeam1D;                 // default beam
    }
    else if (volRenderingTech == "beam1d") {
        volTechnique = EBeamBeam1D;                 // default beam 1D strategy
    }
    else if (volRenderingTech == "beam3d") {
        volTechnique = EBeamBeam3D_Optimized;       // default beam 3D strategy
    }
    else if (volRenderingTech == "beam3d_naive") {
        volTechnique = EBeamBeam3D_Naive;
    }
    else if (volRenderingTech == "beam3d_egsr") {
        volTechnique = EBeamBeam3D_EGSR;
    }
    else if (volRenderingTech == "beam3d_optimized") {
        volTechnique = EBeamBeam3D_Optimized;
    } else if(volRenderingTech == "plane0d") {
        volTechnique = EVolPlane0D;
    } else {
        SLog(EError, "Unknow vol technique: %s", volRenderingTech.c_str());
        volTechnique = EDistance; //< by default
    }
    return volTechnique;
}

enum ELightingEffects {
  ESurf2Surf = 1<<1,
  ESurf2Media = 1<<2,
  EMedia2Surf = 1<<3,
  EMedia2Media = 1<<4,
  EAll2Surf = ESurf2Surf | EMedia2Surf,
  EAll2Media = ESurf2Media | EMedia2Media,
  EAll2All = EAll2Media | EAll2Surf
};

inline bool needRenderSurface(ELightingEffects opt) {
  return opt & EAll2Surf;
}

inline ELightingEffects parseMediaInteractionMode(const std::string& strLightingMode) {
    if(strLightingMode == "all2media") {
        return EAll2Media;
    } else if(strLightingMode == "all2surf") {
        return EAll2Surf;
    } else if(strLightingMode == "surf2surf") {
        return ESurf2Surf;
    } else if(strLightingMode == "media2surf") {
        return EMedia2Surf;
    } else if(strLightingMode == "surf2media") {
        return ESurf2Media;
    } else if(strLightingMode == "media2media") {
        return EMedia2Media;
    } else if(strLightingMode == "all2all") {
        return EAll2All;
    } else {
        SLog(EError, "Invalid media interaction mode: %s", strLightingMode.c_str());
        return EAll2All;
    }
}

// Derived from GPMConfig's getInteractionMode()
inline int parseBsdfInteractionMode(const std::string& mode) {
    std::string lowMode = mode;
    std::transform(lowMode.begin(), lowMode.end(), lowMode.begin(), ::tolower);

    if (lowMode == "all") {
        return BSDF::EAll;
    }
    else if (lowMode == "diffuse") {
        return BSDF::EDiffuse;
    }
    else if (lowMode == "specular") {
        return BSDF::EDelta;
    }
    else if (lowMode == "glossy") {
        return BSDF::EGlossy;
    }
    else {
        SLog(EError, "Invalid BSDF interaction mode: %s", lowMode.c_str());
        return BSDF::ENull;
    }
}

//assume beam is not zero length
inline bool isIntersectedPoint(const Point& sensorPos, const Point& org, const Point& dest, const Float radius){

    Vector beam = dest - org;
    Float beamLengthSq = beam.lengthSquared();

    //add epsilon here ?
    if (beamLengthSq == 0.f){
        //we should not be there
        SLog(EError, "Zero length beam");
    }

    Float t = math::clamp(dot(sensorPos - org, beam) / beamLengthSq, (Float)0.0, (Float)1.0);
    Vector v = (org + t * beam) - sensorPos;

    return (radius*radius) > v.lengthSquared();
}

inline int numberOfBeams(const Path* lt, int minDepth = 0) {
    int nbValidBeams = 0;
    size_t m =  static_cast<size_t>(std::max(1, minDepth));
    for(size_t i = m; i < lt->edgeCount(); i++) {
        if (lt->edge(i)->medium != 0) {
            nbValidBeams += 1;
        }
    }
    return nbValidBeams;
}

inline int numberOfVolumePhotons(const Path* lt, int minDepth = 0) {
    int nbValidVolPhoton = 0;
    size_t m = static_cast<size_t>(std::max(1, minDepth)) + 1;
    for(size_t i = m; i < lt->vertexCount(); i++) {
        if (lt->vertex(i)->isMediumInteraction()) {
            nbValidVolPhoton += 1;
        }
    }
    return nbValidVolPhoton;
}

inline bool bdptDecisionVolume(const std::vector<PathVertex*>& pathVertex, size_t minCameraDepth) {
    size_t lV = 1;

    while(lV < pathVertex.size()) {
        if(pathVertex[lV]->isMediumInteraction()) {
            return minCameraDepth < lV;
        }
        if(pathVertex[lV]->isSurfaceInteraction()) {
            const Intersection& itsL = pathVertex[lV]->getIntersection();
            const BSDF* bsdfL = itsL.getBSDF();
            int16_t comp = pathVertex[lV]->sampledComponentIndex;
            if(bsdfL->getRoughness(itsL, comp == -1 ? 0 : comp) > 0.05) {
                return false;
            }
        }

        lV++;
    }

    return false;
}

inline Float powerOfTwo(Float v) {
    return v*v;
}

const Float POURCENTAGE_BS = 0.01;
inline AABB max_AABB_medium(const Scene *scene) {
    // Just print this information as it can be usefull
    // for setting the scene
    BSphere sceneBSphere = scene->getBSphere();
    SLog(EInfo, "AABB Scene (including sensor): %s", sceneBSphere.toString().c_str());

    // Scene checking: If there is no exterior volume for the smoke, raise an exception
    std::vector<AABB> mediumAABB;
    for(auto media: scene->getMedia()) {
        SLog(EInfo, "Compute AABB for the media: %s", media->toString().c_str());
        AABB currMediaAABB;
        currMediaAABB.reset();
        const Medium* m = media.get();
        bool foundAttachedMesh = false;
        // Check the meshes
        for(auto mesh: scene->getShapes()) {
            if(mesh->getExteriorMedium() != NULL) {
                foundAttachedMesh = foundAttachedMesh || (mesh->getExteriorMedium() == m);
                SLog(EInfo, " - SKIP: Found medium attached as exterior for: %s", mesh->toString().c_str());
            }
            if(mesh->getInteriorMedium() != NULL) {
                foundAttachedMesh = foundAttachedMesh || (mesh->getInteriorMedium() == m);
                SLog(EInfo, " - EXPAND: Found medium attached as interior for: %s", mesh->toString().c_str());
                currMediaAABB.expandBy(mesh->getAABB());
            }
        }
        if(!foundAttachedMesh) {
            SLog(EError, "No shape define interior/exterior limit for the media: %s", m->toString().c_str());
        }
        mediumAABB.emplace_back(currMediaAABB);
    }

    std::sort(mediumAABB.begin(), mediumAABB.end(),
              [](const AABB & a, const AABB & b) -> bool { return a.getVolume() > b.getVolume(); }
    );
    SLog(EInfo, "Max AABB Smoke: %s", mediumAABB.back().toString().c_str());

    return std::move(mediumAABB.back());
}

MTS_NAMESPACE_END

#endif //MITSUBA_VOLUME_UTILS_H
