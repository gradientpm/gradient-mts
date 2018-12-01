
#include <mitsuba/mitsuba.h>

#ifndef MITSUBA_GPT_SHIFTMAPPING_H
#define MITSUBA_GPT_SHIFTMAPPING_H

MTS_NAMESPACE_BEGIN

/// A threshold to use in positive denominators to avoid division by zero.
const Float D_EPSILON = (Float)(1e-19);

// Change some checking which can introduce few errors
// Compare to the orignal code...
#define ORIGINAL_TESTS 0

#define DO_DIRECT_COMPUTATION 1
#define DO_BSDF_COMPUTATION 1

/// Returns whether point1 sees point2.
bool testVisibility(const Scene *scene, const Point3 &point1, const Point3 &point2, Float time) {
    Ray shadowRay;
    shadowRay.setTime(time);
    shadowRay.setOrigin(point1);
    shadowRay.setDirection(point2 - point1);
    shadowRay.mint = Epsilon;
    shadowRay.maxt = (Float) 1.0 - ShadowEpsilon;

    return !scene->rayIntersect(shadowRay);
}

/// Returns whether the given ray sees the environment.
bool testEnvironmentVisibility(const Scene *scene, const Ray &ray, Point3 shiftPos) {
    const Emitter *env = scene->getEnvironmentEmitter();
    if (!env) {
        return false;
    }

    Ray shadowRay(ray);
    shadowRay.setTime(ray.time);
    shadowRay.setOrigin(shiftPos);
    shadowRay.setDirection(ray.d);

    DirectSamplingRecord directSamplingRecord;
    env->fillDirectSamplingRecord(directSamplingRecord, shadowRay);

    shadowRay.mint = Epsilon;
    shadowRay.maxt = ((Float) 1.0 - ShadowEpsilon) * directSamplingRecord.dist;

    return !scene->rayIntersect(shadowRay);
}


/// Classification of vertices into diffuse and glossy.
enum VertexType {
    VERTEX_TYPE_GLOSSY,     ///< "Specular" vertex that requires the half-vector duplication shift.
    VERTEX_TYPE_DIFFUSE     ///< "Non-specular" vertex that is rough enough for the reconnection shift.
};

enum RayConnection {
    RAY_NOT_CONNECTED,      ///< Not yet connected - shifting in progress.
    RAY_RECENTLY_CONNECTED, ///< Connected, but different incoming direction so needs a BSDF evaluation.
    RAY_CONNECTED           ///< Connected, allows using BSDF values from the base path.
};


/// Describes the state of a ray that is being traced in the scene.
struct RayState {
    RayState()
            : radiance(0.0f),
              gradient(0.0f),
              eta(1.0f),
              pdf(1.0f),
              throughput(Spectrum(0.0f)),
              alive(true),
              connection_status(RAY_NOT_CONNECTED) {}

    /// Adds radiance to the ray.
    inline void addRadiance(const Spectrum &contribution, Float weight) {
        Spectrum color = contribution * weight;
        radiance += color;
    }

    inline void scaleValue(Float v) {
        radiance *= v;
        gradient *= v;
    }

    inline const Point &prevPos() const {
        if (lastVolume) {
            return mRec.p;
        } else {
            return rRec.its.p;
        }
    }

    /// Adds gradient to the ray.
    inline void addGradient(const Spectrum &contribution, Float weight) {
        Spectrum color = contribution * weight;
        gradient += color;
    }

    RayDifferential ray;             ///< Current ray.

    Spectrum throughput;             ///< Current throughput of the path.
    Float pdf;                       ///< Current PDF of the path.

    // Note: Instead of storing throughput and pdf, it is possible to store Veach-style weight (throughput divided by pdf), if relative PDF (offset_pdf divided by base_pdf) is also stored. This might be more stable numerically.

    Spectrum radiance;               ///< Radiance accumulated so far.
    Spectrum gradient;               ///< Gradient accumulated so far.

    RadianceQueryRecord rRec;        ///< The radiance query record for this ray.
    Float eta;                       ///< Current refractive index of the ray.
    bool alive;                      ///< Whether the path matching to the ray is still good. Otherwise it's an invalid offset path with zero PDF and throughput.

    RayConnection connection_status; ///< Whether the ray has been connected to the base path, or is in progress.
    MediumSamplingRecord mRec;
    bool lastVolume = false;
};

/// Returns the vertex type of a vertex by its roughness value.
VertexType getVertexTypeByRoughness(Float roughness, const GradientPathTracerConfig &config) {
    if (roughness <= config.m_shiftThreshold) {
        return VERTEX_TYPE_GLOSSY;
    } else {
        return VERTEX_TYPE_DIFFUSE;
    }
}

/// Returns the vertex type (diffuse / glossy) of a vertex, for the purposes of determining
/// the shifting strategy.
///
/// A bare classification by roughness alone is not good for multi-component BSDFs since they
/// may contain a diffuse component and a perfect specular component. If the base path
/// is currently working with a sample from a BSDF's smooth component, we don't want to care
/// about the specular component of the BSDF right now - we want to deal with the smooth component.
///
/// For this reason, we vary the classification a little bit based on the situation.
/// This is perfectly valid, and should be done.
VertexType getVertexType(const BSDF *bsdf, Intersection &its,
                         const GradientPathTracerConfig &config, unsigned int bsdfType) {
    // Return the lowest roughness value of the components of the vertex's BSDF.
    // If 'bsdfType' does not have a delta component, do not take perfect speculars (zero roughness) into account in this.
    Float lowest_roughness = std::numeric_limits<Float>::infinity();

    bool found_smooth = false;
    bool found_dirac = false;
    for (int i = 0, component_count = bsdf->getComponentCount(); i < component_count; ++i) {
        Float component_roughness = bsdf->getRoughness(its, i);

        if (component_roughness == Float(0)) {
            found_dirac = true;
            if (!(bsdfType & BSDF::EDelta)) {
                // Skip Dirac components if a smooth component is requested.
                continue;
            }
        } else {
            found_smooth = true;
        }

        if (component_roughness < lowest_roughness) {
            lowest_roughness = component_roughness;
        }
    }

    // Roughness has to be zero also if there is a delta component but no smooth components.
    if (!found_smooth && found_dirac && !(bsdfType & BSDF::EDelta)) {
        lowest_roughness = Float(0);
    }

    return getVertexTypeByRoughness(lowest_roughness, config);
}

VertexType getVertexType(RayState &ray, const GradientPathTracerConfig &config, unsigned int bsdfType) {
    const BSDF *bsdf = ray.rRec.its.getBSDF(ray.ray);
    return getVertexType(bsdf, ray.rRec.its, config, bsdfType);
}


/// Result of a half-vector duplication shift.
struct HalfVectorShiftResult {
    bool success;   ///< Whether the shift succeeded.
    Float jacobian; ///< Local Jacobian determinant of the shift.
    Vector3 wo;     ///< Tangent space outgoing vector for the shift.
};

HalfVectorShiftResult halfVectorShiftVolume(Vector3 tangentSpaceMainWi,
                                            Vector3 tangentSpaceMainWo,
                                            Vector3 tangentSpaceShiftedWi) {
    HalfVectorShiftResult result;

    // Just copy the H vector (as we are inside the volume)
    Vector tangentSpaceHalfVector = normalize(tangentSpaceMainWi + tangentSpaceMainWo);
    Vector tangentSpaceShiftedWo = reflect(tangentSpaceShiftedWi, tangentSpaceHalfVector);

    Float WoDotH = dot(tangentSpaceShiftedWo, tangentSpaceHalfVector) / dot(tangentSpaceMainWo, tangentSpaceHalfVector);
    Float jacobian = abs(WoDotH);

    result.success = true;
    result.wo = tangentSpaceShiftedWo;
    result.jacobian = jacobian;

    return result;
}


/// Calculates the outgoing direction of a shift by duplicating the local half-vector.
HalfVectorShiftResult
halfVectorShift(Vector3 tangentSpaceMainWi, Vector3 tangentSpaceMainWo, Vector3 tangentSpaceShiftedWi, Float mainEta,
                Float shiftedEta) {
    HalfVectorShiftResult result;

    if (Frame::cosTheta(tangentSpaceMainWi) * Frame::cosTheta(tangentSpaceMainWo) < (Float) 0) {
        // Refraction.

        // Refuse to shift if one of the Etas is exactly 1. This causes degenerate half-vectors.
        if (mainEta == (Float) 1 || shiftedEta == (Float) 1) {
            // This could be trivially handled as a special case if ever needed.
            result.success = false;
            return result;
        }

        // Get the non-normalized half vector.
        Vector3 tangentSpaceHalfVectorNonNormalizedMain;
        if (Frame::cosTheta(tangentSpaceMainWi) < (Float) 0) {
            tangentSpaceHalfVectorNonNormalizedMain = -(tangentSpaceMainWi * mainEta + tangentSpaceMainWo);
        } else {
            tangentSpaceHalfVectorNonNormalizedMain = -(tangentSpaceMainWi + tangentSpaceMainWo * mainEta);
        }

        // Get the normalized half vector.
        Vector3 tangentSpaceHalfVector = normalize(tangentSpaceHalfVectorNonNormalizedMain);

        // Refract to get the outgoing direction.
        Vector3 tangentSpaceShiftedWo = refract(tangentSpaceShiftedWi, tangentSpaceHalfVector, shiftedEta);

        // Refuse to shift between transmission and full internal reflection.
        // This shift would not be invertible: reflections always shift to other reflections.
        if (tangentSpaceShiftedWo.isZero()) {
            result.success = false;
            return result;
        }

        // Calculate the Jacobian.
        Vector3 tangentSpaceHalfVectorNonNormalizedShifted;
        if (Frame::cosTheta(tangentSpaceShiftedWi) < (Float) 0) {
            tangentSpaceHalfVectorNonNormalizedShifted = -(tangentSpaceShiftedWi * shiftedEta + tangentSpaceShiftedWo);
        } else {
            tangentSpaceHalfVectorNonNormalizedShifted = -(tangentSpaceShiftedWi + tangentSpaceShiftedWo * shiftedEta);
        }

        Float hLengthSquared = tangentSpaceHalfVectorNonNormalizedShifted.lengthSquared() /
                               (D_EPSILON + tangentSpaceHalfVectorNonNormalizedMain.lengthSquared());
        Float WoDotH = abs(dot(tangentSpaceMainWo, tangentSpaceHalfVector)) /
                       (D_EPSILON + abs(dot(tangentSpaceShiftedWo, tangentSpaceHalfVector)));

        // Output results.
        result.success = true;
        result.wo = tangentSpaceShiftedWo;
        result.jacobian = hLengthSquared * WoDotH;
    } else {
        // Reflection.
        Vector3 tangentSpaceHalfVector = normalize(tangentSpaceMainWi + tangentSpaceMainWo);
        Vector3 tangentSpaceShiftedWo = reflect(tangentSpaceShiftedWi, tangentSpaceHalfVector);

        Float WoDotH =
                dot(tangentSpaceShiftedWo, tangentSpaceHalfVector) / dot(tangentSpaceMainWo, tangentSpaceHalfVector);
        Float jacobian = abs(WoDotH);

        result.success = true;
        result.wo = tangentSpaceShiftedWo;
        result.jacobian = jacobian;
    }

    return result;
}


/// Result of a reconnection shift.
struct ReconnectionShiftResult {
    bool success;   ///< Whether the shift succeeded.
    Float jacobian; ///< Local Jacobian determinant of the shift.
    Vector3 wo;     ///< World space outgoing vector for the shift.
};

/// Tries to connect the offset path to a specific vertex of the main path.
ReconnectionShiftResult reconnectShift(const Scene *scene,
                                       Point3 mainSourceVertex,
                                       Point3 targetVertex,
                                       Point3 shiftSourceVertex,
                                       Vector3 targetNormal, // FIXME: FIXME FIXME
                                       Float time, bool isVolume = false) {
    ReconnectionShiftResult result;

    // Check visibility of the connection.
    if (!testVisibility(scene, shiftSourceVertex, targetVertex, time)) {
        // Since this is not a light sample, we cannot allow shifts through occlusion.
        result.success = false;
        return result;
    }

    // Calculate the Jacobian.
    Vector3 mainEdge = mainSourceVertex - targetVertex;
    Vector3 shiftedEdge = shiftSourceVertex - targetVertex;

    Float mainEdgeLengthSquared = mainEdge.lengthSquared();
    Float shiftedEdgeLengthSquared = shiftedEdge.lengthSquared();

    Vector3 shiftedWo = -shiftedEdge / sqrt(shiftedEdgeLengthSquared);

    Float mainOpposingCosine = 1.f;
    Float shiftedOpposingCosine = 1.f;

    // Compute the shift connection if needed
    if (!isVolume) {
        mainOpposingCosine = dot(mainEdge, targetNormal) / sqrt(mainEdgeLengthSquared);
        shiftedOpposingCosine = dot(shiftedWo, targetNormal);
    }

    Float jacobian = std::abs(shiftedOpposingCosine * mainEdgeLengthSquared)
                     / (D_EPSILON + std::abs(mainOpposingCosine * shiftedEdgeLengthSquared));

    // Return the results.
    result.success = true;
    result.jacobian = jacobian;
    result.wo = shiftedWo;
    return result;
}

/// Tries to connect the offset path to a the environment emitter.
ReconnectionShiftResult environmentShift(const Scene *scene, const Ray &mainRay, Point3 shiftSourceVertex) {
    const Emitter *env = scene->getEnvironmentEmitter();

    ReconnectionShiftResult result;

    // Check visibility of the environment.
    if (!testEnvironmentVisibility(scene, mainRay, shiftSourceVertex)) {
        // Sampled by BSDF so cannot accept occlusion.
        result.success = false;
        return result;
    }

    // Return the results.
    result.success = true;
    result.jacobian = Float(1);
    result.wo = mainRay.d;

    return result;
}


/// Stores the results of a BSDF sample.
/// Do not confuse with Mitsuba's BSDFSamplingRecord.
struct BSDFSampleResult {
    BSDFSamplingRecord bRec;        ///< The corresponding BSDF sampling record.
    Spectrum weight;///< BSDF weight of the sampled direction.
    Float pdf;                ///< PDF of the BSDF sample.
    Point2 sample;    ///< Random number used
};

struct PhaseSampleResult {
    PhaseFunctionSamplingRecord pRec;
    Float weight;
    Float pdf;
};

class ShiftMapping {
public:
    ShiftMapping(const GradientPathTracerConfig *config) : m_config(config) {};

    virtual ~ShiftMapping() = default;

    virtual void evaluate(RayState &main, RayState *shiftedRays, int secondaryCount, Spectrum &out_veryDirect) = 0;


    struct GradientInfo {
        Spectrum value;
        const bool x_axis;
        Vector2 offset;
        int pixel_a_index; // pixel_a_index - pixel_b_index
        int pixel_b_index;
    };

    virtual void evaluateReuse(RayState *rays, int secondaryCount,
                               Spectrum &out_veryDirect, int id_main,
                               std::vector<GradientInfo> &gradients) = 0;

protected:
    inline BSDFSampleResult sampleBSDF(RayState &rayState) {
        Intersection &its = rayState.rRec.its;
        RadianceQueryRecord &rRec = rayState.rRec;
        RayDifferential &ray = rayState.ray;

        // Note: If the base path's BSDF evaluation uses random numbers,
        // it would be beneficial to use the same random numbers for the offset path's BSDF.
        // This is not done currently.

        const BSDF *bsdf = its.getBSDF(ray);

        // Sample BSDF * cos(theta).
        BSDFSampleResult result = {
                BSDFSamplingRecord(its, rRec.sampler, ERadiance),
                Spectrum(),
                (Float) 0,
                Point2(0.0),
        };

        result.sample = rRec.nextSample2D();
        result.weight = bsdf->sample(result.bRec, result.pdf, result.sample);

        // Variable result.pdf will be 0 if the BSDF sampler failed to produce a valid direction.
        SAssert(result.pdf <= (Float) 0 || fabs(result.bRec.wo.length() - 1.0) < 0.00001);
        return result;
    }

    inline PhaseSampleResult samplePhase(RayState &rayState, const MediumSamplingRecord &mRec) {
        RadianceQueryRecord &rRec = rayState.rRec;
        RayDifferential &ray = rayState.ray;

        // Sample Phase.
        PhaseSampleResult result = {
                PhaseFunctionSamplingRecord(mRec, -ray.d, ERadiance),
                (Float) 0,
                (Float) 0
        };
        result.weight = mRec.getPhaseFunction()->sample(result.pRec, result.pdf, rRec.sampler);

        return result;
    }

    /**
     * This code is responsible of shifting light path
     * using direct reconnection
     */
    struct ShiftRecordLight {
        Float weightDenominator = Float(0.0);
        Spectrum value = Spectrum(0.0);
        bool succeed = false;
    };

    struct MainLightSamplingInfo {
        bool visible = false;
        DirectSamplingRecord dRec;
        Spectrum value = Spectrum(0.0);
        Point2 sample = Point2(0.0);
        Float distSquare = 0.0;
        Float opCos = 0.0;

        MainLightSamplingInfo(const Intersection &its) : dRec(its) {}
    };

    ShiftRecordLight shiftLightSampling(RayState &shifted, /*const*/ RayState &main,
                                        const MainLightSamplingInfo &mainInfo,
                                        const Spectrum &mainBSDFValue, Float mainBsdfPdf, const BSDF *mainBSDF) {
        const Emitter *emitter = static_cast<const Emitter *>(mainInfo.dRec.object);
        SAssert(emitter != nullptr);

        if (!shifted.alive) {
            return ShiftRecordLight{};
        }
        if (shifted.connection_status == RAY_CONNECTED) {
            ShiftRecordLight result;
            result.value = shifted.throughput * (mainBSDFValue * mainInfo.value);
            result.weightDenominator = shifted.pdf * (mainInfo.dRec.pdf + mainBsdfPdf);
            result.succeed = true;
            return result;
        } else if (shifted.connection_status == RAY_RECENTLY_CONNECTED) {
            // Evaluate the BSDF as the incomming direction is different
            Vector3 incomingDirection = normalize(shifted.rRec.its.p - main.rRec.its.p);
            BSDFSamplingRecord bRec(main.rRec.its, main.rRec.its.toLocal(incomingDirection),
                                    main.rRec.its.toLocal(mainInfo.dRec.d), ERadiance);
            Float shiftedBsdfPdf = (emitter->isOnSurface()
                                    && mainInfo.dRec.measure == ESolidAngle && mainInfo.visible) ? mainBSDF->pdf(bRec)
                                                                                                 : 0;
            Spectrum shiftedBsdfValue = mainBSDF->eval(bRec);
            // Return the info
            ShiftRecordLight result;
            result.value = shifted.throughput * (shiftedBsdfValue * mainInfo.value);
            result.weightDenominator = shifted.pdf * (mainInfo.dRec.pdf + shiftedBsdfPdf);
            result.succeed = true;
            return result;
        } else {
            SAssert(shifted.connection_status == RAY_NOT_CONNECTED);
            const BSDF *shiftedBSDF = shifted.rRec.its.getBSDF(shifted.ray);
            // Check the light is not discrete
            // and check that there is a least one component smooth
            VertexType mainVertexType = getVertexType(main, *m_config, BSDF::ESmooth);
            VertexType shiftedVertexType = getVertexType(shifted, *m_config, BSDF::ESmooth);
            bool mainAtPointLight = (mainInfo.dRec.measure == EDiscrete);
            if (!(mainAtPointLight || (mainVertexType == VERTEX_TYPE_DIFFUSE &&
                                       shiftedVertexType == VERTEX_TYPE_DIFFUSE))) {
                // Impossible in this case to do the shift
                // just terminate it
                return ShiftRecordLight{};
            }

            // Sample the light
            DirectSamplingRecord shiftedDRec(shifted.rRec.its);
            std::pair<Spectrum, bool> emitterTuple = main.rRec.scene->sampleEmitterDirectVisible(
                    shiftedDRec, mainInfo.sample);
            bool shiftedEmitterVisible = emitterTuple.second;
            Spectrum shiftedEmitterRadiance = emitterTuple.first * shiftedDRec.pdf;
            Float shiftedDRecPdf = shiftedDRec.pdf;

            // Sample the BSDF.
            Float shiftedDistanceSquared = (mainInfo.dRec.p - shifted.rRec.its.p).lengthSquared();
            Vector emitterDirection =
                    (mainInfo.dRec.p - shifted.rRec.its.p) / sqrt(shiftedDistanceSquared);
            Float shiftedOpposingCosine = -dot(mainInfo.dRec.n, emitterDirection);
            BSDFSamplingRecord bRec(shifted.rRec.its,
                                    shifted.rRec.its.toLocal(emitterDirection), ERadiance);
            if (m_config->m_strictNormals &&
                dot(shifted.rRec.its.geoFrame.n, emitterDirection) *
                Frame::cosTheta(bRec.wo) < 0) {
                return ShiftRecordLight{};
            }

            Spectrum shiftedBsdfValue = shiftedBSDF->eval(bRec);
            Float shiftedBsdfPdf = (emitter->isOnSurface() &&
                                    mainInfo.dRec.measure == ESolidAngle &&
                                    shiftedEmitterVisible) ? shiftedBSDF->pdf(bRec) : 0;
            // TODO: Epsilon can generate error in the jacobian
            Float jacobian = std::abs(shiftedOpposingCosine * mainInfo.distSquare) /
                             (Epsilon + std::abs(mainInfo.opCos * shiftedDistanceSquared));

            ShiftRecordLight result;
            result.value = jacobian * shifted.throughput * (shiftedBsdfValue * shiftedEmitterRadiance);
            result.weightDenominator = (jacobian * shifted.pdf) * (shiftedDRecPdf + shiftedBsdfPdf);
            return result;
        }
    }

    /**
     * This code is reponsible for the shift mapping
     * over the surfaces
     */
    struct ShiftRecordBSDF {
        Float weightDenominator = Float(0.0);
        Spectrum value = Spectrum(0.0);
        bool succeed = false;
        bool exclude_light_sampling = false;
    };

    struct MainBSDFSamplingInfo {
        Intersection previous;
        // --- Emitter
        bool hitEmitter = false;
        Spectrum light_value = Spectrum(0.0);
        Float light_pdf = Float(0.0);
        DirectSamplingRecord dRec;
        // --- Light sampling
        Spectrum bsdf_value = Spectrum(0.0);
        Float bsdf_pdf = Float(0.0);
        Float eta = Float(1.0);
        // --- BSDF classification
        VertexType current = VERTEX_TYPE_DIFFUSE;
        VertexType next = VERTEX_TYPE_DIFFUSE;
        unsigned int sampledType = 0;

        MainBSDFSamplingInfo(const Intersection its) : previous(its), dRec(its) {}
    };

    ShiftRecordBSDF shiftBSDFSampling(RayState &shifted, RayState &main,
                                      const MainBSDFSamplingInfo &mainInfo, const BSDF *mainBSDF, Point2 rndSample) {
        if (!shifted.alive) {
            return ShiftRecordBSDF{};
        }
        Float shiftedPreviousPdf = shifted.pdf;
        if (shifted.connection_status == RAY_CONNECTED) {
            // Update throughput and pdf.
            shifted.throughput *= mainInfo.bsdf_value;
            shifted.pdf *= mainInfo.bsdf_pdf;
            shifted.eta *= mainInfo.eta;
            ShiftRecordBSDF result;
            result.weightDenominator = shiftedPreviousPdf * (mainInfo.bsdf_pdf + mainInfo.light_pdf);
            result.value = shifted.throughput * mainInfo.light_value;
            result.succeed = true;
            return result;
        } else if (shifted.connection_status == RAY_RECENTLY_CONNECTED) {
            // The direct have change, so we need to evaluate the BSDF
            // at the main path
            Vector3 incomingDirection = normalize(shifted.rRec.its.p - main.ray.o);
            BSDFSamplingRecord bRec(mainInfo.previous, mainInfo.previous.toLocal(incomingDirection),
                                    mainInfo.previous.toLocal(main.ray.d), ERadiance);
            // TODO: If it is discrete, I guess the path is dead
            EMeasure measure = (mainInfo.sampledType & BSDF::EDelta) ? EDiscrete : ESolidAngle;
            Spectrum shiftedBsdfValue = mainBSDF->eval(bRec, measure);
            Float shiftedBsdfPdf = mainBSDF->pdf(bRec, measure);
            // Update
            shifted.throughput *= shiftedBsdfValue;
            shifted.pdf *= shiftedBsdfPdf;
            shifted.eta *= mainInfo.eta;
            shifted.connection_status = RAY_CONNECTED;
            ShiftRecordBSDF result;
            result.weightDenominator = shiftedPreviousPdf * (shiftedBsdfPdf + mainInfo.light_pdf);
            result.value = shifted.throughput * mainInfo.light_value;
            result.succeed = true;
            return result;
        } else {
            SAssert(shifted.connection_status == RAY_NOT_CONNECTED);
            const BSDF *shiftedBSDF = shifted.rRec.its.getBSDF(shifted.ray);
            // Base on the decision of the main path...
            VertexType shiftedVertexType = getVertexType(shifted, *m_config, mainInfo.sampledType);
            if (mainInfo.current == VERTEX_TYPE_DIFFUSE && mainInfo.next == VERTEX_TYPE_DIFFUSE &&
                shiftedVertexType == VERTEX_TYPE_DIFFUSE) {
                // We are in the configuration of doing diffuse reconnection
                auto shiftResult = [&]() -> ReconnectionShiftResult {
                    if (main.rRec.its.isValid()) {
                        // This is an actual reconnection shift.
                        return reconnectShift(main.rRec.scene, main.ray.o, main.rRec.its.p,
                                              shifted.rRec.its.p, main.rRec.its.geoFrame.n,
                                              main.ray.time);
                    } else {
                        // This is a reconnection at infinity in environment direction.
                        const Emitter *env = main.rRec.scene->getEnvironmentEmitter();
                        SAssert(env != NULL);
                        return environmentShift(main.rRec.scene, main.ray, shifted.rRec.its.p);
                    }
                }();

                if (!shiftResult.success) {
                    shifted.alive = false;
                    return ShiftRecordBSDF{};
                }

                Vector3 incomingDirection = -shifted.ray.d;
                Vector3 outgoingDirection = shiftResult.wo;
                BSDFSamplingRecord bRec(shifted.rRec.its, shifted.rRec.its.toLocal(incomingDirection),
                                        shifted.rRec.its.toLocal(outgoingDirection), ERadiance);
                // Strict normals check.
                if (m_config->m_strictNormals &&
                    dot(outgoingDirection, shifted.rRec.its.geoFrame.n)
                    * Frame::cosTheta(bRec.wo) <= 0) {
                    shifted.alive = false;
                    return ShiftRecordBSDF{};
                }
                // Evaluate the BRDF to the new direction.
                Spectrum shiftedBsdfValue = shiftedBSDF->eval(bRec);
                Float shiftedBsdfPdf = shiftedBSDF->pdf(bRec);

                // Update throughput and pdf.
                shifted.throughput *= shiftedBsdfValue * shiftResult.jacobian;
                shifted.pdf *= shiftedBsdfPdf * shiftResult.jacobian;
                //shifted.eta *= bRec.eta;
                // ignore eta here
                shifted.connection_status = RAY_RECENTLY_CONNECTED;

                if (!(mainInfo.hitEmitter || main.rRec.its.hasSubsurface())) {
                    // We do not need to gather the emitted radiance,
                    // just return the shift in this case
                    return ShiftRecordBSDF{};
                }
                struct ShiftEmitterIntersection {
                    Spectrum value = Spectrum(0.0);
                    Float pdf = Float(0.0);
                };
                auto shift_emitter = [&]() -> ShiftEmitterIntersection {
                    if (main.rRec.its.isValid()) {
                        // Hit an object.
                        ShiftEmitterIntersection result;
                        if (mainInfo.hitEmitter) {
                            result.value = main.rRec.its.Le(-outgoingDirection);

                            // Evaluate the light sampling PDF of the new segment.
                            DirectSamplingRecord shiftedDRec;
                            shiftedDRec.p = mainInfo.dRec.p;
                            shiftedDRec.n = mainInfo.dRec.n;
                            shiftedDRec.dist = (mainInfo.dRec.p - shifted.rRec.its.p).length();
                            shiftedDRec.d = (mainInfo.dRec.p - shifted.rRec.its.p) / shiftedDRec.dist;
                            shiftedDRec.ref = mainInfo.dRec.ref;
                            shiftedDRec.refN = shifted.rRec.its.shFrame.n;
                            shiftedDRec.object = mainInfo.dRec.object;
                            result.pdf = main.rRec.scene->pdfEmitterDirect(shiftedDRec);
                        }

                        // Sub-surface scattering. Note: Should use the same random numbers as the base path!
                        if (main.rRec.its.hasSubsurface() &&
                            (main.rRec.type & RadianceQueryRecord::ESubsurfaceRadiance)) {
                            SLog(EError, "No support for subsurface scattering");
                        }
                        return result;
                    } else {
                        // Hit the environment.
                        ShiftEmitterIntersection result;
                        result.value = mainInfo.light_value;
                        result.pdf = mainInfo.light_pdf;
                        return result;
                    }
                }();

                ShiftRecordBSDF result;
                result.value = shifted.throughput * shift_emitter.value;
                result.weightDenominator = shiftedPreviousPdf * (shiftedBsdfPdf + shift_emitter.pdf);
                result.succeed = true;
                return result;
            } else {
                // Continue to bounce over the surfaces
                ShiftRecordBSDF result;
                result.exclude_light_sampling = true;

                // Condition on the shift over the surface
                // Just to ensure a consistent classification
                bool bothDelta = (mainInfo.sampledType & BSDF::EDelta) &&
                                 (shiftedBSDF->getType() & BSDF::EDelta);
                bool bothSmooth = (mainInfo.sampledType & BSDF::ESmooth) &&
                                  (shiftedBSDF->getType() & BSDF::ESmooth);
                if (!(bothDelta || bothSmooth)) {
                    shifted.alive = false;
                    return result;
                }
                // Double check that the ray direction is correct
                // If can happens that the ray is wrong. Need to stop in this case.
                if (fabs(shifted.ray.d.lengthSquared() - 1) > 0.01) {
                    SLog(EError, "Bad length?");
                    return result;
                }

                // TODO: The shift mapping need to be agonstic to the shift mapping
                // TODO: Because this code is only for random replay.
                // TODO: Need to change the code
                BSDFSamplingRecord shiftedBRec(shifted.rRec.its, nullptr, ERadiance);
                Float shiftBsdfPdf = Float(0.f);
                Spectrum shiftBsdfValue = shiftedBSDF->sample(shiftedBRec, shiftBsdfPdf, rndSample);
                if (shiftBsdfValue.isZero() || shiftBsdfPdf == 0) {
                    shifted.alive = false;
                    return result;
                }
                Float jacobian = mainInfo.bsdf_pdf / shiftBsdfPdf;
                shiftBsdfValue *= shiftBsdfPdf;
                Vector outgoingDirection = shifted.rRec.its.toWorld(shiftedBRec.wo);
                auto measure = (mainInfo.sampledType & BSDF::EDelta) ? EDiscrete : ESolidAngle;
                if (m_config->m_strictNormals &&
                    dot(outgoingDirection, shifted.rRec.its.geoFrame.n)
                    * Frame::cosTheta(shiftedBRec.wo) <= 0) {
                    shifted.alive = false;
                    return result;
                }

                // FIXME: This line look dubious for me...
                VertexType shiftedVertexType = getVertexType(shifted, *m_config, mainInfo.sampledType);
                Spectrum shiftedEmitterRadiance = Spectrum(0.0);
                shifted.ray = Ray(shifted.rRec.its.p, outgoingDirection, main.ray.time);
                if (!main.rRec.scene->rayIntersect(shifted.ray, shifted.rRec.its)) {
                    // Hit nothing - Evaluate environment radiance.
                    const Emitter *env = main.rRec.scene->getEnvironmentEmitter();
                    if (!env) {
                        shifted.alive = false;
                        return result;
                    }
                    SLog(EError, "Not implemented as we do not support postponed Shift end");
                } else {
                    if (!main.rRec.its.isValid()) {
                        // Deny shifts between env and non-env.
                        shifted.alive = false;
                        return result;
                    }

                    VertexType shiftedNextVertexType = getVertexType(shifted, *m_config,
                                                                     mainInfo.sampledType);
                    if (mainInfo.current == VERTEX_TYPE_DIFFUSE &&
                        shiftedVertexType == VERTEX_TYPE_DIFFUSE &&
                        shiftedNextVertexType == VERTEX_TYPE_DIFFUSE) {
                        // Non-invertible shift: the reverse-shift would use another strategy!
                        shifted.alive = false;
                        return result;
                    }

                    if (shifted.rRec.its.isEmitter()) {
                        // Hit emitter.
                        shiftedEmitterRadiance = shifted.rRec.its.Le(-shifted.ray.d);
                    }
                    // Sub-surface scattering. Note: Should use the same random numbers as the base path!
                    if (shifted.rRec.its.hasSubsurface() &&
                        (shifted.rRec.type & RadianceQueryRecord::ESubsurfaceRadiance)) {
                        SLog(EError, "No support for subscattering");
                    }
                }

                // Finally we can bake the result
                shifted.throughput *= shiftBsdfValue * jacobian;
                shifted.pdf *= shiftBsdfPdf * jacobian;
                // Fill the result
                result.value = shifted.throughput * shiftedEmitterRadiance;
                result.weightDenominator = shiftedPreviousPdf * (shiftBsdfPdf);
                result.succeed = true;
                return result;
            }
        }
    }

protected:
    const GradientPathTracerConfig *m_config;
};

MTS_NAMESPACE_END

#endif //MITSUBA_GPT_SHIFTMAPPING_H
