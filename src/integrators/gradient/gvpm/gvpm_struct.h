#pragma once

#include <algorithm>
#include <string>
#include <utility>

#include <mitsuba/render/gatherproc.h>
#include <mitsuba/bidir/mut_manifold.h>
#include <mitsuba/core/plugin.h>

#include "gvpm_geoOps.h"

// All common lighting effects
#include "../../volume_utils.h"

// To generate the some figure for the paper
// or debugging proposes.
// It does not be use in prod because it makes the code
// a little bit slower
#define HAVE_ADDITIONAL_STATS 0

MTS_NAMESPACE_BEGIN

/// Please don't be confused with ESolidAngle and EArea. 
enum EPdfMeasure {
  EPdfSolidAngle = 0,
  EPdfArea = 1
};

enum ELightShiftType {
  EAllShift = 0,
  EDiffuseShift = 1 << 1,
  EManifoldShift = 1 << 2,
  ENullShift = 1 << 4,
  EMediumShift = 1 << 5,
  EInvalidShift = 1 << 6
};

/// Classification of vertices into diffuse and glossy.
enum EVertexType {
  VERTEX_TYPE_GLOSSY,     ///< "Specular" vertex that requires the half-vector duplication shift.
  VERTEX_TYPE_DIFFUSE,     ///< "Non-specular" vertex that is rough enough for the reconnection shift.
  VERTEX_TYPE_INVALID
};

struct VertexClassifier {
  static Float roughnessThreshold;
  /**
   * \brief Classify a pathvertex into diffuse or specular for gradient tracing use.
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

  static EVertexType type(const PathVertex &vertex, int component) {
    switch (vertex.getType()) {
    case PathVertex::EEmitterSample:return VERTEX_TYPE_DIFFUSE;
    case PathVertex::ESurfaceInteraction:
      return type(vertex.getIntersection().getBSDF(),
                  vertex.getIntersection(),
                  component);
    case PathVertex::EMediumInteraction:
      // FIXME: Double check this code
      return vertex.getMediumSamplingRecord().getPhaseFunction()->getMeanCosine() > 0.5 ? VERTEX_TYPE_GLOSSY
                                                                                        : VERTEX_TYPE_DIFFUSE;
    default:return VERTEX_TYPE_INVALID;
    }
  }

  /// Returns the vertex type of a vertex by its roughness value.
private:
  static EVertexType getVertexTypeByRoughness(Float roughness) {
    if (roughness <= roughnessThreshold) {
      return VERTEX_TYPE_GLOSSY;
    } else {
      return VERTEX_TYPE_DIFFUSE;
    }
  }
};

/// Configuration for the gradient photon mapper.
struct GPMConfig {
  // Photon tracing options
  int maxDepth;
  int minDepth;
  int rrDepth;

  // SPPM parameters
  Float alpha; //< Reduction rate of density estimation kernel
  Float initialScale; //< Initial size for the density estimation kernel used for surface rendering
  Float initialScaleVolume; //< Initial size for the density estimation kernel used for volume rendering
  Float epsilonRadius; //< The minimal kernel size
  int photonCount; //<  Number of photons for surface
  int volumePhotonCount; //< Number of photons for volume
  int maxPasses; //< The number of iteration (at max) before stopping computation

  // Option for the gradient domain
  // to enable or disable certain type of options related to gradient usage (image reconstruction)
  bool estimateGradients;
  bool reconstructL1;
  bool reconstructL2;
  bool reconstructUni;
  bool reconstructWeighted;
  bool reconstructNfor;
  Float reconstructAlpha;
  Float m_shiftThreshold;

  // Additional display options
  int dumpIteration; //< At each iteration number the reconstruction need to be done?

  // Options for debugging.
  // This option is to enable/disable a certain part of the computation
  int lightingInteractionMode;
  int bsdfInteractionMode;
  bool useManifold;
  bool useMIS;
  ELightShiftType debugShift;
  bool noMediumShift; // Replace medium shift by diffuse shift

  // If we use photon beam, do we want to use long beams?
  bool convertLong; // FIXME: Check if this options works properly for photon beams

  // Extra option for the rendering
  int maxManifoldIterations;
  bool forceBlackPixels;
  bool fluxLightDensity;
  bool reusePrimal;
  bool directTracing;
  bool nearSpecularDirectImage;
  bool nanCheck;

  // In the reference mode, the sampler are initialized
  // randomly to be sure to not have the same
  // sequence of gather points generated
  bool referenceMod;

  // To choose the different optins for the volume rendering
  EVolumeTechnique volTechnique;

  // "Advance strategy for the shifts"
  bool useShiftNull;
  bool stratified;
  Float relaxME;  // 0 = ME default behavior

  // Sensor parameter
  int nbCameraSamples;
  size_t minCameraDepth;
  int maxCameraDepth;
  Float cameraSphere;

  // Other options.
  bool adjointCompensation;
  bool deterministic;
  bool pathSet;
  bool powerHeuristic;
  bool newShiftBeam;
  bool use3DKernelReduction;
  std::string forceAPA;

  // Same radius strategy as VCM
  Float initialRadius;

  void load(const Properties &props) {
    // Use the old shift by default
    // this old shift have some problem to keep contrains on the photon beams
    newShiftBeam = props.getBoolean("newShiftBeam", false); // Ok for 3D beams


    initialScale = props.getFloat("initialScale", 1.f);
    initialScaleVolume = props.getFloat("initialScaleVolume", 1.f);
    use3DKernelReduction = props.getBoolean("use3DKernelReduction", false);
    forceAPA = props.getString("forceAPA", "");
    noMediumShift = props.getBoolean("noMediumShift", true);
    powerHeuristic = props.getBoolean("powerHeuristic", false);

    // To exclude single scattering that come from a specular chain
    relaxME = props.getFloat("relaxME", 1.f);
    if (relaxME < 0.0) {
      SLog(EError, "relaxME options need to be 0 or 1.");
    }

    /* The minimum non-zero radius that can be used to start density estimation */
    epsilonRadius = props.getFloat("epsilonRadius", 1e-3f);
    /* Alpha parameter from the paper (influences the speed, at which the photon radius is reduced) */
    alpha = props.getFloat("alpha", .7);
    /* Number of photons to shoot in each iteration */
    photonCount = props.getInteger("photonCount", 250000);
    volumePhotonCount = props.getInteger("volumePhotonCount", 250000);

    /* Longest visualized path length (<tt>-1</tt>=infinite). When a positive value is
    specified, it must be greater or equal to <tt>2</tt>, which corresponds to single-bounce
    (direct-only) illumination */
    maxDepth = props.getInteger("maxDepth", -1);
    minDepth = props.getInteger("minDepth", 0);
    /* Depth to start using russian roulette */
    rrDepth = props.getInteger("rrDepth", 12);
    /* Maximum number of passes to render. -1 renders until the process is stopped. */
    maxPasses = props.getInteger("maxPasses", -1);

    if (maxDepth <= 1 && maxDepth != -1)
      SLog(EError, "Maximum depth must be set to \"2\" or higher!");
    if (maxPasses <= 0 && maxPasses != -1)
      SLog(EError, "Maximum number of passes must either be set to \"-1\" or \"1\" or higher!");

    // Dump a hdr file every N frames
    dumpIteration = props.getInteger("dumpIteration", 1);

    // Reconstruction parameters
    estimateGradients = props.getBoolean("estimateGradients", true);
    reconstructL1 = props.getBoolean("reconstructL1", false);
    reconstructL2 = props.getBoolean("reconstructL2", true);
    reconstructUni = props.getBoolean("reconstructUni", false);
    reconstructWeighted = props.getBoolean("reconstructWeighted", false);
    reconstructNfor = props.getBoolean("reconstructNfor", false);
    reconstructAlpha = props.getFloat("reconstructAlpha", Float(0.2));

    // Vertex classification
    m_shiftThreshold = props.getFloat("shiftThresholdGPM", 0.001);
    if (m_shiftThreshold <= 0.0 || m_shiftThreshold > 1.0) {
      SLog(EError, "Bad roughtness constant: %f", m_shiftThreshold);
    }
    VertexClassifier::roughnessThreshold = m_shiftThreshold;

    // Option for testing our algorithm
    bsdfInteractionMode = parseBsdfInteractionMode(props.getString("bsdfInteractionMode", "all"));
    useManifold = props.getBoolean("useManifold", true);
    initialRadius = props.getFloat("initialRadius", 0.0);

    // MIS selection
    useMIS = [&props]() -> bool {
      std::string mis = props.getString("useMIS", "area");
      std::transform(mis.begin(), mis.end(),
                     mis.begin(), ::tolower);
      if (mis == "area") {
        return true;
      } else if(mis == "none") {
        return false;
      } else {
        SLog(EError, "useMIS: need to be 'none' or 'area'");
        return false;
      }
    }();

    // Shift debugging: To only look at one shift results
    debugShift = [&props]() -> ELightShiftType {
      std::string debugShiftStr = props.getString("debugShift", "0");
      std::transform(debugShiftStr.begin(), debugShiftStr.end(),
                     debugShiftStr.begin(), ::tolower);
      if (debugShiftStr == "me")
        return EManifoldShift;
      else if (debugShiftStr == "diffuse")
        return EDiffuseShift;
      else if (debugShiftStr == "null")
        return ENullShift;
      else
        return EAllShift;
    }();

    maxManifoldIterations = props.getInteger("maxManifoldIterations", 5);
    forceBlackPixels = props.getBoolean("forceBlackPixels", false);
    fluxLightDensity = props.getBoolean("fluxLightDensity", true);
    reusePrimal = props.getBoolean("reusePrimal", true);

    /* Skip gradient computation for pure specular path. Final image = direct image + reconstruction image. */
    directTracing = props.getBoolean("directTracing", false);
    /* Consider near-specular path as pure specular and skip computing gradients for these. */
    nearSpecularDirectImage = props.getBoolean("nearSpecularDirectImage", false);
    nanCheck = props.getBoolean("nanCheck", false);

    referenceMod = props.getBoolean("referenceMod", false);

    /* Option of rendering to avoid certain type of computation */
    std::string strLightingMode = props.getString("lightingInteractionMode", "");
    if (strLightingMode.empty())
      strLightingMode = props.getString("interactionMode", "all2all");
    lightingInteractionMode = (int) parseMediaInteractionMode(strLightingMode);

    // Option on sampled dist
    convertLong = props.getBoolean("convertLong", false);

    // The different options for the volume rendering
    std::string volRenderingTech = props.getString("volTechnique", "distance");
    volTechnique = parseVolumeTechnique(volRenderingTech);
    if (!needSurfaceRendering()) {
      photonCount = 0; // Force to 0 if not needed
    }
    if (!needVolumeRendering()) {
      volumePhotonCount = 0; // Force to 0 if not needed
    }

    // Advance strategies for the shift
    useShiftNull = props.getBoolean("useShiftNull", false);

    // For the camera sample strategies
    stratified = props.getBoolean("stratified", false);
    nbCameraSamples = props.getInteger("nbCameraSamples", 40);
    if (useShiftNull && !EVolumeTechniqueHelper::use3DKernel(volTechnique)
        && !(volTechnique == EBeamBeam1D)) {
      SLog(EError, "Not possible to shift null without using 3D kernel");
    }

    // Camera debug
    minCameraDepth = [&props]() -> size_t {
      auto v = props.getInteger("minCameraDepth", 0);
      if(v < 0) { SLog(EError, "minCamera depth need to be null or positive"); }
      return v;
    }();
    maxCameraDepth = props.getInteger("maxCameraDepth", -1);

    adjointCompensation = true; // Make too much variance.
    deterministic = props.getBoolean("deterministic", false);
    cameraSphere = props.getFloat("cameraSphere", 1.f);

    // Do a global RR
    pathSet = props.getBoolean("pathSet", true);

    if (pathSet && deterministic) {
      SLog(EError, "It is not possible to use pathSet and deterministic option at the same time");
    }
  }

  bool needVolumeRendering() const {
    return (lightingInteractionMode & ESurf2Media) || (lightingInteractionMode & EMedia2Media);
  }
  bool needSurfaceRendering() const {
    return (lightingInteractionMode & ESurf2Surf) || (lightingInteractionMode & EMedia2Surf);
  }
  bool isAPAVolumeEstimator() const {
    return volTechnique == EVolBRE2D || volTechnique == EVolBRE3D ||
        volTechnique == EBeamBeam1D ||
        volTechnique == EBeamBeam3D_Naive ||
        volTechnique == EBeamBeam3D_EGSR ||
        volTechnique == EBeamBeam3D_Optimized ||
        volTechnique == EVolPlane0D;
  }
  bool isPlaneEstimator() const {
    return volTechnique == EVolPlane0D;
  }
};

enum EPixel {
  ELeft = 0,
  ERight = 1,
  ETop = 2,
  EBottom = 3
};

struct SVertexPDF {
  Float pdf = Float(0.f);       // pdf of sampling this path using Monte Carlo, in area measure, for MIS use
  Float jacobian = Float(0.f);  // for adjusting the pdf of the base path due to shift mapping

  Spectrum weight = Spectrum(Float(0.f));  // for a vertex k, store accumulated vertex and edge weight up to vertex (k-1)
  // with Jacobian premultiplied

  Spectrum vertexWeight = Spectrum(Float(0.f));        // for a vertex k, store the BSDF * the pdf that samples the next vertex (k+1)
  // with Jacobian premultiplied
};

#if HAVE_ADDITIONAL_STATS
struct ShiftStats {
  // Light shifts
  size_t nbLightShifts = 0;
  size_t nbSuccessLightShifts = 0;

  // Type of light shifts
  size_t MEShifts = 0;
  size_t DiffuseShifts = 0;
  size_t nullShifts = 0;

  // Just make possible to accumulate the statistic
  void add(const ShiftStats& o) {
    nbLightShifts += o.nbLightShifts;
    nbSuccessLightShifts += o.nbSuccessLightShifts;
    MEShifts += o.MEShifts;
    DiffuseShifts += o.DiffuseShifts;
    nullShifts += o.nullShifts;
  }
};
#endif

struct ShootingThreadData {
  MemoryPool pool;
  ref<Sampler> sampler;

  std::vector<Path *> paths;
  size_t nbPaths = 0;

  void clear() {
    for (Path *p: paths) {
      p->release(pool);
    }
    nbPaths = 0;
  }

  Path *getPath() {
    if (nbPaths < paths.size()) {
      nbPaths += 1;
      return paths[nbPaths - 1];
    } else {
      paths.push_back(new Path);
      nbPaths += 1;
      return paths.back();
    }
  }
};

/// Represents one individual PPM gather point including relevant statistics
class GatherPoint {
public:
  // Radiance information
  Spectrum emission = Spectrum(0.f);
  // total flux of this gather point (eye path weight multiplied)
  Spectrum baseFlux = Spectrum(0.f);
  bool pureSpecular = true;
  Spectrum directEmission = Spectrum(0.f);
  Spectrum mediumFlux = Spectrum(0.f);
  // Light paths where gather point is on emitter is handled separately
  Spectrum currEmission = Spectrum(0.f);

  Spectrum shiftedFlux[4];          // MIS weighted, total flux of the shifted path (eye path weight multiplied)
  Spectrum weightedBaseFlux[4];     // MIS weighted, total flux of the base path (eye path weight multiplied)
  Spectrum shiftedEmitterFlux[4];
  Spectrum weightedEmitterFlux[4];
  Spectrum shiftedMediumFlux[4];
  Spectrum weightedMediumFlux[4];

  // Auxiliary buffers for denoising 
  Spectrum albedo = Spectrum(0.); 
  Vector position = Vector(0.);
  Vector shadingNormal = Vector(0.);

  // Tracing history 
  Path path;

#if HAVE_ADDITIONAL_STATS
  ShiftStats stats;
#endif

  Float radius = 1.f;
  Float N = 0.f;
  Float NVol = 0.f;
  int depth = -1;
  Point2i pixel = Point2i(-1);

  // BSDF information
  // the component index of a multi-lobe BSDF
  int sampledComponent = -1;
  Float pdfComponent = 1.f;
  Float scale = 1.f;
  Float scaleVol = 1.f;
  bool haveSmoke = false;

private:
  std::vector<SVertexPDF> info;

public:
  inline GatherPoint() {
    info.clear();

    for (int i = 0; i < 4; ++i) {
      shiftedFlux[i] = Spectrum(0.);          
      weightedBaseFlux[i] = Spectrum(0.);     
      shiftedEmitterFlux[i] = Spectrum(0.);
      weightedEmitterFlux[i] = Spectrum(0.);
      shiftedMediumFlux[i] = Spectrum(0.);
      weightedMediumFlux[i] = Spectrum(0.);
    }
  }

  virtual ~GatherPoint() = default;

  void clear() {
    depth = -1;
    radius = 1;
    info.clear(); // Clear the current path
    currEmission = Spectrum(0.f);
    pureSpecular = true;
  }

  inline const PathVertex *last() const {
    if (path.vertexCount() < 2)
      SLog(EError, "Not possible to query a last vertex that have no surface interaction");
    return path.vertex(path.vertexCount() - 1);
  }
  inline const Intersection &lastIts() const {
    const PathVertex *l = last();
    if (!l->isSurfaceInteraction())
      SLog(EError, "Not possible to get its with last that is not surface intersection");
    return l->getIntersection();
  }

  inline const SVertexPDF &surfInfo() const {
    // FIXME: double check if this is still correct as vertex info structure is dramatically cleaned up
    // and changed by Feb 21, 2017
    if (!last()->isSurfaceInteraction())
      SLog(EError, "Not possible to get surfInfo if the last vertex is not on surface");
    if (info.empty())
      SLog(EError, "Not possible to get a data to an empty vector");
    return info.back();
  }

  inline SVertexPDF *createInfo() {
    SVertexPDF newInfo;
    if (!info.empty()) {
      newInfo = info.back(); // Copy the previous info
    }
    info.push_back(newInfo);
    return &info.back();
  }

  /**
   * For base path, this function is called to cache vertex info. 
   * For shift path, the cache vertex info is maintained manually as we do not actually store a full path. 
   */
  inline void generateVertexInfo() {
    if (path.vertexCount() <= 1)
      return; // Nothing to create
    SVertexPDF *currInfo = createInfo();

    // The cache at each vertex accounts for the weight up to the current vertex interaction
    // and the subsequent edge.

    currInfo->weight = path.vertex(0)->weight[ERadiance] *
        path.vertex(0)->rrWeight *
        path.edge(0)->weight[ERadiance];
    currInfo->pdf = 1.0f;
    currInfo->jacobian = 1.0f;
    currInfo->vertexWeight = path.vertex(0)->weight[ERadiance] *  path.vertex(0)->rrWeight;

    for (size_t i = 1; (path.vertexCount() - 1) > i; i++) {
      // Allocate new cache vertex
      currInfo = createInfo();

      // Update the weight per vertex
      currInfo->weight *= path.vertex(i)->weight[ERadiance] *
          path.vertex(i)->rrWeight *
          path.edge(i)->weight[ERadiance];
      currInfo->vertexWeight = path.vertex(i)->weight[ERadiance] *  path.vertex(i)->rrWeight;
      currInfo->pdf *= path.vertex(i)->pdf[ERadiance] * path.edge(i)->pdf[ERadiance];

      // Measure transformation
      // to put the pdf in the good measure...
      switch(path.vertex(i)->measure){
      case ELength:
      case ESolidAngle: {
        SLog(EError, "Did not expect to have such measure, abord");
        break;
      }
      case EArea: break; // Already in a good measure
      case EDiscrete: break; // FIXME: Maybe we need to use generalized measure
      default: {
        SLog(EError, "Not covered case of measure transformation. Abort.");
        break;
      }
      };
    }
  }

  inline const SVertexPDF &getVertexInfo(size_t vertexId) const {
    // The vertex ID is the same as the bidirectional layer path.
    if (info.size() <= vertexId)
      SLog(EError, "Try to get a value outside possible value");
    return info[vertexId];
  }

  inline SVertexPDF &getVertexInfo(size_t vertexId) {
    // The vertex ID is the same as the bidirectional layer path.
    if (info.size() <= vertexId)
      SLog(EError, "Try to get a value outside possible value");
    return info[vertexId];
  }

  /**
   * Compute the weight up to the current edge 
   * and the interaction at the last vertex of the beam. 
   */
  inline const Spectrum getWeightBeam(size_t idEdge) const {
    // The edge ID is the same as the bidirectional layer path.
    return getVertexInfo(idEdge).weight;
  }

  inline const Spectrum getWeightVertex(size_t idVertex) const {
    return getVertexInfo(idVertex).vertexWeight;
  }

  /**
   * Compute the Oposite geometry factor from the idVertex to idVertex+1
   * @param idVertex the prev. vertex id
   * @return Opposite geometry factor
   */
  inline Float GOp(size_t idVertex) const {
    return geometryOpposingTerm(path, idVertex, idVertex+1);
  }

  /**
   * For MIS usage:
   * return (pdfOffsetSensor / pdfBaseSensor)*jacobianSensor
   * Note that this doesn't return the Jacobian for the change of "t"
   */
  inline Float sensorMIS(size_t idVertex, const GatherPoint &base, Float sDist, Float bDist) const {
    const SVertexPDF &info = getVertexInfo(idVertex);
    const SVertexPDF &bInfo = base.getVertexInfo(idVertex);
    Float jacobian = info.jacobian;
    Float ratio = info.pdf / bInfo.pdf;
    if (idVertex != 1) {
      Float baseG = base.GOp(idVertex);
      Float currG = GOp(idVertex);

      // Cancel out the embedeed terms
      jacobian *= currG / baseG;
      ratio *= baseG / currG;

      jacobian *= (sDist / bDist) * (sDist / bDist);
      ratio *= (bDist / sDist) * (bDist / sDist);
    }
    if (ratio == 0.f) {
      SLog(EWarn, "PDF ratio for sensor 0");
    }
    if (jacobian == 0.f) {
      SLog(EWarn, "Jacobian is 0");
    }
    return ratio * jacobian;
  }

  inline void rescaleRadii(Float v) {
    baseFlux *= v;

    for (int i = 0; i < 4; ++i) {
      shiftedFlux[i] *= v;
      weightedBaseFlux[i] *= v;
    }
  }
};
typedef std::vector<std::vector<GatherPoint> > GatherBlocks;

struct GPMThreadData {
  ref<ManifoldPerturbation> offsetGenerator;
  MemoryPool pool;

  GPMThreadData(Scene *scene, const GPMConfig &config) {
    offsetGenerator = new ManifoldPerturbation(scene,
                                               nullptr,
                                               pool,
                                               0.f, true, true, 0, 0,
                                               config.m_shiftThreshold,
                                               config.maxManifoldIterations);
  }
};

class GPhoton {
public:
  Intersection its = {};
  const Medium *medium = {};
  Spectrum weight = {};
  int depth = {};
  int sampledComponent = {};
  int prevComponentType = {};     // default is zero, which means accepting all components

  GPhoton() = default;
  GPhoton(const Intersection &_its,
          const Medium *_medium,
          int _sampledComp,
          const Spectrum _w,
          int _d,
          int prevComponentType = 0) :
      its(_its), medium(_medium), weight(_w), depth(_d), sampledComponent(_sampledComp),
      prevComponentType(prevComponentType) {}

  static int getPrevComponentType(const PathVertex *prev) {
    if (prev->isEmitterSample()) {
      // As emitter component type is 0, treat it as diffuse
      return BSDF::EDiffuseReflection;
    } else if (prev->isMediumInteraction()) {
      return BSDF::EDiffuseReflection; //FIXME: Assume diffuse
    } else {
      return prev->componentType;
    }
  }
};

MTS_NAMESPACE_END
