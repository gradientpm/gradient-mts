#include <mitsuba/render/particleproc.h>
#include <mitsuba/bidir/mempool.h>
#include <mitsuba/bidir/path.h>
#include <mitsuba/core/kdtree.h>

#include <list>
#include <tuple>

#include "gvpm_struct.h"
#include "../../volume_utils.h"

#ifndef MITSUBA_GVPM_ACCEL_H
#define MITSUBA_GVPM_ACCEL_H

MTS_NAMESPACE_BEGIN

struct GPhotonNodeData {
  size_t vertexId;         // vertexId - 1 is the number of preceding interactions for both surface and volume.
  const Path *lightPath;
  Spectrum weight;
  unsigned int pathID;

  GPhotonNodeData() :
      vertexId(0), lightPath(nullptr), weight(0.f), pathID(0) {}

  GPhotonNodeData(const Path *_lt, int _vID,
                  Spectrum _weight, unsigned int _pathID) :
      vertexId(_vID), lightPath(_lt), weight(_weight), pathID(_pathID) {}

  //TODO: Optimize this !!!
  void getPhoton(GPhoton &photon) const {
    const size_t vID = vertexId;
    SAssert(vID >= 1);

    const Path *lt = lightPath;
    const PathVertex *v = lt->vertex(vID);

    // Record component type for interaction mode selection
    const PathVertex *prev = lt->vertex(vID - 1);
    const PathEdge *edge = lt->edge(vID - 1);

    size_t depth = vID - 1;
    if (v->type == PathVertex::ESurfaceInteraction) {
      photon = GPhoton(v->getIntersection(), edge->medium,
                       v->sampledComponentIndex, weight, depth, GPhoton::getPrevComponentType(prev));
    } else if (v->type == PathVertex::EMediumInteraction) {
      Intersection its;
      its.p = v->getPosition();
      its.wi = -edge->d; // Global orientation
      photon = GPhoton(its, edge->medium, v->sampledComponentIndex, weight, depth, GPhoton::getPrevComponentType(prev));
    } else {
      SLog(EError, "Bad vertex type !");
      // Default code ...
      Intersection its;
      its.p = v->getPosition();
      its.wi = -edge->d; // Global orientation
      photon = GPhoton(its, edge->medium, v->sampledComponentIndex, weight, depth, GPhoton::getPrevComponentType(prev));
    }
  }

//  const Point& getParentPosition() const {
//      SAssert(vertexId > 1);
//      return lightPath->vertex(vertexId - 1)->getPosition();
//  }
};

// === Acceleration structure
// We need to wrap the photon information inside the Kd-tree node
class GPhotonNodeKD : public SimpleKDNode<Point, GPhotonNodeData> {
public:
  GPhotonNodeKD() = default;

  /// Construct from a photon interaction
  explicit GPhotonNodeKD(GPhotonNodeData _data) :
      SimpleKDNode<Point, GPhotonNodeData>(_data) {

    // Setup the position of the current node
    // using the lightpath information
    GPhoton photon;
    _data.getPhoton(photon);
    this->setPosition(photon.its.p);
  }
};

class SurfaceGradientRecord;
class VolumeGradientRecord;
class GPhotonMap : public SerializableObject {
public:
  typedef PointKDTree<GPhotonNodeKD> PhotonTree;
  typedef PointKDTree<GPhotonNodeKD>::SearchResult SearchResult;

  GPhotonMap(const std::size_t nbPhotons,
             bool storeSurface,
             const Point &sensorPos,
             const Float cameraSphere = 0.f) :
      m_kdtree(0, PhotonTree::ESlidingMidpoint),
      m_storeSurface(storeSurface), m_nbLightPathAdded(0),
      m_cameraSphere(cameraSphere), m_sensorPos(sensorPos), m_photonSkip(0) {
    m_kdtree.reserve(nbPhotons);
  }

  GPhotonMap(Stream *stream, InstanceManager *manager) {
    SLog(EError, "Not allowed to serialize GPhotonMap");
  }

  void serialize(Stream *stream, InstanceManager *manager) const override {
    SLog(EError, "Not allowed to serialize GPhotonMap");
  }

  template<typename T>
  size_t evaluate(T &query, const Point &pos, Float searchRadius) const {
    size_t count = m_kdtree.executeQuery(pos, searchRadius, query);
    return count;
  }

  inline bool outCapacity() const { return m_kdtree.size() >= m_kdtree.capacity(); }
  inline size_t size() const { return m_kdtree.size(); }

  inline int tryAppend(int idWorker,
                       const Path *lightPath,
                       int minDepth) {
    // Check if there is enough space for the volume map
    // Or for the surface map (notify with flags
    if (outCapacity()) {
      return -1;
    }

    const size_t startIndex = static_cast<const size_t>(std::max(2, minDepth + 1));
    if (lightPath->vertexCount() <= startIndex) {
      return 0;
    }

    // Treat each vertex of the path as a photon
    Spectrum importanceWeights(1.f);
    for (size_t i = 0; i < startIndex - 1; i++) {
      importanceWeights *= lightPath->vertex(i)->weight[EImportance] *
          lightPath->vertex(i)->rrWeight *
          lightPath->edge(i)->weight[EImportance];
    }

    int nbAppend = 0;
    // At least start from 2 to avoid super node and light node
    for (size_t i = startIndex; i < lightPath->vertexCount(); i++) {

      // Photon flux
      importanceWeights *= lightPath->vertex(i - 1)->weight[EImportance] *
          lightPath->vertex(i - 1)->rrWeight *
          lightPath->edge(i - 1)->weight[EImportance];

      if (m_storeSurface &&
          lightPath->vertex(i)->isSurfaceInteraction()) {
        // To be consistent with particleproc.cpp and Jensen's two-pass PM in gatherproc.cpp,
        // we should only consider non-specular vertices

        const BSDF *bsdf = lightPath->vertex(i)->getIntersection().getBSDF();
        if (bsdf) {
          int bsdfType = bsdf->getType();
          if (!(bsdfType & BSDF::EDiffuseReflection) && !(bsdfType & BSDF::EGlossyReflection))
            continue;
        }

        if (outCapacity()) {
          continue;
        }

        // Add the surface photon if there is enough
        // space inside the surface photon map
        m_kdtree.push_back(GPhotonNodeKD(GPhotonNodeData(lightPath,
                                                         i,
                                                         importanceWeights,
                                                         m_nbLightPathAdded)));
        nbAppend++;
      } else if ((!m_storeSurface) &&
          lightPath->vertex(i)->isMediumInteraction()) {

        if (outCapacity()) {
          continue;
        }

        // Add the volumic photon if there is enough space inside
        // the volume photon map

        if (!cameraHit(lightPath->vertex(i - 1)->getPosition(),
                       lightPath->vertex(i)->getPosition())) {
          m_kdtree.push_back(GPhotonNodeKD(GPhotonNodeData(lightPath,
                                                           i,
                                                           importanceWeights,
                                                           m_nbLightPathAdded)));
          nbAppend++;
        }
      }
    }

    // Increase the light path ID
    if (nbAppend != 0) {
      m_nbLightPathAdded += 1;
    }
    return nbAppend;
  }

  inline void build(bool recomputeAABB = false) {
    m_kdtree.build(recomputeAABB);
  }

  inline size_t getDepth() const { return m_kdtree.getDepth(); }
  inline GPhotonNodeKD &getPhoton(size_t idx) { return m_kdtree[idx]; }
  inline const GPhotonNodeKD &getPhoton(size_t idx) const { return m_kdtree[idx]; }
  inline size_t nbSkip() const { return m_photonSkip; }

  MTS_DECLARE_CLASS()

public:
  // These functions are for exposing photon data to BRE
  inline size_t nnSearchVolume(const Point &p, Float &sqrSearchRadius,
                               size_t k, SearchResult *results) const {
    return m_kdtree.nnSearch(p, sqrSearchRadius, k, results);
  }
  inline const GPhotonNodeKD &operator[](size_t idx) const { return m_kdtree[idx]; }

private:
  inline bool cameraHit(const Point &prev, const Point &p) {
    if (m_cameraSphere != 0.f) {
      if (isIntersectedPoint(m_sensorPos, prev, p, m_cameraSphere)) {
        m_photonSkip += 1;
        return true;
      }
    }
    return false;
  }

protected:
  // Surface KDTree
  PhotonTree m_kdtree;
  bool m_storeSurface;
  int m_nbLightPathAdded;

  // Configuration
  Float m_cameraSphere;
  Point m_sensorPos;
  size_t m_photonSkip;
};

// --- BRE
class GradientBeamRadianceEstimator : public SerializableObject {
public:
  typedef PhotonMap::IndexType IndexType;

  /**
   * \brief Create a BRE acceleration data structure from
   * an existing volumetric photon map
   */
  GradientBeamRadianceEstimator(const GPhotonMap *photonmap, Float scaleVol = 1.f);

  /**
   * \brief Unserialize a BRE acceleration data structure from
   * a binary data stream
   */
  GradientBeamRadianceEstimator(Stream *stream, InstanceManager *manager) {
    SLog(EError, "Not implemented");
  }

  /// Serialize to a binary data stream
  void serialize(Stream *stream, InstanceManager *manager) const override {
    SLog(EError, "Not implemented");
  }

  /// Compute the beam radiance estimate for the given ray segment and medium
  template<typename T>
  void query(const Ray &ray, const Medium *medium, T &queryRequest, Float randValue) const {
    IndexType *stack = (IndexType *) alloca((m_depth + 1) * sizeof(IndexType));
    IndexType index = 0, stackPos = 1;
    Spectrum result(0.0f);

//      const Spectrum &sigmaT = medium->getSigmaT();
//      const PhaseFunction *phase = medium->getPhaseFunction();

    while (stackPos > 0) {
      const GBRENode &node = m_nodes[index];
      const GPhotonNodeKD &photon = node.photon;

      /* Test against the node's bounding box */
      Float mint, maxt;
      if (!node.aabb.rayIntersect(ray, mint, maxt) || maxt < ray.mint || mint > ray.maxt) {
        index = stack[--stackPos];
        continue;
      }

      /* Recurse on inner photons */
      if (!photon.isLeaf()) {
        if (hasRightChild(index))
          stack[stackPos++] = photon.getRightIndex(index);
        index = photon.getLeftIndex(index);
      } else {
        index = stack[--stackPos];
      }

      Vector originToCenter = node.photon.getPosition() - ray.o;
      Float diskDistance = dot(originToCenter, ray.d), radSqr = node.radius * node.radius;
      Float distSqr = (ray(diskDistance) - node.photon.getPosition()).lengthSquared();

      if (diskDistance > ray.mint && distSqr < radSqr) {

        // Call to compute the contribution of the current photon
        Ray baseRay(ray);
        baseRay.maxt = diskDistance;
        queryRequest.newRayBase(baseRay, medium);

        // Compute all the informations
        queryRequest(photon, node.radius, randValue);
      }
    }
  }

  MTS_DECLARE_CLASS()
protected:
  /// Release all memory
  ~GradientBeamRadianceEstimator() override {
    delete[] m_nodes;
  }

  /// Fit a hierarchy of bounding boxes to the stored photons
  AABB buildHierarchy(IndexType index);

  /// Blurring kernel used by the BRE
  // FIXME: Check if we want to use this kernel or not
  inline Float K2(Float sqrParam) const {
    Float tmp = 1 - sqrParam;
    return (3 / M_PI) * tmp * tmp;
  }

  /**
   * \brief Return whether or not the inner node of the
   * specified index has a right child node.
   *
   * This function is available for convenience and abstracts away some
   * details about the underlying node representation.
   */
  inline bool hasRightChild(IndexType index) const {
    if (Photon::leftBalancedLayout) {
      return 2 * index + 2 < m_photonCount;
    } else {
      return m_nodes[index].photon.getRightIndex(index) != 0;
    }
  }
protected:
  struct GBRENode {
    AABB aabb;
    GPhotonNodeKD photon;
    Float radius;
  };

  GBRENode *m_nodes;
  Float m_scaleFactor;
  size_t m_photonCount;
  size_t m_depth;
};

MTS_NAMESPACE_END

#endif //MITSUBA_GVPM_ACCEL_H
