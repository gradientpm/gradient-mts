#if !defined(__BEAMS_ACCEL_H)
#define __BEAMS_ACCEL_H

#include <mitsuba/render/photonmap.h>
#include <mitsuba/core/kdtree.h> // For subBeams kdTree

#include "beams_struct.h"

MTS_NAMESPACE_BEGIN

struct PhotonSubBeamData {
  const PhotonBeam *beam;
  Float t1, t2; // Start point for the beam
};
struct PhotonSubBeam :
#if MTS_PHOTONMAP_LEFT_BALANCED == 1
    public LeftBalancedKDNode<Point, PhotonSubBeamData> {
#else
    public SimpleKDNode<Point, PhotonSubBeamData> {
#endif
public:
  /// Dummy constructor
  inline PhotonSubBeam() {}

  /// Construct from a photon interaction
  PhotonSubBeam(const Point &pos, const PhotonBeam *beam, Float t1, Float t2);

};

class PhotonSubBeamMap : public SerializableObject {
public:
  typedef PointKDTree<PhotonSubBeam> PhotonSubBeamTree;
  typedef PhotonSubBeamTree::IndexType IndexType;
  typedef PhotonSubBeamTree::SearchResult SearchResult;

  PhotonSubBeamMap(size_t subBeamCount = 0) :
      m_kdtree(0, PhotonSubBeamTree::EBalanced) {
    m_kdtree.reserve(subBeamCount);
  }
  PhotonSubBeamMap(Stream *stream, InstanceManager *manager) {
    SLog(EError, "Not Implemented");
  }

  /// Clear the kd-tree array
  inline void clear() { m_kdtree.clear(); }
  /// Resize the kd-tree array
  inline void resize(size_t size) { m_kdtree.resize(size); }
  /// Reserve a certain amount of memory for the kd-tree array
  inline void reserve(size_t size) { m_kdtree.reserve(size); }
  /// Return the size of the kd-tree
  inline size_t size() const { return m_kdtree.size(); }
  /// Return the capacity of the kd-tree
  inline size_t capacity() const { return m_kdtree.capacity(); }
  /// Append a kd-tree photon to the photon array
  inline void push_back(const PhotonSubBeam &photon) { m_kdtree.push_back(photon); }
  /// Return one of the photons by index
  inline PhotonSubBeam &operator[](size_t idx) { return m_kdtree[idx]; }
  /// Return one of the photons by index (const version)
  inline const PhotonSubBeam &operator[](size_t idx) const { return m_kdtree[idx]; }

  // Build the KD Tree
  inline void build(bool recomputeAABB = false) { m_kdtree.build(recomputeAABB); }
  /// Return the depth of the constructed KD-tree
  inline size_t getDepth() const { return m_kdtree.getDepth(); }

  /// Serialize a photon map to a binary data stream
  void serialize(Stream *stream, InstanceManager *manager) const {
    SLog(EError, "Not Implemented");
  }

  MTS_DECLARE_CLASS()
protected:
  /// Virtual destructor
  virtual ~PhotonSubBeamMap() {
    // FIXME: double check if there is not free memory here
  }
protected:
  PhotonSubBeamTree m_kdtree;
};

template<typename T>
class SubBeamBVH : public SerializableObject {
public:
  typedef PhotonSubBeamMap::IndexType IndexType;

  /**
   * \brief Create a BRE acceleration data structure from
   * an existing volumetric photon map
   */
  SubBeamBVH(const std::vector<std::pair<int, T> > &beams) {
    // First thing, we count how many sub beams we will create.
    // We have two policy:
    // 1) A beam cannot be decimate more than N pieces
    // 2) The sub beams need to be more than N meter length

    // To decide the size of the beam cut
    // We will look at the average size
    Float avgSize = 0.f;
    for (auto b : beams) {
      avgSize += b.second.getLength();
    }
    avgSize /= beams.size();

    Float subbeamSize = avgSize / 10;
    SLog(EInfo, "Number of beams: %i", beams.size());
    SLog(EInfo, "Subbeam size: %f", subbeamSize);

    m_beamCount = 0;
    for (auto b : beams) {
      m_beamCount += b.second.nbSubBeams(subbeamSize);
    }

    if (m_beamCount == 0) {
      SLog(EWarn, "No photon beam generated, abort the map creation");
      m_nodes = 0;
      return;
    }

    SLog(EInfo, "Create %i sub-beams ...", m_beamCount);
    ref<PhotonSubBeamMap> subBeamKDTree = new PhotonSubBeamMap(m_beamCount);
    size_t nbCut = 0;
    for (size_t j = 0; j < beams.size(); j++) {
      int nbSBeams = beams[j].second.nbSubBeams(subbeamSize);
      Float lengthSubBeams = beams[j].second.getLength() / nbSBeams;
      nbCut += nbSBeams;
      const Vector &dPB = beams[j].second.getDir();
      for (int i = 0; i < nbSBeams; i++) {
        PhotonSubBeam s(beams[j].second.getOri() + dPB * lengthSubBeams * (i + 0.5), &beams[j].second,
                        lengthSubBeams * i, lengthSubBeams * (i + 1));
        subBeamKDTree->push_back(s);
      }
    }
    SLog(EInfo, "Avg number sub beam per beams: %f", float(nbCut) / beams.size());
    SLog(EInfo, "Build sub-beam KDTree");
    subBeamKDTree->build();
    m_depth = subBeamKDTree->getDepth();
    SLog(EInfo, "Sub-beam KDTree depth: %i", m_depth);

    SLog(EInfo, "Build sub-beam BVH");
    m_nodes = new SubBeamNode[m_beamCount];
    for (int i = 0; i < (int) m_beamCount; ++i) {
      SubBeamNode &node = m_nodes[i];
      node.subbeam = subBeamKDTree->operator[](i);
    }

    buildHierarchy(0);
  }

  /// Release all memory
  virtual ~SubBeamBVH() {
    if (m_nodes != 0)
      delete[] m_nodes;
  }

  /**
   * \brief Unserialize a SubBeams acceleration data structure from
   * a binary data stream
   */
  SubBeamBVH(Stream *stream, InstanceManager *manager) {
    SLog(EError, "Not implemented");
  }

  /// Serialize to a binary data stream
  void serialize(Stream *stream, InstanceManager *manager) const {
    SLog(EError, "Not Implemented");
  }

  /// Compute the beam radiance estimate for the given ray segment and medium
  template<typename Q>
  Spectrum query(Q &bRadQuery) const {
    if (m_beamCount == 0)
      return Spectrum(0.f); // Case empty subbeams ...
    const Ray &r = bRadQuery.baseCameraRay;
    const Ray ray(r(r.mint), r.d, 0, r.maxt - r.mint, r.time);
    IndexType *stack = (IndexType *) alloca((m_depth + 1) * sizeof(IndexType));
    IndexType index = 0, stackPos = 1;
    Spectrum result(0.0f);

    while (stackPos > 0) {
      const SubBeamNode &node = m_nodes[index];
      const PhotonSubBeam &subbeam = node.subbeam;

      /* Test against the node's bounding box */
      Float mint, maxt;
      if (!node.aabb.rayIntersect(ray, mint, maxt) || maxt < ray.mint || mint > ray.maxt) {
        index = stack[--stackPos];
        continue;
      }

      /* Recurse on inner photons */
      if (!subbeam.isLeaf()) {
        if (hasRightChild(index))
          stack[stackPos++] = subbeam.getRightIndex(index);
        index = subbeam.getLeftIndex(index);
      } else {
        index = stack[--stackPos];
      }

      bRadQuery((T *) subbeam.data.beam, subbeam.data.t1, subbeam.data.t2);
    }

    return result;
  }

protected:

  /// Fit a hierarchy of bounding boxes to the stored photons
  AABB buildHierarchy(IndexType index) {
    SubBeamNode &node = m_nodes[index];

    //Get beam radius
    Float radiusMin = node.subbeam.data.beam->radius;
    Float radiusMax = radiusMin; // Same

    // Get the information of the real beam
    Point p1Beam = node.subbeam.data.beam->getOri();
    Vector dBeam = node.subbeam.data.beam->getDir();

    // Compute the bounding box of the subBeam
    Point p1 = p1Beam + dBeam * node.subbeam.data.t1;
    Point p2 = p1Beam + dBeam * node.subbeam.data.t2;
    node.aabb = AABB(
        p1 - Vector(radiusMin, radiusMin, radiusMin),
        p1 + Vector(radiusMin, radiusMin, radiusMin)
    );
    node.aabb.expandBy(
        AABB(
            p2 - Vector(radiusMax, radiusMax, radiusMax),
            p2 + Vector(radiusMax, radiusMax, radiusMax)
        )
    );

    if (!node.subbeam.isLeaf()) {
      IndexType left = node.subbeam.getLeftIndex(index);
      IndexType right = node.subbeam.getRightIndex(index);
      if (left)
        node.aabb.expandBy(buildHierarchy(left));
      if (right)
        node.aabb.expandBy(buildHierarchy(right));
    }

    return node.aabb;
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
      return 2 * index + 2 < m_beamCount;
    } else {
      return m_nodes[index].subbeam.getRightIndex(index) != 0;
    }
  }
protected:
  struct SubBeamNode {
    AABB aabb;
    PhotonSubBeam subbeam;
  };
  SubBeamNode *m_nodes;
  size_t m_beamCount;
  size_t m_depth;
};

MTS_NAMESPACE_END

#endif
