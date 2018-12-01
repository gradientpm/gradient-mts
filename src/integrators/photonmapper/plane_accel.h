//
// Created by beltegeuse on 8/30/17.
//

#ifndef MITSUBA_PLANE_ACCEL_H
#define MITSUBA_PLANE_ACCEL_H

#include <mitsuba/render/photonmap.h>
#include <mitsuba/core/kdtree.h> // For subBeams kdTree

#include "plane_struct.h"

MTS_NAMESPACE_BEGIN

struct PhotonPlaneKDNode :
#if MTS_PHOTONMAP_LEFT_BALANCED == 1
    public LeftBalancedKDNode<Point, PhotonSubBeamData> {
#else
    public SimpleKDNode<Point, const PhotonPlane *> {
#endif
public:
  /// Dummy constructor
  inline PhotonPlaneKDNode() {}

  /// Construct from a photon interaction
  PhotonPlaneKDNode(const Point &pos, const PhotonPlane *plane) {
    position = pos;
    data = plane;
    flags = 0;
  }

};

class PhotonPlaneMap : public SerializableObject {
public:
  typedef PointKDTree<PhotonPlaneKDNode> PhotonPlaneTree;
  typedef PhotonPlaneTree::IndexType IndexType;
  typedef PhotonPlaneTree::SearchResult SearchResult;

  PhotonPlaneMap(size_t subBeamCount = 0) :
      m_kdtree(0, PhotonPlaneTree::EBalanced) {
    m_kdtree.reserve(subBeamCount);
  }
  PhotonPlaneMap(Stream *stream, InstanceManager *manager) {
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
  inline void push_back(const PhotonPlaneKDNode &photon) { m_kdtree.push_back(photon); }
  /// Return one of the photons by index
  inline PhotonPlaneKDNode &operator[](size_t idx) { return m_kdtree[idx]; }
  /// Return one of the photons by index (const version)
  inline const PhotonPlaneKDNode &operator[](size_t idx) const { return m_kdtree[idx]; }

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
  virtual ~PhotonPlaneMap() {
  }
protected:
  PhotonPlaneTree m_kdtree;
};

template<typename T>
class PhotonPlaneBVH : public SerializableObject {
public:
  typedef PhotonPlaneMap::IndexType IndexType;

  /**
   * \brief Create a BRE acceleration data structure from
   * an existing volumetric photon map
   */
  PhotonPlaneBVH(const std::vector<T> &planes) {
    ref<PhotonPlaneMap> planeKDTree = new PhotonPlaneMap(planes.size());
    for (size_t j = 0; j < planes.size(); j++) {
      PhotonPlaneKDNode p(planes[j].getCenter(), &planes[j]);
      planeKDTree->push_back(p);
    }

    planeKDTree->build();
    m_planeCount = planes.size();
    m_depth = planeKDTree->getDepth();
    m_nodes = new PhotonPlaneNode[planes.size()];
    for (int i = 0; i < (int) planes.size(); ++i) {
      PhotonPlaneNode &node = m_nodes[i];
      node.kdnode = planeKDTree->operator[](i);
    }

    buildHierarchy(0);
  }

  /// Release all memory
  virtual ~PhotonPlaneBVH() {
    if (m_nodes != 0)
      delete[] m_nodes;
  }

  /**
   * \brief Unserialize a Photon Plane acceleration data structure from
   * a binary data stream
   */
  PhotonPlaneBVH(Stream *stream, InstanceManager *manager) {
    SLog(EError, "Not implemented");
  }

  /// Serialize to a binary data stream
  void serialize(Stream *stream, InstanceManager *manager) const {
    SLog(EError, "Not Implemented");
  }

  /// Compute the beam radiance estimate for the given ray segment and medium
  template<typename Q>
  void query(Q &bRadQuery) const {
    if (m_planeCount == 0)
      return; // Early termination
    const Ray &r = bRadQuery.baseCameraRay;
    const Ray ray(r(r.mint), r.d, 0, r.maxt - r.mint, r.time);
    IndexType *stack = (IndexType *) alloca((m_depth + 1) * sizeof(IndexType));
    IndexType index = 0, stackPos = 1;

    while (stackPos > 0) {
      const PhotonPlaneNode &node = m_nodes[index];
      const PhotonPlaneKDNode &nodekd = node.kdnode;

      /* Test against the node's bounding box */
      Float mint, maxt;
      if (!node.aabb.rayIntersect(ray, mint, maxt) || maxt < ray.mint || mint > ray.maxt) {
        index = stack[--stackPos];
        continue;
      }

      /* Recurse on inner photons */
      if (!nodekd.isLeaf()) {
        if (hasRightChild(index))
          stack[stackPos++] = nodekd.getRightIndex(index);
        index = nodekd.getLeftIndex(index);
      } else {
        index = stack[--stackPos];
      }

      bRadQuery((T *) nodekd.getData());
    }
  }

  void topNodeAABBDump() {
    SLog(EInfo, "%s", m_nodes[0].aabb.toString().c_str());
  }

protected:
  /// Fit a hierarchy of bounding boxes to the stored photons
  AABB buildHierarchy(IndexType index) {
    PhotonPlaneNode &node = m_nodes[index];
    node.aabb = node.kdnode.getData()->getAABB();

    if (!node.kdnode.isLeaf()) {
      IndexType left = node.kdnode.getLeftIndex(index);
      IndexType right = node.kdnode.getRightIndex(index);
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
      return 2 * index + 2 < m_planeCount;
    } else {
      return m_nodes[index].kdnode.getRightIndex(index) != 0;
    }
  }

protected:
  struct PhotonPlaneNode {
    AABB aabb;
    PhotonPlaneKDNode kdnode;
  };
  PhotonPlaneNode *m_nodes;
  size_t m_planeCount;
  size_t m_depth;

};

MTS_NAMESPACE_END

#endif //MITSUBA_PLANE_ACCEL_H
