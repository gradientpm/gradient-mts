#pragma once

#include <mitsuba/mitsuba.h>

#include <mitsuba/core/kdtree.h>
#include <mitsuba/core/octree.h>
#include <mitsuba/core/util.h>

#include <fstream>
#include <mitsuba/render/sahkdtree3.h>
#include <mitsuba/render/gkdtreeMedian.h>
#include <mitsuba/core/fstream.h>

#include "pixeldata.h"

MTS_NAMESPACE_BEGIN

class SimpleHeuristic {
public:
  inline SimpleHeuristic(const AABB &aabb) {}

  inline std::pair<Float, Float> operator()(int axis, Float leftWidth, Float rightWidth) const {
    return std::pair<Float, Float>(0.5f, 0.5f);
  }

  inline static Float getQuantity(const AABB &aabb) {
    return 1000.f;
  }
private:
};

template<typename Derived>
class MedianKDTree3D : public GenericKDTreeMedian<AABB, SimpleHeuristic, Derived> {
public:
  typedef GenericKDTreeMedian<AABB, SimpleHeuristic, Derived> Parent;
  typedef typename KDTreeBase<AABB>::SizeType SizeType;
  typedef typename KDTreeBase<AABB>::IndexType IndexType;
  typedef typename KDTreeBase<AABB>::KDNode KDNode;

  using Parent::m_nodes;
  using Parent::m_aabb;
  using Parent::m_indices;

protected:
  void buildInternal() {
    SizeType primCount = cast()->getPrimitiveCount();
    Parent::setRetract(false);
    Parent::setEmptySpaceBonus(0.99999999f);
    KDLog(EInfo, "Constructing a Median kd-tree (%i primitives) ..", primCount);
    Parent::buildInternal();
  }

  /// Cast to the derived class
  inline Derived *cast() {
    return static_cast<Derived *>(this);
  }

  /// Cast to the derived class (const version)
  inline const Derived *cast() const {
    return static_cast<const Derived *>(this);
  }

public:
};

//////////////////////
// ShapeKDTree: Use ShapeKD from MTS
// for storing GPs
//////////////////////
struct GPNodeCopy {
  bool isLeaf;
  float split;
  unsigned char axis;

  int leftChild;
  int rightChild;
  int parent;

  // Left for the user
  int idData;

  void save(FileStream *f) {
    f->writeBool(isLeaf);
    f->writeFloat(split);
    f->writeUChar(axis);
    f->writeInt(leftChild);
    f->writeInt(rightChild);
    f->writeInt(parent);
    f->writeInt(idData);
  }

  void load(FileStream *f) {
    isLeaf = f->readBool();
    split = f->readFloat();
    axis = f->readUChar();
    leftChild = f->readInt();
    rightChild = f->readInt();
    parent = f->readInt();
    idData = f->readInt();
  }
};

template<typename T>
struct GPNodeKDShape {
  AABB aabb;
  T *gp;
};

// Interface for local importance
class LocalImpTree {
public:
  virtual void copyKDTree(std::vector<GPNodeCopy *> &nodes, int maxDepth) = 0;

  virtual ~LocalImpTree() {
  }
};

template<typename T>
class SAHTree : public SAHKDTree3D<SAHTree<T> >, public LocalImpTree {
public:
  typedef SAHKDTree3D<SAHTree<T> > Parent;
  typedef typename SAHKDTree3D<SAHTree<T> >::IndexType IndexType;
  typedef typename SAHKDTree3D<SAHTree<T> >::SizeType SizeType;
  typedef typename SAHKDTree3D<SAHTree<T> >::KDNode KDNode;
  typedef AABB AABBType;

  using Parent::m_nodes;
  using Parent::m_indices;
public:
  SAHTree() :
      m_nodesGP(0), m_nGPS(0) {
    SLog(EError, "No empty constructor");
  }

  SAHTree(std::vector<std::vector<T> > &gps) :
      m_nodesGP(0), m_nGPS(0) {
    // Compute number Gps
    for (int blockIdx = 0; blockIdx < (int) gps.size(); ++blockIdx) {
      std::vector<T> &gatherPoints = gps[blockIdx];
      for (size_t i = 0; i < gatherPoints.size(); ++i) {
        T &gp = gatherPoints[i];
        if (gp.its.isValid() && gp.depth != -1) {
          m_nGPS++;
        }
      }
    }

    // Alloc memory
    m_nodesGP = new GPNodeKDShape<T>[m_nGPS];
    // Fill table
    Float sumRadii = 0; //< To compute the average gatherpoint size
    size_t currId = 0;
    for (int blockIdx = 0; blockIdx < (int) gps.size(); ++blockIdx) {
      std::vector<T> &gatherPoints = gps[blockIdx];
      for (size_t i = 0; i < gatherPoints.size(); ++i) {
        T &gp = gatherPoints[i];
        if (gp.its.isValid() && gp.depth != -1) {
          GPNodeKDShape<T> &node = m_nodesGP[currId];
          currId++;

          node.gp = (&gp);
          Point center = node.gp->its.p;
          Float radius = node.gp->radius;
          sumRadii += radius;
          node.aabb = AABB(
              center - Vector(radius, radius, radius),
              center + Vector(radius, radius, radius)
          );
        }

      }
    }

    SLog(EInfo, "Some stats on gatherpoints pushed in GatherMap: ");
    SLog(EInfo, "  - Number of valid GP: %i", m_nGPS);
    SLog(EInfo, "  - Average GP radius: %f", sumRadii / m_nGPS);

    // build KD-Tree
    Parent::buildInternal();
  }

  virtual ~SAHTree() {
    delete[] m_nodesGP;
  }

  inline SizeType getPrimitiveCount() const {
    return m_nGPS;
  }

  inline AABB getAABB(IndexType primIdx) const {
    return m_nodesGP[primIdx].aabb;
  }

  inline AABB getClippedAABB(IndexType primIdx, const AABBType &aabb) const {
    AABB aabbClip = m_nodesGP[primIdx].aabb;
    aabbClip.clip(aabb);
    return aabbClip;
  }

  bool intersect(const Ray &ray, IndexType idx,
                 Float mint, Float maxt, Float &t, void *tmp) {
    const GPNodeKDShape<T> &node = m_nodesGP[idx];
    Vector3d o = Vector3d(ray.o) - Vector3d(node.gp->its.p);
    Vector3d d(ray.d);

    double A = d.lengthSquared();
    double B = 2 * dot(o, d);
    double C = o.lengthSquared() - node.gp->radius * node.gp->radius;

    double nearT, farT;
    if (!solveQuadraticDouble(A, B, C, nearT, farT))
      return false;

    if (!(nearT <= maxt && farT >= mint)) /* NaN-aware conditionals */
      return false;

    if (nearT < mint) {
      if (farT > maxt)
        return false;
      t = (Float) farT;
    } else {
      t = (Float) nearT;
    }

    return true;
  }

  template<typename Functor>
  void executeQuery(const Point &p, Functor &func) const {

    const KDNode *__restrict currNode = m_nodes;
    while (EXPECT_TAKEN(!currNode->isLeaf())) {
      const Float splitVal = (Float) currNode->getSplit();
      const int axis = currNode->getAxis();

      if (p[axis] <= splitVal) {
        currNode = currNode->getLeft();
        continue;
      } else {
        if (splitVal < p[axis]) {
          currNode = currNode->getRight();
          continue;
        }
      }
    }

    /* Reached a leaf node */
    for (IndexType entry = currNode->getPrimStart(),
             last = currNode->getPrimEnd(); entry != last; entry++) {
      const IndexType primIdx = m_indices[entry];
      func(m_nodesGP[primIdx].gp);
    }
  }

  int copyNode(const KDNode *currNode, std::vector<GPNodeCopy *> &nodes, int leftDepth, int parent) {
    if (leftDepth < 0) {
      return -1;
    } else {

      GPNodeCopy *node = new GPNodeCopy;
      node->parent = parent;
      int idIn = 0;
      if (currNode->isLeaf()) {
        // Create the node and
        // nothing to push into the queue
        node->isLeaf = true;
        node->axis = -1;
        node->split = 0;
        node->leftChild = -1;
        node->rightChild = -1;
        nodes.push_back(node);
        idIn = (int) nodes.size();
      } else {
        nodes.push_back(node);
        idIn = (int) nodes.size();

        node->isLeaf = false;
        node->axis = currNode->getAxis();
        node->split = (Float) currNode->getSplit();
        node->leftChild = copyNode(currNode->getLeft(), nodes, leftDepth - 1, idIn - 1);
        node->rightChild = copyNode(currNode->getRight(), nodes, leftDepth - 1, idIn - 1);
      }

      return idIn - 1;
    }

  }

  virtual void copyKDTree(std::vector<GPNodeCopy *> &nodes, int maxDepth) {
    copyNode(m_nodes, nodes, maxDepth - 1, -1);
  }

  void printAccelStruct(const std::string &outFilename) {
    // === Write down flag
    std::ofstream outFile(outFilename.c_str(), std::ofstream::out | std::ofstream::binary);
    outFile.write("#KDS", sizeof(char) * 4);

    // === Write down AABB
    for (int i = 0; i < 3; i++) {
      outFile.write((const char *) &this->m_aabb.min[i], sizeof(float));
      outFile.write((const char *) &this->m_aabb.max[i], sizeof(float));
    }

    const KDNode **stack = (const KDNode **) alloca((this->getMaxDepth() + 1) * sizeof(KDNode *));
    IndexType index = 0, stackPos = 1;

    const KDNode *__restrict currNode = m_nodes;
    while (stackPos > 0) {
      const Float splitVal = (Float) currNode->getSplit();
      const int axis = currNode->getAxis();

      if (currNode->isLeaf()) {
        int nbGP = currNode->getPrimEnd() - currNode->getPrimStart();
        outFile.write("l", sizeof(char));
        outFile.write((const char *) &nbGP, sizeof(int));
        currNode = stack[--stackPos];
      } else {
        outFile.write("n", sizeof(char));
        outFile.write((const char *) &splitVal, sizeof(Float));
        outFile.write((const char *) &axis, sizeof(int));

        stack[stackPos++] = currNode->getRight();
        currNode = currNode->getLeft();
      }

    }

  }
protected:
  GPNodeKDShape<T> *m_nodesGP;
  SizeType m_nGPS;
};

template<typename T>
class MedianTree : public MedianKDTree3D<MedianTree<T> >, public LocalImpTree {
public:
  typedef MedianKDTree3D<MedianTree<T> > Parent;
  typedef typename MedianKDTree3D<MedianTree<T> >::IndexType IndexType;
  typedef typename MedianKDTree3D<MedianTree<T> >::SizeType SizeType;
  typedef typename MedianKDTree3D<MedianTree<T> >::KDNode KDNode;
  typedef AABB AABBType;

  using Parent::m_nodes;
  using Parent::m_indices;
public:
  MedianTree() :
      m_nodesGP(0), m_nGPS(0) {
    SLog(EError, "No empty constructor");
  }

  MedianTree(std::vector<std::vector<T> > &gps, bool build, Float radiiScale = 1.f) :
      m_nodesGP(0), m_nGPS(0) {
    // Compute number Gps
    for (int blockIdx = 0; blockIdx < (int) gps.size(); ++blockIdx) {
      std::vector<T> &gatherPoints = gps[blockIdx];
      for (size_t i = 0; i < gatherPoints.size(); ++i) {
        T &gp = gatherPoints[i];
        if (gp->its.isValid() && gp->depth != -1) {
          m_nGPS++;
        }
      }
    }

    // Alloc memory
    m_nodesGP = new GPNodeKDShape<T>[m_nGPS];
    // Fill table
    Float sumRadii = 0; //< To compute the average gatherpoint size
    size_t currId = 0;
    for (int blockIdx = 0; blockIdx < (int) gps.size(); ++blockIdx) {
      std::vector<T> &gatherPoints = gps[blockIdx];
      for (size_t i = 0; i < gatherPoints.size(); ++i) {
        T &gp = gatherPoints[i];
        if (gp->its.isValid() && gp->depth != -1) {
          GPNodeKDShape<T> &node = m_nodesGP[currId];
          currId++;

          node.gp = (&(*gp));
          Point center = node.gp->its.p;
          Float radius = node.gp->points->radius * radiiScale;
          sumRadii += radius;
          node.aabb = AABB(
              center - Vector(radius, radius, radius),
              center + Vector(radius, radius, radius)
          );
        }
      }
    }

    SLog(EInfo, "Some stats on gatherpoints pushed in GatherMap: ");
    SLog(EInfo, "  - Number of valid GP: %i", m_nGPS);
    SLog(EInfo, "  - Average GP radius: %f", sumRadii / m_nGPS);

    // build KD-Tree
    if (build) {
      Parent::buildInternal();
    }
  }

  virtual ~MedianTree() {
    delete[] m_nodesGP;
  }

  void forceBuild() {
    Parent::buildInternal();
  }

  void expendGP(std::vector<std::vector<T> > &gps, Float radiiScale = 1.f) {
    int oldNbGPS = m_nGPS;

    // Increase nb GP
    for (int blockIdx = 0; blockIdx < (int) gps.size(); ++blockIdx) {
      std::vector<T> &gatherPoints = gps[blockIdx];
      for (size_t i = 0; i < gatherPoints.size(); ++i) {
        T &gp = gatherPoints[i];
        if (gp->its.isValid() && gp->depth != -1) {
          m_nGPS++;
        }
      }
    }

    SLog(EInfo, "Found extra valid GP: %i", m_nGPS - oldNbGPS);

    // Realloc array
    GPNodeKDShape<T> *oldNodesGP = m_nodesGP;
    m_nodesGP = new GPNodeKDShape<T>[m_nGPS];

    // Copy old data and free memory
    memcpy(m_nodesGP, oldNodesGP, sizeof(GPNodeKDShape<T>) * oldNbGPS);
    delete[] oldNodesGP;

    size_t currId = oldNbGPS; // Begin on old copied GPs
    // And copy all new GPs
    for (int blockIdx = 0; blockIdx < (int) gps.size(); ++blockIdx) {
      std::vector<T> &gatherPoints = gps[blockIdx];
      for (size_t i = 0; i < gatherPoints.size(); ++i) {
        T &gp = gatherPoints[i];
        if (gp->its.isValid() && gp->depth != -1) {
          GPNodeKDShape<T> &node = m_nodesGP[currId];
          currId++;

          node.gp = (&(*gp));
          Point center = node.gp->its.p;
          Float radius = node.gp->points->radius * radiiScale;
          node.aabb = AABB(
              center - Vector(radius, radius, radius),
              center + Vector(radius, radius, radius)
          );
        }
      }
    }

    SLog(EInfo, "Some stats on gatherpoints pushed in GatherMap: ");
    SLog(EInfo, "  - Number of valid GP: %i", m_nGPS);
  }

  inline SizeType getPrimitiveCount() const {
    return m_nGPS;
  }

  inline AABB getAABB(IndexType primIdx) const {
    return m_nodesGP[primIdx].aabb;
  }

  inline AABB getClippedAABB(IndexType primIdx, const AABBType &aabb) const {
    AABB aabbClip = m_nodesGP[primIdx].aabb;
    aabbClip.clip(aabb);
    return aabbClip;
  }

  bool intersect(const Ray &ray, IndexType idx,
                 Float mint, Float maxt, Float &t, void *tmp) {
    const GPNodeKDShape<T> &node = m_nodesGP[idx];
    Vector3d o = Vector3d(ray.o) - Vector3d(node.gp->its.p);
    Vector3d d(ray.d);

    double A = d.lengthSquared();
    double B = 2 * dot(o, d);
    double C = o.lengthSquared() - node.gp->radius * node.gp->radius;

    double nearT, farT;
    if (!solveQuadraticDouble(A, B, C, nearT, farT))
      return false;

    if (!(nearT <= maxt && farT >= mint)) /* NaN-aware conditionals */
      return false;

    if (nearT < mint) {
      if (farT > maxt)
        return false;
      t = (Float) farT;
    } else {
      t = (Float) nearT;
    }

    return true;
  }

  template<typename Functor>
  void executeQuery(const Point &p, Functor &func) const {

    const KDNode *__restrict currNode = m_nodes;
    while (EXPECT_TAKEN(!currNode->isLeaf())) {
      const Float splitVal = (Float) currNode->getSplit();
      const int axis = currNode->getAxis();
      const KDNode *__restrict farChild;

      if (p[axis] <= splitVal) {
        currNode = currNode->getLeft();
        continue;
      } else {
        if (splitVal < p[axis]) {
          currNode = currNode->getRight();
          continue;
        }
      }
    }

    /* Reached a leaf node */
    for (IndexType entry = currNode->getPrimStart(),
             last = currNode->getPrimEnd(); entry != last; entry++) {
      const IndexType primIdx = m_indices[entry];
      func(m_nodesGP[primIdx].gp);
    }
  }

  int copyNode(const KDNode *currNode, std::vector<GPNodeCopy *> &nodes, int leftDepth, int parent) {
    if (leftDepth < 0) {
      return -1;
    } else {

      GPNodeCopy *node = new GPNodeCopy;
      node->parent = parent;
      int idIn = 0;
      if (currNode->isLeaf()) {
        // Create the node and
        // nothing to push into the queue
        node->isLeaf = true;
        node->axis = -1;
        node->split = 0;
        node->leftChild = -1;
        node->rightChild = -1;
        nodes.push_back(node);
        idIn = (int) nodes.size();
      } else {
        nodes.push_back(node);
        idIn = (int) nodes.size();

        node->isLeaf = false;
        node->axis = currNode->getAxis();
        node->split = (Float) currNode->getSplit();
        node->leftChild = copyNode(currNode->getLeft(), nodes, leftDepth - 1, idIn - 1);
        node->rightChild = copyNode(currNode->getRight(), nodes, leftDepth - 1, idIn - 1);
      }

      return idIn - 1;
    }

  }

  virtual void copyKDTree(std::vector<GPNodeCopy *> &nodes, int maxDepth) {
    copyNode(m_nodes, nodes, maxDepth - 1, -1);
  }

  void printAccelStruct(const std::string &outFilename) {
    // === Write down flag
    std::ofstream outFile(outFilename.c_str(), std::ofstream::out | std::ofstream::binary);
    outFile.write("#KDS", sizeof(char) * 4);

    // === Write down AABB
    for (int i = 0; i < 3; i++) {
      outFile.write((const char *) &this->m_aabb.min[i], sizeof(float));
      outFile.write((const char *) &this->m_aabb.max[i], sizeof(float));
    }

    const KDNode **stack = (const KDNode **) alloca((this->getMaxDepth() + 1) * sizeof(KDNode *));
    IndexType index = 0, stackPos = 1;

    const KDNode *__restrict currNode = m_nodes;
    while (stackPos > 0) {
      const Float splitVal = (Float) currNode->getSplit();
      const int axis = currNode->getAxis();

      if (currNode->isLeaf()) {
        int nbGP = currNode->getPrimEnd() - currNode->getPrimStart();
        outFile.write("l", sizeof(char));
        outFile.write((const char *) &nbGP, sizeof(int));
        currNode = stack[--stackPos];
      } else {
        outFile.write("n", sizeof(char));
        outFile.write((const char *) &splitVal, sizeof(Float));
        outFile.write((const char *) &axis, sizeof(int));

        stack[stackPos++] = currNode->getRight();
        currNode = currNode->getLeft();
      }

    }

  }
protected:
  GPNodeKDShape<T> *m_nodesGP;
  SizeType m_nGPS;
};

//////////////////////
// PointKDTree: Use PointKD from MTS
// for storing GPs
//////////////////////
template<typename T>
class GPNodeKD : public SimpleKDNode<Point, T *> {
public:

  GPNodeKD() {
    data = NULL;
  }
  /// Construct from a photon interaction
  GPNodeKD(T *_data) : data(_data) {
    this->setPosition(data->its.p);
  }

  /// Unserialize from a binary data stream
  GPNodeKD(Stream *stream) : data(NULL) {
    SLog(EError, "No serialized possible ...");
  }

public:
  T *data;
};

template<typename T>
class GPAccelPointKD {
public:
public:
  //////////////////////////////////////////////////////////////////////////
  // Public type definition
  //////////////////////////////////////////////////////////////////////////
  typedef PointKDTree<GPNodeKD<T> > GatherTree;
  //typedef GatherTree::IndexType      IndexType;
  typedef typename PointKDTree<GPNodeKD<T> >::SearchResult SearchResult;
  //typedef GPNodeKD<T>						 SearchResult;

public:
  GPAccelPointKD() {
    SLog(EError, "No empty constructor");
  }

  GPAccelPointKD(std::vector<std::vector<T> > &gps) :
      m_kdtree(0, GatherTree::ELeftBalanced) {
    size_t gatherCount = 0;
    m_maxRadius = 0;
    for (int blockIdx = 0; blockIdx < (int) gps.size(); ++blockIdx) {
      std::vector<T> &gatherPoints = gps[blockIdx];
      for (size_t i = 0; i < gatherPoints.size(); ++i) {
        T &gp = gatherPoints[i];
        if (gp->its.isValid() && gp->depth != -1) {
          gatherCount++;
          m_maxRadius = std::max(gp.radius, m_maxRadius); // Data attached to the  pixel
        }
      }
    }

    SLog(EInfo, "Found the max radius as %f", m_maxRadius);

    m_kdtree.reserve(gatherCount);

    // === Create gather map;
    for (int blockIdx = 0; blockIdx < (int) gps.size(); ++blockIdx) {
      std::vector<T> &gatherPoints = gps[blockIdx];
      for (size_t i = 0; i < gatherPoints.size(); ++i) {
        T &gp = gatherPoints[i];
        if (gp->its.isValid() && gp->depth != -1) {// Only gather points on surfaces
          m_kdtree.push_back(GPNodeKD<T>(&(*gp)));
        }
      }
    }
    m_kdtree.build();
  }

  template<typename Functor>
  void executeQuery(const Point &p, Functor &func) const {
    m_kdtree.executeQuery(p, m_maxRadius, func);
  }

  inline size_t nnSearch(const Point &p, size_t k,
                         SearchResult *results) const {
    return m_kdtree.nnSearch(p, k, results);
  }

  size_t getDepth() const { return m_kdtree.getDepth(); }
  inline GPNodeKD<T> *operator[](size_t idx) { return &m_kdtree[idx]; }
  inline size_t size() const { return m_kdtree.size(); }

  int copyNode(GPNodeKD<T> *currNode, int currentIndex,
               std::vector<GPNodeCopy *> &nodes, int leftDepth, int parent) {
    if (leftDepth < 0) {
      return -1;
    } else {
      GPNodeCopy *node = new GPNodeCopy;
      node->parent = parent;
      int idIn = 0;
      if (currNode->isLeaf()) {
        // Create the node and
        // nothing to push into the queue
        node->isLeaf = true;
        node->axis = -1;
        node->split = 0;
        node->leftChild = -1;
        node->rightChild = -1;
        nodes.push_back(node);
        idIn = (int) nodes.size();
      } else {
        nodes.push_back(node);
        idIn = (int) nodes.size();

        node->isLeaf = false;
        node->axis = (unsigned char) currNode->getAxis();
        node->split = (Float) currNode->getPosition()[node->axis];//+Epsilon;
        node->leftChild = copyNode(&m_kdtree[currNode->getLeftIndex(currentIndex)],
                                   currNode->getLeftIndex(currentIndex),
                                   nodes, leftDepth - 1, idIn - 1);
        node->rightChild = copyNode(&m_kdtree[currNode->getRightIndex(currentIndex)],
                                    currNode->getRightIndex(currentIndex),
                                    nodes, leftDepth - 1, idIn - 1);
      }

      return idIn - 1;
    }

  }

  void copyKDTree(std::vector<GPNodeCopy *> &nodes, int maxDepth) {
//    std::vector<CopyNodePointKDTree> vec;
//    m_kdtree.copy(vec, maxDepth);
    copyNode(&m_kdtree[0], 0, nodes, maxDepth - 1, -1);
  }

protected:
  GatherTree m_kdtree;
  Float m_maxRadius;
};

//////////////////////
// BVH: Use PointKDTree to
// build an efficient structure
//////////////////////

template<typename T>
class GPAccelBVH {
protected:
  struct GPNode {
    AABB aabb;
    GPNodeKD<T> *gp;
  };

  GPNode *m_nodes;
  typedef uint32_t IndexType;

public:
  GPAccelBVH() {
    SLog(EError, "No empty constructor");
  }

  GPAccelBVH(std::vector<std::vector<T> > &gps) :
      m_kdtree(gps) {
    m_nodes = new GPNode[m_kdtree.size()];
    for (size_t i = 0; i < m_kdtree.size(); i++) {
      GPNode &n = m_nodes[i];
      n.gp = m_kdtree[i];
//			n.radius = n.gp->data->radius+2*Epsilon;
    }
    buildHierarchy(0);
  }

  virtual ~GPAccelBVH() {
    delete[] m_nodes;
  }

  template<typename Functor>
  void executeQuery(const Point &p, Functor &func) const {
    IndexType *stack = (IndexType *) alloca((m_kdtree.getDepth() + 1) * sizeof(IndexType));
    IndexType index = 0, stackPos = 1;
    Spectrum result(0.0f);

    while (stackPos > 0) {
      const GPNode &node = m_nodes[index];
      const GPNodeKD<T> &gp = (*node.gp);

      /* Test against the node's bounding box */
      if (!node.aabb.contains(p)) {
        index = stack[--stackPos];
        continue;
      }

      /* Recurse on inner photons */
      if (!gp.isLeaf()) {
        if (hasRightChild(index))
          stack[stackPos++] = gp.getRightIndex(index);
        index = gp.getLeftIndex(index);
      } else {
        index = stack[--stackPos];
      }

      // The distance test in inside the functor
      func(gp.data);
    }
  }

  void printAccelStruct(const std::string &outFilename) {
    std::ofstream outFile(outFilename.c_str(), std::ofstream::out | std::ofstream::binary);
    outFile << "BVH";

    IndexType *stack = (IndexType *) alloca((m_kdtree.getDepth() + 1) * sizeof(IndexType));
    IndexType index = 0, stackPos = 1;

    // === Iterate on all hierachy
    while (stackPos > 0) {
      const GPNode &node = m_nodes[index];
      outFile << node.aabb.min[0] << node.aabb.min[1] << node.aabb.min[2];
      outFile << node.aabb.max[0] << node.aabb.max[1] << node.aabb.max[2];
      outFile << (node.gp->isLeaf() ? 1 : 0);

      if (!node.gp->isLeaf()) {
        if (hasRightChild(index))
          stack[stackPos++] = node.gp->getRightIndex(index);
        index = node.gp->getLeftIndex(index);
      } else {
        index = stack[--stackPos];
      }
    }
  }

protected:
  AABB buildHierarchy(IndexType index) {
    GPNode &node = m_nodes[index];

    Point center = node.gp->data->its.p;
    Float radius = node.gp->data->radius;
    node.aabb = AABB(
        center - Vector(radius, radius, radius),
        center + Vector(radius, radius, radius)
    );

    if (!node.gp->isLeaf()) {
      IndexType left = node.gp->getLeftIndex(index);
      IndexType right = node.gp->getRightIndex(index);
      if (left)
        node.aabb.expandBy(buildHierarchy(left));
      if (right)
        node.aabb.expandBy(buildHierarchy(right));
    }

    return node.aabb;
  }

  inline bool hasRightChild(IndexType index) const {
    return m_nodes[index].gp->getRightIndex(index) != 0;
  }

protected:
  GPAccelPointKD<T> m_kdtree;

};

//////////////////////
// Octree: Use Dynamic Octree
// for storing GPs
//////////////////////
template<typename T>
class GPAccelOctree {
public:
  GPAccelOctree() {
    SLog(EError, "No empty constructor");
  }

  GPAccelOctree(std::vector<std::vector<T> > &gps) :
      m_octree(0) {
    // === Compute the AABB
    AABB aabbGPS;
    for (int blockIdx = 0; blockIdx < (int) gps.size(); ++blockIdx) {
      std::vector<T> &gatherPoints = gps[blockIdx];
      for (size_t i = 0; i < gatherPoints.size(); ++i) {
        T &gp = gatherPoints[i];
        if (gp->its.isValid()) {// Only gather points on surfaces
          aabbGPS.expandBy(gp->its.p);
          aabbGPS.expandBy(Point(gp->its.p.x + gp->radius,
                                 gp->its.p.y + gp->radius,
                                 gp->its.p.z + gp->radius));
          aabbGPS.expandBy(Point(gp->its.p.x - gp->radius,
                                 gp->its.p.y - gp->radius,
                                 gp->its.p.z - gp->radius));
        }
      }
    }

    // === Create empty octree
    m_octree = new DynamicOctree<T *>(aabbGPS);

    // === Push into the octree
#if defined(MTS_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
    for (int blockIdx = 0; blockIdx < (int) gps.size(); ++blockIdx) {
      std::vector<T> &gatherPoints = gps[blockIdx];
      for (size_t i = 0; i < gatherPoints.size(); ++i) {
        T &gp = gatherPoints[i];
        if (gp->its.isValid()) {// Only gather points on surfaces
          AABB aabbGP(gp->its.p);
          /*aabbgp->expandBy(Point(gp->its.p.x + gp->radius,
                  gp->its.p.y + gp->radius,
                  gp->its.p.z + gp->radius));
          aabbgp->expandBy(Point(gp->its.p.x - gp->radius,
                  gp->its.p.y - gp->radius,
                  gp->its.p.z - gp->radius));*/
          // aabbgp: Not declared ???
          aabbGP.expandBy(Point(gp->its.p.x + gp->radius,
                                gp->its.p.y + gp->radius,
                                gp->its.p.z + gp->radius));
          aabbGP.expandBy(Point(gp->its.p.x - gp->radius, gp->its.p.y - gp->radius, gp->its.p.z - gp->radius));
          m_octree->insert(&(*gp), aabbGP);
        }

      }
    }
  }

  virtual ~GPAccelOctree() {
    delete m_octree;
  }

  template<typename Functor>
  void executeQuery(const Point &p, Functor &func) const {
    m_octree->lookup(p, func);
  }

protected:
  DynamicOctree<T *> *m_octree;
};

///////////////////////
// Default choice
//////////////////////
#define USE_NODE_FUNC 0
template<class T>
struct GPAccel {
  //typedef GPAccelPointKD<T> Type;
  //typedef GPAccelBVH<T> Type;
  typedef SAHTree<T> Type;
};

MTS_NAMESPACE_END
