#include "gvpm_accel.h"

#include "shift/shift_utilities.h"
#include "shift/shift_surface.h"
#include "shift/shift_volume_photon.h"

MTS_NAMESPACE_BEGIN

// --- BRE
GradientBeamRadianceEstimator::GradientBeamRadianceEstimator(const GPhotonMap *photonmap, Float scaleVol) {

  m_photonCount = photonmap->size();
  m_scaleFactor = 1.f;
  m_depth = photonmap->getDepth();

  Log(EInfo, "Allocating %s of memory for the BRE acceleration data structure",
      memString(sizeof(GBRENode) * m_photonCount).c_str());
  m_nodes = new GBRENode[m_photonCount];

  for (int i = 0; i < (int) m_photonCount; ++i) {
    const GPhotonNodeKD &photon = photonmap->getPhoton(i);
    GBRENode &node = m_nodes[i];
    node.photon = photon;

    /* Compute photon radius based on a locally uniform density assumption */
    node.radius = scaleVol;
  }

  Log(EInfo, "Generating a hierarchy for the beam radiance estimate");
  ref<Timer> timer = new Timer();
  buildHierarchy(0);
  Log(EInfo, "Done (took %i ms)", timer->getMilliseconds());
}

AABB GradientBeamRadianceEstimator::buildHierarchy(IndexType index) {
  GBRENode &node = m_nodes[index];

  Point center = node.photon.getPosition();
  Float radius = node.radius;
  node.aabb = AABB(
      center - Vector(radius, radius, radius),
      center + Vector(radius, radius, radius)
  );

  if (!node.photon.isLeaf()) {
    IndexType left = node.photon.getLeftIndex(index);
    IndexType right = node.photon.getRightIndex(index);
    if (left)
      node.aabb.expandBy(buildHierarchy(left));
    if (right)
      node.aabb.expandBy(buildHierarchy(right));
  }

  return node.aabb;
}

MTS_IMPLEMENT_CLASS_S(GPhotonMap, false, SerializableObject)
MTS_IMPLEMENT_CLASS_S(GradientBeamRadianceEstimator, false, Object)

MTS_NAMESPACE_END
