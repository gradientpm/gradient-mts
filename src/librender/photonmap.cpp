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

#include <mitsuba/render/photonmap.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/phase.h>
#include <fstream>

MTS_NAMESPACE_BEGIN

PhotonMap::PhotonMap(size_t photonCount)
		: m_kdtree(0, PhotonTree::ESlidingMidpoint), m_scale(1.0f) {
	m_kdtree.reserve(photonCount);
	Assert(Photon::m_precompTableReady);
}

PhotonMap::PhotonMap(Stream *stream, InstanceManager *manager)
    : SerializableObject(stream, manager),
	  m_kdtree(0, PhotonTree::ESlidingMidpoint) {
	Assert(Photon::m_precompTableReady);
	m_scale = (Float) stream->readFloat();
	m_kdtree.resize(stream->readSize());
	m_kdtree.setDepth(stream->readSize());
	m_kdtree.setAABB(AABB(stream));
	for (size_t i=0; i<m_kdtree.size(); ++i)
		m_kdtree[i] = Photon(stream);
}

void PhotonMap::serialize(Stream *stream, InstanceManager *manager) const {
	Log(EDebug, "Serializing a photon map (%s)",
		memString(m_kdtree.size() * sizeof(Photon)).c_str());
	stream->writeFloat(m_scale);
	stream->writeSize(m_kdtree.size());
	stream->writeSize(m_kdtree.getDepth());
	m_kdtree.getAABB().serialize(stream);
	for (size_t i=0; i<m_kdtree.size(); ++i)
		m_kdtree[i].serialize(stream);
}

PhotonMap::~PhotonMap() {
}

std::string PhotonMap::toString() const {
	std::ostringstream oss;
	oss << "PhotonMap[" << endl
		<< "  size = " << m_kdtree.size() << "," << endl
		<< "  capacity = " << m_kdtree.capacity() << "," << endl
		<< "  aabb = " << m_kdtree.getAABB().toString() << "," << endl
		<< "  depth = " << m_kdtree.getDepth() << "," << endl
		<< "  scale = " << m_scale << endl
		<< "]";
	return oss.str();
}
void PhotonMap::dumpOBJ(const std::string &filename) {
	std::ofstream os(filename.c_str());
	os << "o Photons" << endl;
	for (size_t i=0; i<m_kdtree.size(); ++i) {
		const Point &p = m_kdtree[i].getPosition();
		os << "v " << p.x << " " << p.y << " " << p.z << endl;
	}

	/// Need to generate some fake geometry so that blender will import the points
	for (size_t i=3; i<=m_kdtree.size(); i++)
		os << "f " << i << " " << i-1 << " " << i-2 << endl;
	os.close();
}

Spectrum PhotonMap::estimateIrradiance(
		const Point &p, const Normal &n,
		Float searchRadius, int maxDepth,
		size_t maxPhotons) const {
	SearchResult *results = static_cast<SearchResult *>(
		alloca((maxPhotons+1) * sizeof(SearchResult)));
	Float squaredRadius = searchRadius*searchRadius;
	size_t resultCount = nnSearch(p, squaredRadius, maxPhotons, results);
	Float invSquaredRadius = 1.0f / squaredRadius;

	/* Sum over all contributions */
	Spectrum result(0.0f);
	for (size_t i=0; i<resultCount; i++) {
		const SearchResult &searchResult = results[i];
		const Photon &photon = m_kdtree[searchResult.index];
		if (photon.getDepth() > maxDepth)
			continue;

		Vector wi = -photon.getDirection();
		Vector photonNormal = photon.getNormal();
		Float wiDotGeoN = dot(photonNormal, wi),
			  wiDotShN  = dot(n, wi);

		/* Only use photons from the top side of the surface */
		if (dot(wi, n) > 0 && dot(photonNormal, n) > 1e-1f && wiDotGeoN > 1e-2f) {
			/* Account for non-symmetry due to shading normals */
			Spectrum power = photon.getPower() * std::abs(wiDotShN / wiDotGeoN);

			/* Weight the samples using Simpson's kernel */
			Float sqrTerm = 1.0f - searchResult.distSquared*invSquaredRadius;

			result += power * (sqrTerm*sqrTerm);
		}
	}

	/* Based on the assumption that the surface is locally flat,
	   the estimate is divided by the area of a disc corresponding to
	   the projected spherical search region */
	return result * (m_scale * 3 * INV_PI * invSquaredRadius);
}

Spectrum PhotonMap::estimateRadiance(const Intersection &its,
		Float searchRadius, size_t maxPhotons) const {
	SearchResult *results = static_cast<SearchResult *>(
		alloca((maxPhotons+1) * sizeof(SearchResult)));
	Float squaredRadius = searchRadius*searchRadius;
	size_t resultCount = nnSearch(its.p, squaredRadius, maxPhotons, results);
	Float invSquaredRadius = 1.0f / squaredRadius;

	/* Sum over all contributions */
	Spectrum result(0.0f);
	const BSDF *bsdf = its.getBSDF();
	for (size_t i=0; i<resultCount; i++) {
		const SearchResult &searchResult = results[i];
		const Photon &photon = m_kdtree[searchResult.index];
		Float sqrTerm = 1.0f - searchResult.distSquared*invSquaredRadius;

		Vector wi = its.toLocal(-photon.getDirection());

		BSDFSamplingRecord bRec(its, wi, its.wi, EImportance);
		result += photon.getPower() * bsdf->eval(bRec) * (sqrTerm*sqrTerm);
	}

	/* Based on the assumption that the surface is locally flat,
	   the estimate is divided by the area of a disc corresponding to
	   the projected spherical search region */
	return result * (m_scale * 3 * INV_PI * invSquaredRadius);
}

struct RawRadianceQuery {
	RawRadianceQuery(const Intersection &its, int maxDepth)
	  : its(its), maxDepth(maxDepth), result(0.0f) {
		bsdf = its.getBSDF();
	}

	inline void operator()(const Photon &photon) {
		Normal photonNormal(photon.getNormal());
		Vector wi = -photon.getDirection();
		Float wiDotGeoN = absDot(photonNormal, wi);

		if (photon.getDepth() > maxDepth
			|| dot(photonNormal, its.shFrame.n) < 1e-1f
			|| wiDotGeoN < 1e-2f)
			return;

		BSDFSamplingRecord bRec(its, its.toLocal(wi), its.wi, EImportance);

		Spectrum value = photon.getPower() * bsdf->eval(bRec);
		if (value.isZero())
			return;

		/* Account for non-symmetry due to shading normals */
		value *= std::abs(Frame::cosTheta(bRec.wi) /
			(wiDotGeoN * Frame::cosTheta(bRec.wo)));

		result += value;
	}

	const Intersection &its;
	const BSDF *bsdf;
	int maxDepth;
	Spectrum result;
};

size_t PhotonMap::estimateRadianceRaw(const Intersection &its,
		Float searchRadius, Spectrum &result, int maxDepth) const {
	RawRadianceQuery query(its, maxDepth);
	size_t count = m_kdtree.executeQuery(its.p, searchRadius, query);
	result = query.result;
	return count;
}

struct RadianceQueryGp {
  RadianceQueryGp(const Intersection &its, int maxDepth, int currSComp, Float radius)
		  : its(its), currSampledComponent(currSComp), maxDepth(maxDepth),
			searchRadius(radius),
			result(0.0f), M(0) {
	  bsdf = its.getBSDF();
  }

  inline void operator()(const Photon &photon) {
	  Normal photonNormal(photon.getNormal());
	  Vector wi = -photon.getDirection();

#ifndef MTS_NOSHADINGNORMAL
	  Float wiDotGeoN = absDot(photonNormal, wi);
#endif

	  // Test the radius
	  Float lengthSqr = (its.p - photon.position).lengthSquared();
	  if ((searchRadius * searchRadius - lengthSqr) < 0)
		  return;

	  if (dot(photonNormal, its.shFrame.n) < 1e-1f
#ifndef MTS_NOSHADINGNORMAL
		  || wiDotGeoN < 1e-2f
#endif
			  )
		  return;

	  // Test the depth of the photon
	  if (maxDepth > 0 && photon.getDepth() > maxDepth)
		  return;

	  // Count the number of photon which we can collect
	  M += 1;

	  // Accumulate the contribution of the photon
	  BSDFSamplingRecord bRec(its, its.toLocal(wi), its.wi, EImportance);
	  bRec.component = currSampledComponent;
#ifdef MTS_NOSHADINGNORMAL
	  Spectrum value = photon.getPower() * (bsdf->eval(bRec)) / std::abs(Frame::cosTheta(bRec.wo));
#else
	  Spectrum value = photon.getPower() * bsdf->eval(bRec);  //Frame::cosTheta(bRec.wo)
#endif

      // pdf of the sampled component already pre-multiplied into photon's weight
	  if (value.isZero())
		  return;

	  // Account for non-symmetry due to shading normals
#ifndef MTS_NOSHADINGNORMAL
	  value *= std::abs(
      Frame::cosTheta(bRec.wi) / (wiDotGeoN * Frame::cosTheta(bRec.wo)));
#endif

	  // Accumulate the results
	  result += value;
  }

  // Information GP
  const Intersection &its;
  const BSDF *bsdf;
  const int currSampledComponent;

  // Other
  int maxDepth;
  Float searchRadius;

  // Results
  Spectrum result;
  int M;
};

size_t PhotonMap::estimateRadianceGP(const Intersection &its,
									 Float searchRadius, Spectrum &result, int maxDepth,
									 int currSComp) const {
	RadianceQueryGp query(its, maxDepth, currSComp, searchRadius);
	m_kdtree.executeQuery(its.p, searchRadius, query);
	result = query.result;
	return query.M;
}

struct RadianceQueryVolume {
  RadianceQueryVolume(const Point &p, const Vector _viewDir, const MediumSamplingRecord &mRec, int maxDepth, Float radius)
	  : pos(p), mRec(mRec), viewDir(_viewDir), maxDepth(maxDepth), searchRadius(radius),
		result(0.0f), M(0) {
  }

  inline void operator()(const Photon &photon) {
	  Vector wi = -photon.getDirection();


	  // Test the radius
	  Float lengthSqr = (pos - photon.position).lengthSquared();
	  if ((searchRadius * searchRadius - lengthSqr) < 0)
		  return;

	  // Test the depth of the photon
	  if (maxDepth > 0 && photon.getDepth() > maxDepth)
		  return;

	  // Count the number of photon which we can collect
	  M += 1;

	  // Accumulate the contribution of the photon
	  PhaseFunctionSamplingRecord pRec(mRec,  wi, -viewDir, EImportance);
	  Spectrum value = photon.getPower() * mRec.getPhaseFunction()->eval(pRec);

	  if (value.isZero())
		  return;

	  // Accumulate the results
	  result += value;
  }

  // Information GP
  const Point& pos;
  const MediumSamplingRecord &mRec;
  const Vector viewDir;
  // Other
  int maxDepth;
  Float searchRadius;

  // Results
  Spectrum result;
  int M;
};


size_t PhotonMap::estimateVolumeRadiance(const MediumSamplingRecord &mRec, const Vector& viewDir,
										 Float searchRadius, Spectrum &result, int maxDepth) const {
  RadianceQueryVolume query(mRec.p, viewDir, mRec, maxDepth, searchRadius);
  m_kdtree.executeQuery(mRec.p, searchRadius, query);
  result = query.result;
  return query.M;
}

MTS_IMPLEMENT_CLASS_S(PhotonMap, false, SerializableObject)
MTS_NAMESPACE_END
