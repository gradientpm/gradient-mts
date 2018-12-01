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

#include <mitsuba/render/scene.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

static StatsCounter avgPathLength("Path tracer", "Average path length", EAverage);

// Primitive path tracer
class MIPathTracer : public MonteCarloIntegrator {
public:
	MIPathTracer(const Properties &props)
		: MonteCarloIntegrator(props) {
		m_roughness = props.getFloat("bounceRoughness", 0.05);
      m_showGlossySurfaces = false;
	}

	/// Unserialize from a binary data stream
	MIPathTracer(Stream *stream, InstanceManager *manager)
		: MonteCarloIntegrator(stream, manager) { }

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
		/* Some aliases and local variables */
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		RayDifferential ray(r);
		Spectrum Li(0.0f);
		bool scattered = false;

		/* Perform the first ray intersection (or ignore if the
		   intersection has already been provided). */
		rRec.rayIntersect(ray);
		ray.mint = Epsilon;

		Spectrum throughput(1.0f);
		Float eta = 1.0f;

		while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
			if (!its.isValid()) {
				/* If no intersection could be found, potentially return
				   radiance from a environment luminaire if it exists */
				if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
					&& (!m_hideEmitters || scattered))
					Li += throughput * scene->evalEnvironment(ray);
				break;
			}

			const BSDF *bsdf = its.getBSDF(ray);

			/* Possibly include emitted radiance if requested */
			if (its.isEmitter())
				Li += throughput * its.Le(-ray.d);

			if ((rRec.depth >= m_maxDepth && m_maxDepth > 0)
				|| (m_strictNormals && dot(ray.d, its.geoFrame.n)
					* Frame::cosTheta(its.wi) >= 0)) {

				/* Only continue if:
				   1. The current path length is below the specifed maximum
				   2. If 'strictNormals'=true, when the geometric and shading
				      normals classify the incident direction to the same side */
				break;
			}

			/* ==================================================================== */
			/*                            BSDF sampling                             */
			/* ==================================================================== */

			/* Sample BSDF * cos(theta) */
			Float bsdfPdf;
			BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);

			// Speedup the process by enumerate the choices
			Float pdfComp = 1.f;
			Point2 randomBSDFChoice = rRec.nextSample2D();
			int compSel = bsdf->sampleComponent(bRec, pdfComp, randomBSDFChoice, m_roughness);
			if(compSel == -1) {
				// In this case, if one of the component is rough (diffuse) all are same
				if(bsdf->getRoughness(its, 0) >= m_roughness) {
					break;
				}
			} else {
				// There is a choice here
				// Count the number of glossy component
				int nbGlossyComp = 0;
				int lastCompGlossy = -1;
				for(int comp = 0; comp < bsdf->getComponentCount(); comp++) {
					if(bsdf->getRoughness(its, comp) < m_roughness) {
						nbGlossyComp++;
						lastCompGlossy = comp;
					}
				}

				bool twoSided = false;
				if(nbGlossyComp > 1) {
					if(bsdf->getClass()->getName() == "TwoSidedBRDF") {
						twoSided = true;
					} else {
						SLog(EError, "Not two sided??? %s", bsdf->toString().c_str());
					}
				}

				// Two cases:
				// - If there is only one, force selection of it
				// - If there is more than one, pick prob of them
				if(nbGlossyComp == 1) {
					compSel = lastCompGlossy;
					bRec.component = compSel;
					pdfComp = 1.f;
				} else {
					if(twoSided) {
						if(nbGlossyComp != 2) {
							SLog(EError, "Too complicated two sided");
						}

						int startID = 0;
						if(Frame::cosTheta(bRec.wi) < 0) {
							// Assumption: twosided => same material
							startID = bsdf->getComponentCount() / 2;
						}

						compSel = -1;
						for(int comp = startID; comp < startID+bsdf->getComponentCount()/2 && compSel == -1; comp++) {
							if(bsdf->getRoughness(its, comp) < m_roughness) {
								compSel = comp;
							}
						}

						if(compSel == -1)
							SLog(EError, "What !?");

						bRec.component = compSel;
						pdfComp = 1.f;
					} else {
						DiscreteDistribution gloDist(nbGlossyComp);
						for(int comp = 0; comp < bsdf->getComponentCount(); comp++) {
							if (bsdf->getRoughness(its, 0) < m_roughness) {
								bRec.component = comp;
								pdfComp = bsdf->pdfComponent(bRec);
								gloDist.append(pdfComp);
							} else {
								gloDist.append(0.f);
							}
						}
						gloDist.normalize();
						compSel = gloDist.sample(rRec.nextSample1D(), pdfComp);
					}
				}
			}

			bRec.component = compSel;
			Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
			if (bsdfWeight.isZero())
				break;
			bsdfWeight /= pdfComp;

            if(m_showGlossySurfaces) {
              Li += bsdfWeight;
              break;
            }

			scattered |= bRec.sampledType != BSDF::ENull;

			/* Prevent light leaks due to the use of shading normals */
			const Vector wo = its.toWorld(bRec.wo);
			Float woDotGeoN = dot(its.geoFrame.n, wo);
			if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
				break;


			/* Trace a ray in this direction */
			ray = Ray(its.p, wo, ray.time);
			if (!scene->rayIntersect(ray, its)) {
				/* Intersected nothing -- perhaps there is an environment map? */
				const Emitter *env = scene->getEnvironmentEmitter();

				if (env) {
					SLog(EError, "Env map not supported yet");
				}

				break;
			}

			/* Keep track of the throughput and relative
			   refractive index along the path */
			throughput *= bsdfWeight;
			eta *= bRec.eta;


			/* ==================================================================== */
			/*                         Indirect illumination                        */
			/* ==================================================================== */

			/* Set the recursive query type. Stop if no surface was hit by the
			   BSDF sample or if indirect illumination was not requested */
			if (!its.isValid() || !(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
				break;
			rRec.type = RadianceQueryRecord::ERadianceNoEmission;

			if (rRec.depth++ >= m_rrDepth) {
				/* Russian roulette: try to keep path weights equal to one,
				   while accounting for the solid angle compression at refractive
				   index boundaries. Stop with at least some probability to avoid
				   getting stuck (e.g. due to total internal reflection) */

				Float q = std::min(throughput.max() * eta * eta, (Float) 0.95f);
				if (rRec.nextSample1D() >= q)
					break;
				throughput /= q;
			}
		}

		/* Store statistics */
		avgPathLength.incrementBase();
		avgPathLength += rRec.depth;

		return Li;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		MonteCarloIntegrator::serialize(stream, manager);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "MIPathTracer[" << endl
			<< "  maxDepth = " << m_maxDepth << "," << endl
			<< "  rrDepth = " << m_rrDepth << "," << endl
			<< "  strictNormals = " << m_strictNormals << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()

	private:
	  Float m_roughness;
      bool m_showGlossySurfaces;
};

MTS_IMPLEMENT_CLASS_S(MIPathTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(MIPathTracer, "MI path tracer");
MTS_NAMESPACE_END
