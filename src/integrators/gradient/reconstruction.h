//
// Created by muliana on 10/4/18.
//

#ifndef MITSUBA_RECONSTRUCTION_H
#define MITSUBA_RECONSTRUCTION_H

#include <map>

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/plugin.h>
// For L2 and L1 reconstruction
#include "poisson_solver/Solver.hpp"

MTS_NAMESPACE_BEGIN
typedef std::vector<Vector3> Pixels;
void run_conjudate_gradient(const Vector2i imgSize,
        const std::vector<float>& dxVar,
        const std::vector<float>& dyVar) {
    // TODO: Implementation of conjugate gradient solver of PG

    size_t n = imgSize.x * imgSize.y;
    auto coord = [&](size_t i) -> Point2i {
        return Point2i(i % imgSize.x, i / imgSize.x);
    };
    auto mult = [](Vector3 v1, Vector3 v2) -> Vector3 {
        return Vector3(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
    };

    // Function for computing...
    auto compute_px = [&](const Pixels& x) -> Pixels {
        std::vector<Vector3> e(n * 2);
        for(size_t i = 0; i < n; i++) {
            auto c = coord(i);
            e[i] = (c.x != imgSize.x - 1) ? x[i+1] - x[i] : Vector3(0.0);
            e[i + n] = (c.y != imgSize.y - 1) ? x[i + imgSize.x] - x[i] : Vector3(0.0);
        }
        return e;
    };
    auto compute_axpy = [&](Pixels& e, const Pixels& b) {
        for(size_t i = 0; i < n; i++) {
            e[i] = b[i] - e[i];
            e[i + n] = b[i + n] - e[i + n];
        }
    };
    auto compute_PTW2x = [&](const std::vector<float>& w2, const Pixels& e) -> Pixels {
        Pixels r(n);
        for(size_t i = 0; i < n; i++) {
            Vector3 val(0.0);
            auto c = coord(i);
            if(c.x != 0) val += e[i - 1] * w2[i - 1];
            if(c.x != imgSize.x - 1) val -= e[i] * w2[i];
            if(c.y != 0) val += e[n + i - imgSize.x] * w2[n + i - imgSize.x];
            if(c.y != imgSize.y - 1) val -= e[n + i] * w2[n + i] ;
            r[i] = val;
        }
    };
    auto compute_Ax_xAx = [&](Pixels& Ap, Pixels& xAp,
            const std::vector<float>& w2, const Pixels& p) {
        for(size_t i = 0; i < n; i++) {
            auto c = coord(i);
            Vector3 Axi(0.0);
            Vector3 xi = p[i];
            if(c.x != 0) Axi += (xi - p[i - 1]) * w2[i-1];
            if(c.x != imgSize.x - 1) Axi += (xi - p[i+1]) * w2[i];
            if(c.y != 0) Axi += (xi - p[i - imgSize.x]) * w2[n + i - imgSize.x];
            if(c.y != imgSize.y - 1) Axi += (xi - p[i + imgSize.x]) * w2[n  + i];
            Ap[i] = Axi;
            xAp[i] = mult(Axi,xi);
        }
    };
    auto compute_r_rz = [&](Pixels& r, const Pixels& Ap, const Vector3& a) {
        for(size_t i = 0; i < n; i++) {
            r[i] = r[i] - mult(Ap[i],a);
        }
    };

    // Compute the weights
    auto w2 = [imgSize, n](const std::vector<float>& dxVar, const std::vector<float>& dyVar) -> std::vector<float> {
        auto compute_weight = [&](float r, float g, float b) -> float {
            return 1.f / (std::max(0.f, r + g + b) + 1e-4f);
        };
        std::vector<float> w2(n  * 2);
        for(size_t i = 0; i < n; i++) {
            w2[i] = compute_weight(dxVar[i * 3],dxVar[i * 3 + 1],dxVar[i * 3 + 2]);
            w2[i + n] = compute_weight(dyVar[i * 3],dyVar[i * 3 + 1],dyVar[i * 3 + 2]);
        }
        return w2;
    }(dxVar, dyVar);

    // Scale the weights
    //TODO

}

struct PostProcessOption {
	bool forceBlackPixels;
	bool clampingValues;

	ref <Bitmap> process(std::vector<float> &rec,
						 std::vector<float> &primal,
						 std::vector<float> &very_direct, const Vector2i imgSize) {
		// Clean the data
		if (forceBlackPixels) {
			size_t len = primal.size();
			for (size_t p = 0; p < len; p += 3) {
				if (primal[p] == 0.f && primal[p + 1] == 0.f && primal[p + 2] == 0.f &&
					very_direct[p] == 0.f && very_direct[p + 1] == 0.f && very_direct[p + 2] == 0.f) {
					rec[p] = rec[p + 1] = rec[p + 2] = 0.f;
				}
			}
		}

		// Export the data inside a bitmap
		ref<Bitmap> bitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, imgSize);
		int x, y;
		Float tmp[3];
		for (int i = 0; i < bitmap->getSize().x * bitmap->getSize().y; i++) {
			y = i / bitmap->getSize().x;
			x = i - y * bitmap->getSize().x;
			if (clampingValues) {
				tmp[0] = std::max(Float(rec[3 * i] + very_direct[3 * i]), Float(0.f));
				tmp[1] = std::max(Float(rec[3 * i + 1] + very_direct[3 * i + 1]), Float(0.f));
				tmp[2] = std::max(Float(rec[3 * i + 2] + very_direct[3 * i + 2]), Float(0.f));
			} else {
				tmp[0] = Float(rec[3 * i] + very_direct[3 * i]);
				tmp[1] = Float(rec[3 * i + 1] + very_direct[3 * i + 1]);
				tmp[2] = Float(rec[3 * i + 2] + very_direct[3 * i + 2]);
			}
			bitmap->setPixel(Point2i(x, y), Spectrum(tmp));
		}
		return bitmap;
	}
};

struct Reconstruction {
	bool reconstructL1 = true;
	bool reconstructL2 = false;
	bool reconstructUni = false;
	bool reconstructWeighted = false;
	bool reconstructL2Weighted = false;

	// Parameters for reconstruction
	float alpha = 0.2;
	float nbIteration = 50;

	struct Result {
		std::string name;
		ref<Bitmap> img;
	};
	struct Variance {
		Variance() = default;
		std::vector<float> primal;
		std::vector<float> dx;
		std::vector<float> dy;
	};

	std::vector<Result> reconstruct(const Vector2i imgSize,
									 std::vector<float> &primal, std::vector<float> &dx, std::vector<float> &dy,
									 std::vector<float> &very_direct,
									 Variance variance,
									 PostProcessOption options) {
		std::vector<Result> results;
		if(reconstructL2) {
			poisson::Solver::Params pL2;
			pL2.setConfigPreset("L2D");
			pL2.alpha = alpha;
			poisson::Solver solverL2(pL2);

			/* apply L2 solver to reconstruct L2 image and store result in recBuff */
			solverL2.importImagesMTS(&dx[0], &dy[0], &primal[0], imgSize.x, imgSize.y);
			solverL2.setupBackend();
			solverL2.solveIndirect();

			auto len = size_t(3 * imgSize.x * imgSize.y);
			std::vector<float> rec(len, 0.f);
			solverL2.exportImagesMTS(&rec[0]);
			results.push_back(Result {
				name: "L2",
				img: options.process(rec, primal, very_direct, imgSize)
			});
		}
		if(reconstructL1) {
			poisson::Solver::Params pL1;
			pL1.setConfigPreset("L1D");
			pL1.alpha = alpha;
			poisson::Solver solverL1(pL1);

			/* apply L1 solver to reconstruct L1 image and store result in recBuff */
			solverL1.importImagesMTS(&dx[0], &dy[0], &primal[0], imgSize.x, imgSize.y);
			solverL1.setupBackend();
			solverL1.solveIndirect();

			auto len = size_t(3 * imgSize.x * imgSize.y);
			std::vector<float> rec(len, 0.f);
			solverL1.exportImagesMTS(&rec[0]);
			results.push_back(Result {
					name: "L1",
					img: options.process(rec, primal, very_direct, imgSize)
			});
		}
		if(reconstructUni) {
			auto len = size_t(3 * imgSize.x * imgSize.y);
			std::vector<std::vector<float>> rec{std::vector<float>(len, 0.f), std::vector<float>(len, 0.f)};
			rec[0].assign(primal.begin(), primal.end());

			int dst;
			for (int iter = 0; iter < nbIteration; iter++) {
				int src = iter % 2;
				dst = 1 - src;

				const int w = imgSize.x;
				const int h = imgSize.y;

#if defined(MTS_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
				for (int i = 0; i < w * h; i++) {
					int x = i % w;
					int y = i / w;
					if (x >= w || y >= h) continue;

					for (int channel = 0; channel < 3; channel++) {
						float color = 0.f;
						float weight = 0.f;
						const float &prim = rec[src][(y * w + x) * 3 + channel];
						color += prim; // GVCM multiply by alpha_sqr
						weight += 1.0; // GVCM multiply by alpha_sqr
						if (x > 0) {
							color += dx[(y * w + x - 1) * 3 + channel];
							color += rec[src][(y * w + x - 1) * 3 + channel];
							weight += 1.f;
						}
						if (x + 1 < w) {
							color -= dx[(y * w + x) * 3 + channel];
							color += rec[src][(y * w + x + 1) * 3 + channel];
							weight += 1.f;
						}
						if (y > 0) {
							color += dy[((y - 1) * w + x) * 3 + channel];
							color += rec[src][((y - 1) * w + x) * 3 + channel];
							weight += 1.f;
						}
						if (y + 1 < h) {
							color -= dy[(y * w + x) * 3 + channel];
							color += rec[src][((y + 1) * w + x) * 3 + channel];
							weight += 1.f;
						}
						rec[dst][(y * w + x) * 3 + channel] = color / weight;
					}
				}
			}
			results.push_back(Result {
				name: "Uni",
				img: options.process(rec[dst], primal, very_direct, imgSize)
			});
		}
		if(reconstructWeighted) {
			auto len = size_t(3 * imgSize.x * imgSize.y);
			std::vector<std::vector<float>> rec{std::vector<float>(len, 0.f), std::vector<float>(len, 0.f)};
			rec[0].assign(primal.begin(), primal.end());


			float variance_reduction_denominator = 1.f;
			float pow_iter = 1.0;
			const float epsilon = 1e-4;

			int dst;
			for (int iter = 0; iter < nbIteration; iter++) {
				variance_reduction_denominator *= (0.01 + 1.f + 4.f * pow_iter);
				pow_iter *= 0.5;
				float variance_reduction_weight = 1.f / variance_reduction_denominator;

				int src = iter % 2;
				dst = 1 - src;

				const int w = imgSize.x;
				const int h = imgSize.y;

#if defined(MTS_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
				for (int i = 0; i < w * h; i++) {
					int x = i % w;
					int y = i / w;
					if (x >= w || y >= h) continue;

					for (int channel = 0; channel < 3; channel++) {
						float color = 0.f;
						float weight = 0.f;
						const float &prim = rec[src][(y * w + x) * 3 + channel];
						const float &prim_var = variance.primal[(y * w + x) * 3 + channel] * variance_reduction_weight;
						{
							float curr_weight = 1.f / (epsilon + prim_var);
							color += prim * curr_weight; // GVCM multiply by alpha_sqr
							weight += curr_weight; // GVCM multiply by alpha_sqr
						}
						if (x > 0) {
							float curr_weight = 1.f / (epsilon + prim_var + variance.dx[(y * w + x - 1) * 3 + channel]);
							color += curr_weight*dx[(y * w + x - 1) * 3 + channel];
							color += curr_weight*rec[src][(y * w + x - 1) * 3 + channel];
							weight += curr_weight;
						}
						if (x + 1 < w) {
							float curr_weight = 1.f / (epsilon + prim_var + variance.dx[(y * w + x) * 3 + channel]);
							color -= curr_weight*dx[(y * w + x) * 3 + channel];
							color += curr_weight*rec[src][(y * w + x + 1) * 3 + channel];
							weight += curr_weight;
						}
						if (y > 0) {
							float curr_weight = 1.f / (epsilon + prim_var + variance.dy[((y - 1) * w + x) * 3 + channel]);
							color += curr_weight*dy[((y - 1) * w + x) * 3 + channel];
							color += curr_weight*rec[src][((y - 1) * w + x) * 3 + channel];
							weight += curr_weight;
						}
						if (y + 1 < h) {
							float curr_weight = 1.f / (epsilon + prim_var + variance.dy[(y * w + x) * 3 + channel]);
							color -= curr_weight*dy[(y * w + x) * 3 + channel];
							color += curr_weight*rec[src][((y + 1) * w + x) * 3 + channel];
							weight += curr_weight;
						}
						rec[dst][(y * w + x) * 3 + channel] = color / weight;
					}
				}
			}
			results.push_back(Result {
					name: "Weighted",
					img: options.process(rec[dst], primal, very_direct, imgSize)
			});
		}
		if(reconstructL2Weighted) {
			poisson::Solver::Params pL2;
			pL2.setConfigPreset("L2D");
			pL2.alpha = alpha;
			poisson::Solver solverL2(pL2);

			/* apply L2 solver to reconstruct L2 image and store result in recBuff */
			solverL2.importImagesMTS(&dx[0], &dy[0], &primal[0], imgSize.x, imgSize.y);
			solverL2.setupBackend();
			if(alpha == 0.0) {
                solverL2.setupWeights(nullptr, &variance.dy[0], &variance.dx[0], alpha);
			} else {
                solverL2.setupWeights(&variance.primal[0], &variance.dy[0], &variance.dx[0], alpha);
			}

			solverL2.solveIndirect();

			auto len = size_t(3 * imgSize.x * imgSize.y);
			std::vector<float> rec(len, 0.f);
			solverL2.exportImagesMTS(&rec[0]);
			results.push_back(Result {
					name: "L2Weighted",
					img: options.process(rec, primal, very_direct, imgSize)
			});
		}

		// Rename the vector if needed
		if(results.size() == 0) {
			SLog(EError, "No reconstruction?");
		} else if(results.size() == 1) {
			results[0].name = "recons";
		}
		return results;
	}
};

std::vector<float> bitmap2vec(Bitmap *bitmap) {
    auto subPixelCount = size_t(3 * bitmap->getSize().x * bitmap->getSize().y);
    std::vector<float> vec(subPixelCount);
    std::transform(bitmap->getFloatData(),
                   bitmap->getFloatData() + subPixelCount,
                   vec.begin(),
                   [](Float x) { return (float) x; });
    return vec;
}

ref<Film> newHDRFilm(Scene *scene) {
    Vector2i size = scene->getFilm()->getSize();
    Properties pHDRFilm("hdrfilm");
    pHDRFilm.setInteger("width", size.x);
    pHDRFilm.setInteger("height", size.y);
    pHDRFilm.setString("fileFormat", "rgbe");
    ref<Film> hdrFilm = dynamic_cast<Film *> (PluginManager::getInstance()->
            createObject(MTS_CLASS(Film), pHDRFilm));
    return hdrFilm;
}

MTS_NAMESPACE_END


#endif //MITSUBA_RECONSTRUCTION_H
