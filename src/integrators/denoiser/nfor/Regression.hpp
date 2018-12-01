#ifndef REGRESSION_HPP_
#define REGRESSION_HPP_

#include "Pixmap.hpp"

#include <vector>

namespace Tungsten {

Pixmap3f collaborativeRegression(const Pixmap3f &image, const Pixmap3f &guide,
        const std::vector<PixmapF> &features, const Pixmap3f &imageVariance,
        int F, int R, float k, float varianceScale = 2.0);

}

#endif /* REGRESSION_HPP_ */
