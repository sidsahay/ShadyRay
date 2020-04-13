//
// Created by walksbynight on 23/3/18.
//

#include "filter.hpp"


Filter::Filter(real radius) : radius(radius) {
}

real Filter::Evaluate(real x, real y) {
    return (real) 1.;
}

MitchellNetravaliFilter::MitchellNetravaliFilter(real b, real c, real radius)
        : Filter(radius), B(b), C(c) {
}

real MitchellNetravaliFilter::Evaluate(real x, real y) {
    return Mitchell1D(x / radius) * Mitchell1D(y / radius);
}

real MitchellNetravaliFilter::Mitchell1D(real x) {
    x = std::abs((real) 2. * x);
    if (x > (real) 1.) {
        return ((-B - 6 * C) * x * x * x + (6 * B + 30 * C) * x * x + (-12 * B - 48 * C) * x + (8 * B + 24 * C)) *
               (1.0f / 6.0f);
    } else {
        return ((12 - 9 * B - 6 * C) * x * x * x + (-18 + 12 * B + 6 * C) * x * x + (6 - 2 * B)) * (1.0f / 6.0f);
    }
}
