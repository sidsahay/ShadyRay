//
// Created by walksbynight on 23/3/18.
//

#ifndef SHADYRAY_FILTER_HPP
#define SHADYRAY_FILTER_HPP

#include "defs.hpp"

#include <cmath>

class Filter {
public:
    explicit Filter(real radius = (real) 2.);

    virtual real Evaluate(real x, real y);

protected:
    real radius;
};

class MitchellNetravaliFilter : public Filter {
public:
    MitchellNetravaliFilter(real b, real c, real radius);

    real Evaluate(real x, real y) override;

protected:
    real Mitchell1D(real x);

    real B;
    real C;
};


#endif //SHADYRAY_FILTER_HPP
