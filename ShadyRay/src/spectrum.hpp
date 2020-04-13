//
// Created by walksbynight on 23/3/18.
//

#ifndef SHADYRAY_SPECTRUM_HPP
#define SHADYRAY_SPECTRUM_HPP

#pragma once

#include "defs.hpp"
#include <algorithm>

class Spectrum {
public:
    real r;
    real g;
    real b;

    explicit Spectrum(real v = (real) 0.);

    Spectrum(real r, real g, real b);

    Spectrum operator+(const Spectrum &s);

    Spectrum operator-(const Spectrum &s);

    Spectrum operator*(real c);

    Spectrum operator*(const Spectrum &s);

    Spectrum MulInverse();

    Spectrum Normalize();

    bool IsBlack();
};

std::ostream &operator<<(std::ostream &out, const Spectrum &pixel);


#endif //SHADYRAY_SPECTRUM_HPP
