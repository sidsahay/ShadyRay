//
// Created by walksbynight on 23/3/18.
//

#include "spectrum.hpp"

std::ostream &operator<<(std::ostream &out, const Spectrum &pixel) {
    out << "[Sp " << pixel.r << ", " << pixel.g << ", " << pixel.b << "]";
    return out;
}

Spectrum::Spectrum(real v) : r(v), g(v), b(v) {
}

Spectrum::Spectrum(real r, real g, real b) : r(r), g(g), b(b) {
}

Spectrum Spectrum::operator+(const Spectrum &s) {
    return {r + s.r, g + s.g, b + s.b};
}

Spectrum Spectrum::operator-(const Spectrum &s) {
    return {r - s.r, g - s.g, b - s.b};
}

Spectrum Spectrum::operator*(real c) {
    return {c * r, c * g, c * b};
}

Spectrum Spectrum::operator*(const Spectrum &s) {
    return {r * s.r, g * s.g, b * s.b};
}

Spectrum Spectrum::MulInverse() {
    auto one = (real) 1.;
    return {one / r, one / g, one / b};
}

Spectrum Spectrum::Normalize() {
    real max = std::max({r, g, b});
    return {r / max, g / max, b / max};
}

bool Spectrum::IsBlack() {
    auto zero = (real) 0.;
    return r == zero && g == zero && b == zero;
}