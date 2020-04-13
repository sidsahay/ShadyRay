//
// Created by walksbynight on 23/3/18.
//

#include "sampler.hpp"


Sampler::Sampler(real scale, int samples_per_pixel)
        : scale(scale), samples_per_pixel(samples_per_pixel) {
    current_sample_count = samples_per_pixel;
}

void Sampler::StartPixel() {
    current_sample_count = 0;
}

bool Sampler::CheckPixelDone() {
    if (current_sample_count >= samples_per_pixel) {
        current_sample_count = 0;
        return false;
    } else {
        return true;
    }
}

real Sampler::GenerateRawSample() {
    return (real) 0.5;
}

real Sampler::Get1D() {
    current_sample_count++;
    return GenerateRawSample() * scale;
}

Point2r Sampler::Get2D() {
    return {Get1D(), Get1D()};
}

UniformSampler::UniformSampler(real scale, int samples_per_pixel)
        : Sampler(scale, samples_per_pixel) {
}

real UniformSampler::GenerateRawSample() {
    return distro(gen) - (real) 0.5;
}

real UniformSampler::Get1D() {
    current_sample_count++;
    return GenerateRawSample() * scale;
}

Point2r UniformSampler::Get2D() {
    current_sample_count += 2;
    return {distro(gen), distro(gen)};
}
