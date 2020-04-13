//
// Created by walksbynight on 23/3/18.
//

#ifndef SHADYRAY_SAMPLER_HPP
#define SHADYRAY_SAMPLER_HPP

#pragma once

#include "linalg.hpp"
#include <random>

//Generates samples from -1 to +1
class Sampler {
public:
    Sampler(real scale = (real) 2., int samples_per_pixel = 2);

    void StartPixel();

    bool CheckPixelDone();

    real GenerateRawSample(); //before scaling
    real Get1D(); //scales and decrements sample count
    Point2r Get2D(); //also scales and decrements sample count

protected:
    real scale;
    int samples_per_pixel;
    int current_sample_count = 0;
};

class UniformSampler : public Sampler {
public:
    UniformSampler(real scale = (real) 2., int samples_per_pixel = 2);

    real GenerateRawSample();

    real Get1D();

    Point2r Get2D();

protected:
    std::minstd_rand gen;
    std::uniform_real_distribution<real> distro;
};

#endif //SHADYRAY_SAMPLER_HPP
