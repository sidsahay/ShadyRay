//
// Created by walksbynight on 25/3/18.
//

#ifndef SHADYRAY_IMAGE_HPP
#define SHADYRAY_IMAGE_HPP

#include "ppm.hpp"
#include "spectrum.hpp"
#include "linalg.hpp"

class Image {
public:
    bool LoadImage(const char *path);

    Spectrum Evaluate(real u, real v);

    void StoreAndBlend(real u, real v, const Spectrum &value);

    void Multiply(real value);

    void Blur();

    ~Image();

    unsigned int width = 0;
    unsigned int height = 0;

    void LoadBlank(unsigned int w, unsigned int h);

protected:
    Spectrum *data = nullptr;
};

#endif //SHADYRAY_IMAGE_HPP
