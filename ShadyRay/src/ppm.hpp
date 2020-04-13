//
// Created by walksbynight on 26/3/18.
//

#ifndef SHADYRAY_PPM_HPP
#define SHADYRAY_PPM_HPP

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cstring>

namespace ppm {

    struct Pixel {
        unsigned char red;
        unsigned char green;
        unsigned char blue;
    };

    struct PPMImage {
        int width;
        int height;
        int maxVal;
        Pixel *data;
    };

    PPMImage* ReadPPM(const char *source);
    bool WritePPM(PPMImage* image, const char *dest);

}
#endif //SHADYRAY_PPM_HPP
