//
// Created by walksbynight on 25/3/18.
//

#include "image.hpp"

bool Image::LoadImage(const char *path) {
    ppm::PPMImage *image = ppm::ReadPPM(path);
    if (!image) {
        return false;
    }

    width = image->width;
    height = image->height;

    data = new Spectrum[width * height];

    for (unsigned int i = 0; i < width * height; ++i) {
        Spectrum image_pixel {(real)image->data[i].red, (real)image->data[i].green, (real)image->data[i].blue};
        data[i] = image_pixel * ((real)1. / 255);
    }

    delete[] image->data;
    delete image;

    return true;
}

Spectrum Image::Evaluate(real u, real v) {
    u = (u + (real)1.) / (real)2.;
    v = (v + (real)1.) / (real)2.;

    unsigned int x_coord = std::floor(u * width);
    unsigned int y_coord = std::floor(v * height);

    x_coord = x_coord <= 0 ? 0 : x_coord - 1;
    y_coord = y_coord <= 0 ? 0 : y_coord - 1;

    return data[y_coord * width + x_coord];
}

Image::~Image() {
    delete[] data;
}

void Image::StoreAndBlend(real u, real v, const Spectrum &value) {
    u = (u + (real)1.) / (real)2.;
    v = (v + (real)1.) / (real)2.;

    unsigned int x_coord = std::floor(u * width);
    unsigned int y_coord = std::floor(v * height);

    x_coord = x_coord <= 0 ? 0 : x_coord - 1;
    y_coord = y_coord <= 0 ? 0 : y_coord - 1;

    int idx = y_coord * width + x_coord;

    /*if (data[idx].IsBlack()) {
        data[idx] = value;
    }
    else {
        data[idx].r = (data[idx].r + value.r) / 2;
        data[idx].g = (data[idx].g + value.g) / 2;
        data[idx].b = (data[idx].b + value.b) / 2;
    }*/

    data[idx] = data[idx] + value;
}

void Image::LoadBlank(unsigned int w, unsigned int h){
    width = w;
    height = h;

    data = new Spectrum[width * height];
    std::memset(data, 0, sizeof(data));
}

void Image::Multiply(real value) {
    for (int y = 0; y < height; ++ y) {
        for (int x = 0; x < width; ++x) {
            int idx = y * width + x;
            data[idx] = data[idx] * value;
        }
    }
}

void Image::Blur() {
    Spectrum *new_data = new Spectrum[width * height];

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            Spectrum new_val;
            int count = 0;

            for (int new_y = y - 2; new_y < y + 2; ++new_y) {
                int actual_new_y = Clamp(new_y, 0, height - 1);

                for (int new_x = x - 2; new_x < x + 2; ++new_x) {
                    int actual_new_x = Clamp(new_x, 0, width - 1);

                    new_val = new_val + data[actual_new_y * width + actual_new_x];
                    count++;
                }
            }

            new_data[y * width + x] = new_val * (1.f / count);
        }
    }

    delete data;
    data = new_data;
}
