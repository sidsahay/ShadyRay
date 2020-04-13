//
// Created by walksbynight on 23/3/18.
//

#ifndef SHADYRAY_LIGHT_HPP
#define SHADYRAY_LIGHT_HPP

#include "linalg.hpp"
#include "ray.hpp"
#include "interaction.hpp"
#include "spectrum.hpp"
#include "texture.hpp"

class VisibilityTester {
public:
    Interaction point_0;
    Interaction point_1;
};

class Light {
public:
    Light(Point3r position, Spectrum color);

    virtual Spectrum Evaluate(const Interaction &hit_point, Ray *light_ray, VisibilityTester *tester);

    real attenuation_constant = (real) 1.;
    real attenuation_linear = (real) 0.;
    real attenuation_quadratic = (real) 0.;

    Spectrum color;
    Point3r position;
};

class TextureProjectionPointLight : public Light{
public:
    TextureProjectionPointLight(const Point3r &position, const Spectrum &color, const char *path = nullptr);

    Spectrum Evaluate(const Interaction &hit_point, Ray *light_ray, VisibilityTester *tester);

    ImageTexture *image_texture;
};


#endif //SHADYRAY_LIGHT_HPP
