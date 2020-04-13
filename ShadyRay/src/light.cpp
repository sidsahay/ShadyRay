//
// Created by walksbynight on 23/3/18.
//

#include "light.hpp"

Light::Light(Point3r position, Spectrum color)
        : position(position), color(color) {
}

Spectrum Light::Evaluate(const Interaction &hit_point, Ray *light_ray, VisibilityTester *tester) {
    light_ray->origin = hit_point.point;
    const Vector3r &v = position - hit_point.point;
    light_ray->direction = Normalize(v);

    real length = v.Length();

    tester->point_0 = hit_point;
    tester->point_1.point = position;

    return color *
           ((real) 1. / (attenuation_constant + attenuation_linear * length + attenuation_quadratic * length * length));
}

TextureProjectionPointLight::TextureProjectionPointLight(const Point3r &position, const Spectrum &color,
                                                         const char *path) : Light(position, color) {
    TextureMap *map = new CubeMap;
    image_texture = new ImageTexture(map, path);
    Bounds3D bounds;
    bounds.centre = position;
    image_texture->SetBounds(bounds);
}

Spectrum TextureProjectionPointLight::Evaluate(const Interaction &hit_point, Ray *light_ray, VisibilityTester *tester) {
    light_ray->origin = hit_point.point;
    const Vector3r &v = position - hit_point.point;
    light_ray->direction = Normalize(v);

    real length = v.Length();

    tester->point_0 = hit_point;
    tester->point_1.point = position;

    return image_texture->Evaluate(hit_point) *
           ((real) 1. / (attenuation_constant + attenuation_linear * length + attenuation_quadratic * length * length));
}


