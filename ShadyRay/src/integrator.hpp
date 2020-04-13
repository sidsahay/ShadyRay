//
// Created by walksbynight on 23/3/18.
//

#ifndef SHADYRAY_INTEGRATOR_HPP
#define SHADYRAY_INTEGRATOR_HPP

#pragma once

#include "linalg.hpp"
#include "ray.hpp"
#include "light.hpp"
#include "interaction.hpp"
#include "sampler.hpp"
#include "filter.hpp"
#include "scene.hpp"
#include "camera.hpp"
#include "film.hpp"
#include "material.hpp"

#define INTEGRATOR_T Integrator<SceneT, CameraT, FilmT, SamplerT, FilterT>
#define LOL_INTEGRATOR_T LolIntegrator<SceneT, CameraT, FilmT, SamplerT, FilterT>

template<class SceneT, class CameraT, class FilmT, class SamplerT, class FilterT>
class Integrator {
public:
    Integrator(SceneT *scene, CameraT *camera, FilmT *film, SamplerT *sampler, FilterT *filter);

    virtual void Render();

    FilmT *film;

protected:
    SceneT *scene;
    SamplerT *sampler;
    FilterT *filter;
    CameraT *camera;
};

template<class SceneT, class CameraT, class FilmT, class SamplerT, class FilterT>
class LolIntegrator : public INTEGRATOR_T {
public:
    LolIntegrator(SceneT *scene, CameraT *camera, FilmT *film, SamplerT *sampler, FilterT *filter);

    void Render() override;

};


template<class SceneT, class CameraT, class FilmT, class SamplerT, class FilterT>
INTEGRATOR_T::Integrator(SceneT *scene, CameraT *camera, FilmT *film, SamplerT *sampler, FilterT *filter)
        : scene(scene), camera(camera), film(film), sampler(sampler), filter(filter) {

}

template<class SceneT, class CameraT, class FilmT, class SamplerT, class FilterT>
void INTEGRATOR_T::Render() {
}

template<class SceneT, class CameraT, class FilmT, class SamplerT, class FilterT>
LOL_INTEGRATOR_T::LolIntegrator(SceneT *scene, CameraT *camera, FilmT *film, SamplerT *sampler, FilterT *filter)
        : INTEGRATOR_T(scene, camera, film, sampler, filter) {
}

template<class SceneT, class CameraT, class FilmT, class SamplerT, class FilterT>
void LOL_INTEGRATOR_T::Render() {
    /*auto film = this->film;
    auto sampler = this->sampler;
    auto camera = this->camera;
    auto filter = this->filter;
    auto scene = this->scene;

    FilmTile tile(0, 0, film->width, film->height);

    real width = film->width;
    real height = film->height;

    for (int j = 0; j < film->height; ++j) {

        for (int i = 0; i < film->width; ++i) {
            Spectrum numerator;
            Spectrum denominator;

            sampler->StartPixel();

            Fresnel fres;//(0.2, 1);
            BeckmannDistribution distro(0.9);
            Spectrum R(3, 1, 0);
            MicrofacetReflection bxdf(R, &distro, &fres);
            // AshikhminShirleyReflection ashk(Spectrum(0, 0, 3),
            //                                 Spectrum(0, 1, 0),
            //                                 &distro);

            do {
                real sample_x = sampler->Get1D();
                real sample_y = sampler->Get1D();

                real filter_factor = filter->Evaluate(sample_x, sample_y);
                denominator = denominator + Spectrum(filter_factor);

                Point3r image_point(i + sample_x - width / (real) 2., j + sample_y - height / (real) 2., (real) 0.);

                Ray camera_ray = camera->GenerateRay(image_point);

                Interaction interaction;

                if (scene->Intersect(camera_ray, &interaction)) {
                    Spectrum irradiance((real) 0.);
                    interaction.normal = Normalize(interaction.normal);

                    for (auto light : scene->lights) {
                        Ray light_ray;
                        VisibilityTester tester;
                        Spectrum light_value = light->Evaluate(interaction, &light_ray, &tester);

                        Ray shadow_ray;
                        shadow_ray.origin = tester.point_0.point;
                        shadow_ray.direction = tester.point_1.point - tester.point_0.point;

                        if (!scene->IntersectCheck(shadow_ray)) {
                            BSDF bsdf(interaction);
                            bsdf.bxdf = &bxdf;
                            real lambert_factor = std::abs(
                                    Dot(Vector3r(Normalize(interaction.normal)), Normalize(shadow_ray.direction)));

                            auto f = bsdf.F(Vector3r(Normalize(camera_ray.direction)) * -1,
                                            Normalize(shadow_ray.direction));
                            irradiance = irradiance + light_value * f * lambert_factor;
                        }
                    }
                    numerator = numerator +
                                (irradiance * filter_factor);
                } else {
                    numerator = numerator + Spectrum((real) 0.3 * filter_factor);
                }
            } while (sampler->CheckPixelDone());

            Spectrum filtered_value = numerator * denominator.MulInverse();
            tile.SetPixel(i, j, filtered_value);
        }
    }

    *film << tile;*/
}

#endif //SHADYRAY_INTEGRATOR_HPP
