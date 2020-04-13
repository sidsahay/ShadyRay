//
// Created by walksbynight on 23/3/18.
//

#ifndef SHADYRAY_PATHINTEGRATOR_HPP
#define SHADYRAY_PATHINTEGRATOR_HPP

#pragma once

#include "integrator.hpp"

#define PATH_INTEGRATOR_T PathIntegrator<SceneT, CameraT, FilmT, SamplerT, FilterT>

template<class SceneT, class CameraT, class FilmT, class SamplerT, class FilterT>
class PathIntegrator : public INTEGRATOR_T {
public:
    PathIntegrator(int max_bounces, SceneT *scene, CameraT *camera, FilmT *film, SamplerT *sampler, FilterT *filter);

    void Render() override;

protected:
    int max_bounces;
};


template<class SceneT, class CameraT, class FilmT, class SamplerT, class FilterT>
PATH_INTEGRATOR_T::PathIntegrator(int max_bounces, SceneT *scene, CameraT *camera, FilmT *film, SamplerT *sampler,
                                  FilterT *filter)
        : max_bounces(max_bounces), INTEGRATOR_T(scene, camera, film, sampler, filter) {
}

template<class SceneT, class CameraT, class FilmT, class SamplerT, class FilterT>
void PATH_INTEGRATOR_T::Render() {
/*    auto film = this->film;
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
            BeckmannDistribution distro_smooth(0.01);
            Spectrum R_smooth(0.5, 0.5, 0.5);
            MicrofacetReflection bxdf_smooth(R_smooth, &distro_smooth, &fres);

            BeckmannDistribution distro_rough(0.9);
            Spectrum R_rough(4, 3, 2);
            MicrofacetReflection bxdf_rough(R_rough, &distro_rough, &fres);

            do {
                real sample_x = sampler->Get1D();
                real sample_y = sampler->Get1D();

                real filter_factor = filter->Evaluate(sample_x, sample_y);
                denominator = denominator + Spectrum(filter_factor);

                Point3r image_point(i + sample_x - width / (real) 2., j + sample_y - height / (real) 2., (real) 0.);

                Ray ray = camera->GenerateRay(image_point);

                Spectrum radiance;
                Spectrum beta((real) 1.);

                for (int bounces = 0;; ++bounces) {
                    Interaction interaction;

                    bool intersected = scene->Intersect(ray, &interaction);

                    //Handle emission on camera ray
                    if (bounces == 0) {
                        if (intersected) {
                            radiance = radiance + beta * Spectrum(0.); // Add object Le if intersection found
                        } else {
                            radiance = radiance + beta * Spectrum(0.3); // Add environmental lighting Le
                        }
                    }

                    //Exit condition
                    if (!intersected || bounces >= this->max_bounces) {
                        break;
                    }

                    //Get light values. Same as LolIntegrator. This is the local direct lighting integral
                    Spectrum direct((real) 0.);
                    interaction.normal = Normalize(interaction.normal);

                    BSDF bsdf(interaction); //Same materials for everything right now

                    bsdf.bxdf = (interaction.geometry_id != 0) ? &bxdf_smooth : &bxdf_rough;

                    for (auto light : scene->lights) {
                        Ray light_ray;
                        VisibilityTester tester;
                        Spectrum light_value = light->Evaluate(interaction, &light_ray, &tester);

                        Ray shadow_ray;
                        shadow_ray.origin = tester.point_0.point;
                        shadow_ray.direction = tester.point_1.point - tester.point_0.point;

                        if (!scene->IntersectCheck(shadow_ray)) {
                            real lambert_factor = std::abs(
                                    Dot(Vector3r(Normalize(interaction.normal)), Normalize(shadow_ray.direction)));

                            auto f = bsdf.F(Vector3r(Normalize(ray.direction)) * -1, Normalize(shadow_ray.direction));
                            direct = direct + light_value * f * lambert_factor;
                        }
                    }

                    radiance = radiance + beta * direct;

                    //Generate new ray for path integral. Also update beta.
                    Vector3r wi, wo = ray.direction * -1;
                    real pdf;

                    Spectrum f = bsdf.SampleF(wo, &wi, sampler->Get2D(), &pdf);

                    if (f.IsBlack() || pdf == 0) {
                        break;
                    }

                    beta = beta * f * std::abs(Dot(wi, Vector3r(Normalize(interaction.normal)))) * (1 / pdf);

                    ray.origin = interaction.point;
                    ray.direction = wi;
                }

                numerator = numerator + radiance * filter_factor;
            } while (sampler->CheckPixelDone());

            Spectrum filtered_value = numerator * denominator.MulInverse();
            tile.SetPixel(i, j, filtered_value);
        }
    }

    *film << tile;*/
}

#endif //SHADYRAY_PATHINTEGRATOR_HPP
