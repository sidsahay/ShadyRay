//
// Created by walksbynight on 23/3/18.
//

#ifndef SHADYRAY_PATHINTEGRATORTH_HPP
#define SHADYRAY_PATHINTEGRATORTH_HPP

#pragma once

#include "pathintegrator.hpp"
#include <tbb/tbb.h>
#include "particletracer.hpp"

#define PATH_INTEGRATOR_TH_T PathIntegratorTh<SceneT, CameraT, FilmT, SamplerT, FilterT>

template<class SceneT, class CameraT, class FilmT, class SamplerT, class FilterT>
class PathIntegratorTh : public PATH_INTEGRATOR_T {
public:
    PathIntegratorTh(int num_x_div, int num_y_div, int max_bounces, SceneT *scene, CameraT *camera, FilmT *film,
                     SamplerT *sampler, FilterT *filter);

    void Render() override;

    void SetGlow(const Spectrum& glow);

    ParticleTracer *tracer = nullptr;

protected:
    int num_x_div;
    int num_y_div;
    Spectrum glow_color;
};

template<class SceneT, class CameraT, class FilmT, class SamplerT, class FilterT>
PATH_INTEGRATOR_TH_T::PathIntegratorTh(int num_x_div, int num_y_div, int max_bounces, SceneT *scene, CameraT *camera,
                                       FilmT *film, SamplerT *sampler, FilterT *filter)
        : num_x_div(num_x_div), num_y_div(num_y_div),
          PATH_INTEGRATOR_T(max_bounces, scene, camera, film, sampler, filter) {
}

template<class SceneT, class CameraT, class FilmT, class SamplerT, class FilterT>
void PATH_INTEGRATOR_TH_T::Render() {
    const auto film = this->film;
    auto sampler = this->sampler;
    const auto camera = this->camera;
    const auto filter = this->filter;
    const auto scene = this->scene;
    const auto lights = this->scene->lights;

    int x_unit = film->width / num_x_div;
    int y_unit = film->height / num_y_div;

    auto *tiles = new FilmTile[num_x_div * num_y_div];
    auto *samplers = new SamplerT[num_x_div * num_y_div];

    for (int j = 0; j < num_y_div; ++j) {
        int offset = j * num_x_div;
        for (int i = 0; i < num_x_div; ++i) {
            int idx = offset + i;

            tiles[idx].Prepare(i * x_unit, j * y_unit, x_unit, y_unit);
            new(samplers + idx) SamplerT(*sampler);
        }
    }

    auto rr_probability = (real) 0.8;

    tbb::parallel_for(0, num_x_div * num_y_div, [&](int idx) {
        auto &tile = tiles[idx];
        auto &tile_sampler = samplers[idx];

        real height_actual = film->height;
        real width_actual = film->width;

        for (int j = 0; j < tile.height; ++j) {

            for (int i = 0; i < tile.width; ++i) {
                real j_actual = j + tile.y_offset;
                real i_actual = i + tile.x_offset;

                Spectrum numerator;
                Spectrum denominator;

                tile_sampler.StartPixel();
                int *path_type_array = new int[this->max_bounces];
                for (int s = 0; s < this->max_bounces; s++)
                    path_type_array[s] = 0;

                do {
                    real sample_x = tile_sampler.Get1D();
                    real sample_y = tile_sampler.Get1D();

                    real filter_factor = filter->Evaluate(sample_x, sample_y);
                    denominator = denominator + Spectrum(filter_factor);

                    Point3r image_point(i_actual + sample_x - width_actual / (real) 2.,
                                        j_actual + sample_y - height_actual / (real) 2., (real) 0.);

                    Ray ray = camera->GenerateRay(image_point);

                    Spectrum radiance;
                    Spectrum beta((real) 1.);

                    bool first_was_lamp = false;
                    unsigned int mask = scene->mask_enable_all;

                    for (int bounces = 0;; ++bounces) {
                        Interaction interaction;

                        bool intersected = scene->Intersect(ray, &interaction, mask);

                        //Handle emission on camera ray
                        if (bounces == 0) {
                            if (intersected) {
                                if (interaction.geometry_id == scene->lamp_geometry_id) {
                                    first_was_lamp = true;
                                    //radiance = radiance + beta * scene->geometries[interaction.geometry_id].material->emission_color; // Add object Le if intersection found
                                }
                                else {
                                    mask = scene->mask_enable_others;
                                }
                            }
                            else {
                                radiance = radiance + beta * Spectrum(0.0); // Add environmental lighting Le
                            }
                        }

                        //Exit condition
                        if (!intersected || bounces >= this->max_bounces) {
                            break;
                        }

                        Material *material = scene->geometries[interaction.geometry_id].material;
                        Spectrum direct((real) 0.);

                        if (first_was_lamp) {
                            if (bounces >= 5) {
                                if (interaction.geometry_id == scene->lamp_geometry_id) {
                                    auto glow = glow_color * tracer->GetBoxValue(interaction.point);
                                    radiance.r = std::max(radiance.r, glow.r);
                                    radiance.g = std::max(radiance.g, glow.g);
                                    radiance.b = std::max(radiance.b, glow.b);

                                    break;
                                }
                            }
                        }
                        //Get light values. Same as LolIntegrator. This is the local direct lighting integral
                        interaction.normal = Normalize(interaction.normal);

                        BSDF *bsdf = material->GetBSDF(interaction, (real) 1.0);

                        Point2r u = tile_sampler.Get2D();
                        bsdf->Choose(u);

                        auto bias_factor = (bounces <= 10) ? (real) 1. : (real) 1. / (1. - rr_probability);

                        //Choose randomlight
                        int num_lights = scene->lights.size();
                        if (num_lights != 0) {
                            int index = num_lights * u.data[1];
                            {
                                const auto &light = scene->lights[index];
                                Ray light_ray;
                                VisibilityTester tester;
                                Spectrum light_value = light->Evaluate(interaction, &light_ray, &tester);

                                Ray shadow_ray;
                                shadow_ray.origin = tester.point_0.point;
                                shadow_ray.direction = tester.point_1.point - tester.point_0.point;

                                if (!scene->IntersectCheck(shadow_ray, mask)) {
                                    real lambert_factor = Clamp(
                                            Dot(Vector3r(Normalize(bsdf->ns)),
                                                Normalize(shadow_ray.direction)), 0.f, 1.f);

                                    auto f = bsdf->F(Vector3r(Normalize(ray.direction)) * -1,
                                                     Normalize(shadow_ray.direction));
                                    direct = direct + light_value * f * lambert_factor;
                                }
                            }
                        }


                        direct = direct * num_lights + /*material->emission_color ;*/ (first_was_lamp && bounces < 3 ? material->emission_color : Spectrum()); //+ (bounces > 1 ? glow_color * tracer->GetBoxValue(interaction.point) * 0.8 : Spectrum());

                        radiance = radiance + beta * direct * bias_factor;

                        //Stop the integration process via Russian Roulette
                        Point2r stop_sample = tile_sampler.Get2D();
                        if (bounces > 10 && stop_sample.data[0] < rr_probability) {
                            break;
                        }

                        //Generate new ray for path integral. Also update beta.
                        Vector3r wi, wo = ray.direction * -1;
                        real pdf;

                        Spectrum f = bsdf->SampleF(wo, &wi, tile_sampler.Get2D(), &pdf);

                        path_type_array[bounces] = bsdf->bxdfs[bsdf->chosen_index]->type;

                        delete bsdf;

                        if (f.IsBlack() || pdf == 0) {
                            break;
                        }

                        beta = beta * f * std::abs(Dot(wi, Vector3r(Normalize(interaction.normal)))) * (1 / pdf);

                        ray.origin = interaction.point;
                        ray.direction = wi;

                    }

                    numerator = numerator + radiance * filter_factor;
                } while (tile_sampler.CheckPixelDone());

                delete path_type_array;

                Spectrum filtered_value = numerator * denominator.MulInverse();
                tile.SetPixel(i, j, filtered_value);
            }
        }
    });


    for (int i = 0; i < num_x_div * num_y_div; ++i) {
        *film << tiles[i];
    }
}

template<class SceneT, class CameraT, class FilmT, class SamplerT, class FilterT>
void PATH_INTEGRATOR_TH_T::SetGlow(const Spectrum &glow) {
    this->glow_color = glow;
}

#endif //SHADYRAY_PATHINTEGRATORTH_HPP
