//
// Created by walksbynight on 23/3/18.
//

#ifndef SHADYRAY_INTEGRATORTH_HPP
#define SHADYRAY_INTEGRATORTH_HPP

#include "integrator.hpp"

#include <tbb/tbb.h>

#define LOL_INTEGRATOR_TH_T LolIntegratorTh<SceneT, CameraT, FilmT, SamplerT, FilterT>

template<class SceneT, class CameraT, class FilmT, class SamplerT, class FilterT>
class LolIntegratorTh : public LOL_INTEGRATOR_T {
public:
    LolIntegratorTh(int x_div, int y_div, SceneT *scene, CameraT *camera, FilmT *film, SamplerT *sampler,
                    FilterT *filter);

    void Render() override;

    int num_x_div;
    int num_y_div;
};


template<class SceneT, class CameraT, class FilmT, class SamplerT, class FilterT>
LOL_INTEGRATOR_TH_T::LolIntegratorTh(int x_div, int y_div, SceneT *scene, CameraT *camera, FilmT *film,
                                     SamplerT *sampler, FilterT *filter)
        : num_x_div(x_div), num_y_div(y_div), LOL_INTEGRATOR_T(scene, camera, film, sampler, filter) {
}

template<class SceneT, class CameraT, class FilmT, class SamplerT, class FilterT>
void LOL_INTEGRATOR_TH_T::Render() {
    auto film = this->film;
    auto sampler = this->sampler;
    auto camera = this->camera;
    auto filter = this->filter;
    auto scene = this->scene;

    int x_unit = film->width / num_x_div;
    int y_unit = film->height / num_y_div;

    FilmTile *tiles = new FilmTile[num_x_div * num_y_div];
    SamplerT *samplers = new SamplerT[num_x_div * num_y_div];

    for (int j = 0; j < num_y_div; ++j) {
        int offset = j * num_x_div;
        for (int i = 0; i < num_x_div; ++i) {
            int idx = offset + i;

            tiles[idx].Prepare(i * x_unit, j * y_unit, x_unit, y_unit);
            new(samplers + idx) SamplerT(*sampler);
        }
    }

    tbb::parallel_for(0, num_x_div * num_y_div, [&](int idx) {
        auto &tile = tiles[idx];
        auto &tile_sampler = samplers[idx];

        real width_actual = film->width;
        real height_actual = film->height;

        for (int j = 0; j < tile.height; ++j) {

            for (int i = 0; i < tile.width; ++i) {
                real j_actual = j + tile.y_offset;
                real i_actual = i + tile.x_offset;

                auto numerator = (real) 0.;
                auto denominator = (real) 0.;

                tile_sampler.StartPixel();

                do {
                    real sample_x = tile_sampler.Get1D();
                    real sample_y = tile_sampler.Get1D();

                    real filter_factor = filter->Evaluate(sample_x, sample_y);
                    denominator += filter_factor;

                    Point3r image_point(i_actual + sample_x - width_actual / (real) 2.,
                                        j_actual + sample_y - height_actual / (real) 2., 0.f);

                    Ray camera_ray = camera->GenerateRay(image_point);

                    Interaction interaction;

                    if (scene->Intersect(camera_ray, &interaction)) {
                        auto irradiance = (real) 1.;
                        numerator += irradiance * filter_factor;
                    } else {
                        numerator += (real) 0.3 * filter_factor;
                    }
                } while (tile_sampler.CheckPixelDone());

                real filtered_value = numerator / denominator;
                //filtered_value *= (real)idx / (num_x_div * num_y_div);

                tile.SetPixel(i, j, Spectrum(filtered_value));
            }
        }
    });

    for (int i = 0; i < num_x_div * num_y_div; ++i) {
        *film << tiles[i];
    }
}

#endif //SHADYRAY_INTEGRATORTH_HPP
