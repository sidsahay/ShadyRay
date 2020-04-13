//
// Created by walksbynight on 16/4/18.
//

#include <memory>
#include "particletracer.hpp"
#include "sampler.hpp"
#include "linalg.hpp"

ParticleTracer::ParticleTracer(Scene *scene) : scene(scene) {

}


void ParticleTracer::PreProcessEmissive() {
    for (auto& emissive_geometry_id : scene->emissive_geometry_ids) {
        //Copy buffers
        RTCScene temp_scene = rtcDeviceNewScene(scene->device_, RTC_SCENE_DYNAMIC, RTC_INTERSECT1);

        auto geometry = scene->geometries[emissive_geometry_id];

        unsigned int temp_geom_id = rtcNewTriangleMesh2(temp_scene, RTC_GEOMETRY_STATIC, geometry.num_faces,
                                                        geometry.num_vertices, 1);

        auto *vertices_dest = (Scene::Vertex_ *) rtcMapBuffer(temp_scene, temp_geom_id, RTC_VERTEX_BUFFER);
        auto *vertices_src = (Scene::Vertex_ *) rtcMapBuffer(scene->scene_, emissive_geometry_id,
                                                             RTC_VERTEX_BUFFER);
        for (int i = 0; i < geometry.num_vertices; ++i) {
            vertices_dest[i] = vertices_src[i];
        }
        rtcUnmapBuffer(scene->scene_, emissive_geometry_id, RTC_VERTEX_BUFFER);

        auto *tindices_dest = (Scene::TriangleIndices_ *) rtcMapBuffer(temp_scene, temp_geom_id, RTC_INDEX_BUFFER);
        auto *tindices_src = (Scene::TriangleIndices_ *) rtcMapBuffer(scene->scene_, emissive_geometry_id,
                                                                      RTC_INDEX_BUFFER);
        for (int i = 0; i < geometry.num_faces; ++i) {
            tindices_dest[i] = tindices_src[i];
        }
        rtcUnmapBuffer(scene->scene_, emissive_geometry_id, RTC_INDEX_BUFFER);

        //rtcCommit(temp_scene);

        //Buffers copied. Now to sample points on the surface via spherical sampling - because I don't have access to a better area sampler.

/*
        const int num_steps_theta = 32;
        const int num_steps_phi = 32;

        real theta_step = 3.141592 / num_steps_theta;

        real phi_step = 3.141592 / num_steps_phi;

        real diameter = (geometry.bounding_box.maximum - geometry.bounding_box.minimum).Length();


        //Sample points via spherical sampling. Will lead to aliasing, can't be helped.

        for (real theta = 0; theta <= 3.141592; theta += theta_step) {
            for (real phi = 0; phi <= 3.141592 * 2; phi += phi_step) {
                Ray sampling_ray;
                sampling_ray.origin = Point3r(SphericalDirection(std::sin(theta), std::cos(theta), phi) * diameter) +
                                      geometry.bounding_box.centre;
                sampling_ray.direction = Normalize(geometry.bounding_box.centre - sampling_ray.origin);

                Interaction interaction;

                if (Scene::IntersectInternal(temp_scene, sampling_ray, &interaction)) {
                    auto normal_vector = Normalize(interaction.point - geometry.bounding_box.centre);

                    sampled_points.push_back({interaction.point, normal_vector});
                }

            }
        }
*/
        const int num_points = 1024;

        auto sampler = std::make_unique<UniformSampler>(UniformSampler(2.0, 1024));

        for (int i = 0; i < num_points; ++i) {
            Point2r samp0 = sampler->Get2D();
            Point2r samp1 = sampler->Get2D();

            int tindex = geometry.num_faces * samp0.data[0];

            auto& triangle = tindices_dest[tindex];

            auto vertex0 = vertices_dest[triangle.v0];
            auto vertex1 = vertices_dest[triangle.v1];
            auto vertex2 = vertices_dest[triangle.v2];

            Point3r v0{vertex0.x, vertex0.y, vertex0.z};
            Point3r v1{vertex1.x, vertex1.y, vertex1.z};
            Point3r v2{vertex2.x, vertex2.y, vertex2.z};

            real u = samp0.data[1];
            real v = samp1.data[0];
            real w = samp1.data[1];

            real sum = u + v + w;

            u /= sum;
            v /= sum;
            w /= sum;

            Point3r sample_point = v0 * u + v1 * v + v2 * w;
            Vector3r normal = Normalize(sample_point - geometry.bounding_box.centre);
            avg_point = avg_point + sample_point;

            sampled_points.push_back({sample_point, normal});
        }

        avg_point = avg_point * (1. / num_points);
        for (auto it = sampled_points.begin(); it != sampled_points.end(); ++it) {
            it->second = Normalize(it->first - avg_point);
        }

        rtcUnmapBuffer(temp_scene, temp_geom_id, RTC_INDEX_BUFFER);
        rtcUnmapBuffer(temp_scene, temp_geom_id, RTC_VERTEX_BUFFER);
        rtcCommit(temp_scene);
        rtcDeleteScene(temp_scene);
    }
}

void ParticleTracer::CalculateRayLengthInVolume() {
    for (auto& geom_id : scene->other_geometry_ids) {
        rtcDisable(scene->scene_, geom_id);
    }
    rtcCommit(scene->scene_);

    const auto max_bounces = 20;
    const auto particles_per_point = 100;

    rays.clear();
    rays.reserve(sampled_points.size() * particles_per_point);

    auto sampler = std::make_unique<UniformSampler>(UniformSampler(2.0, 1024));

    for (auto& sample_point : sampled_points) {
        for (int s = 0; s < particles_per_point; ++s) {
            auto beta = Spectrum{1.f};

            auto bsdf_temp = std::make_unique<BSDF>(BSDF(sample_point.second));

            Vector3r win = CosineSampleHemisphere(sampler->Get2D());
            win = bsdf_temp->LocalToWorld(win);

            Ray ray{sample_point.first, win};
            auto bounces = 0;
            for (; bounces < max_bounces; ++bounces) {

                Interaction interaction;

                if (scene->Intersect(ray, &interaction, scene->mask_enable_all)) {

                    for (auto t = 0; t < interaction.t_max; t += 1.f) {
                        Point3r p = interaction.point + ray.direction * t;
                        StuffIntoBox(p, beta.r);
                    }


                    auto material = scene->geometries[interaction.geometry_id].material;
                    auto bsdf = std::unique_ptr<BSDF>(material->GetBSDF(interaction, 1.f));

                    bsdf->Choose(sampler->Get2D());

                    Vector3r wi, wo = ray.direction * -1;
                    real pdf;

                    auto f = bsdf->SampleF(wo, &wi, sampler->Get2D(), &pdf);

                    if (f.IsBlack() || pdf == 0) {
                        break;
                    }

                    beta = beta * f * std::abs(Dot(wi, Vector3r(Normalize(interaction.normal)))) * (1 / pdf);

                    ray.origin = interaction.point;
                    ray.direction = Normalize(wi);
                }
                else {
                    if (!std::isnan(ray.direction.values.x))
                        rays.push_back(ray);
                    break;
                }
            }
        }
    }

    real max_v = 0.;

    for (int i = 0; i < num_divisions; ++i) {
        for (int j = 0; j < num_divisions; ++j) {
            for (int k = 0; k < num_divisions; ++k) {
                auto value = boxes[((i * num_divisions) + j) * num_divisions + k];
                //std::cout << "[" << i << ", " << j << ", " << k << "]: " << value << "\n";
                max_v = (max_v < value) ? value : max_v;
            }
        }
    }

    max_value = max_v;

    for (int i = 0; i < num_divisions; ++i) {
        for (int j = 0; j < num_divisions; ++j) {
            for (int k = 0; k < num_divisions; ++k) {
                boxes[((i * num_divisions) + j) * num_divisions + k] /= max_value;
            }
        }
    }

    std::cout << std::endl;

    for (auto& geom_id : scene->other_geometry_ids) {
        rtcEnable(scene->scene_, geom_id);
    }
    rtcCommit(scene->scene_);
}

void ParticleTracer::StuffIntoBox(Point3r p, real factor) {
    int id = GetBoxId(p);
    if (id != -1) {
        boxes[id] += factor;
    }
}

real ParticleTracer::GetBoxValue(Point3r p) {
    int id = GetBoxId(p);
    if (id != -1) {
        return boxes[id];
    }
    else {
        return 0;
    }
}

bool IsInside(Point3r value, Vector3r minimum, Vector3r maximum) {
    return value.values.x <= maximum.values.x && value.values.y <= maximum.values.y && value.values.z <= maximum.values.z
            && value.values.x >= minimum.values.x && value.values.y >= minimum.values.y && value.values.z >= minimum.values.z;
}

int ParticleTracer::GetBoxId(Point3r p) {
    auto& geometry = scene->geometries[scene->lamp_geometry_id];
    auto& maxim = geometry.bounding_box.maximum_transformed.values;
    auto& minim = geometry.bounding_box.minimum_transformed.values;

    real extent_minimum = std::min(minim.x, std::min(minim.y, minim.z));
    real extent_maximum = std::max(maxim.x, std::max(maxim.y, maxim.z));

    auto p_transformed = p - geometry.bounding_box.centre;

    if (!IsInside(p, Vector3r{extent_minimum}, Vector3r{extent_maximum})) {
        return -1;
    }

    auto p_canonical = p_transformed * (num_divisions / (extent_maximum - extent_minimum)) + Vector3r{num_divisions/2};

    int x_id = Clamp((int)std::floor(p_canonical.values.x), 0, num_divisions - 1);
    int y_id = Clamp((int)std::floor(p_canonical.values.y), 0, num_divisions - 1);
    int z_id = Clamp((int)std::floor(p_canonical.values.z), 0, num_divisions - 1);

    return ((x_id * num_divisions) + y_id) * num_divisions + z_id;
}

void ParticleTracer::GenerateLights() {
    /*Bounds3D bounds;
    bounds.centre = avg_point;
    auto *tex = new TextureProjectionPointLight(bounds.centre, Spectrum());
    tex->image_texture->SetBounds(bounds);
    tex->image_texture->LoadBlank(256, 256);

    for (auto& ray : rays) {
        Interaction interaction;
        interaction.point = Point3r(ray.direction) + bounds.centre;

        tex->image_texture->StoreAndBlend(interaction, Spectrum(1));
    }

    scene->lights.push_back(tex);
    */
    const int num_clusters = 256;

    std::vector<std::array<real, 3>> data;
    data.reserve(rays.size());

    for (auto& ray : rays) {
        data.push_back({ray.origin.values.x, ray.origin.values.y, ray.origin.values.z});
    }

    auto& geometry = scene->geometries[scene->lamp_geometry_id];

    auto results = dkm::kmeans_lloyd(data, num_clusters);

    data.clear();

    auto& means = std::get<0>(results);
    auto& indices= std::get<1>(results);

    auto lights = std::vector<Light*>{};

    for (auto& mean : means) {
        Point3r centre {mean[0], mean[1], mean[2]};

        auto *tex = new TextureProjectionPointLight(centre, Spectrum());
        tex->image_texture->LoadBlank(64, 64);

        Bounds3D bounds;

        bounds.centre = centre;
        tex->image_texture->SetBounds(bounds);

        lights.push_back(tex);
    }

    for (int i = 0; i < indices.size(); ++i) {
        int index = indices[i];
        auto *tex = dynamic_cast<TextureProjectionPointLight*>(lights[index]);

        Interaction interaction;
        interaction.point = Point3r(rays[i].direction) + tex->position;

        tex->image_texture->StoreAndBlend(interaction, Spectrum(1.));
    }

    int num_particles = 100;

    for (auto& light : lights) {
        auto *tex = dynamic_cast<TextureProjectionPointLight*>(light);
        tex->image_texture->Multiply(1.f / num_particles);
        tex->image_texture->Blur();
        tex->attenuation_quadratic = 1;
        scene->lights.push_back(light);
    }

    rays.clear();

/*    rtcDisable(scene->scene_, scene->lamp_geometry_id);

    for (auto& id : scene->emissive_geometry_ids) {
        rtcDisable(scene->scene_, id);
    }

    rtcCommit(scene->scene_);*/
}
