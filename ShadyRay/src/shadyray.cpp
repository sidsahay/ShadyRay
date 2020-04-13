#include <sstream>

#include "pathintegratorth.hpp"
#include "sdlfilm.hpp"
#include "luaintegration.hpp"
#include "particletracer.hpp"

extern "C" {
    static int render(LUA_state L) {
        LuaContext lc(L);

        auto num_materials = lc.GetTableFieldLocal<int>("materials.num_materials");
        auto material_table = new MaterialTable();

        for (int i = 0; i < num_materials; ++i) {
            std::string prefix = "materials.material_table." + std::to_string(i+1);

            auto name = lc.GetTableFieldLocal<const char*>((prefix + ".mesh_name").c_str());

            auto d_r = (real)lc.GetTableFieldLocal<double>((prefix + ".diffuse_color.r").c_str());
            auto d_g = (real)lc.GetTableFieldLocal<double>((prefix + ".diffuse_color.g").c_str());
            auto d_b = (real)lc.GetTableFieldLocal<double>((prefix + ".diffuse_color.b").c_str());

            auto s_r = (real)lc.GetTableFieldLocal<double>((prefix + ".specular_color.r").c_str());
            auto s_g = (real)lc.GetTableFieldLocal<double>((prefix + ".specular_color.g").c_str());
            auto s_b = (real)lc.GetTableFieldLocal<double>((prefix + ".specular_color.b").c_str());

            auto t_r = (real)lc.GetTableFieldLocal<double>((prefix + ".transmission_color.r").c_str());
            auto t_g = (real)lc.GetTableFieldLocal<double>((prefix + ".transmission_color.g").c_str());
            auto t_b = (real)lc.GetTableFieldLocal<double>((prefix + ".transmission_color.b").c_str());

            auto e_r = (real)lc.GetTableFieldLocal<double>((prefix + ".emission_color.r").c_str());
            auto e_g = (real)lc.GetTableFieldLocal<double>((prefix + ".emission_color.g").c_str());
            auto e_b = (real)lc.GetTableFieldLocal<double>((prefix + ".emission_color.b").c_str());

            Spectrum diffuse_color = {d_r, d_g, d_b};
            Spectrum specular_color = {s_r, s_g, s_b};
            Spectrum transmission_color = {t_r, t_g, t_b};
            Spectrum emission_color = {e_r, e_g, e_b};

            auto material = new Material;

            if (!diffuse_color.IsBlack()) {
                material->material_type |= kMaterialDiffuse;
                material->diffuse_color = diffuse_color;
            }

            if (!specular_color.IsBlack()) {
                auto alpha = (real)lc.GetTableFieldLocal<double>((prefix + ".alpha").c_str());

                material->material_type |= kMaterialSpecular;
                material->specular_color = specular_color;
                material->alpha = alpha;
            }

            if (!transmission_color.IsBlack()) {
                auto alpha = (real)lc.GetTableFieldLocal<double>((prefix + ".alpha").c_str());
                auto ior = (real)lc.GetTableFieldLocal<double>((prefix + ".transmission_color.ior").c_str());

                material->material_type |= kMaterialTransparent;
                material->transmission_color = transmission_color;
                material->alpha = alpha;
                material->ior = ior;
            }

            if (!emission_color.IsBlack()) {
                material->material_type |= kMaterialEmission;
                material->emission_color = emission_color;
            }

            auto mapping = lc.GetTableFieldLocal<const char*>((prefix + ".mapping").c_str());

            auto color_map_type = lc.GetTableFieldLocal<const char*>((prefix + ".color_map.type").c_str());
            if (!strcmp(color_map_type, "Image")) {
                TextureMap *mapping_obj = new PlaneMap(mapping);
                auto color_map_path = lc.GetTableFieldLocal<const char*>((prefix + ".color_map.path").c_str());
                material->color_map = new ImageTexture(mapping_obj, color_map_path);
            }
            else {
                material->color_map = nullptr;
            }

            auto normal_map_type = lc.GetTableFieldLocal<const char*>((prefix + ".normal_map.type").c_str());
            if (!strcmp(normal_map_type, "Image")) {
                TextureMap *mapping_obj = new PlaneMap(mapping);
                auto normal_map_path = lc.GetTableFieldLocal<const char*>((prefix + ".normal_map.path").c_str());
                material->normal_map = new ImageTexture(mapping_obj, normal_map_path);
            }
            else {
                material->normal_map = nullptr;
            }

            auto is_lamp = lc.GetTableFieldLocal<int>((prefix + ".is_lamp").c_str());
            material->is_lamp = is_lamp == 1;

            material_table->material_table[name] = material;
        }

        lc.Pop(1);

        auto scene_path = lc.GetTableFieldLocal<const char*>("render_options.scene.path");
        auto scene_scale = lc.GetTableFieldLocal<double>("render_options.scene.scale");
        auto *scene = new Scene(scene_path, (real)scene_scale, material_table);

        auto *camera = new Camera;

        auto camera_focal_length = lc.GetTableFieldLocal<double>("render_options.camera.focal_length");
        camera->SetFocalLength((real)camera_focal_length);

        auto camera_scale_x = lc.GetTableFieldLocal<double>("render_options.camera.scale.x");
        auto camera_scale_y = lc.GetTableFieldLocal<double>("render_options.camera.scale.y");
        camera->SetScale((real)camera_scale_x, (real)camera_scale_y);

        auto camera_origin_x = lc.GetTableFieldLocal<double>("render_options.camera.origin.x");
        auto camera_origin_y = lc.GetTableFieldLocal<double>("render_options.camera.origin.y");
        auto camera_origin_z = lc.GetTableFieldLocal<double>("render_options.camera.origin.z");
        camera->SetOrigin(Point3r((real)camera_origin_x, (real)camera_origin_y, (real)camera_origin_z));

        auto camera_rotation_x = lc.GetTableFieldLocal<double>("render_options.camera.rotation.x");
        auto camera_rotation_y = lc.GetTableFieldLocal<double>("render_options.camera.rotation.y");
        auto camera_rotation_z = lc.GetTableFieldLocal<double>("render_options.camera.rotation.z");
        camera->SetRotation(Vector3r((real)camera_rotation_x, (real)camera_rotation_y, (real)camera_rotation_z));

        auto scene_width = lc.GetTableFieldLocal<int>("render_options.scene.width");
        auto scene_height = lc.GetTableFieldLocal<int>("render_options.scene.height");
        auto *film = new SDLFilm(scene_width, scene_height);

        auto sampler_multiplier = lc.GetTableFieldLocal<double>("render_options.sampler.multiplier");
        auto sampler_num_samples = lc.GetTableFieldLocal<int>("render_options.sampler.num_samples");
        auto *sampler = new UniformSampler((real)sampler_multiplier, sampler_num_samples);

        auto filter_b = lc.GetTableFieldLocal<double>("render_options.filter.b");
        auto filter_c = lc.GetTableFieldLocal<double>("render_options.filter.c");
        auto filter_radius = lc.GetTableFieldLocal<double>("render_options.filter.radius");
        auto *filter = new MitchellNetravaliFilter((real)filter_b, (real)filter_c, (real)filter_radius);

        auto integrator_x_blocks = lc.GetTableFieldLocal<int>("render_options.integrator.x_blocks");
        auto integrator_y_blocks = lc.GetTableFieldLocal<int>("render_options.integrator.y_blocks");
        auto integrator_num_bounces = lc.GetTableFieldLocal<int>("render_options.integrator.num_bounces");
        PathIntegratorTh<Scene, Camera, SDLFilm, UniformSampler, MitchellNetravaliFilter>
                path_integrator(integrator_x_blocks, integrator_y_blocks, integrator_num_bounces, scene, camera, film, sampler, filter);


        auto tracer = ParticleTracer{scene};
        tracer.PreProcessEmissive();
        tracer.CalculateRayLengthInVolume();
        tracer.GenerateLights();

        std::cout << "Preprocess done." << std::endl;

        path_integrator.tracer = &tracer;
        path_integrator.SetGlow(scene->geometries[scene->emissive_geometry_ids[0]].material->emission_color);
        path_integrator.Render();

        film->FlipY();

        auto window_name = lc.GetTableFieldLocal<const char*>("render_options.window_name");
        film->Prepare(window_name);

        film->BlockingDisplay();

        return 0;
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "Please specify the path to the Lua script. All other settings must go into the script.\n";
    }

    LuaContext lc(argv[1]);
    lc.RegisterFunction("render", render);

    lc.Execute(0, 0);

    return 0;
}