 //
// Created by walksbynight on 23/3/18.
//

#include "scene.hpp"

Scene::Scene(const char *path, real scale, MaterialTable *material_table) : scale(scale), material_table(material_table) {
    device_ = rtcNewDevice();
    scene_ = rtcDeviceNewScene(device_, RTC_SCENE_DYNAMIC, RTC_INTERSECT1);

    if (!LoadScene(path)) {
        std::cout << "Could not load scene in scene ctor. Committing suicide...\n";
        exit(-1);
    }
}

Scene::Scene(real scale) : scale(scale) {
    device_ = rtcNewDevice();
    scene_ = rtcDeviceNewScene(device_, RTC_SCENE_DYNAMIC, RTC_INTERSECT1);
}

bool Scene::LoadScene(const char *path) {
    Assimp::Importer importer;
    const aiScene *scene_data = importer.ReadFile(path, aiProcess_PreTransformVertices);

    if (!scene_data) {
        return false;
    }

    ProcessLights(scene_data);

    ProcessNode(scene_data->mRootNode, scene_data);

    mask_enable_all = 0;
    for (int i = 0; i < geometry_ids.size(); ++i) {
        unsigned int mask = 1 << (geometry_ids[i]);
        mask_enable_all |= mask;
        rtcSetMask(scene_, geometry_ids[i], mask);
    }

    mask_enable_others = 0;
    for (int i = 0; i < other_geometry_ids.size(); ++i) {
        unsigned int mask = 1 << (other_geometry_ids[i]);
        mask_enable_others |= mask;
    }

    rtcCommit(scene_);
    return true;
}

bool Scene::IntersectInternal(RTCScene scene, const Ray &ray, Interaction *interaction, unsigned int mask) {
    RTCRay rtc_ray;

    rtc_ray.mask = mask;

    std::memcpy(rtc_ray.org, ray.origin.data, sizeof(rtc_ray.org));
    std::memcpy(rtc_ray.dir, ray.direction.data, sizeof(rtc_ray.dir));

    rtc_ray.tnear = (real) kNudgeFactor;
    rtc_ray.tfar = (real) 200000.;

    rtc_ray.geomID = RTC_INVALID_GEOMETRY_ID;

    rtcIntersect(scene, rtc_ray);

    if (rtc_ray.geomID == RTC_INVALID_GEOMETRY_ID) {
        return false;
    } else {
        std::memcpy(interaction->normal.data, rtc_ray.Ng, sizeof(rtc_ray.Ng));

        Point3r hit_point = ray.origin + ray.direction * (rtc_ray.tfar);
        std::memcpy(interaction->point.data, hit_point.data, sizeof(hit_point.data));

        interaction->u = rtc_ray.u;
        interaction->v = rtc_ray.v;

        interaction->geometry_id = rtc_ray.geomID;

        interaction->t_max = rtc_ray.tfar;

        return true;
    }
}

bool Scene::Intersect(const Ray &ray, Interaction *interaction, unsigned int mask) {
    return IntersectInternal(scene_, ray, interaction, mask);
}

bool Scene::IntersectCheck(const Ray &ray, unsigned int mask) {
    RTCRay rtc_ray;

    rtc_ray.mask = mask;

    rtc_ray.tfar = ray.direction.Length();

    std::memcpy(rtc_ray.org, ray.origin.data, sizeof(rtc_ray.org));
    std::memcpy(rtc_ray.dir, Normalize(ray.direction).data, sizeof(rtc_ray.dir));

    rtc_ray.tnear = (real) kNudgeFactor;

    rtc_ray.geomID = RTC_INVALID_GEOMETRY_ID;

    rtcOccluded(scene_, rtc_ray);

    return rtc_ray.geomID == 0;
}

void Scene::ProcessNode(const aiNode *node, const aiScene *scene_data) {
    std::cout << "Node name: " << node->mName.C_Str() << std::endl;

    for (unsigned int i = 0; i < node->mNumMeshes; ++i) {
        aiMesh *mesh = scene_data->mMeshes[node->mMeshes[i]];
        ProcessMesh(mesh);
    }

    for (unsigned int i = 0; i < node->mNumChildren; ++i) {
        ProcessNode(node->mChildren[i], scene_data);
    }
}


/*void Scene::ProcessMaterials(const aiScene *scene_data) {
    for (int i = 0; i < scene_data->mNumMaterials; ++i) {
        aiMaterial *mat = scene_data->mMaterials[i];

        Material *material = new Material;

        //Check diffuse material
        aiColor3D diffuse_color;
        if (AI_SUCCESS == mat->Get(AI_MATKEY_COLOR_DIFFUSE, diffuse_color)) {
            material->material_type |= kMaterialDiffuse;

            material->diffuse_color.r = diffuse_color.r;
            material->diffuse_color.g = diffuse_color.g;
            material->diffuse_color.b = diffuse_color.b;
        }

        //Check specular material
        aiColor3D specular_color;
        if (AI_SUCCESS == mat->Get(AI_MATKEY_COLOR_SPECULAR, specular_color)) {
            material->material_type |= kMaterialSpecular;

            material->specular_color.r = specular_color.r;
            material->specular_color.g = specular_color.g;
            material->specular_color.b = specular_color.b;
        }

        //Check emission material
        aiColor3D emission_color;
        if (AI_SUCCESS == mat->Get(AI_MATKEY_COLOR_EMISSIVE, emission_color)) {
            material->material_type |= kMaterialEmission;

            material->emission_color.r = emission_color.r;
            material->emission_color.g = emission_color.g;
            material->emission_color.b = emission_color.b;
        }

        //Check transparent material
        aiColor3D transparent_color;
        if (AI_SUCCESS == mat->Get(AI_MATKEY_COLOR_TRANSPARENT, transparent_color)) {
            material->material_type |= kMaterialTransparent;

            material->transparent_color.r = transparent_color.r;
            material->transparent_color.g = transparent_color.g;
            material->transparent_color.b = transparent_color.b;
        }

        //Check roughness
        float shininess;
        if (AI_SUCCESS == mat->Get(AI_MATKEY_SHININESS, shininess)) {
            material->roughness = (real) 1. / shininess;
        }

        //check IOR
        float ior;
        if (AI_SUCCESS == mat->Get(AI_MATKEY_REFRACTI, ior)) {
            material->eta = ior;
        }

        //No material? Set diffuse white
        if (material->material_type == 0) {
            material->material_type |= kMaterialDiffuse;

            material->diffuse_color.r = diffuse_color.r;
            material->diffuse_color.g = diffuse_color.g;
            material->diffuse_color.b = diffuse_color.b;
        }

        aiString string;
        if (AI_SUCCESS == mat->Get(AI_MATKEY_NAME, string)) {
            std::cout << string.data << " D: " << material->diffuse_color << " S: " << material->specular_color
                      << std::endl;
        }

        materials.push_back(material);
    }
}*/

void Scene::ProcessLights(const aiScene *scene_data) {
    auto num_lights = scene_data->mNumLights;

    for (int i = 0; i < num_lights; ++i) {
        auto light = scene_data->mLights[i];
        if (light->mType == aiLightSource_POINT) {
            auto position = Point3r(light->mPosition.x, light->mPosition.y, light->mPosition.z) * scale;
            auto color = Spectrum(light->mColorDiffuse.r, light->mColorDiffuse.g, light->mColorDiffuse.b);

            Light *temp_light = new Light(position, color);

            temp_light->attenuation_constant = light->mAttenuationConstant;
            temp_light->attenuation_linear = light->mAttenuationLinear / scale;
            temp_light->attenuation_quadratic = light->mAttenuationQuadratic / (scale * scale);

            lights.push_back(temp_light);
        }
    }
}
void Scene::ProcessMesh(const aiMesh *mesh) {
    Bounds3D bounds;

    unsigned int geom_id = rtcNewTriangleMesh2(scene_, RTC_GEOMETRY_STATIC, mesh->mNumFaces, mesh->mNumVertices, 1);

    auto *vertices = (Vertex_ *) rtcMapBuffer(scene_, geom_id, RTC_VERTEX_BUFFER);
    for (unsigned int i = 0; i < mesh->mNumVertices; ++i) {
        auto vertex = mesh->mVertices[i];

        vertices[i].x = vertex.x * scale;
        vertices[i].y = vertex.y * scale;
        vertices[i].z = vertex.z * scale;

        bounds.minimum.values.x = std::min(bounds.minimum.values.x, vertices[i].x);
        bounds.minimum.values.y = std::min(bounds.minimum.values.y, vertices[i].y);
        bounds.minimum.values.z = std::min(bounds.minimum.values.z, vertices[i].z);

        bounds.maximum.values.x = std::max(bounds.maximum.values.x, vertices[i].x);
        bounds.maximum.values.y = std::max(bounds.maximum.values.y, vertices[i].y);
        bounds.maximum.values.z = std::max(bounds.maximum.values.z, vertices[i].z);
    }
    rtcUnmapBuffer(scene_, geom_id, RTC_VERTEX_BUFFER);

    auto *tindices = (TriangleIndices_ *) rtcMapBuffer(scene_, geom_id, RTC_INDEX_BUFFER);
    for (unsigned int i = 0; i < mesh->mNumFaces; ++i) {
        auto face = mesh->mFaces[i];

        if (face.mNumIndices == 3) {
            unsigned int *indices = face.mIndices;

            tindices[i].v0 = indices[0];
            tindices[i].v1 = indices[1];
            tindices[i].v2 = indices[2];
        }
    }
    rtcUnmapBuffer(scene_, geom_id, RTC_INDEX_BUFFER);

    bounds.ComputeCentre();

    auto material = material_table->material_table[mesh->mName.C_Str()];
    std::cout << "Mesh name: " << mesh->mName.C_Str() << std::endl;

    if (material == nullptr) {
        return;
    }
    if (material->color_map != nullptr) {
        material->color_map->SetBounds(bounds);
    }

    if (material->normal_map != nullptr) {
        material->normal_map->SetBounds(bounds);
    }

    Geometry geometry(geom_id, material, bounds);
    geometry.num_vertices = mesh->mNumVertices;
    geometry.num_faces = mesh->mNumFaces;

    geometries[geom_id] = geometry;

    geometry_ids.push_back(geom_id);

    if (TestFlag(material->material_type, kMaterialEmission)) {
        emissive_geometry_ids.push_back(geom_id);
    }
    else if (material->is_lamp) {
        lamp_geometry_id = geom_id;
    }
    else {
        other_geometry_ids.push_back(geom_id);
    }
}

void Scene::DeleteLights() {
    for (Light* light : lights) {
        delete light;
    }

    lights.clear();
}