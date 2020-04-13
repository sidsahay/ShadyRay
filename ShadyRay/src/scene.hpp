//
// Created by walksbynight on 23/3/18.
//

#ifndef SHADYRAY_SCENE_HPP
#define SHADYRAY_SCENE_HPP

#include "embree2/rtcore.h"
#include "embree2/rtcore_ray.h"

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include <iostream>
#include <vector>
#include <map>
#include <cstring>

#include "interaction.hpp"
#include "light.hpp"
#include "geometry.hpp"

static constexpr real kNudgeFactor = (real) 0.00001;
static constexpr unsigned int kMaxGeometries = 32;

struct MaterialTable {
public:
    std::map<std::string, Material*> material_table;
};

class ParticleTracer;

class Scene {
public:

    Scene(const char *path, real scale, MaterialTable *material_table);

    explicit Scene(real scale);

    bool LoadScene(const char *path);

    bool Intersect(const Ray &ray, Interaction *interaction, unsigned int mask);

    bool IntersectCheck(const Ray &ray, unsigned int mask);

    void DeleteLights();

    std::vector<Light *> lights;


    MaterialTable *material_table;

    Geometry geometries[kMaxGeometries];

    real scale;

    friend class ParticleTracer;

    static bool IntersectInternal(RTCScene scene, const Ray &ray, Interaction *interaction, unsigned int mask);

    std::vector<unsigned int> geometry_ids;
    std::vector<unsigned int> emissive_geometry_ids;
    std::vector<unsigned int> other_geometry_ids;

    unsigned int mask_enable_all = 0;
    unsigned int mask_enable_others = 0;

    unsigned int lamp_geometry_id;

protected:
    struct Vertex_ {
        float x, y, z, a;
    };

    struct TriangleIndices_ {
        unsigned int v0, v1, v2;
    };

    void ProcessNode(const aiNode *node, const aiScene *scene_data);

    void ProcessMesh(const aiMesh *mesh);

    //void ProcessMaterials(const aiScene *scene_data);

    void ProcessLights(const aiScene *scene_data);

    RTCDevice device_;
    RTCScene scene_;
};

#endif //SHADYRAY_SCENE_HPP
