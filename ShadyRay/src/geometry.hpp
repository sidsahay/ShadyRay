//
// Created by walksbynight on 23/3/18.
//

#ifndef SHADYRAY_GEOMETRY_HPP
#define SHADYRAY_GEOMETRY_HPP

#include "material.hpp"
#include "linalg.hpp"

class Geometry {
public:
    Geometry() = default;

    Geometry(unsigned int geometry_id, Material *material, const Bounds3D &bounding_box);

    unsigned int num_vertices = 0;
    unsigned int num_faces = 0;

    Material *material;
    Bounds3D bounding_box;
    unsigned int geometry_id;
};

#endif //SHADYRAY_GEOMETRY_HPP
