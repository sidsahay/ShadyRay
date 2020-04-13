//
// Created by walksbynight on 23/3/18.
//

#include "geometry.hpp"

Geometry::Geometry(unsigned int geometry_id, Material *material, const Bounds3D &bounding_box) : geometry_id(geometry_id), material(material), bounding_box(bounding_box) {
}
