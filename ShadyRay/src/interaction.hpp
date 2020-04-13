//
// Created by walksbynight on 23/3/18.
//

#ifndef SHADYRAY_INTERACTION_HPP
#define SHADYRAY_INTERACTION_HPP

#include "linalg.hpp"
#include "spectrum.hpp"

class Interaction {
public:
    Point3r point;
    Normal3r normal;

    real u;
    real v;

    real t_max;

    unsigned int geometry_id;
};

#endif //SHADYRAY_INTERACTION_HPP