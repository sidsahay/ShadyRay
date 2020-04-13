//
// Created by walksbynight on 16/4/18.
//

#ifndef SHADYRAY_PARTICLETRACER_HPP
#define SHADYRAY_PARTICLETRACER_HPP

#include "scene.hpp"
#include "../dkm/dkm.hpp"

class ParticleTracer {
public:

    Point3r avg_point;

    explicit ParticleTracer(Scene *scene);

    void PreProcessEmissive();

    void CalculateRayLengthInVolume();

    void GenerateLights();

    real GetBoxValue(Point3r p);

    static const int num_divisions = 40 ;

    std::array<real, num_divisions * num_divisions * num_divisions> boxes{};
    std::vector<Ray> rays;
    real max_value;

private:
    Scene *scene;
    std::vector<std::pair<Point3r, Vector3r> > sampled_points; //Store the sampled points on the emitter geometry.

    void StuffIntoBox(Point3r p, real factor);
    int GetBoxId(Point3r p);
};

#endif //SHADYRAY_PARTICLETRACER_HPP
