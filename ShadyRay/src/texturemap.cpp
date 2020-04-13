//
// Created by walksbynight on 26/3/18.
//

#include "texturemap.hpp"

void TextureMap::SetBounds(const Bounds3D &bounds) {
    this->bounds = bounds;
}

PlaneMap::PlaneMap(const char *type) {
    if (!strcmp(type, "PlaneX")) {
        primary_axis = 0;
        secondary_axis0 = 1;
        secondary_axis1 = 2;
    }
    else if (!strcmp(type, "PlaneY")) {
        primary_axis = 1;
        secondary_axis0 = 0;
        secondary_axis1 = 2;
    }
    else if (!strcmp(type, "PlaneZ")) {
        primary_axis = 2;
        secondary_axis0 = 0;
        secondary_axis1 = 1;
    }
    else {
        primary_axis = 1;
        secondary_axis0 = 0;
        secondary_axis1 = 2;
    }
}

Point2r PlaneMap::Map(const Interaction &interaction) {
    Point3r transformed_point = Point3r(interaction.point - bounds.centre);
    /*real x_length = bounds.maximum.values.x;
    real z_length = bounds.maximum.values.z;

    real x_extent = transformed_point.values.x;
    real z_extent = transformed_point.values.z;
*/
    real secondary0_length = bounds.maximum.data[secondary_axis0];
    real secondary1_length = bounds.maximum.data[secondary_axis1];

    real secondary0_extent = transformed_point.data[secondary_axis0];
    real secondary1_extent = transformed_point.data[secondary_axis1];

//    return {Clamp(x_extent / x_length, (real)-1., (real)1.), Clamp(z_extent / z_length, (real)-1., (real)1.)};
    return {Clamp(secondary0_extent / secondary0_length, (real)-1., (real)1.), Clamp(secondary1_extent / secondary1_length, (real)-1., (real)1.)};
}

Point2r SphericalMap::Map(const Interaction &interaction) {
    Vector3r direction = Normalize(interaction.point - bounds.centre);

    real theta = Clamp(acos(direction.values.y), 0, 3.1415926);
    real phi = Clamp((real)atan2(direction.values.z, direction.values.x), -3.1515926, 3.1415926);

    return {phi / (real)(3.1415926), 2 * theta / (real)3.1415926 - 1};
}

Point2r CubeMap::Map(const Interaction &interaction) {
    auto direction = Normalize(interaction.point - bounds.centre);

    auto x = direction.values.x;
    auto y = direction.values.y;
    auto z = direction.values.z;

    float absX = fabs(x);
    float absY = fabs(y);
    float absZ = fabs(z);

    int isXPositive = x > 0 ? 1 : 0;
    int isYPositive = y > 0 ? 1 : 0;
    int isZPositive = z > 0 ? 1 : 0;

    float maxAxis, uc, vc;

    int index = 0;

    // POSITIVE X
    if (isXPositive && absX >= absY && absX >= absZ) {
        // u (0 to 1) goes from +z to -z
        // v (0 to 1) goes from -y to +y
        maxAxis = absX;
        uc = -z;
        vc = y;
        index = 0;
    }
    // NEGATIVE X
    if (!isXPositive && absX >= absY && absX >= absZ) {
        // u (0 to 1) goes from -z to +z
        // v (0 to 1) goes from -y to +y
        maxAxis = absX;
        uc = z;
        vc = y;
        index = 1;
    }
    // POSITIVE Y
    if (isYPositive && absY >= absX && absY >= absZ) {
        // u (0 to 1) goes from -x to +x
        // v (0 to 1) goes from +z to -z
        maxAxis = absY;
        uc = x;
        vc = -z;
        index = 2;
    }
    // NEGATIVE Y
    if (!isYPositive && absY >= absX && absY >= absZ) {
        // u (0 to 1) goes from -x to +x
        // v (0 to 1) goes from -z to +z
        maxAxis = absY;
        uc = x;
        vc = z;
        index = 3;
    }
    // POSITIVE Z
    if (isZPositive && absZ >= absX && absZ >= absY) {
        // u (0 to 1) goes from -x to +x
        // v (0 to 1) goes from -y to +y
        maxAxis = absZ;
        uc = x;
        vc = y;
        index = 4;
    }
    // NEGATIVE Z
    if (!isZPositive && absZ >= absX && absZ >= absY) {
        // u (0 to 1) goes from +x to -x
        // v (0 to 1) goes from -y to +y
        maxAxis = absZ;
        uc = -x;
        vc = y;
        index = 5;
    }

    // Convert range from -1 to 1 to 0 to 1
    real u = 0.5f * (uc / maxAxis + 1.0f);
    real v = 0.5f * (vc / maxAxis + 1.0f);

    real u_step = 1.f / 3;
    real v_step = 1.f / 2;

    int u_factor = index % 3;
    int v_factor = index / 3;

    real u_actual = u_factor * u_step + u / 3;
    real v_actual = v_factor * v_step + v / 2;

    return {2 * (u_actual - 0.5f), 2 * (v_actual - 0.5f)};
}
