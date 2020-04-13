//
// Created by walksbynight on 26/3/18.
//

#ifndef SHADYRAY_TEXTUREMAP_HPP
#define SHADYRAY_TEXTUREMAP_HPP

#include "linalg.hpp"
#include "interaction.hpp"

#include <cstring>

class TextureMap {
public:
    explicit TextureMap() = default;

    virtual Point2r Map(const Interaction& interaction) = 0;

    void SetBounds(const Bounds3D &bounds);

    virtual ~TextureMap() = default;

protected:
    Bounds3D bounds;
};

class PlaneMap : public TextureMap {
public:
    explicit PlaneMap(const char *type);

    Point2r Map(const Interaction& interaction) override;

private:
    int primary_axis, secondary_axis0, secondary_axis1;
};


class SphericalMap : public TextureMap {
public:
    Point2r Map(const Interaction& interaction) override;
};


class CubeMap : public TextureMap {
public:
    Point2r Map(const Interaction& interaction) override;
};

#endif //SHADYRAY_TEXTUREMAP_HPP
