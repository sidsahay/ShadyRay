//
// Created by walksbynight on 25/3/18.
//

#ifndef SHADYRAY_TEXTURE_HPP
#define SHADYRAY_TEXTURE_HPP

#include "interaction.hpp"
#include "spectrum.hpp"
#include "image.hpp"
#include "linalg.hpp"
#include "texturemap.hpp"

//T is usually float or Spectrum
class Texture {
public:
    Texture(TextureMap *map);

    void SetBounds (const Bounds3D& bounds);

    virtual Spectrum Evaluate(const Interaction &interaction) = 0;

    virtual ~Texture() = default;

protected:
    TextureMap *map;
};



//Constant texture: useful for initializing material properties

class ConstantTexture : public Texture {
public:
    explicit ConstantTexture(Spectrum constant_value, TextureMap *map);

    Spectrum Evaluate(const Interaction &interaction) override;

    ~ConstantTexture() override = default;

protected:
    Spectrum constant_value;
};


//Image texture: Duh. Should only return Spectrum values, because duh
class ImageTexture : public Texture {
public:
    ImageTexture(TextureMap *map, const char *image_path);

    Spectrum Evaluate(const Interaction &interaction) override;

    void StoreAndBlend(const Interaction &interaction, const Spectrum &value);

    void LoadBlank(unsigned int w, unsigned int h);

    void Multiply(real value);

    void Blur();

protected:
    Image image;
};

#endif //SHADYRAY_TEXTURE_HPP
