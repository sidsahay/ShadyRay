//
// Created by walksbynight on 26/3/18.
//

#include "texture.hpp"


Texture::Texture(TextureMap *map) : map(map){
}

void Texture::SetBounds(const Bounds3D &bounds) {
    map->SetBounds(bounds);
}

ConstantTexture::ConstantTexture(Spectrum constant_value, TextureMap *map):Texture(map), constant_value(constant_value) {
}

Spectrum ConstantTexture::Evaluate(const Interaction &interaction) {
    return constant_value;
}

ImageTexture::ImageTexture(TextureMap *map, const char *image_path) : Texture(map) {
    if (image_path != nullptr) {
        image.LoadImage(image_path);
    }
}

Spectrum ImageTexture::Evaluate(const Interaction &interaction) {
    Point2r uv = map->Map(interaction);
    return image.Evaluate(uv.values.x, uv.values.y);
}

void ImageTexture::StoreAndBlend(const Interaction &interaction, const Spectrum &value) {
    Point2r uv = map->Map(interaction);
    image.StoreAndBlend(uv.values.x, uv.values.y, value);
}

void ImageTexture::LoadBlank(unsigned int w, unsigned int h) {
    image.LoadBlank(w, h);
}

void ImageTexture::Multiply(real value) {
    image.Multiply(value);
}

void ImageTexture::Blur() {
    image.Blur();
}
