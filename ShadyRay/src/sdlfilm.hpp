//
// Created by walksbynight on 23/3/18.
//

#ifndef SHADYRAY_SDLFILM_HPP
#define SHADYRAY_SDLFILM_HPP

#pragma once

#include "linalg.hpp"
#include "film.hpp"

#include <SDL2/SDL.h>

class SDLFilm : public Film {
public:
    explicit SDLFilm(const Film &film);

    SDLFilm(int width, int height);

    void Prepare(const char *window_name);

    void BlockingDisplay();

    ~SDLFilm() override;

private:
    struct SDLRGBAPixel {
        unsigned char r;
        unsigned char g;
        unsigned char b;
        unsigned char a;
    };
    SDLRGBAPixel *sdl_pixel_buffer_;
    SDL_Window *window;
    SDL_Renderer *renderer;
    SDL_Texture *texture;
    SDL_Event event;
};

#endif //SHADYRAY_SDLFILM_HPP
