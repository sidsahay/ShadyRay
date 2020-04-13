//
// Created by walksbynight on 23/3/18.
//

#include "sdlfilm.hpp"

SDLFilm::SDLFilm(const Film &film) : Film(film.width, film.height) {
    std::memcpy(pixel_buffer_, film.GetPixelBuffer(), width * height * sizeof(Spectrum));
}

SDLFilm::SDLFilm(int width, int height) : Film(width, height) {
}

void SDLFilm::Prepare(const char *window_name) {
    SDL_Init(SDL_INIT_VIDEO);

    sdl_pixel_buffer_ = new SDLRGBAPixel[width * height];

    for (std::size_t i = 0; i < width * height; ++i) {
        sdl_pixel_buffer_[i].r = Clamp(powf(pixel_buffer_[i].r, 0.2f), (real) 0., (real) 1.) * 255U;
        sdl_pixel_buffer_[i].g = Clamp(powf(pixel_buffer_[i].g, 0.2f), (real) 0., (real) 1.) * 255U;
        sdl_pixel_buffer_[i].b = Clamp(powf(pixel_buffer_[i].b, 0.2f), (real) 0., (real) 1.) * 255U;
        sdl_pixel_buffer_[i].a = 0U;
    }

    window = SDL_CreateWindow(window_name, SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, width, height,
                              SDL_WINDOW_SHOWN);
    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    SDL_SetRenderDrawColor(renderer, 0U, 0U, 0U, 0U);
    texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA32, SDL_TEXTUREACCESS_STREAMING, width, height);

    int pitch;
    SDLRGBAPixel *pixels;

    SDL_LockTexture(texture, nullptr, (void **) &pixels, &pitch);

    std::memcpy(pixels, sdl_pixel_buffer_, width * height * sizeof(SDLRGBAPixel));

    SDL_UnlockTexture(texture);
}

void SDLFilm::BlockingDisplay() {
    if (window == nullptr) {
        return;
    }

    SDL_RenderClear(renderer);
    SDL_RenderCopy(renderer, texture, nullptr, nullptr);
    SDL_RenderPresent(renderer);

    bool quit = false;

    while (!quit) {
        while (SDL_PollEvent(&event) != 0) {
            if (event.type == SDL_QUIT) {
                quit = true;
            }
        }
    }

    SDL_DestroyTexture(texture);
    texture = nullptr;
    SDL_DestroyRenderer(renderer);
    renderer = nullptr;
    SDL_DestroyWindow(window);
    window = nullptr;
}

SDLFilm::~SDLFilm() {
    delete[] sdl_pixel_buffer_;
}