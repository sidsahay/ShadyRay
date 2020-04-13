//
// Created by walksbynight on 23/3/18.
//

#ifndef SHADYRAY_FILM_HPP
#define SHADYRAY_FILM_HPP

#include <cstring>
#include <iostream>

#include "defs.hpp"
#include "spectrum.hpp"

class Film;

class FilmTile;

class FilmTile {
public:
    FilmTile(int x_offset, int y_offset, int width, int height);

    FilmTile() = default;

    void Prepare(int x_offset, int y_offset, int width, int height);

    void SetPixel(int x, int y, const Spectrum &pixel);

    friend void Merge(Film &film, const FilmTile &tile);

    int x_offset;
    int y_offset;
    int width;
    int height;

    ~FilmTile();

protected:
    Spectrum *pixel_buffer_;
};

class Film {
public:
    Film(int width, int height);

    void FlipY();

    const Spectrum *GetPixelBuffer() const;

    friend void Merge(Film &film, const FilmTile &tile);
    friend void Merge(Film &film, const Film& other_film);

    friend std::ostream &operator<<(std::ostream &out, const Film &film);

    int width;
    int height;

    virtual ~Film();

protected:
    Spectrum *pixel_buffer_;
};

void operator<<(Film &film, const FilmTile &tile);

std::ostream &operator<<(std::ostream &out, const Film &film);

#endif //SHADYRAY_FILM_HPP
