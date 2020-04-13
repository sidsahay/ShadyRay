//
// Created by walksbynight on 23/3/18.
//

#include "film.hpp"


FilmTile::FilmTile(int x_offset, int y_offset, int width, int height)
        : x_offset(x_offset), y_offset(y_offset), width(width), height(height) {
    pixel_buffer_ = new Spectrum[width * height];
}


void FilmTile::Prepare(int x_off, int y_off, int w, int h) {
    x_offset = x_off;
    y_offset = y_off;
    width = w;
    height = h;

    pixel_buffer_ = new Spectrum[width * height];
}

void FilmTile::SetPixel(int x, int y, const Spectrum &pixel) {
    int idx = y * width + x;
    std::memcpy((pixel_buffer_ + idx), &(pixel), sizeof(pixel));
}

FilmTile::~FilmTile() {
    delete[] pixel_buffer_;
}

Film::Film(int width, int height) : width(width), height(height) {
    pixel_buffer_ = new Spectrum[width * height];
}


void Film::FlipY() {
    auto *new_pixel_buffer = new Spectrum[width * height];

    for (int j = 0; j < height; ++j) {
        int offset = j * width;
        int other_offset = (height - j - 1) * width;

        for (int i = 0; i < width; ++i) {
            int idx = offset + i;
            int other_idx = other_offset + i;

            new_pixel_buffer[other_idx] = pixel_buffer_[idx];
        }
    }

    delete[] pixel_buffer_;
    pixel_buffer_ = new_pixel_buffer;
}

const Spectrum *Film::GetPixelBuffer() const {
    return pixel_buffer_;
}

Film::~Film() {
    delete[] pixel_buffer_;
}

void Merge(Film &film, const FilmTile &tile) {
    for (std::size_t j = 0; j < tile.height; ++j) {
        std::size_t film_offset = (j + tile.y_offset) * film.width;
        std::size_t tile_offset = j * tile.width;

        for (std::size_t i = 0; i < tile.width; ++i) {
            std::size_t film_idx = film_offset + tile.x_offset + i;
            std::size_t tile_idx = tile_offset + i;

            std::memcpy(film.pixel_buffer_ + film_idx, tile.pixel_buffer_ + tile_idx, sizeof(Spectrum));
        }
    }
}

void operator<<(Film &film, const FilmTile &tile) {
    Merge(film, tile);
}

std::ostream &operator<<(std::ostream &out, const Film &film) {
    for (std::size_t j = 0; j < film.height; ++j) {
        std::size_t offset = j * film.width;

        for (std::size_t i = 0; i < film.width; ++i) {
            out << "{" << i << ", " << j << "} " << film.pixel_buffer_[offset + i] << std::endl;
        }
    }

    return out;
}

void Merge(Film &film, const Film &other_film) {
    for (std::size_t j = 0; j < film.height; ++j) {
        std::size_t film_offset = j * film.width;

        for (std::size_t i = 0; i < film.width; ++i) {
            std::size_t film_idx = film_offset + i;

            if(film.pixel_buffer_[film_idx].IsBlack()) {
                film.pixel_buffer_[film_idx] = other_film.pixel_buffer_[film_idx];
            }
        }
    }
}
