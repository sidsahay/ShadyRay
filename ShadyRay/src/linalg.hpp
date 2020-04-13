//
// Created by walksbynight on 23/3/18.
//

#ifndef SHADYRAY_LINALG_HPP
#define SHADYRAY_LINALG_HPP
#pragma once

#include "defs.hpp"
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <limits>

struct Vector3r;
struct Vector2r;
struct Point3r;
struct Point2r;
struct Normal3r;

typedef std::pair<Point3r, Point3r> Box;

real Clamp(real value, real lowBound, real highBound);
int Clamp(int value, int lower, int upper);

Vector3r SphericalDirection(real sinTheta, real cosTheta, real phi);

struct Vector3r {
    union {
        struct {
            real x;
            real y;
            real z;
        } values;
        real data[3];
    };

    explicit Vector3r(real val = 0);

    Vector3r(real x, real y, real z);

    explicit Vector3r(const Point3r &point);

    explicit Vector3r(const Normal3r &normal);

    Vector3r operator+(const Vector3r &v) const;

    Vector3r operator-(const Vector3r &v) const;

    Vector3r operator*(real s) const;

    Vector3r operator/(real s) const;

    Vector3r operator-() const;

    void LocalCoordinateSystem(Vector3r *v2, Vector3r *v3);

    real LengthSquared() const;

    real Length() const;
};

Vector3r Normalize(const Vector3r &v);

real Dot(const Vector3r &v0, const Vector3r &v1);

Vector3r Cross(const Vector3r &v0, const Vector3r &v1);

std::ostream &operator<<(std::ostream &os, const Vector3r &v);


struct Vector2r {
    union {
        struct {
            real x;
            real y;
        } values;
        real data[2];
    };

    explicit Vector2r(real val = 0);

    Vector2r(real x, real y);

    explicit Vector2r(const Point2r &point);

    Vector2r operator+(const Vector2r &v) const;

    Vector2r operator-(const Vector2r &v) const;

    Vector2r operator*(real s) const;

    Vector2r operator/(real s) const;

    Vector2r operator-() const;

    real LengthSquared() const;

    real Length() const;
};

Vector2r Normalize(const Vector2r &v);

real Dot(const Vector2r &v0, const Vector2r &v1);

std::ostream &operator<<(std::ostream &os, const Vector2r &v);


struct Point3r {
    union {
        struct {
            real x;
            real y;
            real z;
        } values;
        real data[3];
    };

    explicit Point3r(real val = 0);

    Point3r(real x, real y, real z);

    explicit Point3r(const Vector3r &vector);

    Point3r operator+(const Point3r &p) const;

    Point3r operator*(real s) const;

    Point3r operator/(real s) const;

    Point3r operator-() const;

    Point3r operator+(const Vector3r &v) const;

    Vector3r operator-(const Point3r &p) const;
};

real DistanceSquared(const Point3r &p0, const Point3r &p1);

real Distance(const Point3r &p0, const Point3r &p1);

Point3r Lerp(const Point3r &p0, const Point3r &p1, real t);

std::ostream &operator<<(std::ostream &os, const Point3r &v);

struct Point2r {
    union {
        struct {
            real x;
            real y;
        } values;
        real data[2];
    };

    explicit Point2r(real val = (real) 0.);

    Point2r(real x, real y);

    explicit Point2r(const Vector2r &vector);

    Point2r operator+(const Point2r &p) const;

    Point2r operator*(real s) const;

    Point2r operator/(real s) const;

    Point2r operator-() const;

    Point2r operator+(const Vector2r &v) const;

    Vector2r operator-(const Point2r &p) const;
};

real DistanceSquared(const Point2r &p0, const Point2r &p1);

real Distance(const Point2r &p0, const Point2r &p1);

Point2r Lerp(const Point2r &p0, const Point2r &p1, real t);

std::ostream &operator<<(std::ostream &os, const Point3r &v);

struct Normal3r {
    union {
        struct {
            real x;
            real y;
            real z;
        } values;
        real data[3];
    };

    explicit Normal3r(real val = (real) 0);

    Normal3r(real x, real y, real z);

    explicit Normal3r(const Vector3r &vector);

    Normal3r operator+(const Normal3r &v) const;

    Normal3r operator-(const Normal3r &v) const;

    Normal3r operator*(real s) const;

    Normal3r operator/(real s) const;

    Normal3r operator-() const;

    real LengthSquared() const;

    real Length() const;
};

Normal3r Normalize(const Normal3r &v);

real Dot(const Normal3r &v0, const Normal3r &v1);

std::ostream &operator<<(std::ostream &os, const Normal3r &v);

struct Matrix4r {
    float data[4][4];

    Matrix4r() = default;

    Matrix4r(real m00, real m01, real m02, real m03,
             real m10, real m11, real m12, real m13,
             real m20, real m21, real m22, real m23,
             real m30, real m31, real m32, real m33);

    void LoadIdentity();

};

Matrix4r Transpose(const Matrix4r &m);

//Matrix4r Inverse(const Matrix4r &m);

Matrix4r Multiply(const Matrix4r &m1, const Matrix4r &m2);


std::ostream &operator<<(std::ostream &os, const Matrix4r &m);


struct Transform {
    Matrix4r matrix;
    Matrix4r invMatrix;

    Transform() = default;

    Transform(real m00, real m01, real m02, real m03,
              real m10, real m11, real m12, real m13,
              real m20, real m21, real m22, real m23,
              real m30, real m31, real m32, real m33);

    Vector3r Apply(const Vector3r &v) const;

    Point3r Apply(const Point3r &p) const;

    Normal3r Apply(const Normal3r &n) const;

    Vector3r InvApply(const Vector3r &v) const;

    Point3r InvApply(const Point3r &p) const;

    Normal3r InvApply(const Normal3r &n) const;

    Vector3r operator()(const Vector3r &v) const;

    Point3r operator()(const Point3r &p) const;

    Normal3r operator()(const Normal3r &n) const;

    Transform operator*(const Transform &t) const;

    void Invert();
//    void GenerateInverse();
};

Transform Translate(const Vector3r &v);

Transform Scale(const Vector3r &scales);

Transform Rotate(const Vector3r &axis, real angle);

Transform Rotate(const Vector3r &tbAngles);

std::ostream &operator<<(std::ostream &os, const Transform &m);

class Bounds3D {
public:
    Bounds3D();

    void ComputeCentre();

    bool IsInside(const Point3r &point);

    Point3r centre;
    Point3r minimum;
    Point3r maximum;

    Point3r minimum_transformed;
    Point3r maximum_transformed;
};

#endif //SHADYRAY_LINALG_HPP
