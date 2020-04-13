//
// Created by walksbynight on 23/3/18.
//

#include "linalg.hpp"


Vector3r SphericalDirection(real sinTheta, real cosTheta, real phi) {
    return {sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta};
}

//Vector3r

Vector3r::Vector3r(real val) {
    values.x = values.y = values.z = val;
}

Vector3r::Vector3r(real x, real y, real z) {
    values.x = x;
    values.y = y;
    values.z = z;
}

Vector3r::Vector3r(const Point3r &point) {
    values.x = point.values.x;
    values.y = point.values.y;
    values.z = point.values.z;
}

Vector3r::Vector3r(const Normal3r &normal) {
    values.x = normal.values.x;
    values.y = normal.values.y;
    values.z = normal.values.z;
}

Vector3r Vector3r::operator+(const Vector3r &v) const {
    Vector3r ret;
    for (int i = 0; i < 3; ++i) {
        ret.data[i] = data[i] + v.data[i];
    }

    return ret;
}

Vector3r Vector3r::operator-(const Vector3r &v) const {
    Vector3r ret;
    for (int i = 0; i < 3; ++i) {
        ret.data[i] = data[i] - v.data[i];
    }

    return ret;
}

Vector3r Vector3r::operator*(real s) const {
    Vector3r ret;
    for (int i = 0; i < 3; ++i) {
        ret.data[i] = data[i] * s;
    }

    return ret;
}

Vector3r Vector3r::operator/(real s) const {
    auto f = (real) 1. / s;
    return (*this) * f;
}

Vector3r Vector3r::operator-() const {
    return (*this) * -1.;
}

void Vector3r::LocalCoordinateSystem(Vector3r *v2, Vector3r *v3) {
    if (std::abs(values.x) > std::abs(values.y)) {
        *v2 = Vector3r(-values.z, (real) 0., values.x) / std::sqrt(values.x * values.x + values.z * values.z);
    } else {
        *v2 = Vector3r((real) 0., values.z, -values.y) / std::sqrt(values.y * values.y + values.z * values.z);
    }

    *v3 = Cross(*this, *v2);
}

real Vector3r::LengthSquared() const {
    auto s = (real) 0.;
    for (auto i : data) {
        s += i * i;
    }

    return s;
}

real Vector3r::Length() const {
    return std::sqrt(LengthSquared());
}

real Clamp(real value, real lowBound, real highBound) {
    return (value > highBound) ? highBound : (value < lowBound) ? lowBound : value;
}

int Clamp(int value, int lower, int upper) {
    if (value < lower) {
        value = lower;
    }
    else if (value > upper) {
        value = upper;
    }

    return value;
}

Vector3r Normalize(const Vector3r &v) {
    return v / v.Length();
}

real Dot(const Vector3r &v0, const Vector3r &v1) {
    auto s = (real) 0.;
    for (int i = 0; i < 3; ++i) {
        s += v0.data[i] * v1.data[i];
    }
    return s;
}

Vector3r Cross(const Vector3r &v0, const Vector3r &v1) {
    Vector3r ret;

    ret.values.x = (v0.values.y * v1.values.z) - (v1.values.y * v0.values.z);
    ret.values.y = (v1.values.x * v0.values.z) - (v0.values.x * v1.values.z);
    ret.values.z = (v0.values.x * v1.values.y) - (v1.values.x * v0.values.y);

    return ret;
}

std::ostream &operator<<(std::ostream &os, const Vector3r &v) {
    os << "[Vector3r " << v.values.x << " " << v.values.y << " " << v.values.z << "]";
    return os;
}


//Vector2r

Vector2r::Vector2r(real val) {
    values.x = values.y = val;
}

Vector2r::Vector2r(real x, real y) {
    values.x = x;
    values.y = y;
}

Vector2r::Vector2r(const Point2r &point) {
    values.x = point.values.x;
    values.y = point.values.y;
}

Vector2r Vector2r::operator+(const Vector2r &v) const {
    Vector2r ret;
    for (int i = 0; i < 2; ++i) {
        ret.data[i] = data[i] + v.data[i];
    }

    return ret;
}

Vector2r Vector2r::operator-(const Vector2r &v) const {
    Vector2r ret;
    for (int i = 0; i < 2; ++i) {
        ret.data[i] = data[i] - v.data[i];
    }

    return ret;
}

Vector2r Vector2r::operator*(real s) const {
    Vector2r ret;
    for (int i = 0; i < 2; ++i) {
        ret.data[i] = data[i] * s;
    }

    return ret;
}

Vector2r Vector2r::operator/(real s) const {
    auto f = (real) 1. / s;
    return (*this) * f;
}

Vector2r Vector2r::operator-() const {
    return (*this) * -1.;
}

real Vector2r::LengthSquared() const {
    auto s = (real) 0.;
    for (auto i : data) {
        s += i * i;
    }

    return s;
}

real Vector2r::Length() const {
    return std::sqrt(LengthSquared());
}

Vector2r Normalize(const Vector2r &v) {
    return v / v.Length();
}

real Dot(const Vector2r &v0, const Vector2r &v1) {
    auto s = (real) 0.;
    for (int i = 0; i < 2; ++i) {
        s += v0.data[i] * v1.data[i];
    }
    return s;
}

std::ostream &operator<<(std::ostream &os, const Vector2r &v) {
    os << "[Vector2r " << v.values.x << " " << v.values.y << "]";
    return os;
}

//Point3r

Point3r::Point3r(real val) {
    values.x = values.y = values.z = val;
}

Point3r::Point3r(real x, real y, real z) {
    values.x = x;
    values.y = y;
    values.z = z;
}

Point3r::Point3r(const Vector3r &vector) {
    values.x = vector.values.x;
    values.y = vector.values.y;
    values.z = vector.values.z;
}

Point3r Point3r::operator+(const Point3r &p) const {
    Point3r ret;
    for (int i = 0; i < 3; ++i) {
        ret.data[i] = data[i] + p.data[i];
    }

    return ret;
}

Point3r Point3r::operator*(real s) const {
    Point3r ret;
    for (int i = 0; i < 3; ++i) {
        ret.data[i] = data[i] * s;
    }

    return ret;
}

Point3r Point3r::operator/(real s) const {
    auto f = (real) 1. / s;
    return (*this) * f;
}

Point3r Point3r::operator-() const {
    return (*this) * -1.;
}

Point3r Point3r::operator+(const Vector3r &v) const {
    return {values.x + v.values.x, values.y + v.values.y, values.z + v.values.z};
}

Vector3r Point3r::operator-(const Point3r &p) const {
    return {values.x - p.values.x, values.y - p.values.y, values.z - p.values.z};
}

real DistanceSquared(const Point3r &p0, const Point3r &p1) {
    return (p1 - p0).LengthSquared();
}

real Distance(const Point3r &p0, const Point3r &p1) {
    return std::sqrt(DistanceSquared(p0, p1));
}

Point3r Lerp(const Point3r &p0, const Point3r &p1, real t) {
    return p0 * (1. - t) + p1 * t;
}

std::ostream &operator<<(std::ostream &os, const Point3r &v) {
    os << "[Point3r " << v.values.x << " " << v.values.y << " " << v.values.z << "]";
    return os;
}


//Point2r

Point2r::Point2r(real val) {
    values.x = values.y = val;
}

Point2r::Point2r(real x, real y) {
    values.x = x;
    values.y = y;
}

Point2r::Point2r(const Vector2r &vector) {
    values.x = vector.values.x;
    values.y = vector.values.y;
}

Point2r Point2r::operator+(const Point2r &p) const {
    Point2r ret;
    for (int i = 0; i < 2; ++i) {
        ret.data[i] = data[i] + p.data[i];
    }

    return ret;
}

Point2r Point2r::operator*(real s) const {
    Point2r ret;
    for (int i = 0; i < 2; ++i) {
        ret.data[i] = data[i] * s;
    }

    return ret;
}

Point2r Point2r::operator/(real s) const {
    real f = (real) 1. / s;
    return (*this) * f;
}

Point2r Point2r::operator-() const {
    return (*this) * -1.;
}

Point2r Point2r::operator+(const Vector2r &v) const {
    return {values.x + v.values.x, values.y + v.values.y};
}

Vector2r Point2r::operator-(const Point2r &p) const {
    return {values.x - p.values.x, values.y - p.values.y};
}

real DistanceSquared(const Point2r &p0, const Point2r &p1) {
    return (p1 - p0).LengthSquared();
}

real Distance(const Point2r &p0, const Point2r &p1) {
    return (p1 - p0).Length();
}

Point2r Lerp(const Point2r &p0, const Point2r &p1, real t) {
    return p0 * (1. - t) + p1 * t;
}

std::ostream &operator<<(std::ostream &os, const Point2r &v) {
    os << "[Point2r " << v.values.x << " " << v.values.y << "]";
    return os;
}

//Normal3r

Normal3r::Normal3r(real val) {
    values.x = values.y = values.z = val;
}

Normal3r::Normal3r(real x, real y, real z) {
    values.x = x;
    values.y = y;
    values.z = z;
}

Normal3r::Normal3r(const Vector3r &vector) {
    values.x = vector.values.x;
    values.y = vector.values.y;
    values.z = vector.values.z;
}

Normal3r Normal3r::operator+(const Normal3r &v) const {
    Normal3r ret;
    for (int i = 0; i < 3; ++i) {
        ret.data[i] = data[i] + v.data[i];
    }

    return ret;
}

Normal3r Normal3r::operator-(const Normal3r &v) const {
    Normal3r ret;
    for (int i = 0; i < 3; ++i) {
        ret.data[i] = data[i] - v.data[i];
    }

    return ret;
}

Normal3r Normal3r::operator*(real s) const {
    Normal3r ret;
    for (int i = 0; i < 3; ++i) {
        ret.data[i] = data[i] * s;
    }

    return ret;
}

Normal3r Normal3r::operator/(real s) const {
    real f = (real) 1. / s;
    return (*this) * f;
}

Normal3r Normal3r::operator-() const {
    return (*this) * -1.;
}

real Normal3r::LengthSquared() const {
    auto s = (real) 0.;
    for (float i : data) {
        s += i * i;
    }

    return s;
}

real Normal3r::Length() const {
    return std::sqrt(LengthSquared());
}

Normal3r Normalize(const Normal3r &v) {
    return v / v.Length();
}

real Dot(const Normal3r &v0, const Normal3r &v1) {
    auto s = (real) 0.;
    for (int i = 0; i < 3; ++i) {
        s += v0.data[i] * v1.data[i];
    }
    return s;
}

std::ostream &operator<<(std::ostream &os, const Normal3r &v) {
    os << "[Normal3r " << v.values.x << " " << v.values.y << " " << v.values.z << "]";
    return os;
}


//Matrix4r

/*

Matrix4r::Matrix4r() {

}
*/

Matrix4r::Matrix4r(real m00, real m01, real m02, real m03,
                   real m10, real m11, real m12, real m13,
                   real m20, real m21, real m22, real m23,
                   real m30, real m31, real m32, real m33) {
    data[0][0] = m00;
    data[0][1] = m01;
    data[0][2] = m02;
    data[0][3] = m03;
    data[1][0] = m10;
    data[1][1] = m11;
    data[1][2] = m12;
    data[1][3] = m13;
    data[2][0] = m20;
    data[2][1] = m21;
    data[2][2] = m22;
    data[2][3] = m23;
    data[3][0] = m30;
    data[3][1] = m31;
    data[3][2] = m32;
    data[3][3] = m33;
}


void Matrix4r::LoadIdentity() {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; j++) {
            if (i != j) {
                data[i][j] = (real) 0.;
            } else {
                data[i][j] = (real) 1.;
            }
        }
    }
}

Matrix4r Transpose(const Matrix4r &m) {
    return {m.data[0][0], m.data[1][0], m.data[2][0], m.data[3][0],
            m.data[0][1], m.data[1][1], m.data[2][1], m.data[3][1],
            m.data[0][2], m.data[1][2], m.data[2][2], m.data[3][2],
            m.data[0][3], m.data[1][3], m.data[2][3], m.data[3][3]};
}

/*
Matrix4r Inverse(const Matrix4r &m) {
    Eigen::Matrix4f mat;
    mat << m.data[0][0], m.data[0][1], m.data[0][2], m.data[0][3],
            m.data[1][0], m.data[1][1], m.data[1][2], m.data[1][3],
            m.data[2][0], m.data[2][1], m.data[2][2], m.data[2][3],
            m.data[3][0], m.data[3][1], m.data[3][2], m.data[3][3];

    Eigen::Matrix4f ans = mat.inverse();

    Matrix4r answer;

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; j++) {
            answer.data[i][j] = ans(i, j);
        }
    }

    return answer;
}
*/

Matrix4r Multiply(const Matrix4r &m1, const Matrix4r &m2) {
    Matrix4r m;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; j++) {
            m.data[i][j] =
                    m1.data[i][0] * m2.data[0][j] +
                    m1.data[i][1] * m2.data[1][j] +
                    m1.data[i][2] * m2.data[2][j] +
                    m1.data[i][3] * m2.data[3][j];
        }
    }

    return m;
}

std::ostream &operator<<(std::ostream &os, const Matrix4r &m) {
    os << "[Matrix4r";
    for (auto *i : m.data) {
        os << "\n";
        for (int j = 0; j < 4; j++) {
            os << i[j] << " ";
        }
    }

    os << "]";
    return os;
}


//Transform

/*
Transform::Transform() {
}*/

Transform::Transform(real m00, real m01, real m02, real m03,
                     real m10, real m11, real m12, real m13,
                     real m20, real m21, real m22, real m23,
                     real m30, real m31, real m32, real m33)
        : matrix(m00, m01, m02, m03,
                 m10, m11, m12, m13,
                 m20, m21, m22, m23,
                 m30, m31, m32, m33) {

}

Vector3r Transform::Apply(const Vector3r &v) const {
    real xv = matrix.data[0][0] * v.data[0] + matrix.data[0][1] * v.data[1] + matrix.data[0][2] * v.data[2];
    real yv = matrix.data[1][0] * v.data[0] + matrix.data[1][1] * v.data[1] + matrix.data[1][2] * v.data[2];
    real zv = matrix.data[2][0] * v.data[0] + matrix.data[2][1] * v.data[1] + matrix.data[2][2] * v.data[2];

    return {xv, yv, zv};
}


Normal3r Transform::Apply(const Normal3r &n) const {
    real xn = invMatrix.data[0][0] * n.data[0] + invMatrix.data[0][1] * n.data[1] + invMatrix.data[0][2] * n.data[2];
    real yn = invMatrix.data[1][0] * n.data[0] + invMatrix.data[1][1] * n.data[1] + invMatrix.data[1][2] * n.data[2];
    real zn = invMatrix.data[2][0] * n.data[0] + invMatrix.data[2][1] * n.data[1] + invMatrix.data[2][2] * n.data[2];

    return {xn, yn, zn};
}

Point3r Transform::Apply(const Point3r &p) const {
    real xp = matrix.data[0][0] * p.data[0] + matrix.data[0][1] * p.data[1] + matrix.data[0][2] * p.data[2] +
              matrix.data[0][3];
    real yp = matrix.data[1][0] * p.data[0] + matrix.data[1][1] * p.data[1] + matrix.data[1][2] * p.data[2] +
              matrix.data[1][3];
    real zp = matrix.data[2][0] * p.data[0] + matrix.data[2][1] * p.data[1] + matrix.data[2][2] * p.data[2] +
              matrix.data[2][3];
    real wp = matrix.data[3][0] * p.data[0] + matrix.data[3][1] * p.data[1] + matrix.data[3][2] * p.data[2] +
              matrix.data[3][3];

    if (wp == 1.) {
        return {xp, yp, zp};
    } else {
        return Point3r(xp, yp, zp) / wp;
    }
}

Vector3r Transform::InvApply(const Vector3r &v) const {
    real xv = invMatrix.data[0][0] * v.data[0] + invMatrix.data[0][1] * v.data[1] + invMatrix.data[0][2] * v.data[2];
    real yv = invMatrix.data[1][0] * v.data[0] + invMatrix.data[1][1] * v.data[1] + invMatrix.data[1][2] * v.data[2];
    real zv = invMatrix.data[2][0] * v.data[0] + invMatrix.data[2][1] * v.data[1] + invMatrix.data[2][2] * v.data[2];

    return {xv, yv, zv};
}


Normal3r Transform::InvApply(const Normal3r &n) const {
    real xn = matrix.data[0][0] * n.data[0] + matrix.data[0][1] * n.data[1] + matrix.data[0][2] * n.data[2];
    real yn = matrix.data[1][0] * n.data[0] + matrix.data[1][1] * n.data[1] + matrix.data[1][2] * n.data[2];
    real zn = matrix.data[2][0] * n.data[0] + matrix.data[2][1] * n.data[1] + matrix.data[2][2] * n.data[2];

    return {xn, yn, zn};
}

Point3r Transform::InvApply(const Point3r &p) const {
    real xp = invMatrix.data[0][0] * p.data[0] + invMatrix.data[0][1] * p.data[1] + invMatrix.data[0][2] * p.data[2] +
              invMatrix.data[0][3];
    real yp = invMatrix.data[1][0] * p.data[0] + invMatrix.data[1][1] * p.data[1] + invMatrix.data[1][2] * p.data[2] +
              invMatrix.data[1][3];
    real zp = invMatrix.data[2][0] * p.data[0] + invMatrix.data[2][1] * p.data[1] + invMatrix.data[2][2] * p.data[2] +
              invMatrix.data[2][3];
    real wp = invMatrix.data[3][0] * p.data[0] + invMatrix.data[3][1] * p.data[1] + invMatrix.data[3][2] * p.data[2] +
              invMatrix.data[3][3];

    if (wp == 1.) {
        return {xp, yp, zp};
    } else {
        return Point3r(xp, yp, zp) / wp;
    }
}


Vector3r Transform::operator()(const Vector3r &v) const {
    return Apply(v);
}


Normal3r Transform::operator()(const Normal3r &n) const {
    return Apply(n);
}

Point3r Transform::operator()(const Point3r &p) const {
    return Apply(p);
}

Transform Transform::operator*(const Transform &t) const {
    Transform ret;
    ret.matrix = Multiply(matrix, t.matrix);
//    ret.GenerateInverse();

    return ret;
}

void Transform::Invert() {
    Matrix4r temp = matrix;
    matrix = invMatrix;
    invMatrix = temp;
}

/*
void Transform::GenerateInverse()
{
    invMatrix = Inverse(matrix);
}
*/

Transform Translate(const Vector3r &v) {
    Transform t((real) 1., (real) 0., (real) 0., v.values.x,
                (real) 0., (real) 1., (real) 0., v.values.y,
                (real) 0., (real) 0., (real) 1., v.values.z,
                (real) 0., (real) 0., (real) 0., (real) 1.);

    Matrix4r invMat((real) 1., (real) 0., (real) 0., -v.values.x,
                    (real) 0., (real) 1., (real) 0., -v.values.y,
                    (real) 0., (real) 0., (real) 1., -v.values.z,
                    (real) 0., (real) 0., (real) 0., (real) 1.);

    t.invMatrix = invMat;

    return t;
}

Transform Scale(const Vector3r &s) {
    Transform t(s.values.x, (real) 0., (real) 0., (real) 0.,
                (real) 0., s.values.y, (real) 0., (real) 0.,
                (real) 0., (real) 0., s.values.z, (real) 0.,
                (real) 0., (real) 0., (real) 0., (real) 1.);

    Matrix4r invMat((real) 1. / s.values.x, (real) 0., (real) 0., (real) 0.,
                    (real) 0., (real) 1. / s.values.y, (real) 0., (real) 0.,
                    (real) 0., (real) 0., (real) 1. / s.values.z, (real) 0.,
                    (real) 0., (real) 0., (real) 0., (real) 1.);

    t.invMatrix = invMat;
    return t;
}

Transform Rotate(const Vector3r &axis, real angle) {
    real ux = axis.values.x;
    real uy = axis.values.y;
    real uz = axis.values.z;

    real cosTheta = std::cos(angle);
    real sinTheta = std::sin(angle);

    real oneMinusCosTheta = (real) 1. - cosTheta;

    Transform t(
            cosTheta + ux * ux * oneMinusCosTheta,
            ux * uy * oneMinusCosTheta - uz * sinTheta,
            ux * uz * oneMinusCosTheta + uy * sinTheta,
            (real) 0.,
            uy * ux * oneMinusCosTheta + uz * sinTheta,
            cosTheta + uy * uy * oneMinusCosTheta,
            uy * uz * oneMinusCosTheta - ux * sinTheta,
            (real) 0.,
            uz * ux * oneMinusCosTheta - uy * sinTheta,
            uz * uy * oneMinusCosTheta + ux * sinTheta,
            cosTheta + uz * uz * oneMinusCosTheta,
            (real) 0.,
            (real) 0.,
            (real) 0.,
            (real) 0.,
            (real) 1.);

    t.invMatrix = Transpose(t.matrix);
    return t;
}

Transform Rotate(const Vector3r &tbAngles) {
    real s1 = std::sin(tbAngles.values.x);
    real s2 = std::sin(tbAngles.values.y);
    real s3 = std::sin(tbAngles.values.z);

    real c1 = std::cos(tbAngles.values.x);
    real c2 = std::cos(tbAngles.values.y);
    real c3 = std::cos(tbAngles.values.z);

    Transform t(
            c2 * c3,
            -c2 * s3,
            s2,
            (real) 0.,
            c1 * s3 + c3 * s1 * s2,
            c1 * c3 - s1 * s2 * s3,
            -c2 * s1,
            (real) 0.,
            s1 * s3 - c1 * c3 * s2,
            c3 * s1 + c1 * s2 * s3,
            c1 * c2,
            (real) 0.,
            (real) 0.,
            (real) 0.,
            (real) 0.,
            (real) 1.);

    t.invMatrix = Transpose(t.invMatrix);
    return t;
}

std::ostream &operator<<(std::ostream &os, const Transform &m) {
    os << "[Transform]\n" << m.matrix << "\n" << m.invMatrix;
    return os;
}

Bounds3D::Bounds3D() {
    minimum = {std::numeric_limits<real>::max(), std::numeric_limits<real>::max(), std::numeric_limits<real>::max()};
    maximum = {std::numeric_limits<real>::min(), std::numeric_limits<real>::min(), std::numeric_limits<real>::min()};
}

void Bounds3D::ComputeCentre() {
    centre = (minimum + maximum) / (real)2.;
    minimum_transformed = Point3r(minimum - centre);
    maximum_transformed = Point3r(maximum - centre);
}

bool Bounds3D::IsInside(const Point3r &point) {
    return true;
}
