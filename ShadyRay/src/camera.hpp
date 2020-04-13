//
// Created by walksbynight on 23/3/18.
//

#ifndef SHADYRAY_CAMERA_HPP
#define SHADYRAY_CAMERA_HPP

#include "linalg.hpp"
#include "ray.hpp"

class Camera {
public:
    Camera();

    void UpdateCamera();

    void SetOrigin(const Point3r &org);

    void SetRotation(const Vector3r &axis, real angle);

    void SetRotation(const Vector3r &tb_angles);

    void SetFocalLength(real f_len);

    void SetScale(real x_sc, real y_sc);

    Ray GenerateRay(const Point3r &image_point);

    Point3r origin;
    Point3r focus;
    real focal_length;

    enum RotationMode {
        kAxisAngle,
        kTaitBryan,
    };
    RotationMode rotation_mode;

    union {
        struct {
            Vector3r axis;
            real angle;
        } axis_angle;

        struct {
            Vector3r angles;
        } tait_bryan;
    };

    Transform CameraToWorld;

    real x_scale;
    real y_scale;
};

#endif //SHADYRAY_CAMERA_HPP
