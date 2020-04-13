//
// Created by walksbynight on 23/3/18.
//

#include "camera.hpp"


Camera::Camera() : origin(Vector3r((real) 0.)), focal_length((real) 100.) {
    rotation_mode = kAxisAngle;
    axis_angle.axis = Vector3r((real) 0., (real) 1., (real) 0.);
    axis_angle.angle = (real) 0.;

    CameraToWorld = Translate(Vector3r(origin)) * Rotate(axis_angle.axis, axis_angle.angle);
    focus = Point3r((real) 0., (real) 0., focal_length);
}

void Camera::UpdateCamera() {
    switch (rotation_mode) {
        case kAxisAngle:
            CameraToWorld = Translate(Vector3r(origin)) * Rotate(axis_angle.axis, axis_angle.angle);
            break;

        case kTaitBryan:
            CameraToWorld = Translate(Vector3r(origin)) * Rotate(tait_bryan.angles);
            break;
    }

    focus = CameraToWorld(Point3r((real) 0., (real) 0., focal_length));
}

void Camera::SetOrigin(const Point3r &org) {
    origin = org;
    UpdateCamera();
}

void Camera::SetRotation(const Vector3r &axis, real angle) {
    rotation_mode = kAxisAngle;
    axis_angle.axis = axis;
    axis_angle.angle = angle;

    UpdateCamera();
}

void Camera::SetRotation(const Vector3r &tb_angles) {
    rotation_mode = kTaitBryan;
    tait_bryan.angles = tb_angles;

    UpdateCamera();
}

void Camera::SetFocalLength(real f_len) {
    focal_length = f_len;
    focus = CameraToWorld(Point3r((real) 0., (real) 0., focal_length));
}

Ray Camera::GenerateRay(const Point3r &image_point) {
    Ray ray;
    Point3r scaled_point = {image_point.values.x * x_scale, image_point.values.y * y_scale, (real) 0.};
    ray.origin = CameraToWorld(scaled_point);
    ray.direction = Normalize(ray.origin - focus);
    ray.origin = focus;
    return ray;
}

void Camera::SetScale(real x_sc, real y_sc) {
    x_scale = x_sc;
    y_scale = y_sc;
}
