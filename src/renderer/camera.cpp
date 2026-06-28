#include "camera.hpp"

Camera::Camera() = default;

void Camera::SetPosition(const glm::vec3& position)
{
    this->position = position;
}

void Camera::SetTarget(const glm::vec3& target)
{
    this->target = target;
}

void Camera::SetPerspective(float fovYDegrees, float aspectRatio, float nearPlane, float farPlane)
{
    this->fovYDegrees = fovYDegrees;
    this->aspectRatio = aspectRatio;
    this->nearPlane = nearPlane;
    this->farPlane = farPlane;
}

glm::mat4 Camera::ViewMatrix() const
{
    return glm::lookAt(position, target, up);
}

glm::mat4 Camera::ProjectionMatrix() const
{
    return glm::perspective(glm::radians(fovYDegrees), aspectRatio, nearPlane, farPlane);
}
