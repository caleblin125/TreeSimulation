#include "camera.hpp"

#include <glm/gtc/constants.hpp>

Camera::Camera()
{
    UpdateDirectionVectors();
}

void Camera::SetPosition(const glm::vec3& position)
{
    this->position = position;
}

void Camera::SetTarget(const glm::vec3& target)
{
    glm::vec3 direction = glm::normalize(target - position);
    pitch = glm::degrees(glm::asin(direction.y));
    yaw = glm::degrees(glm::atan(direction.z, direction.x));
    UpdateDirectionVectors();
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
    return glm::lookAt(position, position + front, up);
}

glm::mat4 Camera::ProjectionMatrix() const
{
    return glm::perspective(glm::radians(fovYDegrees), aspectRatio, nearPlane, farPlane);
}

glm::vec3 Camera::GetPosition() const
{
    return position;
}

glm::vec3 Camera::GetFront() const
{
    return front;
}

void Camera::Move(const glm::vec3& offset)
{
    position += offset;
}

void Camera::Rotate(float yawOffset, float pitchOffset, bool constrainPitch)
{
    yaw += yawOffset;
    pitch += pitchOffset;

    if (constrainPitch)
    {
        if (pitch > 89.0f)
            pitch = 89.0f;
        if (pitch < -89.0f)
            pitch = -89.0f;
    }

    UpdateDirectionVectors();
}

void Camera::UpdateDirectionVectors()
{
    glm::vec3 newFront;
    newFront.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
    newFront.y = sin(glm::radians(pitch));
    newFront.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
    front = glm::normalize(newFront);
    right = glm::normalize(glm::cross(front, worldUp));
    up = glm::normalize(glm::cross(right, front));
}
