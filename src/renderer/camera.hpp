#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

class Camera {
public:
    Camera();

    void SetPosition(const glm::vec3& position);
    void SetTarget(const glm::vec3& target);
    void SetPerspective(float fovYDegrees, float aspectRatio, float nearPlane, float farPlane);

    glm::vec3 GetPosition() const;
    glm::vec3 GetFront() const;
    void Move(const glm::vec3& offset);
    void Rotate(float yawOffset, float pitchOffset, bool constrainPitch = true);

    glm::mat4 ViewMatrix() const;
    glm::mat4 ProjectionMatrix() const;

private:
    void UpdateDirectionVectors();
    glm::vec3 position{-14.6793f, 10.6245f, 51.2586f};
    glm::vec3 front{0.0f, 0.0f, -1.0f};
    glm::vec3 up{0.0f, 1.0f, 0.0f};
    glm::vec3 right{1.0f, 0.0f, 0.0f};
    glm::vec3 worldUp{0.0f, 1.0f, 0.0f};

    float yaw = -72.7500f;
    float pitch = -7.6850f;

    float fovYDegrees = 45.0f;
    float aspectRatio = 1.0f;
    float nearPlane = 0.1f;
    float farPlane = 100.0f;
};
