#pragma once

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "limb.hpp"

#include "camera.hpp"
#include "shader.hpp"


class Renderer {
public:
    Renderer();
    ~Renderer();
    bool Update();
    void DrawLimb(Limb* limb);

private:
    void ProcessKeyboardInput(float deltaTime);
    void ProcessMouseInput();

    bool alive = false;
    GLFWwindow* window = nullptr;
    unsigned int planeVAO = 0;
    unsigned int planeVBO = 0;
    Shader shader;
    Camera camera;
    float moveSpeed = 2.5f;
    float mouseSensitivity = 0.15f;
    double lastFrameTime = 0.0;
    double lastMouseX = 400.0;
    double lastMouseY = 300.0;
    bool leftMouseDown = false;
};
