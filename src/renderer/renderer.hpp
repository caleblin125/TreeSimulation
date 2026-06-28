#pragma once

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "camera.hpp"
#include "shader.hpp"

class Renderer {
public:
    Renderer();
    ~Renderer();
    bool Update();

private:
    bool alive = false;
    GLFWwindow* window = nullptr;
    unsigned int planeVAO = 0;
    unsigned int planeVBO = 0;
    Shader shader;
    Camera camera;
};
