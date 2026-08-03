#pragma once

#include <memory>
#include <vector>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "simulation/limb.hpp"
#include "simulation/organism.hpp"

#include "camera.hpp"
#include "shader.hpp"


class Renderer {
public:
    Renderer();
    ~Renderer();
    bool Update();
    void DrawLimb(std::weak_ptr<Limb> limb);
    void DrawBoundingBox(const Organism::BoundingBox& bounds);

private:
    void ProcessKeyboardInput(float deltaTime);
    void ProcessMouseInput();
    void DrawDebugBoundingBoxes();

    bool alive = false;
    GLFWwindow* window = nullptr;

    unsigned int limbvert;
    unsigned int limbVAO = 0;
    unsigned int limbVBO = 0;
    unsigned int limbEBO = 0;
    std::vector<std::weak_ptr<Limb>> limbList;

    unsigned int floorVAO = 0;
    unsigned int floorVBO = 0;

    unsigned int boxVAO = 0;
    unsigned int boxVBO = 0;
    std::vector<Organism::BoundingBox> boundingBoxes;

    Shader shader;
    Camera camera;
    float moveSpeed = 2.5f;
    float mouseSensitivity = 0.15f;
    double lastFrameTime = 0.0;
    double lastMouseX = 400.0;
    double lastMouseY = 300.0;
    bool leftMouseDown = false;
    bool debugBounds = false;
    bool debugKeyWasDown = false;
};
