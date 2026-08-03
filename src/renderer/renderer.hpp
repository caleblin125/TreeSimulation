#pragma once

#include <atomic>
#include <mutex>
#include <memory>
#include <thread>
#include <vector>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "simulation/limb.hpp"
#include "simulation/organism.hpp"

#include "camera.hpp"
#include "shader.hpp"

struct RenderLimb
{
    glm::vec3 root;
    float length = 0.0f;
    float radius = 0.0f;
    glm::quat orientation;
};

struct RenderFrame
{
    std::vector<RenderLimb> limbs;
    std::vector<Organism::BoundingBox> boundingBoxes;
};

class Renderer {
public:
    Renderer();
    ~Renderer();
    bool Update();
    void DrawLimb(std::weak_ptr<Limb> limb);
    void DrawLimb(const RenderLimb& limb);
    void DrawBoundingBox(const Organism::BoundingBox& bounds);
    void SetFrame(RenderFrame frame);
    bool IsAlive() const;

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
    std::vector<RenderLimb> limbList;

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

class RenderPipeline
{
public:
    RenderPipeline();
    ~RenderPipeline();

    RenderPipeline(const RenderPipeline&) = delete;
    RenderPipeline& operator=(const RenderPipeline&) = delete;

    bool IsRunning() const;
    void SubmitFrame(RenderFrame frame);
    void Stop();

private:
    void RenderLoop();

    mutable std::mutex frameMutex;
    RenderFrame pendingFrame;
    bool hasPendingFrame = false;
    std::atomic<bool> running{false};
    std::thread renderThread;
};
