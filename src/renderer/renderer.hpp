#pragma once

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <memory>
#include <optional>
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
    float strain = 0.0f;
};

struct RenderFrame
{
    struct Plant
    {
        unsigned int id = 0;
        double energy = 0.0;
        Organism::Stats stats;
        Organism::BoundingBox bounds{};
        std::vector<RenderLimb> limbs;
    };

    struct SoilCell
    {
        int cellX = 0;
        int cellZ = 0;
        float resource = 1.0f;
    };

    std::vector<Plant> plants;
    std::vector<SoilCell> soilCells;
    std::vector<Organism::NeuralNetworkSample> neuralSamples;
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
    void Zoom(float steps);

private:
    void ProcessKeyboardInput(float deltaTime);
    void ProcessMouseInput();
    void PickPlant(double mouseX, double mouseY);
    bool RayIntersectsBounds(const glm::vec3& origin, const glm::vec3& direction, const Organism::BoundingBox& bounds, float& distance) const;
    bool RayIntersectsLimb(const glm::vec3& origin, const glm::vec3& direction, const RenderLimb& limb, float& distance) const;
    void DrawSoilNutrition();
    void DrawDebugBoundingBoxes();
    void DrawNeuralNetworkWindow();
    void DrawSelectedPlantWindow();
    glm::vec4 PlantColor(const RenderFrame::Plant& plant, double maxEnergy) const;
    glm::vec4 LimbColor(const RenderFrame::Plant& plant, const RenderLimb& limb, double maxEnergy) const;
    const RenderFrame::Plant* FindSelectedPlant() const;

    bool alive = false;
    GLFWwindow* window = nullptr;

    unsigned int limbvert;
    unsigned int limbVAO = 0;
    unsigned int limbVBO = 0;
    unsigned int limbEBO = 0;
    std::vector<RenderFrame::Plant> plants;

    unsigned int floorVAO = 0;
    unsigned int floorVBO = 0;
    unsigned int soilVAO = 0;
    unsigned int soilVBO = 0;

    unsigned int boxVAO = 0;
    unsigned int boxVBO = 0;
    std::vector<Organism::NeuralNetworkSample> neuralSamples;
    std::vector<RenderFrame::SoilCell> soilCells;

    enum class RenderMode
    {
        Default,
        SoilNutrition,
        PlantEnergy,
        PressureStrain
    };

    Shader shader;
    Camera camera;
    float moveSpeed = 2.5f;
    float mouseSensitivity = 0.15f;
    double lastFrameTime = 0.0;
    double lastMouseX = 400.0;
    double lastMouseY = 300.0;
    bool leftMouseDown = false;
    bool leftMouseDragged = false;
    double mouseDownX = 400.0;
    double mouseDownY = 300.0;
    bool debugBounds = false;
    bool debugKeyWasDown = false;
    bool imguiInitialized = false;
    RenderMode renderMode = RenderMode::Default;
    std::optional<unsigned int> selectedPlantId;
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
    std::condition_variable frameCondition;
    RenderFrame pendingFrame;
    bool hasPendingFrame = false;
    bool initialized = false;
    std::atomic<bool> running{false};
    std::thread renderThread;
};
