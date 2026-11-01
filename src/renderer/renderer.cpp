#include "renderer.hpp"

#include "simulation/simulation_constants.hpp"

#include <math.h>
#include <glm/gtc/type_ptr.hpp>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <implot.h>
#include <chrono>
#include <algorithm>
#include <cfloat>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string>
#include <utility>

namespace
{
    std::string ResolveShaderPath(const char *relativePath)
    {
        const std::filesystem::path candidates[] = {
            std::filesystem::path(relativePath),
            std::filesystem::path("..").append(relativePath),
            std::filesystem::path("src").append(relativePath)};

        for (const auto &candidate : candidates)
        {
            if (std::filesystem::exists(candidate))
                return candidate.string();
        }

        return relativePath;
    }
}

void framebuffer_size_callback(GLFWwindow *window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void scroll_callback(GLFWwindow *window, double, double yoffset)
{
    Renderer* renderer = static_cast<Renderer*>(glfwGetWindowUserPointer(window));
    if (renderer)
        renderer->Zoom(static_cast<float>(yoffset));
}

unsigned int MakeLimb(unsigned int N, unsigned int EBO, unsigned int VBO)
{
    // center, radial(N), end -> cone
    float limbVertices[3 * (N + 2)];
    memset(limbVertices, 0, sizeof(limbVertices));
    for (unsigned int i = 0; i < N; i++)
    {
        limbVertices[3 + 3 * i] = cos(M_PI * 2.0f * (float)i / (float)N);     // x
        limbVertices[3 + 3 * i + 2] = sin(M_PI * 2.0f * (float)i / (float)N); // z
    }
    limbVertices[3 + 3 * N + 1] = 1.0f;

    // radial triangles (N), cone triangles (N)
    unsigned int limbIndexes[3 * (N + N)];
    memset(limbIndexes, 0, sizeof(limbIndexes));
    for (unsigned int i = 0; i < N; i++)
    {
        limbIndexes[3 * i + 1] = i + 1;
        limbIndexes[3 * i + 2] = (i + 1) % N + 1;
    }

    for (unsigned int i = 0; i < N; i++)
    {
        limbIndexes[3 * N + 3 * i] = i + 1;
        limbIndexes[3 * N + 3 * i + 1] = (i + 1) % N + 1;
        limbIndexes[3 * N + 3 * i + 2] = 1 + N;
    }

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(limbVertices), limbVertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(limbIndexes), limbIndexes, GL_STATIC_DRAW);
    return sizeof(limbVertices);
}

Renderer::Renderer()
{
    if (!glfwInit())
    {
        std::cerr << "Failed to initialize GLFW\n";
        return;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    window = glfwCreateWindow(800, 600, "TreeSim", nullptr, nullptr);
    if (!window)
    {
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        return;
    }

    glfwMakeContextCurrent(window);
    glfwSetWindowUserPointer(window, this);
    glfwSwapInterval(1);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cerr << "Failed to initialize GLAD\n";
        glfwDestroyWindow(window);
        window = nullptr;
        glfwTerminate();
        return;
    }

    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetScrollCallback(window, scroll_callback);

    glViewport(0, 0, 800, 600);
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");
    imguiInitialized = true;

    const std::string vertexShaderPath = ResolveShaderPath("src/renderer/shaders/plane.vert");
    const std::string fragmentShaderPath = ResolveShaderPath("src/renderer/shaders/plane.frag");

    shader = Shader(vertexShaderPath, fragmentShaderPath);
    if (!shader.IsValid())
    {
        std::cerr << "Failed to load plane shader\n";
        return;
    }

    camera.SetPerspective(45.0f, 800.0f / 600.0f, 0.1f, 200.0f);
    // camera.SetPosition(glm::vec3(0.0f, 1.5f, 3.0f));
    // camera.SetTarget(glm::vec3(0.0f, 0.0f, 0.0f));

    glGenVertexArrays(1, &limbVAO);
    glGenBuffers(1, &limbVBO);
    glGenBuffers(1, &limbEBO);

    glBindVertexArray(limbVAO);
    limbvert = MakeLimb(20, limbEBO, limbVBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);

    glEnableVertexAttribArray(0);
    glBindVertexArray(0);

    float floorVertices[]={
        SimulationConstants::SimulationMaxX, 0.0f, SimulationConstants::SimulationMaxZ,
        SimulationConstants::SimulationMinX, 0.0f, SimulationConstants::SimulationMaxZ,
        SimulationConstants::SimulationMinX, 0.0f, SimulationConstants::SimulationMinZ,
        SimulationConstants::SimulationMinX, 0.0f, SimulationConstants::SimulationMinZ,
        SimulationConstants::SimulationMaxX, 0.0f, SimulationConstants::SimulationMaxZ,
        SimulationConstants::SimulationMaxX, 0.0f, SimulationConstants::SimulationMinZ,
    };

    glGenVertexArrays(1, &floorVAO);
    glGenBuffers(1, &floorVBO);
    glBindVertexArray(floorVAO);
    glBindBuffer(GL_ARRAY_BUFFER, floorVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(floorVertices), floorVertices, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);

    glEnableVertexAttribArray(0);
    glBindVertexArray(0);

    glGenVertexArrays(1, &soilVAO);
    glGenBuffers(1, &soilVBO);
    glBindVertexArray(soilVAO);
    glBindBuffer(GL_ARRAY_BUFFER, soilVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 * 6, nullptr, GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);

    glGenVertexArrays(1, &boxVAO);
    glGenBuffers(1, &boxVBO);

    glBindVertexArray(boxVAO);
    glBindBuffer(GL_ARRAY_BUFFER, boxVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 * 24, nullptr, GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);

    alive = true;
    lastFrameTime = glfwGetTime();
}

void Renderer::ProcessKeyboardInput(float deltaTime)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    const bool debugKeyDown = glfwGetKey(window, GLFW_KEY_B) == GLFW_PRESS;
    if (debugKeyDown && !debugKeyWasDown)
    {
        debugBounds = !debugBounds;
        std::cout << "Debug bounding boxes: " << (debugBounds ? "on" : "off") << std::endl;
    }
    debugKeyWasDown = debugKeyDown;

    const glm::vec3 forward = glm::normalize(camera.GetFront());
    const glm::vec3 right = glm::normalize(glm::cross(forward, glm::vec3(0.0f, 1.0f, 0.0f)));
    const glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);

    glm::vec3 movement(0.0f);
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        movement += forward;
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        movement -= forward;
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        movement -= right;
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        movement += right;
    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
        movement -= up;
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
        movement += up;

    if (glm::length(movement) > 0.0f)
        camera.Move(glm::normalize(movement) * moveSpeed * deltaTime);
}

void Renderer::ProcessMouseInput()
{
    double mouseX, mouseY;
    glfwGetCursorPos(window, &mouseX, &mouseY);
    bool isDown = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;

    if (isDown)
    {
        if (!leftMouseDown)
        {
            leftMouseDown = true;
            leftMouseDragged = false;
            mouseDownX = mouseX;
            mouseDownY = mouseY;
            lastMouseX = mouseX;
            lastMouseY = mouseY;
            return;
        }

        float xoffset = static_cast<float>(mouseX - lastMouseX);
        float yoffset = static_cast<float>(lastMouseY - mouseY);
        if (std::abs(mouseX - mouseDownX) > 3.0 || std::abs(mouseY - mouseDownY) > 3.0)
            leftMouseDragged = true;
        lastMouseX = mouseX;
        lastMouseY = mouseY;

        xoffset *= mouseSensitivity;
        yoffset *= mouseSensitivity;

        camera.Rotate(xoffset, yoffset);
    }
    else
    {
        if (leftMouseDown && !leftMouseDragged && !ImGui::GetIO().WantCaptureMouse)
            PickPlant(mouseX, mouseY);
        leftMouseDown = false;
    }
}

Renderer::~Renderer()
{
    if (imguiInitialized)
    {
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImPlot::DestroyContext();
        ImGui::DestroyContext();
    }

    if (limbVAO)
        glDeleteVertexArrays(1, &limbVAO);
    if (limbVBO)
        glDeleteBuffers(1, &limbVBO);
    if (limbEBO)
        glDeleteBuffers(1, &limbEBO);
    if (boxVAO)
        glDeleteVertexArrays(1, &boxVAO);
    if (boxVBO)
        glDeleteBuffers(1, &boxVBO);
    if (soilVAO)
        glDeleteVertexArrays(1, &soilVAO);
    if (soilVBO)
        glDeleteBuffers(1, &soilVBO);

    if (window)
    {
        glfwDestroyWindow(window);
        window = nullptr;
    }

    shader = Shader();
    glfwTerminate();
}

void Renderer::DrawLimb(std::weak_ptr<Limb> limb)
{
    const std::shared_ptr<Limb> limbPtr = limb.lock();
    if (!limbPtr)
        return;

    const Limb::LimbData data = limbPtr->GetData();
    DrawLimb({data.root, data.length, data.width, data.thickness, data.orientation});
}

void Renderer::DrawLimb(const RenderLimb& limb)
{
    RenderFrame::Plant loosePlant;
    loosePlant.id = 0;
    loosePlant.limbs.push_back(limb);
    plants.push_back(std::move(loosePlant));
}

void Renderer::DrawBoundingBox(const Organism::BoundingBox& bounds)
{
    RenderFrame::Plant loosePlant;
    loosePlant.id = 0;
    loosePlant.bounds = bounds;
    plants.push_back(std::move(loosePlant));
}

void Renderer::SetFrame(RenderFrame frame)
{
    plants = std::move(frame.plants);
    soilCells = std::move(frame.soilCells);
    neuralSamples = std::move(frame.neuralSamples);
    latestMetrics = std::move(frame.metrics);
    if (latestMetrics.timeYears > lastMetricTime)
    {
        lastMetricTime = latestMetrics.timeYears;
        populationHistory.push_back(latestMetrics.population);
        avgEnergyHistory.push_back(latestMetrics.avgEnergy);
        avgGrowthHistory.push_back(latestMetrics.avgGrowth);
        avgReproductionHistory.push_back(latestMetrics.avgReproduction);
        lineageHistory.push_back(latestMetrics.lineages);
        effectiveLineageHistory.push_back(latestMetrics.effectiveLineages);
        constexpr std::size_t maxHistory = 600;
        auto trim = [&](std::vector<float>& values)
        {
            if (values.size() > maxHistory)
                values.erase(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(values.size() - maxHistory));
        };
        trim(populationHistory);
        trim(avgEnergyHistory);
        trim(avgGrowthHistory);
        trim(avgReproductionHistory);
        trim(lineageHistory);
        trim(effectiveLineageHistory);
    }
    if (selectedPlantId && !FindSelectedPlant())
        selectedPlantId.reset();
}

bool Renderer::IsAlive() const
{
    return alive;
}

void Renderer::Zoom(float steps)
{
    const float distance = std::clamp(steps, -5.0f, 5.0f) * moveSpeed * 1.5f;
    camera.Move(glm::normalize(camera.GetFront()) * distance);
}

void Renderer::PickPlant(double mouseX, double mouseY)
{
    int windowWidth = 800;
    int windowHeight = 600;
    glfwGetWindowSize(window, &windowWidth, &windowHeight);
    if (windowWidth <= 0 || windowHeight <= 0)
        return;

    const float x = (2.0f * static_cast<float>(mouseX)) / static_cast<float>(windowWidth) - 1.0f;
    const float y = 1.0f - (2.0f * static_cast<float>(mouseY)) / static_cast<float>(windowHeight);
    const glm::mat4 inverseViewProjection = glm::inverse(camera.ProjectionMatrix() * camera.ViewMatrix());
    const glm::vec4 nearPoint = inverseViewProjection * glm::vec4(x, y, -1.0f, 1.0f);
    const glm::vec4 farPoint = inverseViewProjection * glm::vec4(x, y, 1.0f, 1.0f);
    const glm::vec3 rayOrigin = glm::vec3(nearPoint) / nearPoint.w;
    const glm::vec3 rayFar = glm::vec3(farPoint) / farPoint.w;
    const glm::vec3 rayDirection = glm::normalize(rayFar - rayOrigin);

    std::optional<unsigned int> pickedId;
    float bestDistance = std::numeric_limits<float>::max();
    for (const RenderFrame::Plant& plant : plants)
    {
        float boundsDistance = 0.0f;
        if (!RayIntersectsBounds(rayOrigin, rayDirection, plant.bounds, boundsDistance))
            continue;

        for (const RenderLimb& limb : plant.limbs)
        {
            float limbDistance = 0.0f;
            if (RayIntersectsLimb(rayOrigin, rayDirection, limb, limbDistance) && limbDistance < bestDistance)
            {
                bestDistance = limbDistance;
                pickedId = plant.id;
            }
        }
    }

    selectedPlantId = pickedId;
    if (selectedPlantId)
        std::cout << "Selected plant " << *selectedPlantId << std::endl;
    else
        std::cout << "Selected plant: none" << std::endl;
}

bool Renderer::RayIntersectsBounds(const glm::vec3& origin, const glm::vec3& direction, const Organism::BoundingBox& bounds, float& distance) const
{
    float tMin = 0.0f;
    float tMax = std::numeric_limits<float>::max();
    for (int axis = 0; axis < 3; ++axis)
    {
        const float originAxis = origin[axis];
        const float directionAxis = direction[axis];
        const float minAxis = bounds.min[axis];
        const float maxAxis = bounds.max[axis];

        if (std::abs(directionAxis) < 1.0e-6f)
        {
            if (originAxis < minAxis || originAxis > maxAxis)
                return false;
            continue;
        }

        float t1 = (minAxis - originAxis) / directionAxis;
        float t2 = (maxAxis - originAxis) / directionAxis;
        if (t1 > t2)
            std::swap(t1, t2);

        tMin = std::max(tMin, t1);
        tMax = std::min(tMax, t2);
        if (tMin > tMax)
            return false;
    }

    distance = tMin;
    return true;
}

bool Renderer::RayIntersectsLimb(const glm::vec3& origin, const glm::vec3& direction, const RenderLimb& limb, float& distance) const
{
    const glm::vec3 segmentStart = limb.root;
    const glm::vec3 segmentEnd = limb.root + limb.orientation * glm::vec3(0.0f, limb.length, 0.0f);
    const glm::vec3 segment = segmentEnd - segmentStart;
    const float segmentLengthSq = glm::dot(segment, segment);
    if (segmentLengthSq <= 1.0e-8f)
        return false;

    const glm::vec3 w0 = origin - segmentStart;
    const float a = 1.0f;
    const float b = glm::dot(direction, segment);
    const float c = segmentLengthSq;
    const float d = glm::dot(direction, w0);
    const float e = glm::dot(segment, w0);
    const float denominator = a * c - b * b;

    float rayT = 0.0f;
    float segmentT = 0.0f;
    if (std::abs(denominator) > 1.0e-6f)
        rayT = std::max(0.0f, (b * e - c * d) / denominator);
    segmentT = std::clamp((b * rayT + e) / c, 0.0f, 1.0f);
    rayT = std::max(0.0f, (b * segmentT - d) / a);

    const glm::vec3 rayPoint = origin + direction * rayT;
    const glm::vec3 segmentPoint = segmentStart + segment * segmentT;
    const float hitRadius = std::max(0.5f * std::max(limb.width, limb.thickness), 0.05f);
    if (glm::length(rayPoint - segmentPoint) > hitRadius)
        return false;

    distance = rayT;
    return true;
}

glm::vec4 Renderer::PlantColor(const RenderFrame::Plant& plant, double maxEnergy) const
{
    const float t = static_cast<float>(std::clamp(plant.energy / maxEnergy, 0.0, 1.0));
    return glm::vec4(1.0f - t, 0.18f + 0.75f * t, 0.12f, 1.0f);
}

glm::vec4 Renderer::LimbColor(const RenderFrame::Plant& plant, const RenderLimb& limb, double maxEnergy) const
{
    if (renderMode == RenderMode::PlantEnergy)
        return PlantColor(plant, maxEnergy);

    if (renderMode == RenderMode::PressureStrain)
    {
        const float t = std::clamp(limb.strain, 0.0f, 1.0f);
        return glm::vec4(0.2f + 0.8f * t, 0.95f * (1.0f - t), 0.12f, 1.0f);
    }

    if (renderMode == RenderMode::Lineage)
    {
        // Golden-ratio steps around the colour wheel keep neighbouring ids far apart.
        const float hue = glm::fract(static_cast<float>(plant.lineage) * 0.61803398875f) * 6.0f;
        const float value = 0.65f + 0.35f * glm::fract(static_cast<float>(plant.lineage) * 0.38196601125f);
        const auto channel = [&](float offset)
        {
            const float k = glm::mod(offset + hue, 6.0f);
            return value * (1.0f - 0.8f * std::clamp(std::min(k, 4.0f - k), 0.0f, 1.0f));
        };
        return glm::vec4(channel(5.0f), channel(3.0f), channel(1.0f), 1.0f);
    }

    if (IsLeaf(limb))
        return glm::vec4(0.25f, 0.72f, 0.22f, 1.0f);
    return glm::vec4(1.0f);
}

bool Renderer::IsLeaf(const RenderLimb& limb)
{
    const float narrow = std::max(1.0e-6f, std::min(limb.width, limb.thickness));
    return std::max(limb.width, limb.thickness) / narrow >= SimulationConstants::LeafAspectRatio;
}

void Renderer::DrawSoilNutrition()
{
    if (renderMode != RenderMode::SoilNutrition || soilCells.empty())
        return;

    shader.SetInt("uUseUniformColor", 1);
    shader.SetMat4("uModel", glm::value_ptr(glm::mat4(1.0f)));
    glBindVertexArray(soilVAO);

    constexpr float cellSize = static_cast<float>(SimulationConstants::SoilCellSizeM);
    constexpr float y = 0.015f;
    for (const RenderFrame::SoilCell& cell : soilCells)
    {
        const float x0 = static_cast<float>(cell.cellX) * cellSize;
        const float z0 = static_cast<float>(cell.cellZ) * cellSize;
        const float x1 = x0 + cellSize;
        const float z1 = z0 + cellSize;
        const float vertices[] = {
            x0, y, z0,
            x1, y, z0,
            x1, y, z1,
            x0, y, z0,
            x1, y, z1,
            x0, y, z1,
        };

        const float resource = std::clamp(cell.resource, 0.0f, 1.0f);
        shader.SetVec4("uColor", 0.70f * (1.0f - resource), 0.22f + 0.65f * resource, 0.08f, 0.85f);
        glBindBuffer(GL_ARRAY_BUFFER, soilVBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices);
        glDrawArrays(GL_TRIANGLES, 0, 6);
    }

    glBindVertexArray(0);
    shader.SetInt("uUseUniformColor", 0);
}

const RenderFrame::Plant* Renderer::FindSelectedPlant() const
{
    if (!selectedPlantId)
        return nullptr;

    const auto selected = std::find_if(plants.begin(), plants.end(), [&](const RenderFrame::Plant& plant)
    {
        return plant.id == *selectedPlantId;
    });
    return selected == plants.end() ? nullptr : &*selected;
}

bool Renderer::Update()
{
    if (!window || !alive)
        return false;

    glfwPollEvents();
    if (glfwWindowShouldClose(window))
    {
        alive = false;
        return false;
    }

    const float currentTime = static_cast<float>(glfwGetTime());
    const float deltaTime = static_cast<float>(currentTime - lastFrameTime);
    lastFrameTime = currentTime;

    int framebufferWidth = 800;
    int framebufferHeight = 600;
    glfwGetFramebufferSize(window, &framebufferWidth, &framebufferHeight);
    const float aspect = framebufferHeight > 0 ? static_cast<float>(framebufferWidth) / static_cast<float>(framebufferHeight) : 1.0f;
    camera.SetPerspective(45.0f, aspect, 0.1f, 200.0f);

    ProcessKeyboardInput(deltaTime);
    ProcessMouseInput();

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    glClearColor(0.2f, 0.3f, 0.4f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable(GL_DEPTH_TEST);

    shader.Use();
    glm::mat4 view = camera.ViewMatrix();
    glm::mat4 projection = camera.ProjectionMatrix();

    
    shader.SetMat4("uView", glm::value_ptr(view));
    shader.SetMat4("uProjection", glm::value_ptr(projection));
    shader.SetInt("uUseUniformColor", 0);

    glBindVertexArray(floorVAO);
    glm::mat4 model = glm::mat4(1.0f);
    shader.SetMat4("uModel", glm::value_ptr(model));
    glDrawArrays(GL_TRIANGLES, 0, 6);
    DrawSoilNutrition();

    glBindVertexArray(limbVAO);
    const auto maxEnergyIt = std::max_element(plants.begin(), plants.end(), [](const RenderFrame::Plant& a, const RenderFrame::Plant& b)
    {
        return a.energy < b.energy;
    });
    const double maxEnergy = maxEnergyIt == plants.end() ? 1.0 : std::max(1.0, maxEnergyIt->energy);

    for(const RenderFrame::Plant& plant : plants){
        const bool selected = selectedPlantId && plant.id == *selectedPlantId;

        for (const RenderLimb& limb : plant.limbs)
        {
            const glm::vec4 color = selected
                ? glm::vec4(0.25f, 1.0f, 0.35f, 1.0f)
                : LimbColor(plant, limb, maxEnergy);
            // Default mode keeps the position colouring for wood and paints leaves green.
            const bool useUniformColor = selected || renderMode != RenderMode::Default || IsLeaf(limb);
            shader.SetInt("uUseUniformColor", useUniformColor ? 1 : 0);
            if (useUniformColor)
                shader.SetVec4("uColor", color.r, color.g, color.b, color.a);

            model = glm::mat4(1.0f);
            model = glm::translate(model, limb.root);
            model *= glm::mat4_cast(limb.orientation);
            model = glm::scale(model, {0.5f * limb.width, limb.length, 0.5f * limb.thickness});

            shader.SetMat4("uModel", glm::value_ptr(model));
            glDrawElements(GL_TRIANGLES, limbvert, GL_UNSIGNED_INT, nullptr);
        }
    }
    shader.SetInt("uUseUniformColor", 0);
    glBindVertexArray(0);

    DrawDebugBoundingBoxes();

    DrawNeuralNetworkWindow();
    DrawSelectedPlantWindow();
    DrawMetricsWindow();
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    glfwSwapBuffers(window);
    return true;
}

RenderPipeline::RenderPipeline()
{
    running = true;
    renderThread = std::thread(&RenderPipeline::RenderLoop, this);

    std::unique_lock<std::mutex> lock(frameMutex);
    frameCondition.wait(lock, [this]
    {
        return initialized;
    });
}

RenderPipeline::~RenderPipeline()
{
    Stop();
}

bool RenderPipeline::IsRunning() const
{
    return running.load();
}

void RenderPipeline::SubmitFrame(RenderFrame frame)
{
    {
        std::lock_guard<std::mutex> lock(frameMutex);
        pendingFrame = std::move(frame);
        hasPendingFrame = true;
    }
    frameCondition.notify_one();
}

void RenderPipeline::Stop()
{
    running = false;
    frameCondition.notify_one();
    if (renderThread.joinable())
        renderThread.join();
}

void RenderPipeline::RenderLoop()
{
    Renderer renderer;
    RenderFrame currentFrame;
    bool hasCurrentFrame = false;

    {
        std::lock_guard<std::mutex> lock(frameMutex);
        initialized = true;
        if (!renderer.IsAlive())
            running = false;
    }
    frameCondition.notify_all();

    while (running.load() && renderer.IsAlive())
    {
        {
            std::unique_lock<std::mutex> lock(frameMutex);
            if (!hasCurrentFrame)
            {
                frameCondition.wait(lock, [this]
                {
                    return hasPendingFrame || !running.load();
                });
            }
            else if (!hasPendingFrame)
            {
                frameCondition.wait_for(lock, std::chrono::milliseconds(1), [this]
                {
                    return hasPendingFrame || !running.load();
                });
            }

            if (hasPendingFrame)
            {
                currentFrame = std::move(pendingFrame);
                pendingFrame = {};
                hasPendingFrame = false;
                hasCurrentFrame = true;
            }
        }

        if (!hasCurrentFrame)
            continue;

        renderer.SetFrame(currentFrame);
        if (!renderer.Update())
            break;
    }

    running = false;
}

void Renderer::DrawNeuralNetworkWindow()
{
    ImGui::Begin("Neural network samples");
    int mode = static_cast<int>(renderMode);
    ImGui::RadioButton("Default", &mode, static_cast<int>(RenderMode::Default));
    ImGui::SameLine();
    ImGui::RadioButton("Soil", &mode, static_cast<int>(RenderMode::SoilNutrition));
    ImGui::SameLine();
    ImGui::RadioButton("Energy", &mode, static_cast<int>(RenderMode::PlantEnergy));
    ImGui::SameLine();
    ImGui::RadioButton("Strain", &mode, static_cast<int>(RenderMode::PressureStrain));
    ImGui::SameLine();
    ImGui::RadioButton("Lineage", &mode, static_cast<int>(RenderMode::Lineage));
    renderMode = static_cast<RenderMode>(mode);
    ImGui::Text("Samples: %zu", neuralSamples.size());
    if (selectedPlantId)
        ImGui::Text("Selected plant: %u", *selectedPlantId);
    else
        ImGui::TextUnformatted("Selected plant: none");

    // Direct weights only; the hidden layer is not shown.
    static const char* outputNames[] = {"elongate", "thicken", "root", "split", "branch n", "attach", "reproduce", "widen", "child class", "fine bud"};
    static const char* inputNames[] = {"bias", "energy", "resource", "sun", "underground", "length", "torque", "age", "height", "depth", "tip height", "tip", "class"};

    for (std::size_t sampleIndex = 0; sampleIndex < neuralSamples.size(); ++sampleIndex)
    {
        const Organism::NeuralNetworkSample& sample = neuralSamples[sampleIndex];
        const std::string label = "Organism " + std::to_string(sample.id);
        if (!ImGui::CollapsingHeader(label.c_str(), sampleIndex == 0 ? ImGuiTreeNodeFlags_DefaultOpen : 0))
            continue;

        ImGui::Text("age %.2f years   energy %.1f J", sample.ageYears, sample.energy);
        ImGui::Text("limbs %u   height %.2f m   resource %.2f", sample.stats.limbCount, sample.stats.maxHeight, sample.stats.resourceFactor);
        ImGui::Text("root %.2f   hydraulic %.2f   xylem %.2f   embolism %.2f",
            sample.stats.rootResourceFactor,
            sample.stats.hydraulicResourceFactor,
            sample.stats.xylemConductanceMultiplier,
            sample.stats.embolismLossFraction);
        ImGui::Text("reproduce %.2f   budget %.1f J   split signals %u",
            sample.stats.reproductionPressure,
            sample.stats.reproductionBudget,
            sample.stats.splitSignals);
        ImGui::Text("branches %u   attachment %.2f", sample.stats.branchesSpawned, sample.stats.branchAttachment);
        ImGui::Text("photo %.1f J   transport %.1f J   respiration %.1f J   reserve %.1f J",
            sample.stats.photosynthesisEnergy,
            sample.stats.nutrientTransportEnergy,
            sample.stats.respirationEnergy,
            sample.stats.reserveDecayEnergy);

        if (ImGui::BeginTable(("weights##" + std::to_string(sample.id)).c_str(), static_cast<int>(Organism::GrowthBrainInputCount + 2), ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg))
        {
            ImGui::TableSetupColumn("output");
            for (const char* inputName : inputNames)
                ImGui::TableSetupColumn(inputName);
            ImGui::TableHeadersRow();

            for (std::size_t outputIndex = 0; outputIndex < Organism::GrowthBrainOutputCount; ++outputIndex)
            {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                ImGui::TextUnformatted(outputNames[outputIndex]);
                const std::size_t weightOffset = outputIndex * (Organism::GrowthBrainInputCount + 1);
                for (std::size_t inputIndex = 0; inputIndex < Organism::GrowthBrainInputCount + 1; ++inputIndex)
                {
                    const double weight = sample.weights[weightOffset + inputIndex];
                    ImGui::TableSetColumnIndex(static_cast<int>(inputIndex + 1));
                    ImGui::TextColored(
                        weight >= 0.0 ? ImVec4(0.45f, 0.95f, 0.55f, 1.0f) : ImVec4(0.95f, 0.45f, 0.45f, 1.0f),
                        "%.2f",
                        weight);
                }
            }
            ImGui::EndTable();
        }
    }

    ImGui::End();
}

void Renderer::DrawSelectedPlantWindow()
{
    const RenderFrame::Plant* plant = FindSelectedPlant();
    if (!plant)
        return;

    const Organism::Stats& stats = plant->stats;
    ImGui::Begin("Selected organism");
    ImGui::Text("id %u   lineage %u   age %.2f years   limbs %u", plant->id, plant->lineage, stats.ageYears, stats.limbCount);
    ImGui::Text("energy %.1f J   start %.1f J   alive %s", plant->energy, stats.energyStart, stats.alive ? "yes" : "no");
    ImGui::Separator();
    ImGui::Text("height %.2f m   leaf %.3f m2   root %.3f m2", stats.maxHeight, stats.leafArea, stats.rootSurfaceArea);
    ImGui::Text("leaves %u (%.0f%% of light)   fruiting tips %u   shed %u",
        stats.leafCount,
        100.0 * stats.leafShare,
        stats.fruitingSites,
        stats.limbsShed);
    ImGui::Text("biomass %.3f kg   root %.3f kg in %u limbs   anchorage margin %.2f",
        stats.biomassKg,
        stats.rootBiomassKg,
        stats.rootLimbCount,
        stats.anchorageMargin);
    ImGui::Text("resource %.2f   soil %.2f   root %.2f   hydraulic %.2f",
        stats.resourceFactor,
        stats.soilResourceFactor,
        stats.rootResourceFactor,
        stats.hydraulicResourceFactor);
    ImGui::Text("soil pressure %.3f   sun occlusion %.3f", stats.soilCompetitionPressure, stats.sunlightOcclusionPressure);
    ImGui::Text("xylem %.2f   embolism %.2f   bend %.1f MPa",
        stats.xylemConductanceMultiplier,
        stats.embolismLossFraction,
        stats.effectiveBendingStrength / 1.0e6);
    ImGui::Separator();
    ImGui::Text("photo %.1f J   transport %.1f J", stats.photosynthesisEnergy, stats.nutrientTransportEnergy);
    ImGui::Text("resp %.1f J   reserve decay %.1f J", stats.respirationEnergy, stats.reserveDecayEnergy);
    ImGui::Text("shoot grow %.4f m / %.1f J", stats.shootGrowthLength, stats.shootGrowthCost);
    ImGui::Text("root grow %.4f m / %.1f J", stats.rootGrowthLength, stats.rootGrowthCost);
    ImGui::Text("fracture loss %.1f J", stats.fractureEnergyLoss);
    ImGui::Separator();
    ImGui::Text("reproduce %.2f   budget %.1f J   split signals %u",
        stats.reproductionPressure,
        stats.reproductionBudget,
        stats.splitSignals);
    ImGui::Text("branches %u   attachment %.2f", stats.branchesSpawned, stats.branchAttachment);
    ImGui::Text("bounds min %.2f %.2f %.2f", plant->bounds.min.x, plant->bounds.min.y, plant->bounds.min.z);
    ImGui::Text("bounds max %.2f %.2f %.2f", plant->bounds.max.x, plant->bounds.max.y, plant->bounds.max.z);
    ImGui::End();
}

void Renderer::DrawMetricsWindow()
{
    ImGui::Begin("Simulation metrics");
    ImGui::Text("year %.2f   population %.0f   propagules %.0f",
        latestMetrics.timeYears,
        latestMetrics.population,
        latestMetrics.propagules);
    ImGui::Text("avg energy %.1f J   avg height %.2f m", latestMetrics.avgEnergy, latestMetrics.avgHeight);
    ImGui::Text("avg growth %.4f m/step   avg reproduce %.2f", latestMetrics.avgGrowth, latestMetrics.avgReproduction);
    ImGui::Text("lineages %.0f   effective %.1f", latestMetrics.lineages, latestMetrics.effectiveLineages);
    if (ImGui::TreeNode("trait variation (sd / mean)"))
    {
        for (std::size_t trait = 0; trait < Organism::TraitCount; ++trait)
            ImGui::Text("%-12s %.3f", Organism::TraitNames[trait], latestMetrics.traitVariation[trait]);
        ImGui::TreePop();
    }

    if (ImPlot::BeginPlot("population", ImVec2(-1.0f, 120.0f), ImPlotFlags_CanvasOnly))
    {
        ImPlot::SetupAxes(nullptr, nullptr, ImPlotAxisFlags_NoTickLabels, ImPlotAxisFlags_AutoFit);
        if (!populationHistory.empty())
            ImPlot::PlotLine("trees", populationHistory.data(), static_cast<int>(populationHistory.size()));
        ImPlot::EndPlot();
    }

    if (ImPlot::BeginPlot("lineages", ImVec2(-1.0f, 120.0f)))
    {
        ImPlot::SetupAxes(nullptr, nullptr, ImPlotAxisFlags_NoTickLabels, ImPlotAxisFlags_AutoFit);
        if (!lineageHistory.empty())
            ImPlot::PlotLine("surviving", lineageHistory.data(), static_cast<int>(lineageHistory.size()));
        if (!effectiveLineageHistory.empty())
            ImPlot::PlotLine("effective", effectiveLineageHistory.data(), static_cast<int>(effectiveLineageHistory.size()));
        ImPlot::EndPlot();
    }

    if (ImPlot::BeginPlot("average energy", ImVec2(-1.0f, 120.0f), ImPlotFlags_CanvasOnly))
    {
        ImPlot::SetupAxes(nullptr, nullptr, ImPlotAxisFlags_NoTickLabels, ImPlotAxisFlags_AutoFit);
        if (!avgEnergyHistory.empty())
            ImPlot::PlotLine("J", avgEnergyHistory.data(), static_cast<int>(avgEnergyHistory.size()));
        ImPlot::EndPlot();
    }

    if (ImPlot::BeginPlot("growth and reproduction", ImVec2(-1.0f, 140.0f)))
    {
        ImPlot::SetupAxes("sample", "normalized", ImPlotAxisFlags_NoTickLabels, ImPlotAxisFlags_AutoFit);
        if (!avgGrowthHistory.empty())
            ImPlot::PlotLine("growth", avgGrowthHistory.data(), static_cast<int>(avgGrowthHistory.size()));
        if (!avgReproductionHistory.empty())
            ImPlot::PlotLine("reproduction", avgReproductionHistory.data(), static_cast<int>(avgReproductionHistory.size()));
        ImPlot::EndPlot();
    }

    DrawScatterPlot("growth vs reproduction", latestMetrics.growthVsReproduction, "growth", "reproduce");
    DrawScatterPlot("energy vs height", latestMetrics.energyVsHeight, "energy", "height");
    ImGui::End();
}

void Renderer::DrawScatterPlot(const char* label, const std::vector<RenderFrame::ScatterPoint>& points, const char* xLabel, const char* yLabel)
{
    if (points.empty() || !ImPlot::BeginPlot(label, ImVec2(-1.0f, 160.0f)))
        return;

    std::vector<float> xValues;
    std::vector<float> yValues;
    xValues.reserve(points.size());
    yValues.reserve(points.size());
    for (const RenderFrame::ScatterPoint& point : points)
    {
        xValues.push_back(point.x);
        yValues.push_back(point.y);
    }

    ImPlot::SetupAxes(xLabel, yLabel, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
    ImPlot::PlotScatter(label, xValues.data(), yValues.data(), static_cast<int>(points.size()));
    ImPlot::EndPlot();
}

void Renderer::DrawDebugBoundingBoxes()
{
    if (plants.empty())
        return;

    shader.SetInt("uUseUniformColor", 1);
    shader.SetMat4("uModel", glm::value_ptr(glm::mat4(1.0f)));

    glBindVertexArray(boxVAO);
    glLineWidth(2.0f);

    for (const auto& plant : plants)
    {
        const bool selected = selectedPlantId && plant.id == *selectedPlantId;
        if (!debugBounds && !selected)
            continue;

        shader.SetVec4(
            "uColor",
            selected ? 0.2f : 1.0f,
            selected ? 1.0f : 0.85f,
            selected ? 0.3f : 0.1f,
            1.0f);

        const auto& bounds = plant.bounds;
        const glm::vec3& min = bounds.min;
        const glm::vec3& max = bounds.max;

        const float vertices[] = {
            min.x, min.y, min.z, max.x, min.y, min.z,
            max.x, min.y, min.z, max.x, min.y, max.z,
            max.x, min.y, max.z, min.x, min.y, max.z,
            min.x, min.y, max.z, min.x, min.y, min.z,

            min.x, max.y, min.z, max.x, max.y, min.z,
            max.x, max.y, min.z, max.x, max.y, max.z,
            max.x, max.y, max.z, min.x, max.y, max.z,
            min.x, max.y, max.z, min.x, max.y, min.z,

            min.x, min.y, min.z, min.x, max.y, min.z,
            max.x, min.y, min.z, max.x, max.y, min.z,
            max.x, min.y, max.z, max.x, max.y, max.z,
            min.x, min.y, max.z, min.x, max.y, max.z,
        };

        glBindBuffer(GL_ARRAY_BUFFER, boxVBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices);
        glDrawArrays(GL_LINES, 0, 24);
    }

    glLineWidth(1.0f);
    glBindVertexArray(0);
    shader.SetInt("uUseUniformColor", 0);
}
