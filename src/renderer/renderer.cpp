#include "renderer.hpp"

#include <math.h>
#include <glm/gtc/type_ptr.hpp>
#include <chrono>
#include <filesystem>
#include <iostream>
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

    glViewport(0, 0, 800, 600);

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

    float planeSide = 100.0f;
    float floorVertices[]={
        planeSide, 0.0f, planeSide,
        -planeSide, 0.0f, planeSide,
        -planeSide, 0.0f, -planeSide,
        -planeSide, 0.0f, -planeSide,
        planeSide, 0.0f, planeSide,
        planeSide, 0.0f, -planeSide,
    };

    glGenVertexArrays(1, &floorVAO);
    glGenBuffers(1, &floorVBO);
    glBindVertexArray(floorVAO);
    glBindBuffer(GL_ARRAY_BUFFER, floorVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(floorVertices), floorVertices, GL_STATIC_DRAW);
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
            lastMouseX = mouseX;
            lastMouseY = mouseY;
            return;
        }

        float xoffset = static_cast<float>(mouseX - lastMouseX);
        float yoffset = static_cast<float>(lastMouseY - mouseY);
        lastMouseX = mouseX;
        lastMouseY = mouseY;

        xoffset *= mouseSensitivity;
        yoffset *= mouseSensitivity;

        camera.Rotate(xoffset, yoffset);
    }
    else
    {
        leftMouseDown = false;
    }
}

Renderer::~Renderer()
{
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
    DrawLimb({data.root, data.length, data.radius, data.orientation});
}

void Renderer::DrawLimb(const RenderLimb& limb)
{
    limbList.push_back(limb);
}

void Renderer::DrawBoundingBox(const Organism::BoundingBox& bounds)
{
    boundingBoxes.push_back(bounds);
}

void Renderer::SetFrame(RenderFrame frame)
{
    limbList = std::move(frame.limbs);
    boundingBoxes = std::move(frame.boundingBoxes);
}

bool Renderer::IsAlive() const
{
    return alive;
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
    ProcessKeyboardInput(deltaTime);
    ProcessMouseInput();

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

    glBindVertexArray(limbVAO);
    for(const RenderLimb& limb : limbList){
        model = glm::mat4(1.0f);
        model = glm::translate(model, limb.root);
        model *= glm::mat4_cast(limb.orientation);
        model = glm::scale(model, {limb.radius, limb.length, limb.radius});

        // glm::vec3 dir = data.orientation * glm::vec3(1,0,0);

        //     std::cout
        //         << dir.x << " "
        //         << dir.y << " "
        //         << dir.z << std::endl;
        
        shader.SetMat4("uModel", glm::value_ptr(model));
        glDrawElements(GL_TRIANGLES, limbvert, GL_UNSIGNED_INT, nullptr);
    }
    limbList.clear();
    glBindVertexArray(0);

    DrawDebugBoundingBoxes();
    boundingBoxes.clear();

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

void Renderer::DrawDebugBoundingBoxes()
{
    if (!debugBounds || boundingBoxes.empty())
        return;

    shader.SetInt("uUseUniformColor", 1);
    shader.SetVec4("uColor", 1.0f, 0.85f, 0.1f, 1.0f);
    shader.SetMat4("uModel", glm::value_ptr(glm::mat4(1.0f)));

    glBindVertexArray(boxVAO);
    glLineWidth(2.0f);

    for (const auto& bounds : boundingBoxes)
    {
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
