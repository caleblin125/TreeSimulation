#include "renderer.hpp"

#include <math.h>
#include <glm/gtc/type_ptr.hpp>
#include <filesystem>
#include <iostream>

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

    camera.SetPerspective(45.0f, 800.0f / 600.0f, 0.1f, 100.0f);
    camera.SetPosition(glm::vec3(0.0f, 1.5f, 3.0f));
    camera.SetTarget(glm::vec3(0.0f, 0.0f, 0.0f));

    glGenVertexArrays(1, &limbVAO);
    glGenBuffers(1, &limbVBO);
    glGenBuffers(1, &limbEBO);

    glBindVertexArray(limbVAO);
    limbvert = MakeLimb(20, limbEBO, limbVBO);
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
    limbList.push_back(limb);
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

    glBindVertexArray(limbVAO);
    for(auto& limb : limbList){
        if(limb.expired()){
            continue;
        }
        std::shared_ptr<Limb> limb_ptr = limb.lock();
        Limb::LimbData data = limb_ptr->GetData();
        limb_ptr.reset();

        glm::mat4 model = glm::mat4(1.0f);
        model = glm::translate(model, data.root);
        model *= glm::mat4_cast(data.orientation);
        model = glm::scale(model, {data.radius, data.length, data.radius});

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

    glfwSwapBuffers(window);
    return true;
}