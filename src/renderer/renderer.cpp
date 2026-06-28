#include "renderer.hpp"

#include <glm/gtc/type_ptr.hpp>
#include <filesystem>
#include <iostream>

namespace {
    std::string ResolveShaderPath(const char* relativePath)
    {
        const std::filesystem::path candidates[] = {
            std::filesystem::path(relativePath),
            std::filesystem::path("..").append(relativePath),
            std::filesystem::path("src").append(relativePath)
        };

        for (const auto& candidate : candidates)
        {
            if (std::filesystem::exists(candidate))
                return candidate.string();
        }

        return relativePath;
    }
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

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cerr << "Failed to initialize GLAD\n";
        glfwDestroyWindow(window);
        window = nullptr;
        glfwTerminate();
        return;
    }

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

    float planeVertices[] = {
        -0.5f, -0.5f, 0.0f,
         0.5f, -0.5f, 0.0f,
         0.5f,  0.5f, 0.0f,

        -0.5f, -0.5f, 0.0f,
         0.5f,  0.5f, 0.0f,
        -0.5f,  0.5f, 0.0f
    };

    glGenVertexArrays(1, &planeVAO);
    glGenBuffers(1, &planeVBO);

    glBindVertexArray(planeVAO);
    glBindBuffer(GL_ARRAY_BUFFER, planeVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(planeVertices), planeVertices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);

    alive = true;
}

Renderer::~Renderer()
{
    if (planeVAO)
        glDeleteVertexArrays(1, &planeVAO);
    if (planeVBO)
        glDeleteBuffers(1, &planeVBO);

    if (window)
    {
        glfwDestroyWindow(window);
        window = nullptr;
    }

    glfwTerminate();
}

bool Renderer::Update()
{
    if (!window || glfwWindowShouldClose(window) || !alive)
    {
        alive = false;
        return false;
    }

    glClearColor(0.2f, 0.3f, 0.4f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable(GL_DEPTH_TEST);

    shader.Use();
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 view = camera.ViewMatrix();
    glm::mat4 projection = camera.ProjectionMatrix();
    shader.SetMat4("uModel", glm::value_ptr(model));
    shader.SetMat4("uView", glm::value_ptr(view));
    shader.SetMat4("uProjection", glm::value_ptr(projection));

    glBindVertexArray(planeVAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glBindVertexArray(0);

    glfwSwapBuffers(window);
    glfwPollEvents();
    return true;
}