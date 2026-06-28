#include "shader.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
    std::string ReadFile(const std::string& path)
    {
        std::ifstream file(path);
        if (!file)
            throw std::runtime_error("Unable to open shader file: " + path);

        std::ostringstream stream;
        stream << file.rdbuf();
        return stream.str();
    }

    unsigned int CompileShader(unsigned int type, const char* source)
    {
        unsigned int shader = glCreateShader(type);
        glShaderSource(shader, 1, &source, nullptr);
        glCompileShader(shader);

        int success = 0;
        glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
        if (!success)
        {
            char infoLog[512];
            glGetShaderInfoLog(shader, 512, nullptr, infoLog);
            std::cerr << "Shader compilation failed: " << infoLog << '\n';
            glDeleteShader(shader);
            return 0;
        }

        return shader;
    }
}

Shader::Shader(const std::string& vertexPath, const std::string& fragmentPath)
{
    try
    {
        std::string vertexSource = ReadFile(vertexPath);
        std::string fragmentSource = ReadFile(fragmentPath);

        const char* vertexCode = vertexSource.c_str();
        const char* fragmentCode = fragmentSource.c_str();

        unsigned int vertexShader = CompileShader(GL_VERTEX_SHADER, vertexCode);
        if (!vertexShader)
            return;

        unsigned int fragmentShader = CompileShader(GL_FRAGMENT_SHADER, fragmentCode);
        if (!fragmentShader)
        {
            glDeleteShader(vertexShader);
            return;
        }

        programId = glCreateProgram();
        glAttachShader(programId, vertexShader);
        glAttachShader(programId, fragmentShader);
        glLinkProgram(programId);

        int success = 0;
        glGetProgramiv(programId, GL_LINK_STATUS, &success);
        if (!success)
        {
            char infoLog[512];
            glGetProgramInfoLog(programId, 512, nullptr, infoLog);
            std::cerr << "Shader program linking failed: " << infoLog << '\n';
            glDeleteProgram(programId);
            programId = 0;
        }

        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);
    }
    catch (const std::exception& ex)
    {
        std::cerr << ex.what() << '\n';
    }
}

Shader::~Shader()
{
    if (programId != 0)
        glDeleteProgram(programId);
}

Shader::Shader(Shader&& other) noexcept : programId(other.programId)
{
    other.programId = 0;
}

Shader& Shader::operator=(Shader&& other) noexcept
{
    if (this != &other)
    {
        if (programId != 0)
            glDeleteProgram(programId);

        programId = other.programId;
        other.programId = 0;
    }
    return *this;
}

bool Shader::IsValid() const
{
    return programId != 0;
}

void Shader::Use() const
{
    if (programId != 0)
        glUseProgram(programId);
}

void Shader::SetMat4(const char* name, const float* value) const
{
    if (programId != 0)
    {
        int location = glGetUniformLocation(programId, name);
        if (location >= 0)
            glUniformMatrix4fv(location, 1, GL_FALSE, value);
    }
}

void Shader::SetMat4(const char* name, const void* value) const
{
    SetMat4(name, static_cast<const float*>(value));
}
