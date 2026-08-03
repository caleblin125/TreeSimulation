#pragma once

#include <glad/glad.h>

#include <string>

class Shader {
public:
    Shader() = default;
    Shader(const std::string& vertexPath, const std::string& fragmentPath);
    ~Shader();

    Shader(const Shader&) = delete;
    Shader& operator=(const Shader&) = delete;
    Shader(Shader&& other) noexcept;
    Shader& operator=(Shader&& other) noexcept;

    bool IsValid() const;
    void Use() const;
    void SetMat4(const char* name, const float* value) const;
    void SetMat4(const char* name, const void* value) const;
    void SetVec4(const char* name, float x, float y, float z, float w) const;
    void SetInt(const char* name, int value) const;

private:
    unsigned int programId = 0;
};
