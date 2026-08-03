#version 330 core
layout (location = 0) in vec3 aPos;

uniform mat4 uModel;
uniform mat4 uView;
uniform mat4 uProjection;
uniform vec4 uColor;
uniform bool uUseUniformColor;


out vec4 vertexColor;

void main()
{
    gl_Position = uProjection * uView * uModel * vec4(aPos, 1.0);
    vertexColor = uUseUniformColor ? uColor : vec4(aPos, 1.0f);
}
