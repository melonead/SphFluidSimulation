#version 430 core

layout (location=0) in vec2 position;
layout (location=1) in vec2 quadVertices;


out vec2 fragPosition;
out vec2 center;

uniform mat4 view;
uniform mat4 projection;

void main()
{
    center = position;
    fragPosition = position + quadVertices;
    gl_Position = projection * view * vec4(quadVertices + position, 0.0, 1.0);
}