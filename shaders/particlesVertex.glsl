#version 430 core

layout (location=0) in vec2 position;
layout (location=1) in vec2 quadVertices;
layout (location=2) in float gradientXCoordinate;


out vec2 fragPosition;
out vec2 center;
out float fragRadius;
out float fragXCoord;

uniform float radius;
uniform mat4 view;
uniform mat4 projection;
uniform mat4 model;

vec4 worldPosition;


void main()
{
    worldPosition = vec4(position, 0.0, 0.0) + model * vec4(quadVertices, 0.0, 1.0);
    gl_Position = projection * view * worldPosition;

    center = position;
    fragPosition = worldPosition.xy;
    fragRadius = radius;

    fragXCoord = gradientXCoordinate;
}