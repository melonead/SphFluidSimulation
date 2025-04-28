#version 430 core
layout (location=0) in vec2 quadPosition;

uniform float radius;
uniform float x;
uniform float y;

uniform mat4 view;
uniform mat4 projection;

out float fragRadius;
out vec2 fragPosition;
out vec2 fragCenter;
vec2 position;

void main()
{
    position.x = x + (quadPosition.x / 0.05) * radius;
    position.y = y + (quadPosition.y / 0.05) * radius;
    gl_Position = projection * view * vec4(position, 0.0, 1.0);
    fragRadius = radius;
    fragPosition = position;
    fragCenter.x = x;
    fragCenter.y = y;
}