#version 430 core
layout (location=0) in vec2 position;

uniform float radius;
uniform float x;
uniform float y;

uniform mat4 view;
uniform mat4 projection;
uniform mat4 model;

out float fragRadius;
out vec2 fragPosition;
out vec2 fragCenter;

void main()
{
    vec4 worldPosition = model * vec4(position, 0.0, 1.0);
    gl_Position = projection * view * worldPosition;
    fragPosition = worldPosition.xy;
    fragCenter.x = x;
    fragCenter.y = y;
    fragRadius = radius;
}