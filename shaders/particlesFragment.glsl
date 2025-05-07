#version 430 core

out vec4 color;

float radius = 0.05;

in vec2 center;
in vec2 fragPosition;

void main()
{
    if (length(fragPosition - center) < radius)
    {
        color = vec4(0.0, 0.0, 1.0, 1.0);
    }
    else
    {
        color = vec4(0.0, 0.0, 0.0, 0.0);
    }
}