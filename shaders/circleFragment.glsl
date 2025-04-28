#version 430 core

out vec4 color;

in float fragRadius;
in vec2 fragPosition;
in vec2 fragCenter;

float thickness;

void main()
{   
    thickness = 0.025;
    float len = length(fragPosition - fragCenter);
    if (len > (fragRadius - thickness) && len < fragRadius )
    {
        color = vec4(0.0, 0.0, 1.0, 1.0);
    }
    else
    {
        color = vec4(1.0, 1.0, 0.0, 0.0);
    }
}