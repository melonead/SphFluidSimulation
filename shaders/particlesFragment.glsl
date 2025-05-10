#version 430 core

out vec4 color;

in float fragRadius;

in vec2 center;
in vec2 fragPosition;

void main()
{
    const vec4 cicleColor = vec4(1.0, 0.0, 0.0, 1.0);
    const vec4 transparent = vec4(0.0, 0.0, 0.0, 0.0);

    float dist = length(fragPosition - center);
    float thickness = 0.03;
    float outerRadius = fragRadius;
    float innerRadius = outerRadius - thickness;
    vec4 finalCircleColor = mix(cicleColor, transparent, smoothstep(innerRadius, outerRadius, dist));
    color = finalCircleColor;
}