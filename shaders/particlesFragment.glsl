#version 430 core

out vec4 color;

in float fragRadius;

in vec2 center;
in vec2 fragPosition;
in float fragXCoord;

uniform sampler2D gradientTexture;

void main()
{
    const vec4 transparent = vec4(0.0, 0.0, 0.0, 0.0);

    float dist = length(fragPosition - center);
    float thickness = 0.03;
    float outerRadius = fragRadius;
    float innerRadius = outerRadius - thickness;

    vec2 v = vec2(fragXCoord, 0.5);
    vec4 cicleColor = texture(gradientTexture, v);
    vec4 finalCircleColor = mix(cicleColor, transparent, smoothstep(innerRadius, outerRadius, dist));
    color = finalCircleColor;
    // color = vec4(1.0 - fragXCoord, 0.0, 0.0, 1.0);
}