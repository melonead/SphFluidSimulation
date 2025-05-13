#version 430 core

out vec4 color;

in float fragRadius;

in vec2 center;
in vec2 fragPosition;
in float fragXCoord;

void main()
{


    const vec4 blue = vec4(0.0, 0.0, 1.0, 1.0);
    const vec4 green = vec4(0.0, 1.0, 0.0, 1.0);
    const vec4 yellow = vec4(1.0, 1.0, 0.0, 1.0);
    const vec4 red = vec4(1.0, 0.0, 0.0, 1.0);
    vec4 cicleColor = vec4(1.0, 0.0, 0.0, 1.0);

    if (fragXCoord >= 0.0 && fragXCoord <= 0.25)
    {
        // blue green interpolation
        cicleColor = mix(blue, green, smoothstep(0.0, 0.25, fragXCoord));
    }
    else if (fragXCoord >= 0.25 && fragXCoord <= 0.75)
    {
        // green yellow interpolation
        cicleColor = mix(green, yellow, smoothstep(0.25, 0.75, fragXCoord));
    }
    else if (fragXCoord >= 0.75 && fragXCoord <= 1.0)
    {
        // yellow red interpolation
        cicleColor = mix(yellow, red, smoothstep(0.75, 1.0, fragXCoord));
    }

    const vec4 transparent = vec4(0.0, 0.0, 0.0, 0.0);

    float dist = length(fragPosition - center);
    float thickness = 0.03;
    float outerRadius = fragRadius;
    float innerRadius = outerRadius - thickness;

    vec4 finalCircleColor = mix(cicleColor, transparent, smoothstep(innerRadius, outerRadius, dist));
    color = finalCircleColor;
}