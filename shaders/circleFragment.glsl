#version 430 core

out vec4 color;

in float fragRadius;
in vec2 fragPosition;
in vec2 fragCenter;

float thickness;

void main()
{   
    const vec4 cicleColor = vec4(1.0, 0.0, 0.0, 1.0);
    const vec4 transparent = vec4(1.0, 0.0, 0.0, 0.0);
    thickness = 0.075;
    float dist = length(fragPosition - fragCenter);

    float outerRadius = fragRadius;
    float innerRadius = outerRadius - thickness;
    
    vec4 thicknessColor = mix(transparent, cicleColor, smoothstep(innerRadius, outerRadius, dist));
    color = mix(thicknessColor, transparent, smoothstep(innerRadius, outerRadius, dist));

}