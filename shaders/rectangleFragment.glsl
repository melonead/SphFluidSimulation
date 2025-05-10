#version 430 core

out vec4 color;
in vec2 fragCenter;
in vec2 fragPosition;
in float fragWidth;
in float fragHeight;

void main()
{   
    const vec4 cicleColor = vec4(1.0, 0.0, 0.0, 1.0);
    const vec4 transparent = vec4(2.0, 0.0, 0.0, 0.0);
    float thickness = 0.875;
    float epsilon = 0.001;

    vec2 vecToInner;

    vec2 centerToPosition = fragPosition - fragCenter;
    vec2 dir = normalize(centerToPosition);
    float dist = length(centerToPosition);

    float innerRadius;
    float outerRadius;


    float H = fragHeight / 2.0;
    float X;

    if (centerToPosition.y > epsilon)
        X = H * (centerToPosition.x / centerToPosition.y);
    else
        X = fragWidth / 2.0;

    vecToInner.x = fragCenter.x + dir.x * X;
    vecToInner.y = fragCenter.y + dir.y * H;

    outerRadius = length(vecToInner);
    innerRadius = outerRadius - thickness;



    

    vec4 thicknessColor = mix(transparent, cicleColor, smoothstep(innerRadius, outerRadius, dist));
    color = mix(thicknessColor, transparent, smoothstep(innerRadius, outerRadius, dist));
}