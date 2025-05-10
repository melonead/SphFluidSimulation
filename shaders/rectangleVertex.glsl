#version 430 core
layout (location=0) in vec2 position;

uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;

uniform float cx;
uniform float cy;

uniform float width;
uniform float height;


out vec2 fragCenter;
out vec2 fragPosition;
out float fragWidth;
out float fragHeight;

void main()
{
    
    vec4 worldPosition = model * vec4(position, 0.0, 1.0);
    gl_Position = projection * view  * worldPosition;
    fragCenter.x = cx;
    fragCenter.y = cy;
    fragPosition = worldPosition.xy;
    fragWidth = width;
    fragHeight = height;
}