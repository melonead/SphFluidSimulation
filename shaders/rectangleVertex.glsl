#version 430 core
layout (location=0) in vec2 position;

uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;

void main()
{
    
    vec4 worldPosition = model * vec4(position, 0.0, 1.0);
    gl_Position = projection * view  * worldPosition;
}