#version 430 core

layout (location=0) in vec3 verticalLine;
layout (location=1) in vec3 horizontalLine;

uniform mat4 projection;
uniform mat4 view;


uniform float lineSpacing;
uniform int numberOfLines;

void main()
{	
	int halfOfLinesCount = numberOfLines/2;
	if (gl_InstanceID < halfOfLinesCount)
	{
        // vertical lines
        // Why is my camera left-handed? left axis seems to increase towards the left.
		gl_Position = projection * view * vec4(verticalLine.x - gl_InstanceID * lineSpacing, verticalLine.y, verticalLine.z, 1.0);
	}
	else {
        // horizontal lines
		gl_Position = projection * view * vec4(horizontalLine.x , horizontalLine.y - (gl_InstanceID-halfOfLinesCount) * lineSpacing, horizontalLine.z, 1.0);
	}
	
}