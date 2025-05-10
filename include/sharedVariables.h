#pragma once


const float radiusOfInfluence = 0.50f;
const float cellSize = radiusOfInfluence;
const float containerWidth = 60.0f;
const float containerHeight = 40.0f;
const float PI = 3.141592653589793f;
const unsigned int numParticles {70 * 70};
const float deltaTime = 1.0f / 60.0f;
const float pressureConstant = -4.3f;
const float nearForceConstant = (-1.0f * pressureConstant);
const float viscosityConstant = 10.0f;


static const int SCR_WIDTH =  1200;
static const int SCR_HEIGHT = 600;