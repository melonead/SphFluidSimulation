#pragma once


const float radiusOfInfluence = 0.50f;
const float cellSize = radiusOfInfluence;
const float containerWidth = 40.0f;
const float containerHeight = 20.0f;
const unsigned int numParticles {4032};
const float PI = 3.141592653589793f;
const float deltaTime = 1.0f / 60.0f;
const float pressureConstant = -3.3f;
const float nearForceConstant = (-1.0f * pressureConstant);
const float viscosityConstant = -10.0f;


static const int SCR_WIDTH =  1200;
static const int SCR_HEIGHT = 600;