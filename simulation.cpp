#include <simulation.h>
#include <array>
#include <glm/gtc/matrix_transform.hpp>
#include "Renderer.h"


// taking the particles temporarily so that i don't break the old code completely
// as I reimplement the rendering system
Simulation::Simulation(Shader& circleShaderProgram, Welol::Renderer& glRenderer) 
    : shaderProg{circleShaderProgram}
{

    int size = std::sqrt(numParticles);
    float spacing = 0.2f;


    float startPos[2] = { -size * 0.5f * spacing + 0.05f, size * 0.5f * spacing - 0.05f };
   

    assert(numParticles > 0);
    
    particlesInfo.positions.resize(numParticles);
    particlesInfo.densities.resize(numParticles);
    particlesInfo.velocities = std::vector<std::array<float, 2>>(numParticles);

    for (unsigned int i = 0; i < particlesInfo.velocities.size(); i++)
    {
        particlesInfo.velocities[i][0] = 0.0f;
        particlesInfo.velocities[i][1] = 0.0f;
    }

    for (int y = 0; y < size; y++)
    {
        for (int x = 0; x < size; x++)
        {   
            glm::vec2 pos;
            int index = x + y * size;
            //SimParticles[index].mass = 1;
            pos.x = startPos[0] + (x * spacing);
            pos.y = startPos[1] - (y * spacing);
            particlesInfo.positions[index] = pos;

        }
    }
    
    table.createTable(particlesInfo.positions);
    setUpRendering(glRenderer);

    float containerWidth = table.getHorizontalCellsCount() * table.getCellSize();
    float containerHeight = table.getVerticalCellCount() * table.getCellSize();

    float containerPosX = -(float)table.getHorizontalCellsCount() * table.getCellSize() * 0.5f;
    float containerPosY = (float)table.getVerticalCellCount() * table.getCellSize() * 0.5f;
    vis = Visualization{table.getCellSize(), containerWidth, containerHeight, containerPosX, containerPosY, glRenderer};
}

Simulation::~Simulation() {
} 

void Simulation::update(Welol::Renderer& glRenderer, glm::mat4& view, float dt) {

    computeDensities();
    computeForces();


    
    
    updateRendering(glRenderer, view, perspectiveMatrix);
    vis.drawGrid(glRenderer, view, perspectiveMatrix);
    
    
    /*
    glm::vec2 particlePos = particlesInfo.positions[55];
    vis.drawCircle(particlePos.x, particlePos.y, radiusOfInfluence, view, perspectiveMatrix, glRenderer);
    NeighborQuery nQuery = table.getNeighborIDs(particlePos);
    for (unsigned int i = 0; i < nQuery.size; i++)
    {
        glm::vec2 p = particlesInfo.positions[nQuery.neighborBucket[i]];
        float dist = getDistance(p, particlePos);
        //if (dist < radiusOfInfluence)
        vis.drawCircle(p.x, p.y, 0.05f, view, perspectiveMatrix, glRenderer);
    }
    */

    table.createTable(particlesInfo.positions);
}


void Simulation::setUpRendering(Welol::Renderer& glRenderer) {

    int size = sizeof(glm::vec2);

    std::vector<float> quadVertices = {
        -0.05f, 0.05f, 
         0.05f,-0.05f, 
        -0.05f,-0.05f, 
        -0.05f, 0.05f, 
         0.05f,-0.05f, 
         0.05f, 0.05f
    };

    unsigned int quadVerticesSize = quadVertices.size();
    Welol::VertexAttribute positionAttr{0, Welol::WL_FLOAT2, particlesInfo.positions.data(), numParticles, true};
    Welol::VertexAttribute quadPosAttr{1,  Welol::WL_FLOAT2, quadVertices.data(), quadVerticesSize, false};

    particleRenderOperation.addVertexAttribute(positionAttr);
    particleRenderOperation.addVertexAttribute(quadPosAttr);

    glRenderer.initializeRenderOperation(particleRenderOperation);

}

void Simulation::updateRendering(Welol::Renderer& glRenderer, glm::mat4& view, glm::mat4& projection) {
   
    shaderProg.use();

    glRenderer.render(particleRenderOperation);

    shaderProg.setMatrix4fv("view", view);
    shaderProg.setMatrix4fv("projection", projection);
    glRenderer.updateRenderOperationVertexAttribute(particleRenderOperation, 0, 0, particlesInfo.positions.data());
}


// computeDensities: compute densities for all the particles
void Simulation::computeDensities()
{
    
    for (int i = 0; i < numParticles; i++)
    {
        float dist = 0.0f;
        //int nbSize = particlesInfo.neighbors[i].size();
        float density = 1.0f;

        NeighborQuery nQuery = table.getNeighborIDs(particlesInfo.positions[i]);
        
        unsigned int neighborID;
        // std::cout << "bucket size: " << nQuery.size << std::endl;
        for (unsigned int j = 0; j < nQuery.size; j++)
        {
            neighborID = nQuery.neighborBucket[j];
            if (neighborID == i)
                continue;
            dist = getDistance(particlesInfo.positions[i], particlesInfo.positions[neighborID]);

            /*
                Distance can be zero.
            */
            if (dist <= 0.0f)
                continue;

            if (dist <= radiusOfInfluence) {
                   
                density += mass * cubicSplineKernel(dist);
            }
        }
            
        particlesInfo.densities[i] = density;

    }

}


void Simulation::computeForces()
{
    for (int i = 0; i < numParticles; i++)
    {
        float dist = 0.0f;
        float pressureDensityFieldValue = 0.0f;
        float forceFieldX = 0.0f;
        float forceFieldY = 0.0f;
    
        // Compute the force density field value at this location
        NeighborQuery nQuery = table.getNeighborIDs(particlesInfo.positions[i]);
        // std::cout << "bucket size: " << nQuery.size << std::endl;
        unsigned int neighborID;
        for (unsigned int j = 0; j < nQuery.size; j++)
        {
            neighborID = nQuery.neighborBucket[j];
            if (neighborID == i)
                continue;

            dist = getDistance(particlesInfo.positions[i], particlesInfo.positions[neighborID]);
            /*
                Distance can be zero. Is there a way for C++ say we are dividing by zero?
            */
           
           if (dist <= radiusOfInfluence) 
           {
               glm::vec2 vec = particlesInfo.positions[neighborID] - particlesInfo.positions[i];
                if (dist <= 0.0f)
                {
                    vec.x = (float)std::rand() / (float)std::rand();
                    vec.y = (float)std::rand() / (float)std::rand();
                }
                // Compute pressure field
                float pressNeig = getPressure(particlesInfo.densities[neighborID]);
                float PressCurr = getPressure(particlesInfo.densities[i]);
                float densNeig = particlesInfo.densities[neighborID];
                float densCurr = particlesInfo.densities[i];
                pressureDensityFieldValue += computePressureSingleParticle(PressCurr, pressNeig, densCurr, densNeig, dist);
                //std::cout << "v: " << pressureDensityFieldValue << std::endl;

                forceFieldX += (vec.x/dist) * pressureDensityFieldValue;
                forceFieldY += (vec.y/dist) * pressureDensityFieldValue;
            }
        }

        // Compute acceleration
      
        float accX = -forceFieldX / particlesInfo.densities[i];
        float accY = -forceFieldY / particlesInfo.densities[i];
        // Compute velocity 
        // REVISIT: temporary 
        particlesInfo.velocities[i][0] += accX * deltaTime;
        particlesInfo.velocities[i][1] += accY * deltaTime;
        //particlesInfo.velocities[i][1] += (-9.8f * deltaTime);

        // std::cout << "velx: " << accX << std::endl;
        // std::cout << "vely: " << accY << std::endl;

        // --------------------------------------------------------------------------------------------------------------

        float boundaryX = cellSize * table.getHorizontalCellsCount() * 0.5f;
        float boundaryY = cellSize * table.getVerticalCellCount() * 0.5f;
        
        // clamp velocity
        if (particlesInfo.velocities[i][0] > maxSpeed) {
            particlesInfo.velocities[i][0] = maxSpeed;
        }
        else if (particlesInfo.velocities[i][0] < -maxSpeed) {
            particlesInfo.velocities[i][0] = -maxSpeed;
        }

        if (particlesInfo.velocities[i][1] > maxSpeed) {
            particlesInfo.velocities[i][1] = maxSpeed;
        }
        else if (particlesInfo.velocities[i][1] < -maxSpeed) {
            particlesInfo.velocities[i][1] = -maxSpeed;
        }
        



        particlesInfo.positions[i].x += (particlesInfo.velocities[i][0] * deltaTime);
        //p.predictedPosition[0] = particlesInfo.positions[i].x + particlesInfo.velocities[i][0] * deltaTime;

        float eps = 0.05f;
        if ((particlesInfo.positions[i].x - eps) < -boundaryX)
        {

            particlesInfo.velocities[i][0] *= damp;
            particlesInfo.positions[i].x  += eps;

        }
        else if ((particlesInfo.positions[i].x + eps) > boundaryX) {
            particlesInfo.velocities[i][0] *= damp;
            particlesInfo.positions[i].x -= eps;
        }

        
        particlesInfo.positions[i].y += (particlesInfo.velocities[i][1] * deltaTime);
        //p.predictedPosition[1] = particlesInfo.positions[i].y + particlesInfo.velocities[i][1] * deltaTime;

        if ((particlesInfo.positions[i].y - eps) < -boundaryY)
        {
            particlesInfo.velocities[i][1] *= damp;
            particlesInfo.positions[i].y += eps;

        }
        if ((particlesInfo.positions[i].y + eps) > boundaryY) {
            particlesInfo.velocities[i][1] *= damp;
            particlesInfo.positions[i].y -= eps;
        }
        // --------------------------------------------------------------------------------------------------------------
    }
}

// 
float Simulation::computePressureSingleParticle(
    float currentPressure,
    float neighborPressure, 
    float currentDensity,
    float neighborDensity,
    float dist
)
{
    float symm = (currentPressure + neighborPressure) / 2.0f;
    return mass * symm * cubicSplineKernel(dist);
}

// convert world coordinate to the screen space coordinate
glm::vec2 Simulation::worldToScreen(glm::vec2& position, glm::mat4& view) {
    glm::vec4 clip = perspectiveMatrix * view * glm::vec4(position.x, position.y, 0.0f, 0.0f);
    float clipX = clip.x;
    float clipY = clip.y;
    float screenX = ((1.0f + clipX) * SCR_WIDTH) / 2.0f;
    float screenY = ((1.0f - clipY) * SCR_HEIGHT) / 2.0f;

    return glm::vec2(screenX, screenY);
}

glm::vec4 Simulation::ndcToWorld(glm::vec4& vec, glm::mat4& view, glm::mat4& projection) {
    return glm::inverse(projection * view)  * vec;
}


// getDistance: get distance between two locations
float Simulation::getDistance(glm::vec2& pos1, glm::vec2& pos2)
{	
    float x_dist = (pos2.x - pos1.x);
    float y_dist = (pos2.y - pos1.y);

    return sqrtf(x_dist * x_dist + y_dist * y_dist);
}

// _max: return max of two value
float Simulation::_min(float a, float b)
{
    return a > b ? a : b;
}

// getPressure: return the pressure to be applied to a particle
float Simulation::getPressure(float d)
{	
    float v = pressureConstant * (d - idealDensity);
    return v;
}

float Simulation::poly6Kernel(float dist)
{
    float a = 315.0f / (64.0f * PI * powf(radiusOfInfluence, 9.0f));
    float b = powf(powf(radiusOfInfluence, 2) - powf(dist, 2), 3);
    return a * b;
}

float Simulation::cubicSplineKernel(float dist)
{
    float q = (1.0f / radiusOfInfluence) * dist;
    float a = (40.0f / (7.0f * PI * powf(radiusOfInfluence, 2.0f)));

    if (q <= 0.5f)
        return a * (6.0f * (powf(q, 3.0f) - powf(q, 2.0f)) + 1.0f);
    else 
        return a * (2.0f * powf((1.0f - q), 3.0f));
}

float Simulation::spikyKernel(float dist)
{
    float a = 15.0f / (PI * powf(radiusOfInfluence, 6.0f));
    float b = powf((radiusOfInfluence - dist), 3.0f);
    return a * b;
}