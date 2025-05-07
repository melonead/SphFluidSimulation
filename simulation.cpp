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
    particlesInfo.predictedPositions.resize(numParticles);
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
    
    particlesInfo.predictedPositions = particlesInfo.positions;
    table.createTable(particlesInfo.positions);
    setUpRendering(glRenderer);

    float containerPosX = -containerWidth * 0.5f;
    float containerPosY = containerHeight * 0.5f;
    vis = Visualization{table.getCellSize(), containerWidth, containerHeight, containerPosX, containerPosY, glRenderer};
    vis.setPerspectiveMatrix(perspectiveMatrix);
}

Simulation::~Simulation() {
} 

void Simulation::update(Welol::Renderer& glRenderer, glm::mat4& view, float dt, MouseInfo& mouseInfo) {
    vis.setViewMatrix(view);

    _setMousePosition(mouseInfo.positionX, mouseInfo.positionY, view);
    computeDensities();
    // if (mouseInfo.leftButton || mouseInfo.rightButton)
    //     mouseInteraction(mouseInfo);

    computeForces(mouseInfo); 

    //vis.drawGrid(glRenderer, view, perspectiveMatrix);
    updateRendering(glRenderer, view, perspectiveMatrix);

    vis.drawCircle(particlesInfo.positions[0].x, particlesInfo.positions[0].y, radiusOfInfluence);
    vis.drawCircle(0.0f, 0.0f, radiusOfInfluence * 10.0f);
    vis.drawRectangle(-containerWidth / 2.0f, containerHeight / 2.0f, containerWidth, containerHeight);
    

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

        NeighborQuery& nQuery = table.getNeighborIDs(particlesInfo.predictedPositions[i]);
        
        unsigned int neighborID;
        // std::cout << "bucket size: " << nQuery.size << std::endl;
        for (unsigned int j = 0; j < nQuery.size; j++)
        {
            neighborID = nQuery.neighborBucket[j];
            if (neighborID == i)
                continue;
            dist = getDistance(particlesInfo.predictedPositions[i], particlesInfo.predictedPositions[neighborID]);

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


void Simulation::computeForces(MouseInfo& mouseInfo)
{
    for (int i = 0; i < numParticles; i++)
    {
        float dist = 0.0f;

        float totalForceFieldX = 0.0f;
        float totalForceFieldY = 0.0f;

        float pressureForceFieldX = 0.0f;
        float pressureForceFieldY = 0.0f;
        float viscosityForceFieldX = 0.0f;
        float viscosityForceFieldY = 0.0f;
        /*
            REVISIT: Is it faster to just use the sortedParticlesIDs directly rather than packing them in
            a new container?
        */
        // Compute the force density field value at this location
        NeighborQuery& nQuery = table.getNeighborIDs(particlesInfo.predictedPositions[i]);
        
        unsigned int neighborID;

        for (unsigned int j = 0; j < nQuery.size; j++)
        {
            neighborID = nQuery.neighborBucket[j];
            if (neighborID == i)
                continue;

            dist = getDistance(particlesInfo.predictedPositions[i], particlesInfo.predictedPositions[neighborID]);
            /*
                Distance can be zero. Is there a way for C++ say we are dividing by zero?
            */
           
           if (dist <= radiusOfInfluence) 
           {
               glm::vec2 dir = (particlesInfo.predictedPositions[neighborID] - particlesInfo.predictedPositions[i]) / dist;
                if (dist <= 0.0f)
                {
                    dir.x = (float)std::rand() / (float)std::rand();
                    dir.y = (float)std::rand() / (float)std::rand();
                }

                // Compute pressure field
                float pressNeig = getPressure(particlesInfo.densities[neighborID]);
                float PressCurr = getPressure(particlesInfo.densities[i]);
                float densNeig = particlesInfo.densities[neighborID];
                float densCurr = particlesInfo.densities[i];
                
                float pField = computePressureSingleParticle(PressCurr, pressNeig, densCurr, densNeig, dist);
                pressureForceFieldX += (dir.x) * (pField);
                pressureForceFieldY += (dir.y) * (pField);

                glm::vec2 visField = computeViscosity(particlesInfo.velocities[i], particlesInfo.velocities[neighborID], particlesInfo.densities[neighborID], dist);
                viscosityForceFieldX += (dir.x) * visField.x;
                viscosityForceFieldY += (dir.y) * visField.y;

                // viscosityForceFieldX += (dir.x) * (viscosityDensityFieldAccum);
                // viscosityForceFieldY += (dir.y) * (viscosityDensityFieldAccum);

                // Near Force.
                // Prevent clustering by applying a (near) force inverse propertional to distance between particles.
                // The smaller the distance the higher the force. 
                particlesInfo.velocities[i][0] += (-dir.x) * (1.0f - dist/radiusOfInfluence) * nearForceConstant;
                particlesInfo.velocities[i][1] += (-dir.y) * (1.0f - dist/radiusOfInfluence) * nearForceConstant;
                // move the neighbor in the opposite direction.
                particlesInfo.velocities[neighborID][0] += (dir.x) * (1.0f - dist/radiusOfInfluence) * nearForceConstant;
                particlesInfo.velocities[neighborID][1] += (dir.y) * (1.0f - dist/radiusOfInfluence) * nearForceConstant;

            }
        }

        // Accumulate the forces

        totalForceFieldX += -pressureForceFieldX;
        totalForceFieldY += -pressureForceFieldY;

        // totalForceFieldX += viscosityForceFieldX * viscosityConstant;
        // totalForceFieldY += viscosityForceFieldY * viscosityConstant;

        // Compute acceleration
        float accX = (totalForceFieldX)  / particlesInfo.densities[i];
        float accY = (totalForceFieldY)  / particlesInfo.densities[i];
        // Compute velocity 
        // REVISIT: temporary 
        particlesInfo.velocities[i][0] += accX * deltaTime;
        particlesInfo.velocities[i][1] += accY * deltaTime;
        particlesInfo.velocities[i][1] += (-9.8f * deltaTime);

        if (mouseInfo.leftButton || mouseInfo.rightButton)
            mouseInteraction(mouseInfo, i);
        
        updatePosition(i);
        
    }
}

void Simulation::updatePosition(unsigned int i)
{

    float boundaryX = containerWidth * 0.5f;
    float boundaryY = containerHeight * 0.5f;
    
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
    particlesInfo.predictedPositions[i].x = particlesInfo.positions[i].x;
    
    particlesInfo.positions[i].y += (particlesInfo.velocities[i][1] * deltaTime);
    
    
    if ((particlesInfo.positions[i].y - eps) < -boundaryY)
    {
        particlesInfo.velocities[i][1] *= damp;
        particlesInfo.positions[i].y += eps;

    }

    if ((particlesInfo.positions[i].y + eps) > boundaryY) {
        particlesInfo.velocities[i][1] *= damp;
        particlesInfo.positions[i].y -= eps;
    }
        
    particlesInfo.predictedPositions[i].y = particlesInfo.positions[i].y;
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
    float symm = (currentPressure + neighborPressure) / (2.0f * neighborDensity);
    return mass * symm * cubicSpikyKernel(dist);
}

glm::vec2 Simulation::computeViscosity(std::array<float, 2>& velCurr, std::array<float, 2>& velNeig, float densNeig, float dist)
{

    glm::vec2 res;
    res.x = (velNeig[0] - velCurr[0]) / densNeig;
    res.y = (velNeig[1] - velCurr[1]) / densNeig;

    res = mass * res * viscosityLaplacian(dist);
    return res;
}





void Simulation::mouseInteraction(MouseInfo& mouseInfo, unsigned int ID)
{
    /*
    glm::vec2 pos = mouseWorldPosition;
    float radius = radiusOfInfluence * 30.0f;
    NeighborQuery& nQuery = table.getNeighborIDsForMouse(pos, radius);
 
    glm::vec2 dir;
    float dist;
    float strength = 50.5f;

    if (mouseInfo.rightButton)
        strength *= -1.0f;

    for (unsigned int i = 0; i < nQuery.size; i++)
    {   
        unsigned int nID= nQuery.neighborBucket[i];

        dist = getDistance(particlesInfo.positions[nID], pos);

        if (dist > radius)
            continue;

        dir = (particlesInfo.positions[nID] - pos) / dist;

        //std::cout << "s: " << dir.x * strength << "s1: " << dir.y * strength << std::endl;
        particlesInfo.velocities[nID][0] += 1.0f;
        particlesInfo.velocities[nID][1] += 1.0f;
    }
    */

    float radius = radiusOfInfluence * 5.0f;
    float dist = getDistance(particlesInfo.positions[ID], mouseWorldPosition);
    glm::vec2 dir;
    float strength = 2.0f;

    if (mouseInfo.rightButton)
        strength *= -1.0f;

    if (dist < radius)
    {
        dir = (particlesInfo.positions[ID] - mouseWorldPosition) / dist;
        particlesInfo.velocities[ID][0] += dir.x * strength;
        particlesInfo.velocities[ID][1] += dir.y * strength;
    }
}

glm::vec2 Simulation::screenToNdc(float x, float y)
{
    glm::vec2 result;
    result.x = ((2.0f * x) / SCR_WIDTH -  1.0f);
    result.y = (1.0f - (2.0f * y) / SCR_HEIGHT);

    return result;
}

glm::vec2 Simulation::ndcToWorld(float x, float y, glm::mat4& view)
{   
    
    glm::vec4 interm(x, y, -1.0f, 1.0f);

    interm = glm::inverse(perspectiveMatrix * view) * interm;

    
    glm::vec2 result;
    result.x  = interm.x * 30.0f;
    result.y = interm.y * 30.0f;
    
    return result;
}

void Simulation::_setMousePosition(float x, float y, glm::mat4& view) 
{
    glm::vec2 mouseNdc = screenToNdc(x, y);    
    mouseWorldPosition = ndcToWorld(mouseNdc.x, mouseNdc.y, view);
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

float Simulation::viscosityLaplacian(float dist)
{
    float a = 45.0f / (PI * powf(radiusOfInfluence, 6.0f));
    float b = (radiusOfInfluence - dist);
    return a * b;
}

float Simulation::quadraticSpkikyKernel(float dist)
{
    return powf(1.0f - (dist / radiusOfInfluence), 2.0f);
}

float Simulation::cubicSpikyKernel(float dist)
{
    return powf(1.0f - (dist / radiusOfInfluence), 3.0f);
}