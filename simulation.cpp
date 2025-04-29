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
    
    //particlesInfo.neighbors = std::vector<std::vector<unsigned int>>(numParticles);
    particlesInfo.positions.resize(numParticles);
    particlesInfo.densities.resize(numParticles);
    particlesInfo.pressures.resize(numParticles);
    particlesInfo.forces = std::vector<std::array<float, 2>>(numParticles);
    // particlesInfo.accelerations = std::vector<std::array<float, 2>>(numParticles);
    particlesInfo.velocities = std::vector<std::array<float, 2>>(numParticles);

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


    // 
    std::cout << table.getCellSize() << std::endl;
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
    // computePressures();
    computeForces();
    //LeapFrogIntegration(dt);
    // table.clearTable();

    glm::vec2 particlePos = particlesInfo.positions[55];
    std::vector<unsigned int> neighbors = table.getNeighborIDs(particlePos);


    updateRendering(glRenderer, view, perspectiveMatrix);
    vis.drawGrid(glRenderer, view, perspectiveMatrix);

    vis.drawCircle(particlePos.x, particlePos.y, radiusOfInfluence, view, perspectiveMatrix, glRenderer);

    for (unsigned int i = 0; i < neighbors.size(); i++)
    {
        glm::vec2 p = particlesInfo.positions[neighbors[i]];
        float dist = getDistance(p, particlePos);
        //if (dist < radiusOfInfluence)
        vis.drawCircle(p.x, p.y, 0.05f, view, perspectiveMatrix, glRenderer);
    }

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

        std::vector<unsigned int> neighborIDs = table.getNeighborIDs(particlesInfo.positions[i]);
    
        unsigned int size = neighborIDs.size();
        for (unsigned int id: neighborIDs)
        {
            if (id == i)
                continue;
            dist = getDistance(particlesInfo.positions[i], particlesInfo.positions[id]);

            if (dist <= radiusOfInfluence) {
                   
                density += mass * poly6(dist) * poly6Const;
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
        std::vector<unsigned int> neighborIDs = table.getNeighborIDs(particlesInfo.positions[i]);
        float forceFieldX = 0.0f;
        float forceFieldY = 0.0f;
    
        unsigned int size = neighborIDs.size();

        // Compute the force density field value at this location
        for (unsigned int id: neighborIDs)
        {
            if (id == i)
                continue;
            dist = getDistance(particlesInfo.positions[i], particlesInfo.positions[id]);

            if (dist <= radiusOfInfluence) 
            {
                // Compute pressure field
                float pressNeig = getPressure(particlesInfo.densities[id]);
                float PressCurr = getPressure(particlesInfo.densities[i]);
                float densNeig = particlesInfo.densities[id];
                pressureDensityFieldValue += computePressureSingleParticle(PressCurr, pressNeig, densNeig, dist);
                //std::cout << "v: " << pressureDensityFieldValue << std::endl;
                glm::vec2 vec = particlesInfo.positions[id] - particlesInfo.positions[i];
                forceFieldX += (vec.x/dist) * pressureDensityFieldValue;
                forceFieldY += (vec.y/dist) * pressureDensityFieldValue;
            }
        }

        float pressureMultiplier = 1.0f;
        // Compute acceleration
        float accX = (-forceFieldX * pressureMultiplier) / particlesInfo.densities[i];
        float accY = (-forceFieldY * pressureMultiplier) / particlesInfo.densities[i];
        // Compute velocity 
        // REVISIT: temporary 
        float deltaTime = (1.0f/60.0f);
        particlesInfo.velocities[i][0] += accX * deltaTime;
        particlesInfo.velocities[i][1] += accY * deltaTime;
        //particlesInfo.velocities[i][1] += (-9.8f * deltaTime) ;






        // --------------------------------------------------------------------------------------------------------------
        float ax = 0.0f;
        float ay = 0.0f;
        float boundaryX = cellSize * table.getHorizontalCellsCount() * 0.5f;
        float boundaryY = cellSize * table.getVerticalCellCount() * 0.5f;
    
        if (particlesInfo.velocities[i][0] > maxSpeed) {
            particlesInfo.velocities[i][0] = std::min(particlesInfo.velocities[i][0], maxSpeed);
        }
        else if (particlesInfo.velocities[i][0] < -maxSpeed) {
            particlesInfo.velocities[i][0] = std::max(particlesInfo.velocities[i][0], -maxSpeed);
        }
        particlesInfo.positions[i].x += (particlesInfo.velocities[i][0] * deltaTime + (ax * pow(deltaTime, 2)) / 2.0);
        //p.predictedPosition[0] = particlesInfo.positions[i].x + particlesInfo.velocities[i][0] * deltaTime;

        float eps = 0.01f;
        if ((particlesInfo.positions[i].x + eps) < -boundaryX)
        {

            particlesInfo.velocities[i][0] *= damp;
            particlesInfo.positions[i].x = -boundaryX + eps;

        }

        if ((particlesInfo.positions[i].x + eps) > boundaryX) {
            particlesInfo.velocities[i][0] *= damp;
            particlesInfo.positions[i].x = boundaryX - eps;
        }

        particlesInfo.velocities[i][1] += gravity * deltaTime;
        if (particlesInfo.velocities[i][1] > maxSpeed) {
            particlesInfo.velocities[i][1] = std::min(particlesInfo.velocities[i][1], maxSpeed);
        }
        else if (particlesInfo.velocities[i][1] < -maxSpeed) {
            particlesInfo.velocities[i][1] = std::max(particlesInfo.velocities[i][1], -maxSpeed);
        }
        particlesInfo.positions[i].y += (particlesInfo.velocities[i][1] * deltaTime + (ay * pow(deltaTime, 2.0)) / 2.0);
        //p.predictedPosition[1] = particlesInfo.positions[i].y + particlesInfo.velocities[i][1] * deltaTime;

        if ((particlesInfo.positions[i].y - eps) < -boundaryY)
        {
            particlesInfo.velocities[i][1] *= damp;
            particlesInfo.positions[i].y = -boundaryY + eps;

        }
        if ((particlesInfo.positions[i].y + eps) > boundaryY) {
            particlesInfo.velocities[i][1] *= damp;
            particlesInfo.positions[i].y = boundaryY - eps;
        }
        // --------------------------------------------------------------------------------------------------------------
    }
}

// 
float Simulation::computePressureSingleParticle(
    float currentPressure,
    float neighborPressure, 
    float neighborDensity, 
    float dist
)
{
    float symm = (currentPressure + neighborPressure) / (2.0f);
    float grad = (15.0f / (PI * powf(radiusOfInfluence, 6))) * powf((radiusOfInfluence - dist), 3);

    return mass * symm * grad;
}


/*
// computePressures: compute the pressure values of the particles
 void Simulation::computePressures() {

     for (unsigned int i = 0; i < numParticles; i++) {
        std::vector<unsigned int> neighborCellsKeys = table.computeNeighboringCellsKeys(particlesInfo.positions[i]);
        unsigned int size = neighborCellsKeys.size();

        float dist = 0.0f;
        float pressure = 0.0f;
        for (unsigned int key : neighborCellsKeys)
        {

            for (unsigned int neighborIndex : table.getCell(key)) {

                if (neighborIndex == i)
                    continue;

                dist = getDistance(particlesInfo.positions[i], particlesInfo.positions[neighborIndex]);

                if (dist < radiusOfInfluence) {
                    dist = getDistance(particlesInfo.positions[i], particlesInfo.positions[neighborIndex]);
                    pressure += computePressureSingleParticle(
                        particlesInfo.pressures[i], 
                        particlesInfo.pressures[neighborIndex],
                        particlesInfo.densities[neighborIndex],
                        dist);
                }

            }

            particlesInfo.pressures[i] = pressure * pressureMultiplier;
        }
        
    }

}

// computeFroces: calculate the total force acting on a particles
// this force will be used to compute acceleration -> velocity -> position
void Simulation::computeForces() {
    for (unsigned int i = 0; i < numParticles; i++) {
        std::array<float, 2> directionVector = { 0.0f, 0.0f };
        std::array<float, 2> pressureForce = { 0.0f, 0.0f };

        std::vector<unsigned int> neighborCellsKeys = table.computeNeighboringCellsKeys(particlesInfo.positions[i]);
        unsigned int size = neighborCellsKeys.size();

        float dist = 0.0f;
        float pressure = 0.0f;
        for (unsigned int key : neighborCellsKeys)
        {

            for (unsigned int neighborIndex : table.getCell(key)) {

                if (neighborIndex == i)
                    continue;

                dist = getDistance(particlesInfo.positions[i], particlesInfo.positions[neighborIndex]);

                if (dist < radiusOfInfluence) {
                    dist = getDistance(particlesInfo.positions[i], particlesInfo.positions[neighborIndex]);

                    glm::vec2& pos1 = particlesInfo.positions[i];
                    glm::vec2& pos2 = particlesInfo.positions[neighborIndex];
                    
                    directionVector[0] = (pos2.x - pos1.x) / dist;
                    directionVector[1] = (pos2.y - pos1.y) / dist;

                    pressureForce[0] += particlesInfo.pressures[i] * directionVector[0];
                    pressureForce[1] += particlesInfo.pressures[i] * directionVector[1];
                }

            }

            particlesInfo.forces[i] = pressureForce;


        }
    }
}


// simulateParticle: complete run on a single particle's simulation
void simulateParticle(Particle& p, std::vector<Particle> &nbs, int N)
{
    float pressureForce[2] = {0.0, 0.0};
    float viscForce[2] = {0.0, 0.0};
    float dirVec[2] = { 0.0, 0.0 };


    computeDensities();
    computePressures();
    computeForces();
    for (int i = 0; i < N; i++)
    {
        Particle& nb = nbs.at(i);
        if (nb.id == p.id)
            continue;
        float dist = getDistance(p.pos, nb.pos);
        dirVec[0] = (nb.pos[0] - p.pos[0]) / dist;
        dirVec[1] = (nb.pos[1] - p.pos[1]) / dist;
        // compute pressure
        float pInf = computePressure(p, nb, dist);
        pressureForce[0] += pInf * dirVec[0];
        pressureForce[1] += pInf * dirVec[1];
        // compute viscosity
        float vInf = computeViscosity(p, nb, dist);
        viscForce[0] += (nb.velocity[0] - p.velocity[0]) * vInf * dirVec[0];
        viscForce[1] += (nb.velocity[1] - p.velocity[1]) * vInf * dirVec[1];
    }

    p.totalForce[0] = pressureForce[0] + viscForce[0] * VISCMULTIPLIER;
    p.totalForce[1] = pressureForce[1] + viscForce[1] * VISCMULTIPLIER;

}

*/

// LeapFrogIntegration: integrate particle's position and velocity
void Simulation::LeapFrogIntegration(float deltaTime)
{

    // std::cout << "position.x: " << particlesInfo.positions[0].x << std::endl;
    float ax = 0.0f;
    float ay = 0.0f;
    float boundaryX = cellSize * table.getHorizontalCellsCount() * 0.5f;
    float boundaryY = cellSize * table.getVerticalCellCount() * 0.5f;
    for(unsigned int i = 0; i < numParticles; i++ ) {
        if (particlesInfo.velocities[i][0] > maxSpeed) {
            particlesInfo.velocities[i][0] = std::min(particlesInfo.velocities[i][0], maxSpeed);
        }
        else if (particlesInfo.velocities[i][0] < -maxSpeed) {
            particlesInfo.velocities[i][0] = std::max(particlesInfo.velocities[i][0], -maxSpeed);
        }
        particlesInfo.positions[i].x += (particlesInfo.velocities[i][0] * deltaTime + (ax * pow(deltaTime, 2)) / 2.0);
        //p.predictedPosition[0] = particlesInfo.positions[i].x + particlesInfo.velocities[i][0] * deltaTime;

        float eps = 0.01f;
        if ((particlesInfo.positions[i].x + eps) < -boundaryX)
        {

            particlesInfo.velocities[i][0] *= damp;
            particlesInfo.positions[i].x = -boundaryX + eps;

        }

        if ((particlesInfo.positions[i].x + eps) > boundaryX) {
            particlesInfo.velocities[i][0] *= damp;
            particlesInfo.positions[i].x = boundaryX - eps;
        }

        particlesInfo.velocities[i][1] += gravity * deltaTime;
        if (particlesInfo.velocities[i][1] > maxSpeed) {
            particlesInfo.velocities[i][1] = std::min(particlesInfo.velocities[i][1], maxSpeed);
        }
        else if (particlesInfo.velocities[i][1] < -maxSpeed) {
            particlesInfo.velocities[i][1] = std::max(particlesInfo.velocities[i][1], -maxSpeed);
        }
        particlesInfo.positions[i].y += (particlesInfo.velocities[i][1] * deltaTime + (ay * pow(deltaTime, 2.0)) / 2.0);
        //p.predictedPosition[1] = particlesInfo.positions[i].y + particlesInfo.velocities[i][1] * deltaTime;

        if ((particlesInfo.positions[i].y - eps) < -boundaryY)
        {
            particlesInfo.velocities[i][1] *= damp;
            particlesInfo.positions[i].y = -boundaryY + eps;

        }
        if ((particlesInfo.positions[i].y + eps) > boundaryY) {
            particlesInfo.velocities[i][1] *= damp;
            particlesInfo.positions[i].y = boundaryY - eps;
        }

    }
    
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

// poly6a: poly6 kernel
float Simulation::poly6(float dist)
{	
    return powf((radiusOfInfluence * radiusOfInfluence - dist * dist), 3);
}

// _max: return max of two value
float Simulation::_min(float a, float b)
{
    return a > b ? a : b;
}

// getPressure: return the pressure to be applied to a particle
float Simulation::getPressure(float d)
{	
    float v = pressureMultiplier * (idealDensity - d);
    return std::max(v, 0.0f);
}


// getPressureGradient: the pressure gradient kernel
float Simulation::getPressureGradient(float dist)
{
    float s = 15.0f / (pi * pow(radiusOfInfluence, 6.0f));
    return s * (pow((radiusOfInfluence - dist), 3.0));
}

// getViscGradient: the viscosity gradient kernel
float Simulation::getViscGradient(float dist)
{
    return spikyLap * (radiusOfInfluence - dist);
}

// computeViscosity: calculate viscous force between two particles
// float Simulation::computeViscosity(Particle& p, Particle& nb, float dist)
// {
//     return (mass / nb.density) * getViscGradient(dist);
// }

/*
void Simulation::interactWithMouse(
    Particle&p,
    glm::vec4 mousePos,
    bool mouseClick[],
    std::vector<Particle>(&Table)[NUMCELLS]
)
{
    // get distance from mouse position
    float xDist = mousePos.x - p.pos[0];
    float yDist = mousePos.y - p.pos[1];
    float speed = 0.85;
    float infDist = 2.0;
    float distSqr = xDist * xDist + yDist * yDist;

    if (distSqr < (infDist * infDist) )
    {
        float dist = sqrt(distSqr);
        float dir[2] = { xDist / dist, yDist / dist };
        // repel
        if (dist < infDist && mouseClick[0])
        {
            p.velocity[0] += (-dir[0] * speed);
            p.velocity[1] += (-dir[1] * speed);
        }

        // attract
        if (dist < infDist && mouseClick[1])
        {
            p.velocity[0] += (dir[0] * speed);
            p.velocity[1] += (dir[1] * speed);
        }
    }
}
*/