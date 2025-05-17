#include <simulation.h>
#include <array>
#include <glm/gtc/matrix_transform.hpp>
#include "Renderer.h"
#include "Camera.h"
#include "Texture.h"
#include "omp.h"

// taking the particles temporarily so that i don't break the old code completely
// as I reimplement the rendering system
Simulation::Simulation(Shader& circleShaderProgram, Welol::Renderer& glRenderer) 
    : shaderProg{circleShaderProgram}
{

    
    float containerPosX = -settings->containerWidth * 0.5f;
    float containerPosY = settings->containerHeight * 0.5f;
    vis = Visualization{table.getCellSize(), settings->containerWidth, settings->containerHeight, containerPosX, containerPosY, glRenderer};
    vis.setRenderer(glRenderer);
    vis.setPerspectiveMatrix(perspectiveMatrix);

    setParticles(settings->numParticles, glRenderer);
    
    setUpRendering(glRenderer);
}

Simulation::~Simulation() {
}

void Simulation::setParticles(unsigned int numParticles, Welol::Renderer& glRenderer)
{
    //assert(settings->numParticles > 0);

    
    // Add particles if any have been added.
    if (settings->numParticles > settings->prevNumParticles)
    {
        unsigned int particleDifference = settings->numParticles - settings->prevNumParticles;
        int size = std::sqrt(settings->numParticles);
        int height = (int) (settings->numParticles / size);
        float spacing = 0.4f;
        
        float startPos[2] = { -size * 0.5f * spacing + 0.05f, size * 0.5f * spacing - 0.05f };
        
        // REVISIT(PERFORMANCE): Do we have to resize the whole vector? could we
        // just add the additional size instead of resizing the whole thing? Perhaps
        // it's not bad because this only happens when use add particles.
        particlesInfo.positions.resize(settings->numParticles);
        particlesInfo.predictedPositions.resize(settings->numParticles);
        particlesInfo.densities.resize(settings->numParticles);
        particlesInfo.gradientTextureCoordinates.resize(settings->numParticles);
        particlesInfo.velocities = std::vector<std::array<float, 2>>(settings->numParticles);
        
        for (unsigned int i = 0; i < settings->numParticles; i++)
        {
            particlesInfo.velocities[i][0] = 0.0f;
            particlesInfo.velocities[i][1] = 0.0f;
        
        }
        
        
        int extraparticleCount = settings->numParticles - height * size;
        int index = 0;
        for (int y = 0; y < height; y++)
        {
            for (int x = 0; x < size; x++)
            {   
                index = x + y * size;
                glm::vec2 pos;
                //SimParticles[index].settings->mass = 1;
                pos.x = startPos[0] + (x * spacing);
                pos.y = startPos[1] - (y * spacing);
                particlesInfo.positions[index] = pos;
            }
        
        
        }
        
        // Append remaining particles.
        for (int ex = 0; ex < extraparticleCount; ex++)
        {
            glm::vec2 pos;
            index++;
            //SimParticles[index].settings->mass = 1;
            pos.x = startPos[0] + (ex * spacing);
            pos.y = startPos[1] - ((size + 1) * spacing);
            particlesInfo.positions[index] = pos;
        }
        
        
        particlesInfo.predictedPositions = particlesInfo.positions;
        
        // 
        unsigned int offset = sizeof(glm::vec2) * settings->prevNumParticles;
        unsigned int numberOfBytes = sizeof(glm::vec2) * particleDifference;
    
        unsigned int offsetIntoBucket = sizeof(glm::vec2) * (settings->prevNumParticles - 1);

        glRenderer.addSizeBytesToBuffer(particleRenderOperation, 0, offset, numberOfBytes, particlesInfo.positions.data());
        // Generate indices before creating the table.
        table.generateParticlesIDs(settings->numParticles);
        table.createTable(particlesInfo.positions);

        particleRenderOperation.setNumberOfInstancesToRender(particlesInfo.positions.size());
    } 
    else // Remove particles if any have been removed.
    {
        unsigned int diff = settings->prevNumParticles - settings->numParticles;

        if (diff <= 0)
            return;
        // Remove diff particles from positions list
        for (unsigned int i = 0; i < diff; i++)
        {
            particlesInfo.positions.pop_back();
        }

        particleRenderOperation.setNumberOfInstancesToRender(particlesInfo.positions.size());
    }
}

void Simulation::update(Welol::Renderer& glRenderer, Welol::Camera& camera, float dt, MouseInfo& mouseInfo) {

    vis.setViewMatrix(camera.getViewMatrix());

    if (settings->updateParticlesCount)
    {
        setParticles(settings->numParticles, glRenderer);
    }

    if (settings->startSimulation)
    {
        _setMousePosition(mouseInfo.positionX, mouseInfo.positionY, camera);
        computeDensities();
        computeForces(mouseInfo);
    }
    
    
    updateRendering(glRenderer, camera.getViewMatrix(), perspectiveMatrix);
    vis.drawGrid(glRenderer);

    vis.drawRectangle(-settings->containerWidth * 0.5f, settings->containerHeight * 0.5f, settings->containerWidth, settings->containerHeight);
    

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
    Welol::VertexAttribute positionAttr{0, Welol::WL_FLOAT2, particlesInfo.positions.data(), (unsigned int) settings->maxParticles, true};
    Welol::VertexAttribute quadPosAttr{1,  Welol::WL_FLOAT2, quadVertices.data(), quadVerticesSize, false};
    Welol::VertexAttribute gradientCoordinates{2, Welol::WL_FLOAT, particlesInfo.gradientTextureCoordinates.data(), (unsigned int) settings->maxParticles, true};

    particleRenderOperation.addVertexAttribute(positionAttr);
    particleRenderOperation.addVertexAttribute(quadPosAttr);       
    particleRenderOperation.addVertexAttribute(gradientCoordinates);       

    glRenderer.initializeRenderOperation(particleRenderOperation);

    gradientTexture.attachImageData(cameraTexturePath);
}

void Simulation::updateRendering(Welol::Renderer& glRenderer, glm::mat4& view, glm::mat4& projection) {
   
    shaderProg.use();

    
    shaderProg.setMatrix4fv("view", view);
    shaderProg.setMatrix4fv("projection", projection);
    shaderProg.setMatrix4fv("model", particleModelMatrix);
    shaderProg.setFloat("radius", particleRadius);
    gradientTexture.update(shaderProg);
    glRenderer.render(particleRenderOperation);



    //glRenderer.updateRenderOperationVertexAttribute(particleRenderOperation, 0, 0, particlesInfo.positions.data());
    unsigned int numberOfBytes = sizeof(glm::vec2) * particlesInfo.positions.size();
    glRenderer.addSizeBytesToBuffer(particleRenderOperation, 0, 0, numberOfBytes, particlesInfo.positions.data());
    unsigned int gradientBytes = sizeof(float) * particlesInfo.gradientTextureCoordinates.size();
    glRenderer.addSizeBytesToBuffer(particleRenderOperation, 2, 0, gradientBytes, particlesInfo.gradientTextureCoordinates.data());
    
}


// computeDensities: compute densities for all the particles
void Simulation::computeDensities()
{
    #pragma omp parallel for
    for (int i = 0; i < settings->numParticles; i++)
    {
        // compute predicted positions.
        particlesInfo.predictedPositions[i].x = particlesInfo.positions[i].x + particlesInfo.velocities[i][0] * settings->deltaTime;
        particlesInfo.predictedPositions[i].y = particlesInfo.positions[i].y + particlesInfo.velocities[i][1] * settings->deltaTime;

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
            dist = getDistance(particlesInfo.predictedPositions[i], particlesInfo.predictedPositions[neighborID]);

            /*
                Distance can be zero.
            */
            if (dist <= 0.0f)
                continue;

            if (dist <= settings->radiusOfInfluence) {
                   
                density += settings->mass * poly6Kernel(dist);
            }
        }
        
        particlesInfo.densities[i] = density;

    }

}


void Simulation::computeForces(MouseInfo& mouseInfo)
{
    #pragma omp parallel for
    for (int i = 0; i < settings->numParticles; i++)
    {
        float dist = 0.0f;

        float totalForceFieldX = 0.0f;
        float totalForceFieldY = 0.0f;

        float pressureForceFieldX = 0.0f;
        float pressureForceFieldY = 0.0f;
        float viscosityForceFieldX = 0.0f;
        float viscosityForceFieldY = 0.0f;
        float surfaceTensionX = 0.0f;
        float surfaceTensionY = 0.0f;

        float nFx = 0.0f;
        float nFy = 0.0f;

        /*
            REVISIT: Is it faster to just use the sortedParticlesIDs directly rather than packing them in
            a new container?
        */
        // Compute the force density field value at this location
        NeighborQuery nQuery = table.getNeighborIDs(particlesInfo.predictedPositions[i]);  
        
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
           
           if (dist <= settings->radiusOfInfluence) 
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
                pressureForceFieldX += dir.x * pField;
                pressureForceFieldY += dir.y * pField;

                glm::vec2 visField = computeViscosity(particlesInfo.velocities[i], particlesInfo.velocities[neighborID], particlesInfo.densities[neighborID], dist);
                viscosityForceFieldX += visField.x;
                viscosityForceFieldY += visField.y;

                // Near Force.
                // Prevent clustering by applying a (near) force inverse propertional to distance between particles.
                // The smaller the distance the higher the force. 
                // move the neighbor in the opposite direction.
                nFx += (-dir.x) * (1.0f - dist/settings->radiusOfInfluence) * settings->nearForceConstant;
                nFy += (-dir.y) * (1.0f - dist/settings->radiusOfInfluence) * settings->nearForceConstant;



                // Surface Tension
                float sTGrad = surfaceTensionField(dist, densNeig);
                float colField = colorField(dist, densNeig);

                if (abs(colField) < 0.005f)
                    continue;

                float lap = surfaceTensionLaplacian(dist, densNeig);

                surfaceTensionX = (1.0f - sTGrad) * dir.x;
                surfaceTensionY = (1.0f - sTGrad) * dir.y;

                // surfaceTensionX = surfaceTensionField(dist, densNeig);
            }
        }

        // Accumulate the forces

        totalForceFieldX += pressureForceFieldX * -1.0f;
        totalForceFieldY += pressureForceFieldY * -1.0f;

        totalForceFieldX += nFx;
        totalForceFieldY += nFy;

        totalForceFieldX += viscosityForceFieldX * settings->viscosityConstant;
        totalForceFieldY += viscosityForceFieldY * settings->viscosityConstant;

        // totalForceFieldX += surfaceTensionX;
        // totalForceFieldY += surfaceTensionY;

        // Compute acceleration
        float accX = (totalForceFieldX)  / particlesInfo.densities[i];
        float accY = (totalForceFieldY)  / particlesInfo.densities[i];
        // Compute velocity 
        // REVISIT: temporary 
        float halfDeltaVelX = particlesInfo.velocities[i][0] - 0.5f * settings->deltaTime * accX;
        float halfDeltaVelY = particlesInfo.velocities[i][1] - 0.5f * settings->deltaTime * accY;
        particlesInfo.velocities[i][0] = halfDeltaVelX + accX * settings->deltaTime;
        particlesInfo.velocities[i][1] = halfDeltaVelY + accY * settings->deltaTime;
        particlesInfo.velocities[i][1] += (settings->gravity * settings->deltaTime);

        // compute the magnitude of the velocity and normalize it.
        // This value will be used to determine the color of the particle in the color gradient.
        float gradientImageWidth = 1920.0f;
        float mag = (particlesInfo.velocities[i][0] * particlesInfo.velocities[i][0] + particlesInfo.velocities[i][1] * particlesInfo.velocities[i][1]) / (settings->maxSpeed * settings->maxSpeed);
        //particlesInfo.gradientTextureCoordinates[i] = mag;
        particlesInfo.gradientTextureCoordinates[i] = mag;
        

        if (mouseInfo.leftButton || mouseInfo.rightButton)
            mouseInteraction(mouseInfo, i);
        
        updatePosition(i);
        
    }
}

void Simulation::updatePosition(unsigned int i)
{

    float boundaryX = settings->containerWidth * 0.5f;
    float boundaryY = settings->containerHeight * 0.5f;
    
    // clamp velocity
    if (particlesInfo.velocities[i][0] > settings->maxSpeed) {
        particlesInfo.velocities[i][0] = settings->maxSpeed;
    }
    else if (particlesInfo.velocities[i][0] < -settings->maxSpeed) {
        particlesInfo.velocities[i][0] = -settings->maxSpeed;
    }

    if (particlesInfo.velocities[i][1] > settings->maxSpeed) {
        particlesInfo.velocities[i][1] = settings->maxSpeed;
    }
    else if (particlesInfo.velocities[i][1] < -settings->maxSpeed) {
        particlesInfo.velocities[i][1] = -settings->maxSpeed;
    }

    particlesInfo.positions[i].x += (particlesInfo.velocities[i][0] * settings->deltaTime);
    

    float eps = 0.1f;
    if ((particlesInfo.positions[i].x - eps) < -boundaryX)
    {

        particlesInfo.velocities[i][0] *= damp;
        particlesInfo.positions[i].x  += eps;

    }
    else if ((particlesInfo.positions[i].x + eps) > boundaryX) {
        particlesInfo.velocities[i][0] *= damp;
        particlesInfo.positions[i].x -= eps;
    }
    
    particlesInfo.positions[i].y += (particlesInfo.velocities[i][1] * settings->deltaTime);
    
    
    if ((particlesInfo.positions[i].y - eps) < -boundaryY)
    {
        particlesInfo.velocities[i][1] *= damp;
        particlesInfo.positions[i].y += eps;

    }

    if ((particlesInfo.positions[i].y + eps) > boundaryY) {
        particlesInfo.velocities[i][1] *= damp;
        particlesInfo.positions[i].y -= eps;
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
    float symm = (currentPressure + neighborPressure) / (2.0f * neighborDensity);
    return settings->mass * symm * spikyKernel(dist);
}

glm::vec2 Simulation::computeViscosity(std::array<float, 2>& velCurr, std::array<float, 2>& velNeig, float densNeig, float dist)
{

    glm::vec2 res;
    res.x = (velNeig[0] - velCurr[0]) / densNeig;
    res.y = (velNeig[1] - velCurr[1]) / densNeig;

    res = settings->mass * res * viscosityLaplacian(dist);
    return res;
}





void Simulation::mouseInteraction(MouseInfo& mouseInfo, unsigned int ID)
{

    float radius = settings->radiusOfInfluence * 10.0f;
    float dist = getDistance(particlesInfo.positions[ID], mouseWorldPosition);
    glm::vec2 dir;
    float strength = settings->mouseStrength;

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

glm::vec2 Simulation::ndcToWorld(float x, float y, Welol::Camera& camera)
{   
    
    glm::vec4 interm(x, y, -1.0f, 1.0f);

    interm = glm::inverse(perspectiveMatrix) * interm;

    // eye space ray
    interm = glm::vec4(interm.x, interm.y, -1.0f, 0.0f);

    // world space ray
    interm = glm::inverse(camera.getViewMatrix()) * interm;

    glm::vec3 unitRay = glm::normalize(glm::vec3(interm.x, interm.y, interm.z));
    
    glm::vec3 planeNormal(0.0f, 0.0f, 1.0f);
    float planeDistanceFromOrigin = 0.0f;
    float distance = -(glm::dot(camera.getPosition(), planeNormal) + planeDistanceFromOrigin) / glm::dot(unitRay, planeNormal);
    
    glm::vec2 result(unitRay.x * distance, unitRay.y * distance);

    return result;
}

void Simulation::_setMousePosition(float x, float y, Welol::Camera& camera) 
{
    glm::vec2 mouseNdc = screenToNdc(x, y);    
    mouseWorldPosition = ndcToWorld(mouseNdc.x, mouseNdc.y, camera);

}

// getDistance: get distance between two locations
float Simulation::getDistance(glm::vec2& pos1, glm::vec2& pos2)
{	
    return glm::length(pos2 - pos1);
}

// _max: return max of two value
float Simulation::_min(float a, float b)
{
    return a > b ? a : b;
}

// getPressure: return the pressure to be applied to a particle
float Simulation::getPressure(float d)
{	
    float v = (settings->pressureConstant) * (settings->idealDensity / 7.0f) * (powf((d/settings->idealDensity), 7.0f) - 1.0f);
    return v;
}

float Simulation::poly6Kernel(float dist)
{
    float a = 315.0f / (64.0f * settings->PI * powf(settings->radiusOfInfluence, 9.0f));
    float b = powf(powf(settings->radiusOfInfluence, 2) - powf(dist, 2), 3);
    return a * b;
}

float Simulation::poly6Gradient(float dist)
{
    float res = (29.53125 / (settings->PI * powf(settings->radiusOfInfluence, 9.0f)));
    res *= powf(powf(settings->radiusOfInfluence, 2) - powf(dist, 2), 2) * dist;
    return res;
}

float Simulation::spikyKernel(float dist)
{
    float a = -45.0f / (settings->PI * powf(settings->radiusOfInfluence, 6.0f));
    float b = powf((settings->radiusOfInfluence - dist), 2.0f);
    return a * b;
}

float Simulation::viscosityLaplacian(float dist)
{
    float a = 45.0f / (settings->PI * powf(settings->radiusOfInfluence, 6.0f));
    float b = (settings->radiusOfInfluence - dist);
    return a * b;
}

float Simulation::colorField(float dist, float density)
{
    return settings->mass * (1.0f / density) * poly6Kernel(dist);
}

float Simulation::surfaceTensionField(float dist, float density)
{
    // I omitted the mass which give much sharper results.
    return (1.0f / density) * poly6Gradient(dist);
}

float Simulation::surfaceTensionLaplacian(float dist, float density)
{
    float res = powf(settings->radiusOfInfluence, 2) - powf(dist, 2);
    res *= (powf(settings->radiusOfInfluence, 2) - 5.0f * powf(dist, 2));
    res *= (29.53125 / (settings->PI * powf(settings->radiusOfInfluence, 9.0f)));
    return res;
}

// REVISIT: compute surface tension here.
float Simulation::computeSurfaceTension(float dist, float density)
{
    float gradientField = surfaceTensionField(dist, density);
    float laplacian = surfaceTensionLaplacian(dist, density);
    return 0;
}