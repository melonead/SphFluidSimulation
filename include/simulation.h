#pragma once
#include "Hashtable.h"
#include <glad/glad.h>
#include <vector>
#include "shader.h"
#include "Hashtable.h"
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "Renderer.h"
#include "Camera.h"
#include "visualization.h"
#include "sharedVariables.h"

/*
    Simulation brings manages the whole simulation
*/
class Simulation {
public:
    Simulation(Shader& circleShaderProgram, Welol::Renderer& glRenderer);
    ~Simulation();
    void update(Welol::Renderer& glRenderer, glm::mat4& view, float dt);
private:
    /*
        Information used by OpenGl to render the particles.
    */
    struct DisplayInformation {
        glm::mat4 projection;
        unsigned int instancingVbo;
        unsigned int quadVao;
        unsigned int quadVbo;
    } displayInfo;


    /*
        ParticleInformation contains information about the particles: positions, velocities
        Velocities are going to be sent to the gpu for rendering color based on the speed of
        the particle.
        Positions are used for rendering too.
    */
    struct ParticlesInformation {
        std::vector<glm::vec2> positions;
        std::vector<float> densities;
        std::vector<float> pressures;
        std::vector<std::array<float, 2>> forces;
        //std::vector<std::array<float, 2>> accelerations;
        std::vector<std::array<float, 2>> velocities;
        //std::vector<std::vector<unsigned int>> neighbors;

    } particlesInfo;

    // The number of particles to be simulated
    unsigned int numParticles {50 * 50};
    float gravity {0.0f};
    float SCR_WIDTH{ 600.0f };
    float SCR_HEIGHT{ 600.0f };

    float aspect = SCR_WIDTH / SCR_HEIGHT;

    //glm::mat4 orthoMatrix = glm::ortho(0.0f, 1200.0f, 600.0f, 0.0f, -1.0f, 1.0f);
    glm::mat4 perspectiveMatrix = glm::perspective(glm::radians(45.0f), aspect, 0.1f, 100.0f);

    float pi = 3.145;
    const float poly6Const = 315.0f / (64.0f * pi * pow(radiusOfInfluence, 9.0f));
    const float spikyGrad = -45.0f / (pi * pow(radiusOfInfluence, 6.0f));
    const float spikyLap = 45.0f / (pi * pow(radiusOfInfluence, 6.0f));

    float maxSpeed{ 3.0f };
    float damp{ -0.8f };
    float idealDensity{ 1000.0f };
    float pressureMultiplier{ 1.0f };
    float mass{ 1.0f };


    unsigned int numberOfVertices = 6;
    Welol::RenderOperation particleRenderOperation{Welol::WL_TRIANGLES, numberOfVertices, 0, numParticles, true, false};

    // container boundaries
    glm::vec2 horizontalBoundaries = glm::vec2(-1.25f, 1.25f);
    glm::vec2 verticalBoundaries = glm::vec2(-1.25f, 1.25f);

    /*
        The hash table enables fast search for particles withing the radius of influence
        of a particles.
    */
    HashTable table{ numParticles};


    /*
        setUpRendering prepares for the rendering of the particles.
    */

    Shader shaderProg;
    void setUpRendering(Welol::Renderer& glRenderer);
    void updateRendering(Welol::Renderer& glRenderer, glm::mat4& view, glm::mat4& projection);
    glm::vec2 worldToScreen(glm::vec2& position, glm::mat4& view);
    glm::vec4 ndcToWorld(glm::vec4& vec, glm::mat4& view, glm::mat4& projection);


    /*
        Visualization of the hash table
    */

    Visualization vis;


    /*
        Smoothed Particle Hydrodynamics related functions
    */

    // getDistance: get distance between two locations
    float getDistance(glm::vec2& pos1, glm::vec2& pos2);

    // poly6a: poly6 kernel
    float poly6(float dist);

    // _max: return max of two value
    float _min(float a, float b);

    // getPressure: return the pressure to be applied to a particle
    float getPressure(float d);

    // avPressure: return the average pressure between 
    // two particles
    float avPressure(float a, float b);

    // getPressureGradient: the pressure gradient kernel
    float getPressureGradient(float dist);

    // getViscGradient: the viscosity gradient kernel
    float getViscGradient(float dist);

    // computeDensities: compute densities for all the particles
    void computeDensities();

    // computeViscosity: calculate viscous force between two particles
    /*
    float computeViscosity(Particle& p, Particle& nb, float dist);

    void interactWithMouse(
        Particle&p,
        glm::vec4 mousePos,
        bool mouseClick[],
        std::vector<Particle>(&Table)[NUMCELLS]
    );
    */
    // computePressures: compute the pressure values of the particles
    void computePressures();

    // compute pressure of a two neighboring particles
    float computePressureSingleParticle(
        float currentPressure,
        float neighborPressure,
        float neighborDensity,
        float dist
    );

    // computeFroces: calculate the total force acting on a particles
    // this force will be used to compute acceleration -> velocity -> position
    void computeForces();

    // LeapFrogIntegration: integrate particle's position and velocity
    void LeapFrogIntegration(float deltaTime);



};