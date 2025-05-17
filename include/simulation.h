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
#include "Settings.h"
#include "Input.h"
#include "Texture.h"

/*
    Simulation brings manages the whole simulation
*/
class Simulation {
public:
    Simulation(Shader& circleShaderProgram, Welol::Renderer& glRenderer);
    ~Simulation();
    void update(Welol::Renderer& glRenderer, Welol::Camera& camera, float dt, MouseInfo& mouseInfo);
private:

    SettingsSingleton* settings = SettingsSingleton::instance();
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
        std::vector<glm::vec2> predictedPositions;
        std::vector<float> gradientTextureCoordinates;
        std::vector<float> densities;
        std::vector<std::array<float, 2>> velocities;

    } particlesInfo;

    // The number of particles to be simulated
    
    float aspect = SCR_WIDTH / SCR_HEIGHT;

    //glm::mat4 orthoMatrix = glm::ortho(0.0f, 1200.0f, 600.0f, 0.0f, -1.0f, 1.0f);
    glm::mat4 perspectiveMatrix = glm::perspective(glm::radians(45.0f), aspect, 0.1f, 100.0f);

    float particleRadius{0.05f * 4.0f};
    glm::mat4x4 particleModelMatrix = glm::scale(glm::mat4(1.0f), glm::vec3(particleRadius/0.05f, particleRadius/0.05f, 1.0f));

    float damp{ -0.8f };

    unsigned int numberOfVertices = 6;

    Welol::RenderOperation particleRenderOperation{Welol::WL_TRIANGLES, numberOfVertices, 0, settings->maxParticles, true, false};

    /*
        The hash table enables fast search for particles withing the radius of influence
        of a particles.
    */
    HashTable table{};


    glm::vec2 mouseWorldPosition;

    //Welol::Texture gradientTexture;

    // Gradient Texture
    //std::string gTexPath = "C:\\Users\\brian\\programming_projects\\WelolRenderer\\WelolRenderer\\FluidSim\\particleGradient.png";
    // std::string gTexPath = "C:\\Users\\brian\\programming_projects\\WelolRenderer\\WelolRenderer\resource\\skybox\\cubemap\\cubemap_negy.png";
    std::string cameraTexturePath = "C:\\Users\\brian\\programming_projects\\WelolRenderer\\WelolRenderer\\FluidSim\\particleGradient.png";
    std::string cameraTextureName = "gradientTexture";
    unsigned int texUnit = 0;
    Welol::Texture gradientTexture{Welol::WL_TEX_2D, Welol::WL_RGB, cameraTexturePath, 0, cameraTextureName, texUnit};

    /*
        setUpRendering prepares for the rendering of the particles.
    */
    void setUpRendering(Welol::Renderer& glRenderer);

    Shader shaderProg;
    void updateRendering(Welol::Renderer& glRenderer, glm::mat4& view, glm::mat4& projection);
    void setParticles(unsigned int numParticles, Welol::Renderer& glRenderer);
    /*
        Visualization of the hash table
    */

    Visualization vis;
    


    void updatePosition(unsigned int i);
    /*
        Smoothed Particle Hydrodynamics related functions
    */

    // getDistance: get distance between two locations
    float getDistance(glm::vec2& pos1, glm::vec2& pos2);

    // poly6a: poly6 kernel
    float poly6Kernel(float dist);
    float poly6Gradient(float dist);
    // Color field has nothing to do with colors of the particles
    // this is a name for this quantity in the literature.
    float colorField(float dist, float density);
    float surfaceTensionLaplacian(float dist, float density);
    float spikyKernel(float dist);
    float viscosityLaplacian(float dist);

    glm::vec2 getMouseWorldCoords(float mouseX, float mouseY);

    // _max: return max of two value
    float _min(float a, float b);

    // getPressure: return the pressure to be applied to a particle
    float getPressure(float d);

    // computeDensities: compute densities for all the particles
    void computeDensities();


    // compute pressure of a two neighboring particles
    float computePressureSingleParticle(
        float currentPressure,
        float neighborPressure,
        float currentDensity,
        float neighborDensity,
        float dist
    );

    float surfaceTensionField(float dist, float density);
    float computeSurfaceTension(float dist, float density);

    glm::vec2 computeViscosity(std::array<float, 2>& velCurr, std::array<float, 2>& velNeig, float densNeig, float dist);

    // computeFroces: calculate the total force acting on a particles
    // this force will be used to compute acceleration -> velocity -> position
    void computeForces(MouseInfo& mouseInfo);

    void _setMousePosition(float x, float y, Welol::Camera& camera);

    glm::vec2 screenToNdc(float X, float Y);

    glm::vec2 ndcToWorld(float X, float Y, Welol::Camera& camera);

    void mouseInteraction(MouseInfo& mouseInfo, unsigned int ID);
};