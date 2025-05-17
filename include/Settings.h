#pragma once

const int SCR_WIDTH =  1200;
const int SCR_HEIGHT = 600;

class SettingsSingleton
{

public:
    static SettingsSingleton* instance();
    const unsigned int maxParticles{10000};
    float test{0.0f};
    float radiusOfInfluence = 0.50f;
    float cellSize = radiusOfInfluence;
    float containerWidth = 60.0f;
    float containerHeight = 40.0f;
    float PI = 3.141592653589793f;
    int numParticles {0};
    int prevNumParticles{numParticles};
    float deltaTime = 1.0f / 60.0f;
    // float pressureConstant = -4.3f;
    float pressureConstant = 10.0f;
    float nearForceConstant{1000.0f};
    float viscosityConstant = 0.5f;
    float gravity{-9.803f};
    float mass{0.07f};
    float maxSpeed{20.0f};
    float idealDensity{ 50.0f };
    float mouseStrength{1.0f};
    bool startSimulation{false};
    bool updateParticlesCount{false};

    void displayUi();
protected:
    SettingsSingleton();
private:
    static SettingsSingleton* _instance;


};
