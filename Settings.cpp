#include "Settings.h"
#include "imgui.h"
#include <iostream>

SettingsSingleton* SettingsSingleton::instance()
{
    if (_instance == 0)
    {
        _instance = new SettingsSingleton;
    }
    return _instance;
}

SettingsSingleton::SettingsSingleton()
{

}

void SettingsSingleton::displayUi()
{
    ImGui::SeparatorText("Simulation state");
    {
        ImGui::DragFloat("Pressure Constant", &pressureConstant, 0.05f);
        ImGui::SliderFloat("Radius of Influence", &radiusOfInfluence, 0.0f, 2.5f);
        ImGui::DragFloat("viscosity Multiplier", &viscosityConstant, 0.01f);
        ImGui::SliderFloat("Near Force Multiplier", &nearForceConstant, -1000.0f, 2000.0f);
        ImGui::DragFloat("Particle Mass", &mass, 0.01f);
        ImGui::DragFloat("Max Speed", &maxSpeed, 0.1f);

        prevNumParticles = numParticles;
   
        if (ImGui::SliderInt("Particles", &numParticles, 1, maxParticles))
        {
            updateParticlesCount = true;
        }
        else
        {
            updateParticlesCount = false;
        }

        ImGui::SliderFloat("gravity", &gravity, -9.803f, 9.803f);
        ImGui::SliderFloat("Ideal Density", &idealDensity, -1000.0f, 1000.0f);
        ImGui::DragFloat("Mouse Strength", &mouseStrength, 0.01f);        

        // ImGui::PushID(0);
        // ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor::HSV(0.01f, 0.6f, 0.6f));
        if (ImGui::Button("Start Simulation"))
        {
            startSimulation = true;
        }
        // ImGui::PopStyleColor(3);
        // ImGui::PopID();
    }
}

SettingsSingleton* SettingsSingleton::_instance = 0;