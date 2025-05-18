#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <stb_image.h>
#include <chrono>
#include <vector>
#include "Hashtable.h"
#include "shader.h"
#include "Text.h"
#include "simulation.h"
#include "Renderer.h"
#include "Camera.h"
#include "Settings.h"
#include "Input.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"



void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window, bool& rotateCamera);
void window_reshape_callback(GLFWwindow* window, int newWidth, int newHeight);
float randomFloat(float Min, float Max);

void mouse_callback(GLFWwindow* window,double xpos,double ypos);
void mouse_scroll_callback(GLFWwindow*window, double x, double y);

double lastX = 0.0;
double lastY = 0.0;
float pitch = 0.0f;
float yaw = 0.0f;
bool firstMouse = true;
bool forwardCamera = false;
float cameraScrollDirection = 1.0f;

MouseInfo mouseInfo;



int main()
{
    // glfw: initialize and configure
    // ------------------------------
    if (!glfwInit())
    {
        std::cout << "Failed to initialize GLFW" << std::endl;
        return -1;
    }
    
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 4);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // glfw window creation
    // --------------------
    GLFWwindow* window; 
    window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Fluid Simulation", NULL, NULL);
    //glfwSetWindowSizeLimits(window, 600, 300, 1200, 600);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, mouse_scroll_callback);
    glfwSetFramebufferSizeCallback(window, window_reshape_callback);


    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // 

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;         // Enable Docking
    io.ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;       // Enable Multi-Viewport / Platform Windows
    //io.ConfigViewportsNoAutoMerge = true;
    //io.ConfigViewportsNoTaskBarIcon = true;

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsLight();

    // When viewports are enabled we tweak WindowRounding/WindowBg so platform windows can look identical to regular ones.
    ImGuiStyle& style = ImGui::GetStyle();
    if (io.ConfigFlags & ImGuiConfigFlags_ViewportsEnable)
    {
        style.WindowRounding = 0.0f;
        style.Colors[ImGuiCol_WindowBg].w = 1.0f;
    }

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 430");

    // Our state
    bool show_demo_window = true;
    bool show_another_window = false;
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

    //glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    //glEnable(GL_DEPTH_TEST);

    /*
        shader path
    */

    const std::string vertex_shader_path = "C:\\Users\\brian\\programming_projects\\WelolRenderer\\WelolRenderer\\FluidSim\\shaders\\particlesVertex.glsl";
    const std::string fragment_shader_path = "C:\\Users\\brian\\programming_projects\\WelolRenderer\\WelolRenderer\\FluidSim\\shaders\\particlesFragment.glsl";

    Welol::Renderer glRenderer;


     // new width&height provided by the callback
    glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT); // sets screen region associated with framebuffer 

    // time
    double lastTime = glfwGetTime();

    // initialize

    glfwSetWindowSizeCallback(window, window_reshape_callback);

    // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    
    Shader particleShaderProgram{vertex_shader_path, fragment_shader_path};


     /*
        simulation stuff
    */

    Simulation simulation{particleShaderProgram, glRenderer};


    Welol::Camera camera{glm::vec3(0.0f, 0.0f, 30.0f), glm::vec3(0.0f, 0.0f, 0.0f)};
    bool rotateCamera = false;

    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    Text text{glRenderer};

    SettingsSingleton* settings = SettingsSingleton::instance();

    const GLubyte* version = glGetString(GL_VERSION);

    while (!glfwWindowShouldClose(window))
    {


        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        processInput(window, rotateCamera);

        //glEnable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        double mousePosX;
        double mousePosY;
        glfwGetCursorPos(window, &mousePosX, &mousePosY);
        
        //vMat = glm::translate(vMat, glm::vec3(0.0f, 0.0f, 0.5f));

        if (!rotateCamera)
        {
            yaw = 0.0f;
            pitch = 0.0f;
        }
        
        if (rotateCamera)
        {
            camera.update(cameraScrollDirection, yaw, pitch, mouseInfo.positionX, mouseInfo.positionY, forwardCamera);
        } else{
            camera.setLastMousePos(mouseInfo.positionX, mouseInfo.positionY);
        }
        if (forwardCamera)
        camera.moveForward(cameraScrollDirection);
        forwardCamera = false;

        simulation.update(glRenderer, camera, (float) settings->deltaTime, mouseInfo);

        //text.render("text", 20.0f, 20.0f, 2.0f, glm::vec3(1.0f, 0.0f, 0.0f), glRenderer);

        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        settings->displayUi();

        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);

        // Rendering
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        // Update and Render additional Platform Windows
        // (Platform functions may change the current OpenGL context, so we save/restore it to make it easier to paste this code elsewhere.
        //  For this specific demo app we could also call glfwMakeContextCurrent(window) directly)
        if (io.ConfigFlags & ImGuiConfigFlags_ViewportsEnable)
        {
            GLFWwindow* backup_current_context = glfwGetCurrentContext();
            ImGui::UpdatePlatformWindows();
            ImGui::RenderPlatformWindowsDefault();
            glfwMakeContextCurrent(backup_current_context);
        }

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    
    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}




float randomFloat(float Min, float Max)
{
    return ((float(rand()) / float(RAND_MAX)) * (Max - Min)) + Min;
}


void window_reshape_callback(GLFWwindow* window, int newWidth, int newHeight/*, float aspect, glm::mat4& pMat*/ ) {
    //aspect = (float)newWidth / (float)newHeight; // new width&height provided by the callback
    glViewport(0, 0, newWidth, newHeight); // sets screen region associated with framebuffer 
    //pMat = glm::perspective(glm::radians(45.0f), aspect, 0.1f, 100.0f);
}

void processInput(GLFWwindow* window, bool& rotateCamera)
{

    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_3) == GLFW_PRESS)
    {
        rotateCamera = true;
    }

    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_3) == GLFW_RELEASE)
    {
        rotateCamera = false;
    }

    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_1) == GLFW_PRESS)
    {
        mouseInfo.leftButton = true;
    }

    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_1) == GLFW_RELEASE)
    {
        mouseInfo.leftButton = false;
    }

    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_2) == GLFW_PRESS)
    {
        mouseInfo.rightButton = true;
    }

    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_2) == GLFW_RELEASE)
    {
        mouseInfo.rightButton = false;
    }

}


void mouse_callback(GLFWwindow* window,double xpos,double ypos)
 {
    if(firstMouse)
    {
        lastX=xpos;
        lastY=ypos;
        firstMouse=false;
    }

    yaw = xpos - lastX;
    pitch = lastY - ypos;

    if (abs(yaw) < 2.5f)
        yaw = 0.0f;
    if (abs(pitch) < 2.5f)
        pitch = 0.0f;

    mouseInfo.positionX = xpos;
    mouseInfo.positionY = ypos;

    float yawSensitivity = 0.005f;
    float pitchSensitivity = 0.005f;
    yaw   *= yawSensitivity;
    pitch *= pitchSensitivity;

 }

 void mouse_scroll_callback(GLFWwindow*window, double x, double y)
 {
    forwardCamera = true;
    cameraScrollDirection = y;
 }