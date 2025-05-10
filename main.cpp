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
#include "sharedVariables.h"
#include "Input.h"
// #include <SOIL2/SOIL2.h>


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


    Welol::Camera* camera = new Welol::Camera(glm::vec3(0.0f, 0.0f, 30.0f), glm::vec3(0.0f, 0.0f, 0.0f));
    bool rotateCamera = false;

    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    //Text text{glRenderer};

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
            camera->update(cameraScrollDirection, yaw, pitch, mouseInfo.positionX, mouseInfo.positionY, forwardCamera);
        } else{
            camera->setLastMousePos(mouseInfo.positionX, mouseInfo.positionY);
        }
        if (forwardCamera)
        camera->moveForward(cameraScrollDirection);
        forwardCamera = false;

        simulation.update(glRenderer, camera->getViewMatrix(), (float) deltaTime, mouseInfo);

        //text.render("text", 20.0f, 20.0f, 2.0f, glm::vec3(1.0f, 0.0f, 0.0f), glRenderer);


        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    delete camera;

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
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