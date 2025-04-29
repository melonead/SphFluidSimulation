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
// #include <SOIL2/SOIL2.h>


void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window, bool& rotateCamera);
void window_reshape_callback(GLFWwindow* window, int newWidth, int newHeight);
float randomFloat(float Min, float Max);

void mouse_callback(GLFWwindow* window,double xpos,double ypos);
void mouse_scroll_callback(GLFWwindow*window, double x, double y);

static int SCR_WIDTH =  800;
static int SCR_HEIGHT = 800;


double lastX = 0.0;
double lastY = 0.0;
float pitch = 0.0f;
float yaw = 0.0f;
bool firstMouse = true;
bool forwardCamera = false;
float cameraScrollDirection = 1.0f;


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
    double deltaTime = 0.0;

    // initialize

    glfwSetWindowSizeCallback(window, window_reshape_callback);

    // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    
    Shader particleShaderProgram{vertex_shader_path, fragment_shader_path};


     /*
        simulation stuff
    */

    Simulation simulation{particleShaderProgram, glRenderer};


    Welol::Camera* camera = new Welol::Camera(glm::vec3(0.0f, 0.0f, -3.0f), glm::vec3(0.0f, 0.0f, 0.0f));
    bool rotateCamera = false;

    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    while (!glfwWindowShouldClose(window))
    {
        
        deltaTime = glfwGetTime() - lastTime;
        lastTime = glfwGetTime();


        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        processInput(window, rotateCamera);

        //glEnable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        double mousePosX;
        double mousePosY;
        glfwGetCursorPos(window, &mousePosX, &mousePosY);
        
        //vMat = glm::translate(vMat, glm::vec3(0.0f, 0.0f, 0.5f));
        deltaTime = (1.0 / 60.0);

        if (!rotateCamera)
        {
            yaw = 0.0f;
            pitch = 0.0f;
        }
        
        if (rotateCamera)
        {
            camera->update(cameraScrollDirection, yaw, pitch, lastX, lastY, forwardCamera);
        } else{
            camera->setLastMousePos(lastX, lastY);
        }
        if (forwardCamera)
        camera->moveForward(cameraScrollDirection);
        forwardCamera = false;


        simulation.update(glRenderer, camera->getViewMatrix(), (float) deltaTime);


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
    lastX = xpos;
    lastY = ypos;
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