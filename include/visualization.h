#pragma once
#include "Renderer.h"


class Visualization
{
public:

    Visualization() = default;
    Visualization(
        float cellSize,
        float containerWidth,
        float containerHeight,
        float containerPosX,
        float containerPosY,
        Welol::Renderer& rend
    )
        : cSize(cellSize)
    {
        setRenderer(rend);
        /*
            Grid render operation set up.
        */

        std::vector<float> verticalLine = {
            -containerPosX, containerPosY, 0.0f,
            -containerPosX, containerPosY - containerHeight, 0.0f,
        };

        std::vector<float> horizontalLine = {
            containerPosX, containerPosY, 0.0f,
            containerPosX + containerWidth, containerPosY, 0.0f,
        };

        unsigned int numberOfVerticesPerLine = 2;

        gridInfo.numberOfLines = (int) (containerWidth / cellSize) + (containerHeight / cellSize);
        gridInfo.verticalLines = (int) (containerWidth / cellSize);

        
        Grid = Welol::RenderOperation(Welol::WL_LINES, numberOfVerticesPerLine, 0, gridInfo.numberOfLines + 1, true, false);
        Welol::VertexAttribute verticalLinesAtt{0, Welol::WL_FLOAT3, verticalLine.data(), numberOfVerticesPerLine, false};
        Welol::VertexAttribute horizontalLineAtt{ 1, Welol::WL_FLOAT3, horizontalLine.data(), numberOfVerticesPerLine, false };
        Grid.addVertexAttribute(verticalLinesAtt);
        Grid.addVertexAttribute(horizontalLineAtt);
        renderer.initializeRenderOperation(Grid);


        // grid shader
        std::string gridVertexShaderPath = "C:\\Users\\brian\\programming_projects\\WelolRenderer\\WelolRenderer\\FluidSim\\shaders\\gridVertex.glsl";
        std::string gridFragmentShaderPath = "C:\\Users\\brian\\programming_projects\\WelolRenderer\\WelolRenderer\\FluidSim\\shaders\\gridFragment.glsl";

        gridShader = Shader{gridVertexShaderPath, gridFragmentShaderPath};

        std::vector<float> quadVertices = {
            -0.05f, 0.05f, 
            0.05f,-0.05f, 
            -0.05f,-0.05f, 
            -0.05f, 0.05f, 
            0.05f,-0.05f, 
            0.05f, 0.05f
        };

        std::vector<float> quadLoopVertices = {
            -0.05f, 0.05f, 
            0.05f, 0.05f,
            0.05f,-0.05f, 
            -0.05f,-0.05f
        };

        Welol::VertexAttribute quadPositions{0, Welol::WL_FLOAT2, quadVertices.data(), 6, false};
        Circle.addVertexAttribute(quadPositions);
        renderer.initializeRenderOperation(Circle);

        
        // circle shader
        std::string circleVertexShaderPath = "C:\\Users\\brian\\programming_projects\\WelolRenderer\\WelolRenderer\\FluidSim\\shaders\\circleVertex.glsl";
        std::string circleFragmentShaderPath = "C:\\Users\\brian\\programming_projects\\WelolRenderer\\WelolRenderer\\FluidSim\\shaders\\circleFragment.glsl";
        
        circleShader = Shader{circleVertexShaderPath, circleFragmentShaderPath};

        std::vector<float> quadVertices1 = {
            -0.05f, 0.05f, 
            0.05f,-0.05f, 
            -0.05f,-0.05f, 
            -0.05f, 0.05f, 
            0.05f,-0.05f, 
            0.05f, 0.05f
        };
        
        Welol::VertexAttribute quadLoopPositions{0, Welol::WL_FLOAT2, quadVertices.data(), 6, false};
        Rectangle.addVertexAttribute(quadLoopPositions);
        renderer.initializeRenderOperation(Rectangle);

        // circle shader
        std::string rectVertexShaderPath = "C:\\Users\\brian\\programming_projects\\WelolRenderer\\WelolRenderer\\FluidSim\\shaders\\rectangleVertex.glsl";
        std::string rectFragmentShaderPath = "C:\\Users\\brian\\programming_projects\\WelolRenderer\\WelolRenderer\\FluidSim\\shaders\\rectangleFragment.glsl";

        rectangleShader = Shader{rectVertexShaderPath, rectFragmentShaderPath};
    }

    /*
        The grid represents the hash table.
    */
    void drawGrid(Welol::Renderer& renderer)
    {
        gridShader.use();
        gridShader.setMatrix4fv("view", viewMatrix);
        gridShader.setMatrix4fv("projection", perspectiveMatrix);
        gridShader.setFloat("lineSpacing", cSize);
        gridShader.setInt("numberOfLines", gridInfo.numberOfLines);
        gridShader.setInt("numberOfVerticalLines", gridInfo.verticalLines);

        renderer.render(Grid);
    }

    void drawRectangle(float x, float y, float width, float height)
    {
        // float lineThickness = 0.175f;

        // float quadWidth = 0.05f;
        // float scaleX = (width * 0.5f) / quadWidth;
        // float scaleY = (height * 0.5f) / quadWidth;
        
        // rectModel = glm::translate(glm::mat4(1.0f), glm::vec3(x + (width * 0.5f), y - (height * 0.5f), 0.0f));
        // rectModel = glm::scale(rectModel, glm::vec3(scaleX, scaleY, 1.0f));

        // rectangleShader.use();
        // rectangleShader.setMatrix4fv("view", viewMatrix);
        // rectangleShader.setMatrix4fv("projection", perspectiveMatrix);
        // rectangleShader.setMatrix4fv("model", rectModel);

        // rectangleShader.setInt("cx", x);
        // rectangleShader.setInt("cy", y);

        // rectangleShader.setInt("width", width);
        // rectangleShader.setInt("height", height);

        float lineThickness = 0.075f;

        float quadWidth = 0.05f;
        float scaleX = (width * 0.5f) / quadWidth;
        float scaleY = (height * 0.5f) / quadWidth;
        
        rectModel = glm::translate(glm::mat4(1.0f), glm::vec3(x, y, 0.0f));
        rectModel = glm::scale(rectModel, glm::vec3(scaleX, scaleY, 1.0f));

        rectangleShader.use();
        rectangleShader.setMatrix4fv("view", viewMatrix);
        rectangleShader.setMatrix4fv("projection", perspectiveMatrix);
        rectangleShader.setMatrix4fv("model", rectModel);

        rectangleShader.setFloat("cx", x);
        rectangleShader.setFloat("cy", y);

        rectangleShader.setFloat("width", width);
        rectangleShader.setFloat("height", height);
  
        renderer.render(Rectangle);
    }

    void drawCircle(float x, float y, float radius)
    {
        float quadWidth = 0.05f;
        float scale = radius / quadWidth;

        circleModel = glm::translate(glm::mat4(1.0f), glm::vec3(x, y, 0.0f));
        circleModel = glm::scale(circleModel, glm::vec3(scale, scale, 1.0f));
        circleShader.use();
        circleShader.setMatrix4fv("view", viewMatrix);
        circleShader.setMatrix4fv("projection", perspectiveMatrix);
        circleShader.setMatrix4fv("model", circleModel);
        circleShader.setFloat("radius", radius);
        circleShader.setFloat("x", x);
        circleShader.setFloat("y", y);
        renderer.render(Circle);
    }

    void setPerspectiveMatrix(glm::mat4& matrix)
    {
        perspectiveMatrix = matrix;
    }

    void setViewMatrix(glm::mat4& matrix)
    {
        viewMatrix = matrix;
    }

    void setRenderer(Welol::Renderer& rend)
    {
        renderer = rend;
    }

private:
    float cSize;
    unsigned int maxCircles {20};
    Welol::RenderOperation Circle{Welol::WL_TRIANGLES, 6, 0, maxCircles, true, false};
    Welol::RenderOperation Rectangle{Welol::WL_TRIANGLES, 6, 0, maxCircles, false, false};
    Welol::RenderOperation Grid;
    Shader gridShader;
    Shader circleShader;
    Shader rectangleShader;
    glm::mat4 rectModel{1.0f};
    glm::mat4 circleModel{1.0f};
    glm::mat4 perspectiveMatrix;
    glm::mat4 viewMatrix;
    Welol::Renderer renderer;


    struct {
        int numberOfLines{0};
        int verticalLines{0};
    } gridInfo;
};