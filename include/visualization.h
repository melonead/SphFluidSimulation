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
        Welol::Renderer& renderer
    )
        : cSize(cellSize)
    {

        /*
            Grid render operation set up.
        */

        std::vector<float> verticalLine = {
            -containerPosX, containerPosY, 0.0f,
            -containerPosX, containerPosY - containerWidth, 0.0f,
        };

        std::vector<float> horizontalLine = {
            containerPosX, containerPosY, 0.0f,
            containerPosX + containerWidth, containerPosY, 0.0f,
        };

        unsigned int numberOfVerticesPerLine = 2;
        gridInfo.numberOfLines = (int) (containerWidth / cellSize) + (containerHeight / cellSize);
        
        Grid = Welol::RenderOperation(Welol::WL_LINES, numberOfVerticesPerLine, 0, gridInfo.numberOfLines, true, false);
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

        Welol::VertexAttribute quadPositions{Welol::WL_TRIANGLES, Welol::WL_FLOAT2, quadVertices.data(), 6, false};
        Circle.addVertexAttribute(quadPositions);
        renderer.initializeRenderOperation(Circle);

        // circle shader
        std::string circleVertexShaderPath = "C:\\Users\\brian\\programming_projects\\WelolRenderer\\WelolRenderer\\FluidSim\\shaders\\circleVertex.glsl";
        std::string circleFragmentShaderPath = "C:\\Users\\brian\\programming_projects\\WelolRenderer\\WelolRenderer\\FluidSim\\shaders\\circleFragment.glsl";

        circleShader = Shader{circleVertexShaderPath, circleFragmentShaderPath};
    }

    /*
        The grid represents the hash table.
    */
    void drawGrid(Welol::Renderer& renderer, glm::mat4& view, glm::mat4& projection)
    {
        gridShader.use();
        gridShader.setMatrix4fv("view", view);
        gridShader.setMatrix4fv("projection", projection);
        gridShader.setFloat("lineSpacing", cSize);
        gridShader.setInt("numberOfLines", gridInfo.numberOfLines);

        renderer.render(Grid);
    }

    /*
        Draw the container containing the fluid.
    */
    void drawContainer()
    {

    }

    void drawCircle(float x, float y, float radius, glm::mat4& view, glm::mat4& projection, Welol::Renderer& renderer)
    {
        circleShader.use();
        circleShader.setMatrix4fv("view", view);
        circleShader.setMatrix4fv("projection", projection);
        circleShader.setFloat("radius", radius);
        circleShader.setFloat("x", x);
        circleShader.setFloat("y", y);
        renderer.render(Circle);
    }

private:
    float cSize;
    unsigned int maxCircles {20};
    Welol::RenderOperation Circle{Welol::WL_TRIANGLES, 6, 0, maxCircles, true, false};
    Welol::RenderOperation Grid;
    Shader gridShader;
    Shader circleShader;

    struct {
        int numberOfLines{0};
    } gridInfo;
};