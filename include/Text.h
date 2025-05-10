#include <glm/glm.hpp>
#include "ft2build.h"
#include FT_FREETYPE_H
#include "shader.h"
#include "Renderer.h"

struct Character
{
    unsigned int Texture_ID;
    glm::ivec2 size;
    glm::ivec2 bearing;
    unsigned int advance;
};

class Text
{
public:
    Text(Welol::Renderer& renderer);
    void load_glyphs();
    float render(const std::string& text, float x, float y, float scale, glm::vec3 color, Welol::Renderer& renderer);
private:
    Shader shader;
    glm::mat4 projection;
    std::map<char, Character> Characters;
    FT_Face face;
    FT_Library ft;
    Welol::RenderOperation fontRenderOperation;

};
