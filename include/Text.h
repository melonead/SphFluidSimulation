#include <glm/glm.hpp>
#include "ft2build.h"
#include "shader.h"

class Text
{
    Text();
    void load_glyphs();
    void render(const std::string& text, float x, float y, float scale, glm::vec3 color);
};
