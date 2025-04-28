#include "Text.h"
#include "glm/gtc/matrix_transform.hpp"
#include "shader.h"

/*
Text::Text() : shader{ "shaders/textVert.glsl", "shaders/textFrag.glsl" } {

    //freetype loading
    if(FT_Init_FreeType(&ft)) {
        std::cout << "COULD NOT INITILALIZE FT LIBRARY" << std::endl;
    }

    if(FT_New_Face(ft, "data/fonts/Antonio-Regular.ttf", 0, &face)) {
        std::cout << "FAILED TO LOAD THE FONT" << std::endl;
    }

    FT_Set_Pixel_Sizes(face, 0, 24);

    // load the ascii characters
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    load_glyphs();

    projection = glm::ortho(0.0f, 600.0f, 0.0f, 600.0f);

    vao.bind();
    vbo.bind();
    vbo.attach_none(24);
    vbo.unbind();
    vao.unbind();

}

void Text::load_glyphs() {
    for(unsigned char c = 0; c < 128; c++) {
        if (FT_Load_Char(face, c, FT_LOAD_RENDER)) {
            std::cout << "failed to load " << c << " glyph" << std::endl;
            continue;
        }
        // texture
        unsigned int texture_id;
        glGenTextures(1, &texture_id);
        glBindTexture(GL_TEXTURE_2D, texture_id);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, face->glyph->bitmap.width, face->glyph->bitmap.rows, 0, GL_RED, GL_UNSIGNED_BYTE, face->glyph->bitmap.buffer);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        Character ch = {
            texture_id,
            glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
            glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
            face->glyph->advance.x
        };
        Characters.insert(std::pair<char, Character>(c, ch));
    }

    FT_Done_Face(face);
    FT_Done_FreeType(ft);
}

float Text::render(const std::string& text, float x, float y, float scale, glm::vec3 color) {
    shader.use();

    float original_x = x;
    shader.setVec3f("text_color", color);
    shader.setMatrix4fv("projection", projection);
    glActiveTexture(GL_TEXTURE0);
    vao.bind();

    // iterate through all characters
    float length = 0.0f;
    std::string::const_iterator c;
    for (c = text.begin(); c != text.end(); c++) 
    {
        Character ch = Characters[*c];

        float xpos = x + ch.bearing.x * scale;
        float ypos = y - (ch.size.y - ch.bearing.y) * scale;

        float w = ch.size.x * scale;
        float h = ch.size.y * scale;
        // update VBO for each character
        float vertices[6][4] = {
            { xpos,     ypos + h,   0.0f, 0.0f },            
            { xpos,     ypos,       0.0f, 1.0f },
            { xpos + w, ypos,       1.0f, 1.0f },

            { xpos,     ypos + h,   0.0f, 0.0f },
            { xpos + w, ypos,       1.0f, 1.0f },
            { xpos + w, ypos + h,   1.0f, 0.0f }           
        };
        // render glyph texture over quad
        //setFloat("text", ch.Texture_ID, shader_id);
        glBindTexture(GL_TEXTURE_2D, ch.Texture_ID);
        // update content of VBO memory
        vbo.bind();
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices);
        vbo.unbind();
        // render quad
        glDrawArrays(GL_TRIANGLES, 0, 6);
        // now advance cursors for next glyph (note that advance is number of 1/64 pixels)
        x += (ch.advance >> 6) * scale; // bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))

        length += 6.0f;
    }
    vao.unbind();
    glBindTexture(GL_TEXTURE_2D, 0);
    return length;
}

*/