/** @file
    @brief Header

    @date 2014

    @author
    Russell Taylor
    <russ@reliasolve.com>
    <http://sensics.com>
*/

// Copyright 2014 Sensics, Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "undistort_shader.h"

/// Must have a blank line between glew and gl
#include <GL/glew.h>

#include <GL/gl.h>
#include <cstdio>
#include <cstdlib>

class Undistort_Shader_Private {
  public:
    Undistort_Shader_Private() : d_shader_id(Undistort_Shader::NO_SHADER){};

    // Read a shader string from a file into a string.  Returns an empty
    // string on failure.
    std::string readShaderFromFile(std::string filename);

    // Load, compile, and link the shaders.
    static int loadShaders(const char* vertexShader, const char* fragmentShader);

    // TODO: Figure out which parameters can be uniform
    GLuint d_shader_id;   //< The index of our shader program
    GLint d_k1RedParam;   //< The location of the K1 param for Red
    GLint d_k1GreenParam; //< The location of the K1 param for Green
    GLint d_k1BlueParam;  //< The location of the K1 param for Blue
    GLint d_centerParam;
    GLint d_radiusParam;

    GLfloat d_k1Red;     // K1 red value to use in shader
    GLfloat d_k1Green;   // K1 green value to use in shader
    GLfloat d_k1Blue;    // K1 blue value to use in shader
    GLfloat d_center[2]; // Center value to use in shader
    GLfloat d_radius;    // Center value to use in shader
};

// XXX List of shader attributes
struct shader_bind_attribute_list {
    GLuint index;
    const char* name; // Was GLchar
};

const shader_bind_attribute_list sbal_undistort[] = {
    {0, "vVertex"}, {1, "vNormal"}, {2, "vTexture0"}, {3, "vVaryingColor"}, {0, NULL}};

// Handle a shader error
// Input
//  shader:   Shader handle
//  name:     name of shader to use for pretty printing errors
static int shaderError(GLuint shader, const char* name) {
    GLint status;

    glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
    if (status == GL_FALSE) {
        char message[2000];
        glGetShaderInfoLog(shader, sizeof(message), NULL, message);
        fprintf(stderr, "Error in %s shader:\n%s\n", name, message);
        return -1;
    }
    return 0;
}

// Read a shader program string from a file.
// Returns an empty string on failure.
std::string Undistort_Shader::readShaderFromFile(std::string filename) {
    std::string ret;

    // TODO: Convert this to using only standard library calls.
    FILE* f = fopen(filename.c_str(), "r");
    if (f == NULL) {
        printf("No shader file with name %s found;", filename.c_str());
        return ret;
    }

    // Read each line of the file and append it to the string.
    char line[4096];
    while (fgets(line, sizeof(line), f) != NULL) {
        ret += line;
    }

    // Done; send the result after closing the file.
    fclose(f);
    return ret;
}

// Load a pair of vertex and fragment shaders
// Inputs
//  vertexShader:   vertex shader source code
//  fragmentShader: fragment shader source code
//  sbal:       list of generic vertex attributes to bind
// Outputs
//  handle of the linked shader program, or NO_SHADER if a problem
int Undistort_Shader::loadShaders(const char* vertexShader, const char* fragmentShader) {
    GLuint program;
    GLint temp;

    GLuint vertexShaderHandle = glCreateShader(GL_VERTEX_SHADER);
    GLuint fragmentShaderHandle = glCreateShader(GL_FRAGMENT_SHADER);

    glShaderSource(vertexShaderHandle, 1, &vertexShader, NULL);
    glShaderSource(fragmentShaderHandle, 1, &fragmentShader, NULL);

    glCompileShader(vertexShaderHandle);
    glCompileShader(fragmentShaderHandle);

    if (shaderError(vertexShaderHandle, "vertex") || shaderError(fragmentShaderHandle, "fragment")) {
        glDeleteShader(vertexShaderHandle);
        glDeleteShader(fragmentShaderHandle);
        return NO_SHADER;
    }

    program = glCreateProgram();
    glAttachShader(program, vertexShaderHandle);
    glAttachShader(program, fragmentShaderHandle);

    const shader_bind_attribute_list* sbal = sbal_undistort;
    if (sbal != NULL) {
        while (sbal->name != NULL) {
            glBindAttribLocation(program, sbal->index, sbal->name);
            sbal++;
        }
    }

    glLinkProgram(program);

    glDeleteShader(vertexShaderHandle);
    glDeleteShader(fragmentShaderHandle);

    glGetProgramiv(program, GL_LINK_STATUS, &temp);
    if (temp == GL_FALSE) {
        glDeleteProgram(program);
        return NO_SHADER;
    }
    return program;
}

// Undistort_Shader Constructor
Undistort_Shader::Undistort_Shader(std::string vert_shader_file_name, std::string frag_shader_file_name)
    : d_p(new Undistort_Shader_Private) {

    // Read the vertex shader and the fragment shader from the specified
    // files.  Bail if we can't get them.
    std::string vertexProgram = readShaderFromFile(vert_shader_file_name);
    std::string fragmentProgram = readShaderFromFile(frag_shader_file_name);
    if ((vertexProgram.size() == 0) || (fragmentProgram.size() == 0)) {
        return;
    }
    glewInit();
    // Load, compile, and link the shaders
    if ((d_p->d_shader_id = loadShaders(vertexProgram.c_str(), fragmentProgram.c_str())) == NO_SHADER) {
        return;
    }

    // Get the uniform variable locations
    d_p->d_k1RedParam = glGetUniformLocation(d_p->d_shader_id, "k1Red");
    d_p->d_k1GreenParam = glGetUniformLocation(d_p->d_shader_id, "k1Green");
    d_p->d_k1BlueParam = glGetUniformLocation(d_p->d_shader_id, "k1Blue");

    // Set the default values
    SetDefaultValues();
}

// Set the Default Values for the Color Processing
void Undistort_Shader::SetDefaultValues(void) {
    setK1Red(0.0f);
    setK1Green(0.0f);
    setK1Blue(0.0f);
}

// Set the K1 parameter for the Red channel
// Input
//  val:    The value to use in the shader
void Undistort_Shader::setK1Red(float val) {
    d_p->d_k1Red = val;
    glUseProgram(d_p->d_shader_id);
    glUniform1f(d_p->d_k1RedParam, d_p->d_k1Red);
    // TODO: Push and pop to go back to the original.
    glUseProgram(0);
}

// Set the K1 parameter for the Green channel
// Input
//  val:    The value to use in the shader
void Undistort_Shader::setK1Green(float val) {
    d_p->d_k1Green = val;
    glUseProgram(d_p->d_shader_id);
    glUniform1f(d_p->d_k1GreenParam, d_p->d_k1Green);
    // TODO: Push and pop to go back to the original.
    glUseProgram(0);
}

// Set the K1 parameter for the Blue channel
// Input
//  val:    The value to use in the shader
void Undistort_Shader::setK1Blue(float val) {
    d_p->d_k1Blue = val;
    glUseProgram(d_p->d_shader_id);
    glUniform1f(d_p->d_k1BlueParam, d_p->d_k1Blue);
    // TODO: Push and pop to go back to the original.
    glUseProgram(0);
}

// Use the shader for rendering (set up the program)
void Undistort_Shader::useShader() {

    glUseProgram(d_p->d_shader_id);

    /* TODO Uniform parameters
        glUniform1i(locCM, 0);
        glUniform1i(locAL, 1);
    */
}

// Destructor
Undistort_Shader::~Undistort_Shader() {
    if (d_p->d_shader_id != NO_SHADER) {
        glDeleteProgram(d_p->d_shader_id);
    }
}
