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

#pragma once
#include <string>

// Forward declaration to portions declared in the .cpp file, to avoid
// #include conflicts with other source code.
class Undistort_Shader_Private;

class Undistort_Shader {
  public:
    Undistort_Shader(const std::string& vert_shader_file_name = "./quadratic_tri_color_vert.glsl",
                     const std::string& frag_shader_file_name = "./quadratic_tri_color_frag.glsl");
    ~Undistort_Shader();

    // Use the shader for rendering.
    void useShader();

    // Set individual shader parameters
    // XXX window/texture size?  COP.  K params for each color.

    // Set the Parameters for the shader
    void setK1Red(float val);
    void setK1Green(float val);
    void setK1Blue(float val);

    // Sets the values back to their defaults.
    void SetDefaultValues(void);

    // No shader yet loaded.
    static const int NO_SHADER = 9999;

  protected:
    // Read a shader string from a file into a string.  Returns an empty
    // string on failure.
    std::string readShaderFromFile(const std::string& filename);

    // Load, compile, and link the shaders.
    static int loadShaders(const char* vertexShader, const char* fragmentShader);

  private:
    Undistort_Shader_Private* d_p; //< Private objects requiring GL/GLEW
};