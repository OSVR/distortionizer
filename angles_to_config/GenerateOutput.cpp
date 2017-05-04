/** @file
    @brief Implementation

    @date 2015-2017

    @author
    Sensics, Inc.
    <http://sensics.com/osvr>
*/

// Copyright 2015-2017 Sensics, Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//        http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#define SENSICS_ENABLE_MESH_DEBUG

// Internal Includes
#include "GenerateOutput.h"

// Library/third-party includes
// - none

// Standard includes
#include <iomanip>
#include <iostream>

#ifdef SENSICS_ENABLE_MESH_DEBUG
#include <fstream>
#include <sstream>
#endif // SENSICS_ENABLE_MESH_DEBUG

static const auto PRECISION = 4;

void writeMesh(std::ostream& s, MeshDescription const& mesh) {
#ifdef SENSICS_ENABLE_MESH_DEBUG
    static int num = 0;
    num++;
    std::ostringstream ss;
    ss << "debugmesh" << num << ".csv";
    std::cerr << "Additionally writing mesh data to " << ss.str() << std::endl;
    std::ofstream debugCsv(ss.str());
    debugCsv << "u1,v1,u2,v2" << std::endl;
#endif // SENSICS_ENABLE_MESH_DEBUG
    s << "[" << std::endl;
    for (size_t i = 0; i < mesh.size(); i++) {
        if (i == 0) {
            s << " ";
        } else {
            s << ",";
        }
        s << std::setprecision(PRECISION) << "[ [" << mesh[i][0][0] << "," << mesh[i][0][1] << "], [" << mesh[i][1][0]
          << "," << mesh[i][1][1] << "] ]" << std::endl;
#ifdef SENSICS_ENABLE_MESH_DEBUG
        debugCsv << std::setprecision(PRECISION) << mesh[i][0][0] << "," << mesh[i][0][1] << "," << mesh[i][1][0] << ","
                 << mesh[i][1][1] << std::endl;
#endif // SENSICS_ENABLE_MESH_DEBUG
    }
    s << "]" << std::endl;
}
