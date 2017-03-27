/** @file
    @brief Produces an example table mapping angles to locations, which
           can be fed into AnglesToConfig to produce a mesh.

    @date 2016

    @author
    Russ Taylor working through ReliaSolve.com for Sensics, Inc.
    <http://sensics.com/osvr>
*/

// Copyright 2016 Sensics, Inc.
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

// Internal Includes

// Standard includes
#include <cmath>
#include <cstdlib> // For exit()
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// Global constants and variables
static bool g_verbose = false;
#define MY_PI (4.0 * atan(1.0))

void Usage(const std::string& name) {
    std::cerr << "Usage: " << name << std::endl
              << "  This program produces on standard output an example angles"
              << " to display location file with a 90 degree field of view and."
              << " no distortion." << std::endl;
    exit(1);
}

int main(int argc, char* argv[]) {
    // Set defaults
    double xFOVDeg = 90.0;
    double yFOVDeg = 90.0;
    int xCount = 11;
    int yCount = 11;

    // Parse the command line
    int realParams = 0;
    for (int i = 1; i < argc; i++) {
        if (std::string("-verbose") == argv[i]) {
            g_verbose = true;
        } else if ((argv[i][0] == '-') && (atof(argv[i]) == 0.0)) {
            Usage(argv[0]);
        } else {
            switch (++realParams) {
            case 1:
            default:
                Usage(argv[0]);
            }
        }
    }
    if (realParams != 0) {
        Usage(argv[0]);
    }

    // Produce a mesh with the specified number of elements in X and
    // Y that cover the specified total FOV in X and Y.
    double xMin = -xFOVDeg / 2;
    double xStep = xFOVDeg / (xCount - 1);
    double yMin = -yFOVDeg / 2;
    double yStep = yFOVDeg / (yCount - 1);
    for (int x = 0; x < xCount; x++) {
        double xDeg = xMin + x * xStep;
        double xRad = xDeg * MY_PI / 180;
        for (int y = 0; y < yCount; y++) {
            double yDeg = yMin + y * yStep;
            double yRad = yDeg * MY_PI / 180;
            std::cout << xDeg << " " << yDeg << " " << tan(xRad) << " " << tan(yRad) << std::endl;
        }
    }

    return 0;
}
