/** @file
    @brief Produces a distortion mesh and partial display description
           from a table of display locations to angle inputs.

    @date 2015

    @author
    Russ Taylor working through ReliaSolve.com for Sensics, Inc.
    <http://sensics.com/osvr>
*/

// Copyright 2015 Sensics, Inc.
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
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <stdlib.h> // For exit()

// Global constants and variables
static bool g_verbose = false;

// Forward declarations
static int testAlgorithms();

#include "types.h"
#include "helper.h"

static void writeMesh(std::ostream &s, MeshDescription const &mesh)
{
  s << "[" << std::endl;
  for (size_t i = 0; i < mesh.size(); i++) {
    if (i == 0) { s << " "; }
    else { s << ","; }
    s << std::setprecision(4)
      << "[ [" << mesh[i][0][0] << "," << mesh[i][0][1] << "], ["
      << mesh[i][1][0] << "," << mesh[i][1][1] << "] ]"
      << std::endl;
  }
  s << "]" << std::endl;
}

// Produce a mapping that is reflected around X=0 in both angles and
// screen coordinates.
static std::vector<Mapping> reflect_mapping(std::vector<Mapping> mapping)
{
  std::vector<Mapping> ret;
  for (size_t i = 0; i < mapping.size(); i++) {
    ret.push_back(mapping[i]);
    ret[i].xyLatLong.longitude *= -1;
    ret[i].xyLatLong.x *= -1;
  }

  return ret;
}

void Usage(std::string name)
{
  std::cerr << "Usage: " << name
    << " [-eye right|left] (default is right)"
    << " [-depth_meters D] (default is 2.0)"
    << " [-mm] (screen distance units in the config file, default is meters)"
    << " [-verbose] (default is not)"
    << " [-screen screen_left_meters screen_bottom_meters screen_right_meters screen_top_meters]"
    <<   " (default auto-compute based on ranges seen)"
    << std::endl
    << "  This program reads from standard input a configuration that has a list of" << std::endl
    << "x,y screen coordinates in meters followed by long,lat angles in" << std::endl
    << "degrees where (0,0) is straight ahead from the eye, positive" << std::endl
    << "longitude is left and positive latitude is up." << std::endl
    << "  It produces on standard output a partial OSVR display configuration file." << std::endl
    << std::endl;
  exit(1);
}

int main(int argc, char *argv[])
{
  // Parse the command line
  bool useRightEye = true;
  bool computeBounds = true;
  double left, right, bottom, top;
  double depth = 2.0;
  double toMeters = 1.0;
  int realParams = 0;
  for (int i = 1; i < argc; i++) {
    if (std::string("-mm") == argv[i]) {
      toMeters = 1e-3;  // Convert input in millimeters to meters
    } else if (std::string("-verbose") == argv[i]) {
      g_verbose = true;
    } else if (std::string("-depth_meters") == argv[i]) {
      if (++i >= argc) { Usage(argv[0]); }
      depth = atof(argv[i]);
    } else if (std::string("-screen") == argv[i]) {
      computeBounds = false;
      if (++i >= argc) { Usage(argv[0]); }
      left = atof(argv[i]);
      if (++i >= argc) { Usage(argv[0]); }
      bottom = atof(argv[i]);
      if (++i >= argc) { Usage(argv[0]); }
      right = atof(argv[i]);
      if (++i >= argc) { Usage(argv[0]); }
      top = atof(argv[i]);
    }
    else if (std::string("-eye") == argv[i]) {
      if (++i >= argc) { Usage(argv[0]); }
      std::string eye = argv[i];
      if (eye == "left") {
        useRightEye = false;
      } else if (eye == "right") {
        useRightEye = true;
      } else {
        std::cerr << "Bad value for -eye: " << eye << ", expected left or right" << std::endl;
        Usage(argv[0]);
      }
    } else if ((argv[i][0] == '-') && (atof(argv[i]) == 0.0) ) {
      Usage(argv[0]);
    }
    else switch (++realParams) {
    case 1:
    default:
      Usage(argv[0]);
    }
  }
  if (realParams != 0) { Usage(argv[0]); }

  //====================================================================
  // Run our algorithm test to make sure things are working properly.
  int ret;
  if ((ret = testAlgorithms()) != 0) {
    std::cerr << "Error testing basic algorithms, code " << ret << std::endl;
    return 100;
  }

  //====================================================================
  // Parse the angle-configuration information from standard input.  Expect white-space
  // separation between numbers and also between entries (which may be on separate
  // lines).
  std::vector<Mapping> mapping;
  while (!std::cin.eof()) {
    // Read the mapping info from the input file.
    Mapping map;
    std::cin >> map.xyLatLong.longitude >> map.xyLatLong.latitude >> map.xyLatLong.x >> map.xyLatLong.y;
    mapping.push_back(map);
  }
  // There will have been one extra added, when running into EOF.
  mapping.pop_back();
  if (g_verbose) {
    std::cerr << "Found " << mapping.size() << " points" << std::endl;
  }
  if (mapping.size() == 0) {
    std::cerr << "Error: No input points found" << std::endl;
    return 2;
  }

  //====================================================================
  // If we've been asked to auto-range the screen coordinates, compute
  // them here.
  if (computeBounds) {
    left = right = mapping[0].xyLatLong.x;
    bottom = top = mapping[0].xyLatLong.y;
    for (size_t i = 1; i < mapping.size(); i++) {
      double x = mapping[i].xyLatLong.x;
      double y = mapping[i].xyLatLong.y;
      if (x < left) { left = x; }
      if (x > right) { right = x; }
      if (y < bottom) { bottom = y; }
      if (y > top) { top = y; }
    }
    left *= toMeters;
    right *= toMeters;
    bottom *= toMeters;
    top *= toMeters;
  }
  if (g_verbose) {
    std::cerr << "Left, bottom, right, top = " << left << ", "
      << bottom << ", " << right << ", " << top << std::endl;
  }

  //====================================================================
  // Make an inverse mapping for the opposite eye.  Invert around X in
  // angle and viewing direction.  Depending on whether we are using the
  // left or right eye, set the eyes appropriately.
  //  Also make a different set of screen boundaries for each, again
  // inverting around X = 0.
  std::vector<Mapping> leftMapping;
  std::vector<Mapping> rightMapping;
  double leftScreenLeft, leftScreenRight, leftScreenBottom, leftScreenTop;
  double rightScreenLeft, rightScreenRight, rightScreenBottom, rightScreenTop;
  rightScreenBottom = leftScreenBottom = bottom;
  rightScreenTop = leftScreenTop = top;
  if (useRightEye) {
    rightMapping = mapping;
    rightScreenLeft = left;
    rightScreenRight = right;

    leftMapping = reflect_mapping(mapping);
    leftScreenLeft = -right;
    leftScreenRight = -left;
  } else {
    leftMapping = mapping;
    leftScreenLeft = left;
    leftScreenRight = right;

    rightMapping = reflect_mapping(mapping);
    rightScreenLeft = -right;
    rightScreenBottom = -left;
  }

  //====================================================================
  // Convert the input values into normalized coordinates and into 3D
  // locations.
  convert_to_normalized_and_meters(leftMapping, toMeters, depth,
    leftScreenLeft, leftScreenBottom, leftScreenRight, leftScreenTop);
  convert_to_normalized_and_meters(rightMapping, toMeters, depth,
    rightScreenLeft, rightScreenBottom, rightScreenRight, rightScreenTop);

  //====================================================================
  // Determine the screen description and distortion mesh based on the
  // input points and screen parameters.
  ScreenDescription leftScreen, rightScreen;
  MeshDescription leftMesh, rightMesh;
  if (!findScreenAndMesh(leftMapping, leftScreenLeft, leftScreenBottom,
    leftScreenRight, leftScreenTop, leftScreen, leftMesh, g_verbose)) {
    std::cerr << "Error: Could not find left screen or mesh" << std::endl;
    return 3;
  }
  if (leftMesh.size() != mapping.size()) {
    std::cerr << "Error: Left mesh size " << leftMesh.size()
      << " does not match mapping size" << mapping.size() << std::endl;
    return 4;
  }
  if (!findScreenAndMesh(rightMapping, rightScreenLeft, rightScreenBottom,
    rightScreenRight, rightScreenTop, rightScreen, rightMesh, g_verbose)) {
    std::cerr << "Error: Could not find right screen or mesh" << std::endl;
    return 5;
  }
  if (rightMesh.size() != mapping.size()) {
    std::cerr << "Error: Right mesh size " << rightMesh.size()
      << " does not match mapping size" << mapping.size() << std::endl;
    return 6;
  }

  //====================================================================
  // Construct Json screen description.
  // We do this by hand rather than using JsonCPP because we need
  // to control the printed precision of the numbers to avoid making
  // a huge file.
  std::cout << "{" << std::endl;
  std::cout << " \"display\": {" << std::endl;
  std::cout << "  \"hmd\": {" << std::endl;

  std::cout << "   \"field_of_view\": {" << std::endl;
  std::cout << "    \"monocular_horizontal\": "
    << rightScreen.hFOVDegrees
    << "," << std::endl;
  std::cout << "    \"monocular_vertical\": "
    << rightScreen.vFOVDegrees
    << "," << std::endl;
  std::cout << "    \"overlap_percent\": "
    << rightScreen.overlapPercent
    << "," << std::endl;
  std::cout << "    \"pitch_tilt\": 0" << std::endl;
  std::cout << "   }," << std::endl; // field_of_view

  std::cout << "   \"distortion\": {" << std::endl;
  std::cout << "    \"type\": \"mono_point_samples\"," << std::endl;
  std::cout << "    \"mono_point_samples\": [" << std::endl;
  writeMesh(std::cout, leftMesh);
  std::cout << "," << std::endl;
  // @todo Figure out which mesh to invert (left or right) and write
  // two different meshes, rather than the same one twice.
  writeMesh(std::cout, rightMesh);
  std::cout << "    ]" << std::endl; // mono_point_samples
  std::cout << "   }," << std::endl; // distortion

  std::cout << "   \"eyes\": [" << std::endl;
  std::cout << "    {" << std::endl;
  std::cout << "     \"center_proj_x\": "
    << leftScreen.xCOP
    << "," << std::endl;
  std::cout << "     \"center_proj_y\": "
    << leftScreen.yCOP
    << "," << std::endl;
  std::cout << "     \"rotate_180\": 0" << std::endl;
  std::cout << "    }," << std::endl;
  std::cout << "    {" << std::endl;
  std::cout << "     \"center_proj_x\": "
    << rightScreen.xCOP
    << "," << std::endl;
  std::cout << "     \"center_proj_y\": "
    << rightScreen.yCOP
    << "," << std::endl;
  std::cout << "     \"rotate_180\": 0" << std::endl;
  std::cout << "    }" << std::endl;
  std::cout << "   ]" << std::endl; // eyes

  std::cout << "  }" << std::endl;  // hmd
  std::cout << " }" << std::endl;   // display
  std::cout << "}" << std::endl;    // Closes outer object

  return 0;
}

static bool small(double d)
{
  return fabs(d) <= 1e-5;
}

static int testAlgorithms()
{
  if (g_verbose) {
    std::cerr << "=================== Starting testAlgorithms()" << std::endl;
  }

  // Construct a set of points that should make a square screen with no distortion
  // that has a 90 degree horizontal and vertical field of view and the identity
  // distortion correction and that is oriented straight ahead.
  Mapping p1(XYLatLong(0, 0, -45, -45), XYZ(-1, -1, -1));
  Mapping p2(XYLatLong(1, 0,  45, -45), XYZ( 1, -1, -1));
  Mapping p3(XYLatLong(1, 1,  45,  45), XYZ( 1,  1, -1));
  Mapping p4(XYLatLong(0, 1, -45,  45), XYZ(-1,  1, -1));
  std::vector<Mapping> mapping;
  mapping.push_back(p1);
  mapping.push_back(p2);
  mapping.push_back(p3);
  mapping.push_back(p4);

  // Find the screen associated with this mapping.
  ScreenDescription screen;
  MeshDescription mesh;
  if (!findScreenAndMesh(mapping, 0, 0, 1, 1, screen, mesh, g_verbose)) {
    std::cerr << "testAlgorithms(): Could not find screen" << std::endl;
    return 100;
  }

  // Make sure the screen has the expected behavior.
  if (!small(screen.hFOVDegrees - 90)) {
    std::cerr << "testAlgorithms(): hFOV not near 90: " << screen.hFOVDegrees << std::endl;
    return 201;
  }
  if (!small(screen.vFOVDegrees - 90)) {
    std::cerr << "testAlgorithms(): vFOV not near 90: " << screen.vFOVDegrees << std::endl;
    return 202;
  }
  if (!small(screen.overlapPercent - 100)) {
    std::cerr << "testAlgorithms(): Overlap percent not near 100: " << screen.overlapPercent << std::endl;
    return 203;
  }
  if (!small(screen.xCOP - 0.5)) {
    std::cerr << "testAlgorithms(): xCOP not near 0.5: " << screen.xCOP << std::endl;
    return 204;
  }
  if (!small(screen.yCOP - 0.5)) {
    std::cerr << "testAlgorithms(): yCOP not near 0.5: " << screen.yCOP << std::endl;
    return 205;
  }

  // Make sure the mesh has the expected behavior.
  if (mesh.size() != mapping.size()) {
    std::cerr << "testAlgorithms(): Mesh size does not match mapping size: " << mesh.size() << std::endl;
    return 301;
  }
  for (int entry = 0; entry < mesh.size(); entry++) {
    size_t outIndex = 1;
    if (!small(mesh[entry][outIndex][0] - mapping[entry].xyLatLong.x)) {
      std::cerr << "testAlgorithms(): X normalized mesh mismatch for element: " << entry << std::endl;
      return 400 + entry;
    }
    if (!small(mesh[entry][outIndex][1] - mapping[entry].xyLatLong.y)) {
      std::cerr << "testAlgorithms(): Y normalized mesh mismatch for element: " << entry
        << " (got " << mesh[entry][outIndex][1]
        << ", expected " << mapping[entry].xyLatLong.y << ")" << std::endl;
      return 500 + entry;
    }
  }

  // Construct a set of points that should make a square screen with no distortion
  // that has a 90 degree horizontal and vertical field of view and the identity
  // distortion correction, but which is rotated 45 degrees with respect to straight
  // ahead.
  Mapping rp1(XYLatLong(0, 0, 0, -45), XYZ(0, -1, -sqrt(2)));
  Mapping rp2(XYLatLong(1, 0, 90, -45), XYZ(sqrt(2), -1, 0));
  Mapping rp3(XYLatLong(1, 1, 90, 45), XYZ(sqrt(2), 1, 0));
  Mapping rp4(XYLatLong(0, 1, 0, 45), XYZ(0, 1, -sqrt(2)));
  std::vector<Mapping> rmapping;
  rmapping.push_back(rp1);
  rmapping.push_back(rp2);
  rmapping.push_back(rp3);
  rmapping.push_back(rp4);

  // Find the screen associated with this mapping.
  ScreenDescription rscreen;
  MeshDescription rmesh;
  if (!findScreenAndMesh(rmapping, 0, 0, 1, 1, rscreen, rmesh, g_verbose)) {
    std::cerr << "testAlgorithms(): Could not find rotated screen" << std::endl;
    return 600;
  }

  // Make sure the screen has the expected behavior.
  if (!small(rscreen.hFOVDegrees - 90)) {
    std::cerr << "testAlgorithms(): rotated hFOV not near 90: " << rscreen.hFOVDegrees << std::endl;
    return 701;
  }
  if (!small(rscreen.vFOVDegrees - 90)) {
    std::cerr << "testAlgorithms(): rotated vFOV not near 90: " << rscreen.vFOVDegrees << std::endl;
    return 702;
  }
  if (!small(rscreen.overlapPercent - 0)) {
    std::cerr << "testAlgorithms(): Rotated overlap percent not near 0: " << rscreen.overlapPercent << std::endl;
    return 703;
  }
  if (!small(rscreen.xCOP - 0.5)) {
    std::cerr << "testAlgorithms(): Rotated xCOP not near 0.5: " << rscreen.xCOP << std::endl;
    return 704;
  }
  if (!small(rscreen.yCOP - 0.5)) {
    std::cerr << "testAlgorithms(): Rotated yCOP not near 0.5: " << rscreen.yCOP << std::endl;
    return 705;
  }

  // Make sure the mesh has the expected behavior.
  if (rmesh.size() != rmapping.size()) {
    std::cerr << "testAlgorithms(): Rotated mesh size does not match mapping size: " << mesh.size() << std::endl;
    return 801;
  }
  for (int entry = 0; entry < rmesh.size(); entry++) {
    size_t outIndex = 1;
    if (!small(rmesh[entry][outIndex][0] - rmapping[entry].xyLatLong.x)) {
      std::cerr << "testAlgorithms(): Rotated X normalized mesh mismatch for element: " << entry << std::endl;
      return 800 + entry;
    }
    if (!small(rmesh[entry][outIndex][1] - rmapping[entry].xyLatLong.y)) {
      std::cerr << "testAlgorithms(): Rotated Y normalized mesh mismatch for element: " << entry
        << " (got " << rmesh[entry][outIndex][1]
        << ", expected " << rmapping[entry].xyLatLong.y << ")" << std::endl;
      return 900 + entry;
    }
  }

  if (g_verbose) {
    std::cerr << "=================== Successfully finished testAlgorithms()" << std::endl;
  }
  return 0;
}
