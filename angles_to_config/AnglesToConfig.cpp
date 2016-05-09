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
#include <fstream>
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
    << " [-latlong] (use latitude/longitude angles, default is field angles)"
    << " [-mm] (screen distance units in the config file, default is meters)"
    << " [-verbose] (default is not)"
    << " [-screen screen_left_meters screen_bottom_meters screen_right_meters screen_top_meters]"
    <<   " (default auto-compute based on ranges seen)"
    << " [-verify_angles xx xy yx yy max_degrees]"
    << "   The vector (xx, xy) points in screen space in the direction of +longitude (left)"
    << "   The vector (yx, yy) points in screen space in the direction of +latitude (up)"
    << "   The max_degrees tells how far the screen-space neighbor vector can differ from it corresponding angle-space vector"
    << " [-mono in_config_mono_file_name ] (default standard input)"
    << " [-rgb in_config_red_file_name in_config_green_file_name in_config_blue_file_name]"
    << std::endl
    << "  This program reads one or three configurations with lists of" << std::endl
    << "x,y screen coordinates in meters followed by long,lat angles in" << std::endl
    << "degrees where (0,0) is straight ahead from the eye, positive" << std::endl
    << "longitude is left and positive latitude is up." << std::endl
    << "  A single color from standard input is the default, input files" << std::endl
    << "can be optionally specified." << std::endl
    << "  It produces on standard output a partial OSVR display configuration file." << std::endl
    << std::endl;
  exit(1);
}

int main(int argc, char *argv[])
{
  // Parse the command line
  std::vector<std::string> inputFileNames;
  bool useRightEye = true;
  bool computeBounds = true;
  bool useFieldAngles = true;
  bool verifyAngles = false;
  double xx, xy, yx, yy;
  double maxAngleDiffDegrees;
  double left, right, bottom, top;
  double depth = 2.0;
  double toMeters = 1.0;
  int realParams = 0;
  for (int i = 1; i < argc; i++) {
    if (std::string("-mm") == argv[i]) {
      toMeters = 1e-3;  // Convert input in millimeters to meters
    } else if (std::string("-verbose") == argv[i]) {
      g_verbose = true;
    } else if (std::string("-latlong") == argv[i]) {
      useFieldAngles = false;
    } else if (std::string("-depth_meters") == argv[i]) {
      if (++i >= argc) { Usage(argv[0]); }
      depth = atof(argv[i]);
    } else if (std::string("-mono") == argv[i]) {
      if (++i >= argc) { Usage(argv[0]); }
      inputFileNames.push_back(argv[i]);
    } else if (std::string("-rgb") == argv[i]) {
      if (++i >= argc) { Usage(argv[0]); }
      inputFileNames.push_back(argv[i]);
      if (++i >= argc) { Usage(argv[0]); }
      inputFileNames.push_back(argv[i]);
      if (++i >= argc) { Usage(argv[0]); }
      inputFileNames.push_back(argv[i]);
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
    } else if (std::string("-eye") == argv[i]) {
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
    } else if (std::string("-verify_angles") == argv[i]) {
      verifyAngles = true;
      if (++i >= argc) { Usage(argv[0]); }
      xx = atof(argv[i]);
      if (++i >= argc) { Usage(argv[0]); }
      xy = atof(argv[i]);
      if (++i >= argc) { Usage(argv[0]); }
      yx = atof(argv[i]);
      if (++i >= argc) { Usage(argv[0]); }
      yy = atof(argv[i]);
      if (++i >= argc) { Usage(argv[0]); }
      maxAngleDiffDegrees = atof(argv[i]);
    }
    else if ((argv[i][0] == '-') && (atof(argv[i]) == 0.0)) {
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
  // The output screens and meshes.  There is one mesh per color, so one
  // for mono and three for RGB.
  ScreenDescription leftScreen, rightScreen;
  std::vector<MeshDescription> leftMeshes, rightMeshes;

  //====================================================================
  // Parse the angle-configuration information from standard or from the set
  // of input files specified.  Expect white-space separation between numbers
  // and also between entries (which may be on separate lines).
  std::vector<std::istream *> inFiles;
  if (inputFileNames.size() == 0) {
    inFiles.push_back(&std::cin);
    inputFileNames.push_back("standard input");
  } else {
    for (size_t i = 0; i < inputFileNames.size(); i++) {
      if (g_verbose) {
        std::cerr << "Opening file " << inputFileNames[i] << std::endl;
      }
      std::ifstream *inFile = new std::ifstream;
      inFile->open(inputFileNames[i].c_str(), std::ifstream::in);
      if (!inFile->good()) {
        std::cerr << "Error: Could not open " << inputFileNames[i] << std::endl;
        return 1;
      }
      inFiles.push_back(inFile);
    }
  }
  std::vector< std::vector<Mapping> > mappings;
  for (size_t i = 0; i < inFiles.size(); i++) {
    std::vector<Mapping> mapping = read_from_infile(*inFiles[i]);
    if (g_verbose) {
      std::cerr << "Found " << mapping.size() << " points in "
        << inputFileNames[i] << std::endl;
    }
    if (mapping.size() == 0) {
      std::cerr << "Error: No input points found in " << inputFileNames[i]
        << std::endl;
      return 2;
    }
    mappings.push_back(mapping);
  }

  //====================================================================
  // If we've been asked to verify the angles on the meshes, do so now.
  // This makes sure that the direction between neighbors in angle space
  // is consistent with their direction in screen space, removing points
  // that don't satisfy the criterion.  This removes inconsistent points
  // from the simulation (caused by multiple ray bounces or other
  // singularities in the simulation).
  if (verifyAngles) {
    for (size_t m = 0; m < mappings.size(); m++) {
      int ret = remove_invalid_points_based_on_angle(
        mappings[m], xx, xy, yx, yy, maxAngleDiffDegrees);
      if (ret < 0) {
        std::cerr << "Error verifying angles for mesh "
          << m << std::endl;
        return 60;
      }
      if (g_verbose) {
        std::cerr << "Removed " << ret
          << " points from mesh " << m << std::endl;
      }
    }
  }

  //====================================================================
  // If we've been asked to auto-range the screen coordinates, compute
  // them here.  Look at all of the points from all of the colors and
  // make a bound on all of them.
  if (computeBounds) {
    left = right = mappings[0][0].xyLatLong.x;
    bottom = top = mappings[0][0].xyLatLong.y;
    for (size_t m = 0; m < mappings.size(); m++) {
      for (size_t i = 1; i < mappings[m].size(); i++) {
        double x = mappings[m][i].xyLatLong.x;
        double y = mappings[m][i].xyLatLong.y;
        if (x < left) { left = x; }
        if (x > right) { right = x; }
        if (y < bottom) { bottom = y; }
        if (y > top) { top = y; }
      }
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
  // Compute left and right screen boundaries that are mirror images
  // of each other.
  double leftScreenLeft, leftScreenRight, leftScreenBottom, leftScreenTop;
  double rightScreenLeft, rightScreenRight, rightScreenBottom, rightScreenTop;
  rightScreenBottom = leftScreenBottom = bottom;
  rightScreenTop = leftScreenTop = top;
  if (useRightEye) {
    rightScreenLeft = left;
    rightScreenRight = right;

    leftScreenLeft = -right;
    leftScreenRight = -left;
  } else {
    leftScreenLeft = left;
    leftScreenRight = right;

    rightScreenLeft = -right;
    rightScreenBottom = -left;
  }

  //====================================================================
  // Compute a left- and right-eye mappings that are mirrors of each
  // other, so that we can produce distortion maps for both eyes.
  std::vector< std::vector<Mapping> > leftMappings;
  std::vector< std::vector<Mapping> > rightMappings;
  for (size_t i = 0; i < mappings.size(); i++) {
    // Do each mapping in turn, one per color.
    std::vector<Mapping> &mapping = mappings[i];

    //====================================================================
    // Make an inverse mapping for the opposite eye.  Invert around X in
    // angle and viewing direction.  Depending on whether we are using the
    // left or right eye, set the eyes appropriately.
    //  Also make a different set of screen boundaries for each, again
    // inverting around X = 0.
    std::vector<Mapping> leftMapping;
    std::vector<Mapping> rightMapping;
    rightScreenBottom = leftScreenBottom = bottom;
    rightScreenTop = leftScreenTop = top;
    if (useRightEye) {
      rightMapping = mapping;
      leftMapping = reflect_mapping(mapping);
    } else {
      leftMapping = mapping;
      rightMapping = reflect_mapping(mapping);
    }

    //====================================================================
    // Convert the input values into normalized coordinates and into 3D
    // locations.
    convert_to_normalized_and_meters(leftMapping, toMeters, depth,
      leftScreenLeft, leftScreenBottom, leftScreenRight, leftScreenTop,
      useFieldAngles);
    convert_to_normalized_and_meters(rightMapping, toMeters, depth,
      rightScreenLeft, rightScreenBottom, rightScreenRight, rightScreenTop,
      useFieldAngles);

    //====================================================================
    // Store the mappings.
    leftMappings.push_back(leftMapping);
    rightMappings.push_back(rightMapping);
  }

  //====================================================================
  // Construct mappings that include all points for all colors that
  // we will use to determine the screen boundaries in a manner that
  // encompasses all of them.
  std::vector<Mapping> leftFullMapping, rightFullMapping;
  for (size_t i = 0; i < mappings.size(); i++) {
    for (size_t j = 0; j < leftMappings[i].size(); j++) {
      leftFullMapping.push_back(leftMappings[i][j]);
    }
    for (size_t j = 0; j < rightMappings[i].size(); j++) {
      rightFullMapping.push_back(rightMappings[i][j]);
    }
  }
  if (!findScreen(leftFullMapping, leftScreenLeft, leftScreenBottom,
    leftScreenRight, leftScreenTop, leftScreen, g_verbose)) {
    std::cerr << "Error: Could not find left screen" << std::endl;
    return 3;
  }
  if (!findScreen(rightFullMapping, rightScreenLeft, rightScreenBottom,
    rightScreenRight, rightScreenTop, rightScreen, g_verbose)) {
    std::cerr << "Error: Could not find right screen" << std::endl;
    return 5;
  }

  //====================================================================
  // Compute the three colored mappings based on the screen boundaries
  // we found above.
  for (size_t i = 0; i < mappings.size(); i++) {
    // Do each pair of mappings in turn, one per color.
    std::vector<Mapping> &leftMapping = leftMappings[i];
    std::vector<Mapping> &rightMapping = rightMappings[i];

    //====================================================================
    // Determine the screen description and distortion mesh based on the
    // input points and screen parameters.
    // This will re-compute the screen each time, but it will get the same
    // answer because we're using the same bounds for each of them.
    MeshDescription leftMesh, rightMesh;
    if (!findMesh(leftMapping, leftScreenLeft, leftScreenBottom,
      leftScreenRight, leftScreenTop, leftScreen, leftMesh, g_verbose)) {
      std::cerr << "Error: Could not find left mesh" << std::endl;
      return 30;
    }
    if (leftMesh.size() != leftMapping.size()) {
      std::cerr << "Error: Left mesh size " << leftMesh.size()
        << " does not match mapping size" << leftMapping.size() << std::endl;
      return 4;
    }
    leftMeshes.push_back(leftMesh);

    if (!findMesh(rightMapping, rightScreenLeft, rightScreenBottom,
      rightScreenRight, rightScreenTop, rightScreen, rightMesh, g_verbose)) {
      std::cerr << "Error: Could not find right mesh" << std::endl;
      return 50;
    }
    if (rightMesh.size() != rightMapping.size()) {
      std::cerr << "Error: Right mesh size " << rightMesh.size()
        << " does not match mapping size" << rightMapping.size() << std::endl;
      return 6;
    }
    rightMeshes.push_back(rightMesh);
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
  switch (leftMeshes.size()) {
  case 1:
    std::cout << "    \"type\": \"mono_point_samples\"," << std::endl;
    std::cout << "    \"mono_point_samples\": [" << std::endl;
    writeMesh(std::cout, leftMeshes[0]);
    std::cout << "," << std::endl;
    writeMesh(std::cout, rightMeshes[0]);
    std::cout << "    ]" << std::endl; // mono_point_samples
    std::cout << "   }," << std::endl; // distortion
    break;
  case 3:
    std::cout << "    \"type\": \"rgb_point_samples\"," << std::endl;
    std::cout << "    \"red_point_samples\": [" << std::endl;
      writeMesh(std::cout, leftMeshes[0]);
      std::cout << "," << std::endl;
      writeMesh(std::cout, rightMeshes[0]);
    std::cout << "    ]," << std::endl; // red_point_samples
    std::cout << "    \"green_point_samples\": [" << std::endl;
      writeMesh(std::cout, leftMeshes[1]);
      std::cout << "," << std::endl;
      writeMesh(std::cout, rightMeshes[1]);
    std::cout << "    ]," << std::endl; // green_point_samples
    std::cout << "    \"blue_point_samples\": [" << std::endl;
      writeMesh(std::cout, leftMeshes[2]);
      std::cout << "," << std::endl;
      writeMesh(std::cout, rightMeshes[2]);
    std::cout << "    ]" << std::endl; // blue_point_samples
    std::cout << "   }," << std::endl; // distortion
    break;
  default:
    std::cerr << "Error: Unexpected number of meshes: " << leftMeshes.size()
      << std::endl;
    return 3;
  }

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

  // The latitude to use is less than 45 degrees at the corners of the
  // cube.  Wolfram's online calculator reports it as the arcsine of
  // the square root of 2/3, but that is the angle down from Pi/2, so
  // we convert it here to that coordinate system and to degrees.
  double phi = 90 - (180 / MY_PI * asin(sqrt(2.0 / 3.0)));

  // The longitude of the corners of the cube is 45 degrees.
  double theta = 45;

  // The distance to the corners of the unit cube is the square root
  // of the sum of the squares of the distance.
  double depth = sqrt(3);

  // Construct a set of points that should make a square screen with no distortion
  // that has a 90 degree horizontal and vertical field of view and the identity
  // distortion correction and that is oriented straight ahead.
  Mapping p1(XYLatLong(0, 0, -phi, -theta), XYZ(-1, -1, -1));
  Mapping p2(XYLatLong(1, 0, -phi,  theta), XYZ( 1, -1, -1));
  Mapping p3(XYLatLong(1, 1,  phi,  theta), XYZ( 1,  1, -1));
  Mapping p4(XYLatLong(0, 1,  phi, -theta), XYZ(-1,  1, -1));
  std::vector<Mapping> mapping;
  mapping.push_back(p1);
  mapping.push_back(p2);
  mapping.push_back(p3);
  mapping.push_back(p4);

  // Find the screen associated with this mapping.
  if (!convert_to_normalized_and_meters(mapping, 1, depth,
    0, 0, 1, 1)) {
    std::cerr << "testAlgorithms(): Could not normalize points" << std::endl;
    return 100;
  }
  ScreenDescription screen;
  MeshDescription mesh;
  if (!findScreen(mapping, 0, 0, 1, 1, screen, g_verbose)) {
    std::cerr << "testAlgorithms(): Could not find screen" << std::endl;
    return 101;
  }
  if (!findMesh(mapping, 0, 0, 1, 1, screen, mesh, g_verbose)) {
    std::cerr << "testAlgorithms(): Could not find mesh" << std::endl;
    return 102;
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
      std::cerr << "testAlgorithms(): X normalized mesh mismatch for element: " << entry
        << " (got " << mesh[entry][outIndex][0]
        << ", expected " << mapping[entry].xyLatLong.x << ")" << std::endl;
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
  // In this case, the distance to the points is sqrt(2).
  Mapping rp1(XYLatLong(0, 0, 0, -theta), XYZ(0, -1, -sqrt(2)));
  Mapping rp2(XYLatLong(1, 0, 90, -theta), XYZ(sqrt(2), -1, 0));
  Mapping rp3(XYLatLong(1, 1, 90, theta), XYZ(sqrt(2), 1, 0));
  Mapping rp4(XYLatLong(0, 1, 0, theta), XYZ(0, 1, -sqrt(2)));
  std::vector<Mapping> rmapping;
  rmapping.push_back(rp1);
  rmapping.push_back(rp2);
  rmapping.push_back(rp3);
  rmapping.push_back(rp4);

  // Find the screen associated with this mapping.
  ScreenDescription rscreen;
  MeshDescription rmesh;
  if (!findScreen(rmapping, 0, 0, 1, 1, rscreen, g_verbose)) {
    std::cerr << "testAlgorithms(): Could not find rotated screen" << std::endl;
    return 600;
  }
  if (!findMesh(rmapping, 0, 0, 1, 1, rscreen, rmesh, g_verbose)) {
    std::cerr << "testAlgorithms(): Could not find rotated mesh" << std::endl;
    return 601;
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

  // Construct a set of points that should make a square screen with distortion
  // that has a 90 degree horizontal and vertical field of view 
  // that is oriented straight ahead.  Make the points not cover the entire
  // screen so that we can test the impact of the oversize parameters.

  // The expected results, based on a screen that goes from 0-1 in X and Y.
  // The mapping goes from physical display-normalized coordinates to
  // normalized coordinates in the canonical display.
  // This places the normalized output texture coordinates in the range from
  // 0.1 to 0.9 and the normalized input texture coordinates into the range
  // 0 to 1 (covering the whole virtual viewport).  The center point has an
  // input coordinate in the center of the input space but not in the center
  // of the output space.
  std::vector<double> dExpectedXIn, dExpectedYIn;
  std::vector<double> dExpectedXOut, dExpectedYOut;

  Mapping dp1(XYLatLong(0.1, 0.1, -phi, -theta), XYZ(-1, -1, -1));
  dExpectedXIn.push_back(0.1); dExpectedYIn.push_back(0.1);
  dExpectedXOut.push_back(0.0); dExpectedYOut.push_back(0.0);
  Mapping dp2(XYLatLong(0.9, 0.1, -phi, theta), XYZ(1, -1, -1));
  dExpectedXIn.push_back(0.9); dExpectedYIn.push_back(0.1);
  dExpectedXOut.push_back(1.0); dExpectedYOut.push_back(0.0);
  Mapping dp3(XYLatLong(0.9, 0.9, phi, theta), XYZ(1, 1, -1));
  dExpectedXIn.push_back(0.9); dExpectedYIn.push_back(0.9);
  dExpectedXOut.push_back(1.0); dExpectedYOut.push_back(1.0);
  Mapping dp4(XYLatLong(0.1, 0.9, phi, -theta), XYZ(-1, 1, -1));
  dExpectedXIn.push_back(0.1); dExpectedYIn.push_back(0.9);
  dExpectedXOut.push_back(0.0); dExpectedYOut.push_back(1.0);
  Mapping dp5(XYLatLong(0.2, 0.3, 0,  0), XYZ(0, 0, -depth));
  dExpectedXIn.push_back(0.2); dExpectedYIn.push_back(0.3);
  dExpectedXOut.push_back(0.5); dExpectedYOut.push_back(0.5);
  std::vector<Mapping> dmapping;
  dmapping.push_back(dp1);
  dmapping.push_back(dp2);
  dmapping.push_back(dp3);
  dmapping.push_back(dp4);
  dmapping.push_back(dp5);

  // Normalize the mapping and make sure the results are as expected.
  std::vector<Mapping> ndmapping = dmapping;
  if (!convert_to_normalized_and_meters(ndmapping, 1, depth,
      0, 0, 1, 1)) {
    std::cerr << "testAlgorithms(): Could not normalize distorted points" << std::endl;
    return 1100;
  }
  for (int entry = 0; entry < ndmapping.size(); entry++) {
    size_t outIndex = 1;

    if (!small(ndmapping[entry].xyz.x - dmapping[entry].xyz.x)) {
      std::cerr << "testAlgorithms(): Distorted X coord mismatch for element: " << entry
      << " (found " << ndmapping[entry].xyz.x
      << ", expected " << dmapping[entry].xyz.x << ")"
      << std::endl;
      return 1100 + 3*entry;
    }
    if (!small(ndmapping[entry].xyz.y - dmapping[entry].xyz.y)) {
      std::cerr << "testAlgorithms(): Distorted Y coord mismatch for element: " << entry
      << " (found " << ndmapping[entry].xyz.y
      << ", expected " << dmapping[entry].xyz.y << ")"
      << std::endl;
      return 1101 + 3*entry;
    }
    if (!small(ndmapping[entry].xyz.z - dmapping[entry].xyz.z)) {
      std::cerr << "testAlgorithms(): Distorted Z coord mismatch for element: " << entry
      << " (found " << ndmapping[entry].xyz.z
      << ", expected " << dmapping[entry].xyz.z << ")"
      << std::endl;
      return 1102 + 3*entry;
    }

    if (!small(ndmapping[entry].xyLatLong.x - dExpectedXIn[entry])) {
      std::cerr << "testAlgorithms(): Distorted X coord mismatch for element: " << entry
        << " (found " << ndmapping[entry].xyLatLong.x
        << ", expected " << dExpectedXIn[entry] << ")"
        << std::endl;
      return 1103 + 3*entry;
    }
    if (!small(ndmapping[entry].xyLatLong.y - dExpectedYIn[entry])) {
      std::cerr << "testAlgorithms(): Distorted Y coord mismatch for element: " << entry
        << " (found " << ndmapping[entry].xyLatLong.y
        << ", expected " << dExpectedYIn[entry] << ")"
        << std::endl;
      return 1104 + 3 * entry;
    }
  }

  // Find the screen associated with this mapping.
  ScreenDescription dscreen;
  MeshDescription dmesh;
  if (!findScreen(ndmapping, 0, 0, 1, 1, dscreen, g_verbose)) {
    std::cerr << "testAlgorithms(): Could not find distorted screen" << std::endl;
    return 1200;
  }
  if (!findMesh(ndmapping, 0, 0, 1, 1, dscreen, dmesh, g_verbose)) {
    std::cerr << "testAlgorithms(): Could not find distorted mesh" << std::endl;
    return 1201;
  }

  // Make sure the screen has the expected behavior.
  if (!small(dscreen.hFOVDegrees - 90)) {
    std::cerr << "testAlgorithms(): Distorted hFOV not near 90: " << dscreen.hFOVDegrees << std::endl;
    return 1201;
  }
  if (!small(dscreen.vFOVDegrees - 90)) {
    std::cerr << "testAlgorithms(): Distorted vFOV not near 90: " << dscreen.vFOVDegrees << std::endl;
    return 1202;
  }
  if (!small(dscreen.overlapPercent - 100)) {
    std::cerr << "testAlgorithms(): Distorted overlap percent not near 100: " << dscreen.overlapPercent << std::endl;
    return 1203;
  }
  if (!small(dscreen.xCOP - 0.5)) {
    std::cerr << "testAlgorithms(): Distorted xCOP not near 0.5: " << dscreen.xCOP << std::endl;
    return 1204;
  }
  if (!small(dscreen.yCOP - 0.5)) {
    std::cerr << "testAlgorithms(): Distorted yCOP not near 0.5: " << dscreen.yCOP << std::endl;
    return 1205;
  }

  // Make sure the mesh has the expected behavior.
  if (dmesh.size() != dmapping.size()) {
    std::cerr << "testAlgorithms(): Destorted mesh size does not match mapping size: " << mesh.size() << std::endl;
    return 1301;
  }
// @todo Why do we get coordinates (0,1) for the fifth element, when we expect (1/3, 2/3)?
  for (int entry = 0; entry < dmesh.size(); entry++) {
    size_t outIndex = 1;
    if (!small(dmesh[entry][outIndex][0] - dExpectedXOut[entry])) {
      std::cerr << "testAlgorithms(): Distorted X normalized mesh mismatch for element: " << entry
      << " (found " << dmesh[entry][outIndex][0]
      << ", expected " << dExpectedXOut[entry] << ")"
      << std::endl;
      return 1400 + entry;
    }
    if (!small(dmesh[entry][outIndex][1] - dExpectedYOut[entry])) {
      std::cerr << "testAlgorithms(): Distorted Y normalized mesh mismatch for element: " << entry
        << " (found " << dmesh[entry][outIndex][1]
        << ", expected " << dExpectedYOut[entry] << ")" << std::endl;
      return 1500 + entry;
    }
  }

  // @todo Insert a test for field angles in normalization.

  if (g_verbose) {
    std::cerr << "=================== Successfully finished testAlgorithms()" << std::endl;
  }
  return 0;
}
