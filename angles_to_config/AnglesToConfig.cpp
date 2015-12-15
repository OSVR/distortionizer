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
#include <json/json.h>

// Standard includes
#include <string>
#include <iostream>
#include <cmath>
#include <stdlib.h> // For exit()

#define MY_PI (4.0*atan(1.0))

// Screen-space to/from angle-space map entry
typedef struct {
  double x;
  double y;
  double latitude;
  double longitude;
} XYLatLong;

// 3D coordinate
typedef struct {
  double x;
  double y;
  double z;

  /// Return the rotation about the Y axis, where 0 rotation points along
  // the -Z axis and positive rotation heads towards the -X axis.
  double rotationAboutY() const {
    return atan2(-x, -z);
  }
} XYZ;

// Mapping entry, along with its associated 3D coordinate
typedef struct {
  XYLatLong xyLatLong;
  XYZ xyz;
} Mapping;


void Usage(std::string name)
{
  std::cerr << "Usage: " << name
    << " [-eye right|left] (default is right)"
    << " [-depth_meters D] (default is 2.0)"
    << " [-mm] (default is meters)"
    << " screen_left_meters screen_bottom_meters screen_right_meters screen_top_meters"
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
  bool rightEye = true;
  double left, right, bottom, top;
  double depth = 2.0;
  double toMeters = 1.0;
  int realParams = 0;
  for (int i = 1; i < argc; i++) {
    if (std::string("-mm") == argv[i]) {
      toMeters = 1e-3;  // Convert input in millimeters to meters
    } else if (std::string("-depth_meters") == argv[i]) {
      if (++i >= argc) { Usage(argv[0]); }
      depth = atof(argv[i]);
    } else if (std::string("-eye") == argv[i]) {
      if (++i >= argc) { Usage(argv[0]); }
      std::string eye = argv[i];
      if (eye == "left") {
        rightEye = false;
      } else if (eye == "right") {
        rightEye = true;
      } else {
        std::cerr << "Bad value for -eye: " << eye << ", expected left or right" << std::endl;
        Usage(argv[0]);
      }
    } else if ((argv[i][0] == '-') && (atof(argv[i]) == 0.0) ) {
      Usage(argv[0]);
    }
    else switch (++realParams) {
    case 1:
      left = atof(argv[i]);
      break;
    case 2:
      bottom = atof(argv[i]);
      break;
    case 3:
      right = atof(argv[i]);
      break;
    case 4:
      top = atof(argv[i]);
      break;
    default:
      Usage(argv[0]);
    }
  }
  if (realParams != 4) { Usage(argv[0]); }

  // Parse the angle-configuration information from standard input.  Expect white-space
  // separation between numbers and also between entries (which may be on separate
  // lines).
  std::vector<Mapping> mapping;
  while (!std::cin.eof()) {
    // Read the mapping info from the input file.
    Mapping map;
    std::cin >> map.xyLatLong.x >> map.xyLatLong.y >> map.xyLatLong.latitude >> map.xyLatLong.longitude;

    //  Convert the input coordinates from its input space into meters
    // and then convert (using the screen dimensions) into normalized screen units.
    map.xyLatLong.x *= toMeters;
    map.xyLatLong.x = (map.xyLatLong.x - left) / (right - left);
    map.xyLatLong.y *= toMeters;
    map.xyLatLong.y = (map.xyLatLong.y - bottom) / (top - bottom);

    // Convert the input latitude and longitude from degrees to radians.
    map.xyLatLong.latitude *= MY_PI / 180;
    map.xyLatLong.longitude *= MY_PI / 180;

    // Compute the 3D coordinate of the point w.r.t. the eye at the origin.
    // longitude = 0, lattitude = 0 points along the -Z axis in eye space.
    map.xyz.y = depth * sin(map.xyLatLong.latitude);
    map.xyz.z = -depth * cos(map.xyLatLong.longitude) * cos(map.xyLatLong.latitude);
    map.xyz.x = -depth * sin(map.xyLatLong.longitude) * cos(map.xyLatLong.latitude);

    mapping.push_back(map);
  }
  if (mapping.size() == 0) {
    std::cerr << "Error: No input points found" << std::endl;
    return 2;
  }

  // Figure out the X screen-space extents.
  // The X screen-space extents are defined by the lines perpendicular to the Y axis passing through:
  //  left: the point location whose reprojection into the Y = 0 plane has the most -
  //        positive angle(note that this may not be the point with the largest
  //        longitudinal coordinate, because of the impact of changing latitude on
  //        X - Z position).
  //  right : the point location whose reprojection into the Y = 0 plane has the most -
  //        negative angle(note that this may not be the point with the smallest
  //        longitudinal coordinate, because of the impact of changing latitude on
  //        X - Z position).
  XYZ screenLeft, screenRight;
  screenLeft = screenRight = mapping[0].xyz;
  for (size_t i = 0; i < mapping.size(); i++) {
    if (mapping[i].xyz.rotationAboutY() > screenLeft.rotationAboutY()) {
      screenLeft = mapping[i].xyz;
    }
    if (mapping[i].xyz.rotationAboutY() < screenRight.rotationAboutY()) {
      screenRight = mapping[i].xyz;
    }
  }
  if (screenLeft.rotationAboutY() - screenRight.rotationAboutY() >= MY_PI) {
    std::cerr << "Error: Field of view > 180 degrees: found " <<
      180 / MY_PI * (screenLeft.rotationAboutY() - screenRight.rotationAboutY())
      << std::endl;
    return 3;
  }

  // Figure out the Y screen-space extents.
  // The Y screen-space extents are symmetric and correspond to the lines parallel
  //  to the screen X axis that are within the plane of the X line specifying the
  //  axis extents at the largest magnitude angle up or down from the horizontal.
  XYZ screenTop, screenBottom;
  screenTop = screenBottom = mapping[0].xyz;

  return 0;
}

