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
#include "types.h"
#include "helper.h"

// Standard includes
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

bool convert_to_normalized_and_meters(
  std::vector<Mapping> &mapping, double toMeters, double depth,
  double left, double bottom, double right, double top)
{
  for (size_t i = 0; i < mapping.size(); i++) {
    //  Convert the input coordinates from its input space into meters
    // and then convert (using the screen dimensions) into normalized screen units.
    mapping[i].xyLatLong.x *= toMeters;
    mapping[i].xyLatLong.x = (mapping[i].xyLatLong.x - left) / (right - left);
    mapping[i].xyLatLong.y *= toMeters;
    mapping[i].xyLatLong.y = (mapping[i].xyLatLong.y - bottom) / (top - bottom);

    // Convert the input latitude and longitude from degrees to radians.
    mapping[i].xyLatLong.latitude *= MY_PI / 180;
    mapping[i].xyLatLong.longitude *= MY_PI / 180;

    // Compute the 3D coordinate of the point w.r.t. the eye at the origin.
    // longitude = 0, lattitude = 0 points along the -Z axis in eye space.
    // Positive rotation in longitude is towards -X and positive rotation in
    // latitude points towards +Y.
    mapping[i].xyz.y = depth * sin(mapping[i].xyLatLong.latitude);
    mapping[i].xyz.z = -depth * cos(mapping[i].xyLatLong.longitude) * cos(mapping[i].xyLatLong.latitude);
    mapping[i].xyz.x = -depth * (-sin(mapping[i].xyLatLong.longitude)) * cos(mapping[i].xyLatLong.latitude);
  }

  // Make sure that the normalized screen coordinates are all within 0 and 1.
  for (size_t i = 0; i < mapping.size(); i++) {
    if ((mapping[i].xyLatLong.x < 0) || (mapping[i].xyLatLong.x > 1)) {
      std::cerr << "Error: Point " << i << " x out of range [0,1]: "
        << mapping[i].xyLatLong.x << " (increase bounds on -screen or don't specify it)"
        << std::endl;
    }
    if ((mapping[i].xyLatLong.y < 0) || (mapping[i].xyLatLong.y > 1)) {
      std::cerr << "Error: Point " << i << " y out of range [0,1]: "
        << mapping[i].xyLatLong.y << " (increase bounds on -screen or don't specify it)"
        << std::endl;
    }
  }

  return true;
}
