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

#pragma once

// Internal Includes

// Standard includes
#include <string>
#include <cmath>
#include <array>
#include <vector>

// Global constants and variables
#define MY_PI (4.0*atan(1.0))

// Screen-space to/from angle-space map entry
class XYLatLong {
public:
  double x;
  double y;
  double latitude;
  double longitude;

  XYLatLong(double px, double py, double plat, double plong)
  {
    x = px; y = py; latitude = plat; longitude = plong;
  }
  XYLatLong() { x = y = latitude = longitude = 0; }
};

// 3D coordinate
class XYZ {
public:
  double x;
  double y;
  double z;

  XYZ(double px, double py, double pz)
  {
    x = px; y = py; z = pz;
  }
  XYZ() { x = y = z = 0; }

  /// Return the rotation about the Y axis, where 0 rotation points along
  // the -Z axis and positive rotation heads towards the -X axis.
  // The X axis in atan space corresponds to the -z axis in head space,
  // and the Y axis in atan space corresponds to the -x axis in head space.
  double rotationAboutY() const {
    return atan2(-x, -z);
  }

  /// Project from the origin through our point onto a plane whose
  // equation is specified.
  XYZ projectOntoPlane(double A, double B, double C, double D) const {
    XYZ ret;

    // Solve for the value of S that satisfies:
    //    Asx + Bsy + Csz + D = 0,
    //    s = -D / (Ax + By + Cz)
    // Then for the location sx, sy, sz.

    double s = -D / (A*x + B*y + C*z);
    ret.x = s*x;
    ret.y = s*y;
    ret.z = s*z;

    return ret;
  }

  /// Return the rotation distance from another point.
  double distanceFrom(const XYZ &p) const {
    return sqrt( (x-p.x)*(x-p.x) + (y-p.y)*(y-p.y) + (z-p.z)*(z-p.z) );
  }
};

// Mapping entry, along with its associated 3D coordinate
class Mapping{
public:
  XYLatLong xyLatLong;
  XYZ xyz;

  Mapping(XYLatLong const &ll, XYZ const &x)
  {
    xyLatLong = ll;
    xyz = x;
  }
  Mapping() {};
};

// Description of a screen
typedef struct {
  double hFOVDegrees;
  double vFOVDegrees;
  double overlapPercent;
  double xCOP;
  double yCOP;

  // These are quantities computed along the way to getting the
  // screen that are needed by the mesh calculations, so they
  // are stored in the screen to pass from the findScreen to
  // the findMesh functions.
  double A, B, C, D;  //!< Ax + By + Cz + D = 0 screen plane
  XYZ screenLeft, screenRight;  //!< Left-most and right-most points on screen
  double maxY;  //!< Maximum absolute value of Y for points on screen
} ScreenDescription;

/// Holds a list of mappings from physical-display normalized
/// coordinates to canonical-display normalized coordinates.
typedef std::vector<        //!< Vector of mappings
    std::array<            //!< 2-vector of from, to coordinates
      std::array<double,2>   //!< 2-vector of unit coordinates (x,y)
    ,2 >
> MeshDescription;
