/** @file
    @brief Produces a distortion polynomial and partial display description
           from a table of display locations to angle inputs.

    @date 2019

    @author
    Russ Taylor working through ReliaSolve.com for ValitXR.
*/

// Copyright 2019 ValityXR.
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

// Standard includes
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdlib.h> // For exit()

// Local includes
#include "types.h"

// Global constants and variables
static bool g_verbose = false;
static double MY_PI = (4.0*atan(1.0));

// Forward declarations
static int testAlgorithms();

static void writePolynomial(std::ostream &s, std::vector<double> const &poly)
{
  s << "[";
  for (size_t i = 0; i < poly.size(); i++) {
    if (i == 0) { s << " "; }
    else { s << ","; }
    s << std::setprecision(4) << poly[i];
  }
  std::cout << " ]";
}

/// @brief Evaluate the coefficients with both radius and vector in meters.
/// @return Angle in Radians associated with radial point.
static double evaluateAngleRadians(double r, std::vector<double> const &coefs)
{
  double ret = 0;
  for (size_t i = 0; i < coefs.size(); i++) {
    ret += coefs[i] * pow(r, 1 + 2 * i);
  }
  return atan(ret);
}

/// @brief Determine the OSVR equivalent FOV for an axis.
/// @param left Left side of the display (or Bottom) (at depth)
/// @param right Right side of the display (or Top) (at depth)
/// @param depth Perpendicular distance to the screen
static double computeFOVDegrees(double left, double right, double depth)
{
  // Find half of the width of the axis, because we need to compute
  // the FOV as if we were centered.
  double halfWidth = fabs((right - left) / 2);

  // Find twice the angle to that point at the specified depth
  double rad = 2 * atan2(halfWidth, depth);

  // Convert to Radians
  return 180 / MY_PI * rad;
}

void Usage(std::string name)
{
  std::cerr << "Usage: " << name
    << " [-eye right|left] (default is right)"
    << " [-depth_meters D] (default is 2.0)"
    << " [-verbose] (default is not)"
    << " display_size"
    << " C1 [C3 [C5 [...]]"
    << std::endl
    << "  This program solves for the polynomial distortion correction that will correct for" << std::endl
    << "the screen-to-angle mapping described in the coefficients.  The first-order term C1" << std::endl
    << "must be specified, and other higher odd-numbered terms may also be specified (C3, C5, ...)" << std::endl
    << "The coefficents are from the equation tan(theta) = C1*r + C3*r^3 + C5*r^5 + ..." << std::endl
    << "  The input r for the mapping is in millimeters." << std::endl
    << "  The output of the mapping is the tangent of the angle towards the visible point" << std::endl
    << "  The size of the display is specified in millimeters and must be the." << std::endl
    << "same for both width and height until the code is generalized." << std::endl
    << std::endl;
  exit(1);
}

int main(int argc, char *argv[])
{
  // Parse the command line
  std::vector<std::string> inputFileNames;
  bool useRightEye = true;
  double left, right, bottom, top;
  std::vector<double> coeffs;
  double depth = 2.0;
  int realParams = 0;
  for (int i = 1; i < argc; i++) {
    if (std::string("-verbose") == argv[i]) {
      g_verbose = true;
    } else if (std::string("-depth_meters") == argv[i]) {
      if (++i >= argc) { Usage(argv[0]); }
      depth = atof(argv[i]);
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
    } else if ((argv[i][0] == '-') && (atof(argv[i]) == 0.0)) {
      Usage(argv[0]);
    } else switch (++realParams) {
    case 1:
      right = top = atof(argv[i])/2;
      left = bottom = -right;
      break;
    default:
      // We got another coefficient, so put it onto the list after scaling it.
      coeffs.push_back(atof(argv[i]));
    }
  }
  if (realParams < 5) { Usage(argv[0]); }

  //====================================================================
  // Ensure that our parameters are not more general than we can handle.

  //====================================================================
  // Run our algorithm test to make sure things are working properly.
  /* @todo Put this back in when the test is implemented
  int ret;
  if ((ret = testAlgorithms()) != 0) {
    std::cerr << "Error testing basic algorithms, code " << ret << std::endl;
    return 100;
  }
  */

  //====================================================================
  // The output screens and polynomial.
  ScreenDescription leftScreen, rightScreen;
  std::vector<double> polynomial;

  //====================================================================
  // Compute left and right screen boundaries that are mirror images
  // of each other.
  // Compute the boundaries of the screen based on the depth and the maximum
  // angles in +/- X and Y projected to that distance.  These need to be
  // projected from the display space into the distance that the screen
  // lies at.
  double leftAngle = evaluateAngleRadians(left, coeffs);
  double rightAngle = evaluateAngleRadians(right, coeffs);
  double topAngle = evaluateAngleRadians(top, coeffs);
  double bottomAngle = evaluateAngleRadians(bottom, coeffs);
  if (g_verbose) {
    std::cout << "leftAngle (degrees): " << leftAngle * 180 / MY_PI << std::endl;
    std::cout << "rightAngle (degrees): " << rightAngle * 180 / MY_PI << std::endl;
    std::cout << "bottomAngle (degrees): " << bottomAngle * 180 / MY_PI << std::endl;
    std::cout << "topAngle (degrees): " << topAngle * 180 / MY_PI << std::endl;
  }
  double leftScreenLeft, leftScreenRight, leftScreenBottom, leftScreenTop;
  double rightScreenLeft, rightScreenRight, rightScreenBottom, rightScreenTop;
  rightScreenBottom = leftScreenBottom = depth * tan(bottomAngle);
  rightScreenTop = leftScreenTop = depth * tan(topAngle);
  if (useRightEye) {
    rightScreenLeft = depth * tan(leftAngle);
    rightScreenRight = depth * tan(rightAngle);

    leftScreenLeft = -depth * tan(rightAngle);
    leftScreenRight = -depth * tan(leftAngle);
  } else {
    leftScreenLeft = depth * tan(leftAngle);
    leftScreenRight = depth * tan(rightAngle);

    rightScreenLeft = -depth * tan(rightAngle);
    rightScreenRight = -depth * tan(leftAngle);
  }
  if (g_verbose) {
    std::cout << "rightScreenLeft: " << rightScreenLeft << std::endl;
    std::cout << "rightScreenRight: " << rightScreenRight << std::endl;
    std::cout << "rightScreenBottom: " << rightScreenBottom << std::endl;
    std::cout << "rightScreenTop: " << rightScreenTop << std::endl;
  }

  //====================================================================
  // Convert the screen boundaries into the OSVR description format.
  leftScreen.xCOP = (0 - leftScreenLeft) / (leftScreenRight - leftScreenLeft);
  leftScreen.yCOP = (0 - leftScreenBottom) / (leftScreenTop - leftScreenBottom);
  leftScreen.overlapPercent = 100; /// @todo Consider generalizing
  leftScreen.hFOVDegrees = computeFOVDegrees(leftScreenLeft, leftScreenRight, depth);
  leftScreen.vFOVDegrees = computeFOVDegrees(leftScreenBottom, leftScreenTop, depth);

  rightScreen.xCOP = (0 - rightScreenLeft) / (rightScreenRight - rightScreenLeft);
  rightScreen.yCOP = (0 - rightScreenBottom) / (rightScreenTop - rightScreenBottom);
  rightScreen.overlapPercent = 100; /// @todo Consider generalizing
  rightScreen.hFOVDegrees = computeFOVDegrees(rightScreenLeft, rightScreenRight, depth);
  rightScreen.vFOVDegrees = computeFOVDegrees(rightScreenBottom, rightScreenTop, depth);

  if (g_verbose) {
    std::cout << "hFOV (degrees): " << rightScreen.hFOVDegrees << std::endl;
    std::cout << "vFOV (degrees): " << rightScreen.vFOVDegrees << std::endl;
  }

  //====================================================================
  // Compute a polynomial that will map from the linear, rendered canonical
  // image to the correct location that matches that viewing angle when seen
  // through HMD lens.  It will leave points at the center of projection
  // at the same location and should map a point at the furthest location on
  // the same vertical or horizontal screen to the same location (if the
  // screen is centered, this will be points on the horizontal and vertical
  // centers of the edges).

  // The coefficients provide a mapping from radial distance in millimeters
  // away from the center of projection to the tangent of the angle at which
  // the point will appear.  OSVR wants a polynomial that moves points from
  // an initial location in a space to an offset location in that same space.
  //   The tangent space is scaled differently than the pixel space, so we
  // need to transform the coefficients so that they are not.
  //   We are free to choose either space (or another one entirely), so long
  // as the distance metrics match.
  //   We choose to scale the input space to match the output tangent space
  // in a manner that will map the furthest input pixel to the furthest tangent,
  // which means that we need to convert the input units to tangent space by
  // multiplying them by tangent_range/mm_range.
  //   This basically means that we need to scale the polynomial coefficients
  // by this inverse of this factor raised to the power that they are applying
  // so that they will do the conversion for us.
  double units = (left - right) / (tan(leftAngle) - tan(rightAngle));

  // This means that the coefficients are providing us the correct mapping,
  // but we need to convert them into the appropriate distance scale an put
  // them into the appropriate coefficients.  The appropriate coefficients
  // are the 1st, 3rd, and so forth.  We skip all even coefficients.
  // Always push back 0 for the 0th-order term.
  polynomial.push_back(0);
  for (size_t i = 0; i < coeffs.size(); i++) {
    // Zero all even terms.  We already zeroed the 0th-order term
    if (i > 0) { polynomial.push_back(0); }
    // Scale by the units conversion.
    polynomial.push_back(coeffs[i] * pow(units,1 + 2 * i));
  }
  // Set the distance units to twice the distance from the center of
  // the display to the farthest direction so that the point at the center
  // of the furthest edge will map to itself.  This is in the tangent
  // space.
  // (For the square case, this is just twice the right edge.)
  double distance_scale = 2*tan(rightAngle);

  /// @todo Generalize for non-square, non-centered displays.  In that
  // case, we want the farthest edge from the center to map to itself and
  // the other edges will overfill to some extent.

  // Determine the longest distance from the center of projection to an
  // edge in millimeters.
  /*
  double maxEdge = fabs(left);
  maxEdge = std::max(maxEdge, fabs(right));
  maxEdge = std::max(maxEdge, fabs(bottom));
  maxEdge = std::max(maxEdge, fabs(top));
  */

  // Report the scale factor needed to offset the pixel slope at the
  // center of the image (where it is presumably the largest) so that we
  // preserve the full resolution after scaling the texture by the
  // distortion map.
  if (g_verbose && polynomial.size() > 2) {
    // Derivative of the polynomial at 0 depends only on the first-order coefficient
    std::cout << "Display pixel slope at center: " << polynomial[1] << std::endl;
    if (polynomial[1] > 0) {
      std::cout << "renderManagerConfig:renderManagerConfig:renderOversampleFactor needed to compensate: " << 1/polynomial[1] << std::endl;
    }
  }

  //====================================================================
  // Construct Json screen description.
  // We do this by hand rather than using JsonCPP because we want
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
  std::cout << "     \"distance_scale_x\": " << distance_scale << "," << std::endl;
  std::cout << "     \"distance_scale_y\": " << distance_scale << "," << std::endl;
  std::cout << "     \"polynomial_coeffs_red\": ";
    writePolynomial(std::cout, polynomial);
  std::cout << "," << std::endl;
  std::cout << "     \"polynomial_coeffs_green\": ";
  writePolynomial(std::cout, polynomial);
  std::cout << "," << std::endl;
  std::cout << "     \"polynomial_coeffs_blue\": ";
  writePolynomial(std::cout, polynomial);
  std::cout << std::endl;
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

  // @todo
  std::cerr << "Warning: Test not yet implemented" << std::endl;

  if (g_verbose) {
    std::cerr << "=================== Successfully finished testAlgorithms()" << std::endl;
  }
  return 0;
}
