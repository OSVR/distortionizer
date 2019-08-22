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
#include <array>
#include <stdlib.h> // For exit()

// Local includes
#include "types.h"

// Global constants and variables
static bool g_verbose = false;
static double MY_PI = (4.0*atan(1.0));

// Forward declarations
static int testAlgorithms();

std::vector<std::string> split(const std::string& str, const std::string& delim)
{
  std::vector<std::string> tokens;
  size_t prev = 0, pos = 0;
  do {
    pos = str.find(delim, prev);
    if (pos == std::string::npos) { pos = str.length(); }
    std::string token = str.substr(prev, pos - prev);
    if (!token.empty()) { tokens.push_back(token); }
    prev = pos + delim.length();
  } while (pos < str.length() && prev < str.length());
  return tokens;
}

std::vector<double> parseColorLine(std::vector<std::string> const& words) {
  std::vector<double> ret;
  for (size_t i = 1; i < words.size(); i++) {
    ret.push_back(atof(words[i].c_str()));
  }
  return ret;
}

static void writePolynomial(std::ostream &s, std::vector<double> const &poly)
{
  s << "[";
  for (size_t i = 0; i < poly.size(); i++) {
    if (i == 0) { s << " "; }
    else { s << ","; }
    s << poly[i];
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
    << " [-verbose] (default is not)" << std::endl
    << " Then one of the following: " << std::endl
    << "  display_size C1 [C3 [C5 [...]]" << std::endl
    << "  config_file_name"
    << std::endl
    << "  This program solves for the polynomial distortion correction that will correct for" << std::endl
    << "the screen-to-angle mapping described in the coefficients.  The first-order term C1" << std::endl
    << "must be specified, and other higher odd-numbered terms may also be specified (C3, C5, ...)" << std::endl
    << "The coefficents are from the equation tan(theta) = C1*r + C3*r^3 + C5*r^5 + ..." << std::endl
    << "  The display_size is in millimeters, as is the space the coefficients map to." << std::endl
    << "  The output of the mapping is the tangent of the angle towards the visible point" << std::endl
    << "  The size of the display is specified in millimeters and must be the" << std::endl
    << "same for both width and height until the code is generalized." << std::endl
    << "  If the config file is specified, it has four lines.  The first has the display size in mm." << std::endl
    << "The second starts with the string 'leftred' (no apostrophes) and then has red coefficients" << std::endl
    << "for the left eye. The next two start with 'leftgreen' and 'leftblue' and have these coefficients." << std::endl
    << "The next three lines are 'rightred', 'rightgreen', and 'rightblue', with corresponding descriptors."
    << std::endl;
  exit(1);
}

int main(int argc, char *argv[])
{
  //====================================================================
  // Information about each color channel for each eye.
  class Color {
  public:
    std::vector<double> coeffs; //< Coefficients for the input polynomial
    std::vector<double> polynomial;   //< Output polynomial
  };
  std::array< std::array<Color, 3>, 2> eyeColors;

  // Parse the command line
  std::string inputFileName;
  double left, right, bottom, top;  //< Screen corners in coefficient space
  bool useRightEye = true;
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
      {
        // @todo Generalize to non-square displays
        double val = atof(argv[i]);
        if (val > 0) {
          right = top = val / 2;
          left = bottom = -right;
        } else {
          // They have specified a configuration file to be read.
          inputFileName = argv[i];
        }
      }
      break;
    default:
      // We got another coefficient, so put it onto the list after scaling it.
      if (inputFileName.size() > 0) {
        std::cerr << "Cannot add coefficients to command line when using config file" << std::endl;
        return 1;
      }
      // Monochrome; all colors are the same for both eyes.
      double val = atof(argv[i]);
      for (auto &eye : eyeColors) {
        for (auto &color : eye) {
          color.coeffs.push_back(val);
        }
      }
    }
  }
  if (realParams < 1) { Usage(argv[0]); }
  if ((realParams < 2) && (inputFileName.size() == 0)) { Usage(argv[0]); }

  //====================================================================
  // Parse the input configuration file if we have one specified.
  if (inputFileName.size() > 0) {
    std::ifstream s;
    s.open(inputFileName, std::ifstream::in | std::ifstream::binary);
    if (!s.good()) {
      std::cerr << "Could not open " << inputFileName << " for reading" << std::endl;
      return 2;
    }

    std::string line;
    std::vector<std::string> words;

    // Parse the display_size line
    // @todo Generalize to non-square displays
    std::getline(s, line);
    double val = atof(line.c_str());
    if (val <= 0) {
      std::cerr << "Bad display size line in " << inputFileName << ": " << line << std::endl;
      return 3;
    }
    right = top = val / 2;
    left = bottom = -right;

    // Read and parse the red, green, and blue lines for the left eye
    std::getline(s, line);
    words = split(line, " ");
    if (words[0] != "leftred") {
      std::cerr << "Bad leftred line in " << inputFileName << ": " << line << std::endl;
      return 4;
    }
    eyeColors[0][0].coeffs = parseColorLine(words);

    std::getline(s, line);
    words = split(line, " ");
    if (words[0] != "leftgreen") {
      std::cerr << "Bad leftgreen line in " << inputFileName << ": " << line << std::endl;
      return 5;
    }
    eyeColors[0][1].coeffs = parseColorLine(words);

    std::getline(s, line);
    words = split(line, " ");
    if (words[0] != "leftblue") {
      std::cerr << "Bad leftblue line in " << inputFileName << ": " << line << std::endl;
      return 6;
    }
    eyeColors[0][2].coeffs = parseColorLine(words);

    // Read and parse the red, green, and blue lines for the right eye
    std::getline(s, line);
    words = split(line, " ");
    if (words[0] != "rightred") {
      std::cerr << "Bad righttred line in " << inputFileName << ": " << line << std::endl;
      return 4;
    }
    eyeColors[1][0].coeffs = parseColorLine(words);

    std::getline(s, line);
    words = split(line, " ");
    if (words[0] != "rightgreen") {
      std::cerr << "Bad rightgreen line in " << inputFileName << ": " << line << std::endl;
      return 5;
    }
    eyeColors[1][1].coeffs = parseColorLine(words);

    std::getline(s, line);
    words = split(line, " ");
    if (words[0] != "rightblue") {
      std::cerr << "Bad rightblue line in " << inputFileName << ": " << line << std::endl;
      return 6;
    }
    eyeColors[1][2].coeffs = parseColorLine(words);
  }

  //====================================================================
  // Ensure that our parameters are not more general than we can handle.
  /// @todo

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
  // Compute left and right screen boundaries that are mirror images
  // of each other.
  // Compute the boundaries of the screen based on the depth and the maximum
  // angles in +/- X and Y projected to that distance.  These need to be
  // projected from the display space into the distance that the screen
  // lies at.
  ScreenDescription leftScreen, rightScreen;  //< Screen descriptions for each eye

  // Find the maximum angle for each direction across all colors in all eyes.
  double maxLeftAngle = 0, maxRightAngle = 0, maxTopAngle = 0, maxBottomAngle = 0;

  for (auto eye : eyeColors) {
    for (auto color : eye) {
      double leftAngle = evaluateAngleRadians(left, color.coeffs);
      if (fabs(leftAngle) > fabs(maxLeftAngle)) { maxLeftAngle = leftAngle; }
      double rightAngle = evaluateAngleRadians(right, color.coeffs);
      if (fabs(rightAngle) > fabs(maxRightAngle)) { maxRightAngle = rightAngle; }
      double topAngle = evaluateAngleRadians(top, color.coeffs);
      if (fabs(topAngle) > fabs(maxTopAngle)) { maxTopAngle = topAngle; }
      double bottomAngle = evaluateAngleRadians(bottom, color.coeffs);
      if (fabs(bottomAngle) > fabs(maxBottomAngle)) { maxBottomAngle = bottomAngle; }
    }
  }
  if (g_verbose) {
    std::cout << "leftAngle (degrees): " << maxLeftAngle * 180 / MY_PI << std::endl;
    std::cout << "rightAngle (degrees): " << maxRightAngle * 180 / MY_PI << std::endl;
    std::cout << "bottomAngle (degrees): " << maxBottomAngle * 180 / MY_PI << std::endl;
    std::cout << "topAngle (degrees): " << maxTopAngle * 180 / MY_PI << std::endl;
  }
  double leftScreenLeft, leftScreenRight, leftScreenBottom, leftScreenTop;
  double rightScreenLeft, rightScreenRight, rightScreenBottom, rightScreenTop;
  rightScreenBottom = leftScreenBottom = depth * tan(maxBottomAngle);
  rightScreenTop = leftScreenTop = depth * tan(maxTopAngle);
  if (useRightEye) {
    rightScreenLeft = depth * tan(maxLeftAngle);
    rightScreenRight = depth * tan(maxRightAngle);

    leftScreenLeft = -depth * tan(maxRightAngle);
    leftScreenRight = -depth * tan(maxLeftAngle);
  } else {
    leftScreenLeft = depth * tan(maxLeftAngle);
    leftScreenRight = depth * tan(maxRightAngle);

    rightScreenLeft = -depth * tan(maxRightAngle);
    rightScreenRight = -depth * tan(maxLeftAngle);
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
  double units = (left - right) / (tan(maxLeftAngle) - tan(maxRightAngle));

  // This means that the coefficients are providing us the correct mapping,
  // but we need to convert them into the appropriate distance scale and put
  // them into the appropriate coefficients.  The appropriate coefficients
  // are the 1st, 3rd, and so forth.  We skip all even coefficients.
  // Always push  back 0 for the 0th-order term.
  for (auto &eye : eyeColors) {
    for (auto &color : eye) {
      color.polynomial.push_back(0);
      for (size_t i = 0; i < color.coeffs.size(); i++) {
        // Zero all even terms.  We already zeroed the 0th-order term
        if (i > 0) { color.polynomial.push_back(0); }
        // Scale by the units conversion.
        color.polynomial.push_back(color.coeffs[i] * pow(units, 1 + 2 * i));
      }
    }
  }

  // Set the distance units to twice the distance from the center of
  // the display to the farthest direction so that the point at the center
  // of the furthest edge will map to itself.  This is in the tangent
  // space.
  // (For the square case, this is just twice the right edge.)
  // NOTE: We take the absolute value of this result here, in case the
  // distortion calibration inverted the image when reporting its
  // results.  This leaves the resulting polynomial parameters positive
  // while making sure the distance scale is also positive.
  double distance_scale = fabs(2*tan(maxRightAngle));

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
  if (g_verbose && eyeColors[0][0].polynomial.size() > 2) {
    // Derivative of the polynomial at 0 depends only on the first-order coefficient
    std::cout << "Display pixel slope at center: " << eyeColors[0][0].polynomial[1] << std::endl;
    if (eyeColors[0][0].polynomial[1] > 0) {
      std::cout << "renderManagerConfig:renderManagerConfig:renderOversampleFactor needed to compensate: "
        << 1/eyeColors[0][0].polynomial[1] << std::endl;
    }
  }

  //====================================================================
  // Construct Json screen description.
  // We do this by hand rather than using JsonCPP because we want
  // to control the printed precision of the numbers to avoid making
  // a huge file.
  std::cout << std::setprecision(4);
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

  std::cout << "   \"eyes\": [" << std::endl;
  std::cout << "    {" << std::endl;
  std::cout << "     \"distortion\": {" << std::endl;
  std::cout << "       \"distance_scale_x\": " << distance_scale << "," << std::endl;
  std::cout << "       \"distance_scale_y\": " << distance_scale << "," << std::endl;
  std::cout << "       \"polynomial_coeffs_red\": ";
  writePolynomial(std::cout, eyeColors[0][0].polynomial);
  std::cout << "," << std::endl;
  std::cout << "       \"polynomial_coeffs_green\": ";
  writePolynomial(std::cout, eyeColors[0][1].polynomial);
  std::cout << "," << std::endl;
  std::cout << "       \"polynomial_coeffs_blue\": ";
  writePolynomial(std::cout, eyeColors[0][2].polynomial);
  std::cout << std::endl;
  std::cout << "     }," << std::endl; // distortion
  std::cout << "     \"center_proj_x\": "
    << leftScreen.xCOP
    << "," << std::endl;
  std::cout << "     \"center_proj_y\": "
    << leftScreen.yCOP
    << "," << std::endl;
  std::cout << "     \"rotate_180\": 0" << std::endl;
  std::cout << "    }," << std::endl;
  std::cout << "    {" << std::endl;
  std::cout << "     \"distortion\": {" << std::endl;
  std::cout << "       \"distance_scale_x\": " << distance_scale << "," << std::endl;
  std::cout << "       \"distance_scale_y\": " << distance_scale << "," << std::endl;
  std::cout << "       \"polynomial_coeffs_red\": ";
  writePolynomial(std::cout, eyeColors[1][0].polynomial);
  std::cout << "," << std::endl;
  std::cout << "       \"polynomial_coeffs_green\": ";
  writePolynomial(std::cout, eyeColors[1][1].polynomial);
  std::cout << "," << std::endl;
  std::cout << "       \"polynomial_coeffs_blue\": ";
  writePolynomial(std::cout, eyeColors[1][2].polynomial);
  std::cout << std::endl;
  std::cout << "     }," << std::endl; // distortion
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
