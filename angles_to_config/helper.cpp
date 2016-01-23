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
#include <iomanip>
#include <cmath>
#include <map>

std::vector<Mapping> read_from_infile(std::istream &in)
{
  std::vector<Mapping> mapping;

  while (!in.eof()) {
    // Read the mapping info from the input file.
    Mapping map;
    in >> map.xyLatLong.longitude >> map.xyLatLong.latitude >> map.xyLatLong.x >> map.xyLatLong.y;
    mapping.push_back(map);
  }
  // There will have been one extra added, when running into EOF.
  mapping.pop_back();

  return mapping;
}


bool convert_to_normalized_and_meters(
  std::vector<Mapping> &mapping, double toMeters, double depth,
  double left, double bottom, double right, double top,
  bool useFieldAngles)
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

    if (useFieldAngles) {
      // These are expressed as angles with respect to a screen that is
      // straight ahead, independent in X and Y.  The -Z axis is straight
      // ahead.  Positive rotation in longitude points towards +X,
      // positive rotation in latitude points towards +Y.
      mapping[i].xyz.x =  depth * tan(mapping[i].xyLatLong.longitude);
      mapping[i].xyz.y =  depth * tan(mapping[i].xyLatLong.latitude);
      mapping[i].xyz.z = -depth;
    } else {
      // Compute the 3D coordinate of the point w.r.t. the eye at the origin.
      // longitude = 0, latitude = 0 points along the -Z axis in eye space.
      // Positive rotation in longitude is towards -X and positive rotation in
      // latitude points towards +Y.
      double theta = mapping[i].xyLatLong.longitude;
      double phi = MY_PI / 2 - mapping[i].xyLatLong.latitude;
      mapping[i].xyz.y = depth * cos(phi);
      mapping[i].xyz.z = -depth * cos(theta) * sin(phi);
      mapping[i].xyz.x = -depth * (-sin(theta)) * sin(phi);
    }
  }

  // Make sure that the normalized screen coordinates are all within the range 0 to 1.
  for (size_t i = 0; i < mapping.size(); i++) {
    if ((mapping[i].xyLatLong.x < 0) || (mapping[i].xyLatLong.x > 1)) {
      std::cerr << "Warning: Point " << i << " (line " << i+1 << " in the file):"
        << " x out of range [0,1]: "
        << mapping[i].xyLatLong.x << " (increase bounds on -screen or don't specify it)"
        << std::endl;
    }
    if ((mapping[i].xyLatLong.y < 0) || (mapping[i].xyLatLong.y > 1)) {
      std::cerr << "Warning: Point " << i << " (line " << i + 1 << " in the file):"
        << " y out of range [0,1]: "
        << mapping[i].xyLatLong.y << " (increase bounds on -screen or don't specify it)"
        << std::endl;
    }
  }

  return true;
}

bool findScreen(const std::vector<Mapping> &mapping,
  double left, double bottom, double right, double top,
  ScreenDescription &screen, bool verbose)
{
  if (mapping.size() == 0) {
    std::cerr << "findScreen(): Error: No points in mapping" 
      << std::endl;
    return false;
  }

  //====================================================================
  // Figure out the X screen-space extents.
  // The X screen-space extents are defined by the lines perpendicular to the
  // Y axis passing through:
  //  left: the point location whose reprojection into the Y = 0 plane has the most -
  //        positive angle(note that this may not be the point with the largest
  //        longitudinal coordinate, because of the impact of changing latitude on
  //        X - Z position).
  //  right : the point location whose reprojection into the Y = 0 plane has the most -
  //        negative angle(note that this may not be the point with the smallest
  //        longitudinal coordinate, because of the impact of changing latitude on
  //        X - Z position).
  XYZ &screenLeft = screen.screenLeft;
  XYZ &screenRight = screen.screenRight;;
  screenLeft = screenRight = mapping[0].xyz;
  if (verbose) {
    std::cerr << "First point rotation about Y (degrees): "
      << screenLeft.rotationAboutY() * 180 / MY_PI << std::endl;
  }
  for (size_t i = 0; i < mapping.size(); i++) {
    if (mapping[i].xyz.rotationAboutY() > screenLeft.rotationAboutY()) {
      screenLeft = mapping[i].xyz;
    }
    if (mapping[i].xyz.rotationAboutY() < screenRight.rotationAboutY()) {
      screenRight = mapping[i].xyz;
    }
  }
  if (verbose) {
    std::cerr << "Horizontal angular range: "
      << 180 / MY_PI * (screenLeft.rotationAboutY() - screenRight.rotationAboutY())
      << std::endl;
  }
  if (screenLeft.rotationAboutY() - screenRight.rotationAboutY() >= MY_PI) {
    std::cerr << "findScreen(): Error: Field of view > 180 degrees: found " <<
      180 / MY_PI * (screenLeft.rotationAboutY() - screenRight.rotationAboutY())
      << std::endl;
    return false;
  }

  //====================================================================
  // Find the plane of the screen, using the equation that has the normal
  // pointing towards the origin.  This is AX + BY + CZ + D = 0, where the
  // normal is in A, B, C and the offset is in D.
  //   Two points on the plane are given above.  Two more are the projection
  // of each of these points into the Y=0 plane.  We take the cross
  // product of the line from the left-most projected point to the right-
  // most projected point with the vertical line to get the normal to that
  // plane that points towards the origin.  Then we normalize this and plug
  // it back into the equation to solve for D.
  //   We're crossing with the vector (0, 1, 0), so we get:
  //   x = -dz
  //   y = 0
  //   z = dx
  double dx = screenRight.x - screenLeft.x;
  double dz = screenRight.z - screenLeft.z;
  double &A = screen.A;
  double &B = screen.B;
  double &C = screen.C;
  double &D = screen.D;
  A = -dz;
  B = 0;
  C = dx;
  double len = sqrt(A*A + B*B + C*C);
  A /= len;
  B /= len;
  C /= len;
  D = -(A*screenRight.x + B*screenRight.y + C*screenRight.z);
  if (verbose) {
    std::cerr << "Plane of the screen A,B,C, D: "
      << A << "," << B << "," << C << ", " << D
      << std::endl;
  }

  //====================================================================
  // Figure out the Y screen-space extents.
  // The Y screen-space extents are symmetric and correspond to the lines parallel
  //  to the screen X axis that are within the plane of the X line specifying the
  //  axis extents at the largest magnitude angle up or down from the horizontal.
  // Find the highest-magnitude Y value of all points when they are
  // projected into the plane of the screen.
  double &maxY = screen.maxY;
  maxY = fabs(mapping[0].xyz.projectOntoPlane(A, B, C, D).y);
  for (size_t i = 1; i < mapping.size(); i++) {
    double Y = fabs(mapping[i].xyz.projectOntoPlane(A, B, C, D).y);
    if (Y > maxY) { maxY = Y; }
  }
  if (verbose) {
    std::cerr << "Maximum-magnitude Y projection: " << maxY << std::endl;
  }

  //====================================================================
  // Figure out the monocular horizontal field of view for the screen.
  // Find the distance between the left and right points projected
  // into the Y=0 plane.  The FOV is twice the arctangent of half of this
  // distance divided by the distance to the screen.  In the configuration
  // file, this angle assumes that the center of projection is the center
  // of the screen, and an eye at the specified distance from the plane
  // of the screen for any COP will have the same distance for the COP
  // shifted to be centered.
  XYZ leftProj = screenLeft;
  XYZ rightProj = screenRight;
  leftProj.y = 0;
  rightProj.y = 0;
  double screenWidth = leftProj.distanceFrom(rightProj);
  if (verbose) {
    std::cerr << "Screen width: " << screenWidth << std::endl;
  }
  double hFOVRadians = 2 * atan((screenWidth / 2) / fabs(D));
  double hFOVDegrees = hFOVRadians * 180 / MY_PI;
  if (verbose) {
    std::cerr << "Horizontal field of view (degrees): " << hFOVDegrees << std::endl;
  }

  //====================================================================
  // Figure out the monocular vertical field of view for the screen.
  // The FOV is twice the arctangent of half of the Y
  // distance divided by the distance to the screen.
  double vFOVRadians = 2 * atan(maxY / fabs(D));
  double vFOVDegrees = vFOVRadians * 180 / MY_PI;
  if (verbose) {
    std::cerr << "Vertical field of view (degrees): " << vFOVDegrees << std::endl;
  }

  //====================================================================
  // Figure out the overlap percent for the screen that corresponds to
  // the angle between straight ahead and the normal to the plane.  First
  // find the angle itself, and then the associated overlap percent.
  // The percent overlap assumes that the center of projection is at the
  // center of the screen and that both eyes are at the origin (discounts
  // the IPD).
  // The angle is determined by finding the location at which the vector
  // along the -Z axis from the eye location (which is at the origin)
  // pierces the screen.  The angle between this vector and the normal
  // to the plane of the screen determines how much to rotate the screen.
  // Because the angle between forward direction and the -Z axis is 0,
  // the angle depends only on the unit normal to the plane,
  // which is (A,B,C), but B = 0 and we only care about rotation
  // around the Y axis.  For the purpose of the atan function, the part
  // of X is played by the -Z axis and the part of Y is played by the
  // -X axis.  A is associated with the X axis and C with the Z axis.
  // Here is the code we are inverting to go from angle to overlap...
  //  double overlapFrac = m_params.m_displayConfiguration.getOverlapPercent();
  //  const auto hfov = m_params.m_displayConfiguration.getHorizontalFOV();
  //  const auto angularOverlap = hfov * overlapFrac;
  //  rotateEyesApart = (hfov - angularOverlap) / 2;
  // Here is the inversion:
  //  rotateEyesApart = (hfov - (hfov * overlapFrac)) / 2;
  //  2 * rotateEyesApart = (hfov - (hfov * overlapFrac));
  //  2 * rotateEyesApart - hfov = - hfov * overlapFrac;
  //  1 - 2*rotateEyesApart/hfov = overlapFrac
  double angleRadians = fabs(atan2(A, C));
  if (verbose) {
    std::cerr << "Angle degrees: " << angleRadians * 180 / MY_PI << std::endl;
  }
  double overlapFrac = 1 - 2 * angleRadians / hFOVRadians;
  double overlapPercent = overlapFrac * 100;
  if (verbose) {
    std::cerr << "Overlap percent: " << overlapPercent << std::endl;
  }

  //====================================================================
  // Figure out the center of projection for the screen.  This is the
  // location where a line from the origin perpendicular to the screen
  // pierces the screen.
  // Then figure out the normalized coordinates of this point in screen
  // space, which is the fraction of the way from the left to the right
  // of the screen.  It is always (by construction above) in the center
  // of the screen in Y.
  // Also, by construction it is at a distance D along the (A,B,C) unit
  // vector from the origin.
  double yCOP = 0.5;
  XYZ projection;
  projection.x = -D * A;
  projection.y = -D * B;
  projection.z = -D * C;
  double xCOP = leftProj.distanceFrom(projection) / leftProj.distanceFrom(rightProj);
  if (verbose) {
    std::cerr << "Center of projection x,y: " << xCOP << ", " << yCOP << std::endl;
  }

  //====================================================================
  // Fill in the screen parameters.
  screen.hFOVDegrees = hFOVDegrees;
  screen.vFOVDegrees = vFOVDegrees;
  screen.overlapPercent = overlapPercent;
  screen.xCOP = xCOP;
  screen.yCOP = yCOP;

  return true;
}

bool findMesh(const std::vector<Mapping> &mapping,
  double left, double bottom, double right, double top,
  ScreenDescription const &screen, MeshDescription &mesh, bool verbose)
{
  if (mapping.size() == 0) {
    std::cerr << "findMesh(): Error: No points in mapping"
      << std::endl;
    return false;
  }

  const double &A = screen.A;
  const double &B = screen.B;
  const double &C = screen.C;
  const double &D = screen.D;

  //====================================================================
  // Map each incoming mesh coordinate into the corresponding output
  // coordinate, storing them into the output mesh.
  mesh.clear();

  // Scale and offset to apply to the points projected onto the plane
  // to convert their values into normalized screen coordinates.  This
  // checks the assumption that the screen has some width in the X
  // coordinate (not rotated 90 degrees) and uses only the X value to
  // scale X.
  const XYZ &leftProj = screen.screenLeft;
  const XYZ &rightProj = screen.screenRight;
  if (leftProj.x == rightProj.x) {
    std::cerr << "Error computing mesh: screen has no X extent" << std::endl;
    return false;
  }
  double xOutOffset = -leftProj.x;
  double xOutScale = 1 / (rightProj.x - leftProj.x);
  double yOutOffset = screen.maxY; // Negative of negative maxY is maxY
  double yOutScale = 1 / (2 * screen.maxY);

  for (size_t i = 0; i < mapping.size(); i++) {

    // Input point coordinates are already normalized.
    double xNormIn = mapping[i].xyLatLong.x;
    double yNormIn = mapping[i].xyLatLong.y;
    std::array<double, 2> in;
    in[0] = xNormIn;
    in[1] = yNormIn;

    // Project the 3D points back into the plane of the screen and determine
    // the normalized coordinates in the coordinate system with the lower left
    // corner at (0,0) and the upper right at (1,1).  Because we oversized the
    // screen, these will all be in this range.  Otherwise, they might not be.
    XYZ onScreen = mapping[i].xyz.projectOntoPlane(A, B, C, D);
    double xNormOut = (onScreen.x + xOutOffset) * xOutScale;
    double yNormOut = (onScreen.y + yOutOffset) * yOutScale;
    std::array<double, 2> out;
    out[0] = xNormOut;
    out[1] = yNormOut;

    // Add a new entry onto the mesh
    std::array< std::array<double, 2>, 2 > element;
    element[0] = in;
    element[1] = out;

    mesh.push_back(element);
  }

  return true;
}

static void normalize(std::array<double, 2> &v)
{
  double len = sqrt(v[0] * v[0] + v[1] * v[1]);
  if (len > 0) {
    v[0] /= len;
    v[1] /= len;
  }
}

/// Finds out whether the neighbor of the specified index
/// in the mapping violates the strictures of the
/// remove_invalid_points_based_on_angle function.
/// @return True if the neighbor angle difference is too large.
static bool neighbor_error(
  const std::vector<Mapping> &mapping,
  size_t index0, size_t index1,
  double xx, double xy,
  double yx, double yy, double minDotProduct)
{
  // Map the difference in angle space between the point and
  // its neighbor into screen space.
  std::array<double, 2> angleVec;
  angleVec[0] = mapping[index1].xyLatLong.longitude -
    mapping[index0].xyLatLong.longitude;
  angleVec[1] = mapping[index1].xyLatLong.latitude -
    mapping[index0].xyLatLong.latitude;
  std::array<double, 2> screenMappedVec;
  screenMappedVec[0] = angleVec[0] * xx + angleVec[1] * yx;
  screenMappedVec[1] = angleVec[0] * xy + angleVec[1] * yy;

  // Find the screen-space difference between the point and its
  // neighbor.
  std::array<double, 2> screenVec;
  screenVec[0] = mapping[index1].xyLatLong.x -
    mapping[index0].xyLatLong.x;
  screenVec[1] = mapping[index1].xyLatLong.y -
    mapping[index0].xyLatLong.y;

  // Normalize the two vectors
  normalize(screenMappedVec);
  normalize(screenVec);

  // Find the dot product between the two vectors.
  double dotProduct = screenMappedVec[0] * screenVec[0] +
    screenMappedVec[1] * screenVec[1];

  // If the dot product is too small, the angle is too large
  // and so we fail.
  return dotProduct < minDotProduct;
}

static double point_distance(double x1, double y1, double x2, double y2) {
  return std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

/// Finds out how many neighbors of the specified index
/// in the mapping violate the strictures of the
/// remove_invalid_points_based_on_angle function.
/// Checks up to 8 nearest neighbors in angle space.
/// @return How many neighbor angle differences are too large.
static size_t neighbor_errors(
  const std::vector<Mapping> &mapping, size_t index,
  double xx, double xy,
  double yx, double yy, double minDotProduct)
{
  // Sort all points other than the one we're looking for
  // based on distance from the one we are looking for.
  typedef std::multimap<double, size_t> PointDistanceIndexMap;
  PointDistanceIndexMap map;
  for (size_t i = 0; (i < mapping.size()) && (i != index); i++) {
    map.insert(std::make_pair(
      point_distance(
        mapping[index].xyLatLong.longitude,
        mapping[index].xyLatLong.latitude,
        mapping[i].xyLatLong.longitude,
        mapping[i].xyLatLong.latitude
      ), i));
  }

  // Check the first 8 points on the list (if there are that many)
  // and count up how many violate the angle condition.
  size_t ret = 0;
  PointDistanceIndexMap::const_iterator it = map.begin();
  for (size_t i = 0; (i < map.size()) && (i < 8); i++) {
    if (neighbor_error(mapping, index, it->second,
        xx, xy, yx, yy, minDotProduct)) {
      ret++;
    }
    it++;
  }
  return ret;
}

/// Looks through the mapping to see if there is a point whose
/// neighbor angles violate the strictures of the
/// remove_invalid_points_based_on_angle function.
/// @return index of the point that has the largest number
/// of invalid neighbor angles if one is found, size of mapping
/// if not.
static size_t find_index_of_angle_worst_offender(
  std::vector<Mapping> &mapping, double xx, double xy,
  double yx, double yy, double minDotProduct)
{
  size_t worstIndex = 0;
  size_t worstCount = neighbor_errors(mapping, 0,
    xx, xy, yx, yy, minDotProduct);
  for (size_t i = 1; i < mapping.size(); i++) {
    size_t count = neighbor_errors(mapping, i,
      xx, xy, yx, yy, minDotProduct);
    if (count > worstCount) {
      worstCount = count;
      worstIndex = i;
    }
  }
  if (worstCount == 0) { return mapping.size(); }
  return worstIndex;
}

int remove_invalid_points_based_on_angle(
  std::vector<Mapping> &mapping, double xx, double xy,
  double yx, double yy, double maxAngleDegrees)
{
  int ret = 0;  // No points yet removed from the mesh.

  // Find the dot product associated with two unit vectors
  // separated by the angle specified.  This is the cosine
  // of the angle.
  double minDotProduct = cos(maxAngleDegrees / 180.0 * MY_PI);

  // We remove the worst offender from the list each time,
  // then re-start.  Assuming that we get the actual outlier,
  // as opposed to one of its neighbors, this avoids trimming
  // too many points from the vector.
  bool foundOutlier;
  do {
    foundOutlier = false;
    size_t off = find_index_of_angle_worst_offender(
      mapping, xx, xy, yx, yy, minDotProduct);
    if (off < mapping.size()) {
      mapping.erase(mapping.begin() + off);
      foundOutlier = true;
      ret++;
    }
  } while (foundOutlier);

  return ret;
}
