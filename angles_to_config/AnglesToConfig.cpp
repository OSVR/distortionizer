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
#include <cmath>
#include <vector>
#include <stdlib.h> // For exit()

// Global constants and variables
#define MY_PI (4.0*atan(1.0))
static bool g_verbose = false;

// Forward declarations
static int testAlgorithms();

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
} ScreenDescription;

/// Holds a list of mappings from physical-display normalized
/// coordinates to canonical-display normalized coordinates.
typedef std::vector<        //!< Vector of mappings
    std::vector<            //!< 2-vector of from, to coordinates
      std::vector<double>   //!< 2-vector of unit coordinates (x,y)
  >
> MeshDescription;

bool findScreenAndMesh(const std::vector<Mapping> &mapping,
  double left, double bottom, double right, double top,
  ScreenDescription &screen, MeshDescription &mesh)
{
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
  XYZ screenLeft, screenRight;
  screenLeft = screenRight = mapping[0].xyz;
  if (g_verbose) {
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
  if (g_verbose) {
    std::cerr << "Horizontal angular range: "
      << 180 / MY_PI * (screenLeft.rotationAboutY() - screenRight.rotationAboutY())
      << std::endl;
  }
  if (screenLeft.rotationAboutY() - screenRight.rotationAboutY() >= MY_PI) {
    std::cerr << "Error: Field of view > 180 degrees: found " <<
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
  double A = -dz;
  double B = 0;
  double C = dx;
  double len = sqrt(A*A + B*B + C*C);
  A /= len;
  B /= len;
  C /= len;
  double D = -(A*screenRight.x + B*screenRight.y + C*screenRight.z);
  if (g_verbose) {
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
  double maxY = fabs(mapping[0].xyz.projectOntoPlane(A, B, C, D).y);
  for (size_t i = 1; i < mapping.size(); i++) {
    double Y = fabs(mapping[i].xyz.projectOntoPlane(A, B, C, D).y);
    if (Y > maxY) { maxY = Y; }
  }
  if (g_verbose) {
    std::cerr << "Maximum-magnitude Y projection: " << maxY << std::endl;
  }

  //====================================================================
  // Figure out the monocular horizontal field of view for the screen.
  // Find the distance between the left and right points projected
  // into the Y=0 plane.  The FOV is twice the arctangent of half of this
  // distance divided by the distance to the screen.
  XYZ leftProj = screenLeft;
  XYZ rightProj = screenRight;
  leftProj.y = 0;
  rightProj.y = 0;
  double screenWidth = leftProj.distanceFrom(rightProj);
  if (g_verbose) {
    std::cerr << "Screen width: " << screenWidth << std::endl;
  }
  double hFOVRadians = 2 * atan((screenWidth / 2) / fabs(D));
  double hFOVDegrees = hFOVRadians * 180 / MY_PI;
  if (g_verbose) {
    std::cerr << "Horizontal field of view (degrees): " << hFOVDegrees<< std::endl;
  }

  //====================================================================
  // Figure out the monocular vertical field of view for the screen.
  // The FOV is twice the arctangent of half of the Y
  // distance divided by the distance to the screen.
  double vFOVRadians = 2 * atan(maxY / fabs(D));
  double vFOVDegrees = vFOVRadians * 180 / MY_PI;
  if (g_verbose) {
    std::cerr << "Vertical field of view (degrees): " << vFOVDegrees << std::endl;
  }

  //====================================================================
  // Figure out the overlap percent for the screen that corresponds to
  // the angle between straight ahead and the normal to the plane.  First
  // find the angle itself, and then the associated overlap percent.
  // The angle is determined based on the unit normal to the plane,
  // which is (A,B,C), but B = 0 and we only care about rotation
  // around the Y axis.  For the purpose of the atan function, the part
  // of X is played by the -Z axis and the part of Y is played by the
  // -X axis.  A is associated with the X axis and C with the Z axis.
  // Here is the code we are inverting...
  //  double overlapFrac = m_params.m_displayConfiguration.getOverlapPercent();
  //  const auto hfov = m_params.m_displayConfiguration.getHorizontalFOV();
  //  const auto angularOverlap = hfov * overlapFrac;
  //  rotateEyesApart = (hfov - angularOverlap) / 2.;
  // Here is the inversion:
  //  rotateEyesApart = (hfov - (hfov * overlapFrac));
  //  rotateEyesApart - hfov = - hfov * overlapFrac;
  //  1 - rotateEyesApart/hfov = overlapFrac
  double angleRadians = fabs(atan2(A, C));
  double overlapFrac = 1 - angleRadians / hFOVRadians;
  double overlapPercent = overlapFrac * 100;
  if (g_verbose) {
    std::cerr << "Overlap percent: " << overlapPercent<< std::endl;
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
  if (g_verbose) {
    std::cerr << "Center of projection x,y: " << xCOP << ", " << yCOP << std::endl;
  }

  //====================================================================
  // Fill in the screen parameters.
  screen.hFOVDegrees = hFOVDegrees;
  screen.vFOVDegrees = vFOVDegrees;
  screen.overlapPercent = overlapPercent;
  screen.xCOP = xCOP;
  screen.yCOP = yCOP;

  //====================================================================
  // Map each incoming mesh coordinate into the corresponding output
  // coordinate, storing them into the output mesh.
  mesh.clear();

  // Scale and offset to apply to the points projected onto the plane
  // to convert their values into normalized screen coordinates.  This
  // checks the assumption that the screen has some width in the X
  // coordinate (not rotated 90 degrees) and uses only the X value to
  // scale X.
  if (leftProj.x == rightProj.x) {
    std::cerr << "Error computing mesh: screen has no X extent" << std::endl;
    return false;
  }
  double xOutOffset = -leftProj.x;
  double xOutScale = 1 / (rightProj.x - leftProj.x);
  double yOutOffset = maxY; // Negative of negative maxY
  double yOutScale = 1 / (2 * maxY);

  for (size_t i = 0; i < mapping.size(); i++) {

    // Input point coordinates are already normalized.
    double xNormIn = mapping[i].xyLatLong.x;
    double yNormIn = mapping[i].xyLatLong.y;
    std::vector<double> in;
    in.push_back(xNormIn);
    in.push_back(yNormIn);

    // Project the 3D points back into the plane of the screen and determine
    // the normalized coordinates in the coordinate system with the lower left
    // corner at (0,0) and the upper right at (1,1).  Because we oversized the
    // screen, these will all be in this range.  Otherwise, they might not be.
    XYZ onScreen = mapping[i].xyz.projectOntoPlane(A, B, C, D);
    double xNormOut = (onScreen.x + xOutOffset) * xOutScale;
    double yNormOut = (onScreen.y + yOutOffset) * yOutScale;
    std::vector<double> out;
    out.push_back(xNormOut);
    out.push_back(yNormOut);

    // Add a new entry onto the mesh
    std::vector< std::vector<double> > element;
    element.push_back(in);
    element.push_back(out);

    mesh.push_back(element);
  }

  return true;
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
  // Convert the input values into normalized coordinates and into 3D
  // locations.
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

    if (g_verbose) {
      if (i == 0) {
        std::cerr << "First point:" << std::endl
          << " Lat/long rad: " << mapping[i].xyLatLong.latitude
          << "/" << mapping[i].xyLatLong.longitude
          << ", norm x,y: " << mapping[i].xyLatLong.x << "," << mapping[i].xyLatLong.y
          << ", meters x,y,z: " << mapping[i].xyz.x << "," << mapping[i].xyz.y << "," << mapping[i].xyz.z
          << std::endl;
      }
    }
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

  //====================================================================
  // Determine the screen description and distortion mesh based on the
  // input points and screen parameters.
  ScreenDescription screen;
  MeshDescription mesh;
  if (!findScreenAndMesh(mapping, left, bottom, right, top, screen, mesh)) {
    std::cerr << "Error: Could not find screen or mesh" << std::endl;
    return 3;
  }
  if (mesh.size() != mapping.size()) {
    std::cerr << "Error: Mesh size " << mesh.size()
      << " does not match mapping size" << mapping.size() << std::endl;
    return 4;
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
  // @todo
  std::cout << "   }," << std::endl; // field_of_view

  std::cout << "   \"eyes\": {" << std::endl;
  // @todo
  std::cout << "   }" << std::endl; // eyes

  std::cout << "  }" << std::endl;  // hmd
  std::cout << " }" << std::endl;   // display
  std::cout << "}" << std::endl;    // Closes outer object

  /*
  Json::Value jRoot;
  Json::Value jDisplay;
  Json::Value jHmd;
  Json::Value jFOV;
  Json::Value jHFOV = screen.hFOVDegrees;
  Json::Value jVFOV = screen.vFOVDegrees;
  Json::Value jOverlap = screen.overlapPercent;
  Json::Value jPitch = 0.0;
  jFOV["monocular_horizontal"] = jHFOV;
  jFOV["monocular_vertical"] = jVFOV;
  jFOV["overlap_percent"] = jOverlap;
  jFOV["pitch_tilt"] = jPitch;
  jHmd["field_of_view"] = jFOV;

  // Construct Json eye description
  Json::Value eyes;
  Json::Value leftEye, rightEye;
  double leftEyeXCOP, rightEyeXCOP;
  double invertXCOP = 1.0 - screen.xCOP;
  if (useRightEye) {
    leftEyeXCOP = invertXCOP;
    rightEyeXCOP = screen.xCOP;
  } else {
    leftEyeXCOP = screen.xCOP;
    rightEyeXCOP = invertXCOP;
  }
  leftEye["center_proj_x"] = leftEyeXCOP;
  leftEye["center_proj_y"] = screen.yCOP;
  leftEye["rotate_180"] = 0;
  rightEye["center_proj_x"] = rightEyeXCOP;
  rightEye["center_proj_y"] = screen.yCOP;
  rightEye["rotate_180"] = 0;
  eyes[0] = leftEye;
  eyes[1] = rightEye;
  jHmd["eyes"] = eyes;

  // Construct the Json distortion mesh description and add it to
  // the existing HMD description.
  Json::Value jDistortionType = "mono_point_samples";
  Json::Value jDistortion;
  Json::Value jMesh;
  for (int i = 0; i < mesh.size(); i++) {
    Json::Value jIn;
    jIn[0] = mesh[i][0][0];
    jIn[1] = mesh[i][0][1];
    Json::Value jOut;
    jOut[0] = mesh[i][1][0];
    jOut[1] = mesh[i][1][1];
    Json::Value jEntry;
    jEntry[0] = jIn;
    jEntry[1] = jOut;
    jMesh[i] = jEntry;
  }
  jDistortion["type"] = jDistortionType;
  jDistortion["mono_point_samples"] = jMesh;
  jHmd["distortion"] = jDistortion;

  // Store all of this into the hierarchy.
  jDisplay["hmd"] = jHmd;
  jRoot["display"] = jDisplay;

  //====================================================================
  // Write the complete description.
  Json::StyledWriter jWriter;
  std::cout << jWriter.write(jRoot) << std::endl;
  */

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
  // distortion correction.
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
  if (!findScreenAndMesh(mapping, 0, 0, 1, 1, screen, mesh)) {
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
  // @todo

  if (g_verbose) {
    std::cerr << "=================== Successfully finished testAlgorithms()" << std::endl;
  }
  return 0;
}
