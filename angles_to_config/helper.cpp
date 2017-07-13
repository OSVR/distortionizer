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
#include "helper.h"
#include "types.h"
#include "EigenStdArrayInterop.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

// Standard includes
#include <cmath>
#include <iomanip>
#include <map>
#include <string>

template <typename T> inline T radToDegree(T radians) {
    return static_cast<T>(static_cast<double>(radians) * 180.0 / MY_PI);
}

std::vector<Mapping> read_from_infile(std::istream& in) {
    std::vector<Mapping> mapping;

    while (!in.eof()) {
        // Read the mapping info from the input file.
        Mapping map;
        in >> map.xyLatLong.longitude >> map.xyLatLong.latitude >> map.xyLatLong.x >> map.xyLatLong.y;
        mapping.push_back(map);
    }
    if (!mapping.empty()) {
        // There will have been one extra added, when running into EOF.
        mapping.pop_back();
    }

    return mapping;
}

InputMeasurements readInputMeasurements(std::string const& inputSource, std::istream& in) {
    InputMeasurements ret;
    ret.inputSource = inputSource;
    std::size_t line = 1;
    while (!in.eof()) {
        // Read the mapping info from the input file.
        InputMeasurement meas;
        in >> meas.viewAnglesDegrees.longitude() >> meas.viewAnglesDegrees.latitude() >> meas.screen[0] >>
            meas.screen[1];
        meas.lineNumber = line;
        ret.measurements.push_back(meas);
        line++;
    }
    if (!ret.empty()) {
        // There will have been one extra added, when running into EOF.
        ret.measurements.pop_back();
    }

    return ret;
}

std::vector<LongLat> readAdditionalAngles(std::istream& in) {
    std::vector<LongLat> ret;
    try {
        while (in.good()) {
            // Read the data in from the file.
            LongLat elt;
            in >> elt.longitude() >> elt.latitude();
            ret.push_back(elt);
        }
        if (!ret.empty()) {
            // There will have been one extra added, when running into EOF.
            ret.pop_back();
        }
    } catch (std::bad_alloc const&) {
        std::cerr << "Out of memory loading additional angles file - skipping." << std::endl;
        return std::vector<LongLat>{};
    }

    return ret;
}
using Plane = Eigen::Hyperplane<double, 3>;

/// takes in angles (long/lat or field angles) in degrees, returns
inline XYZ longLatToWorldSpace(LongLat longLat, bool useFieldAngles, double depth) {
    using std::sin;
    using std::cos;
    using std::tan;
    // Convert the input latitude and longitude from degrees to radians.
    ei::map(longLat.longLat) *= MY_PI / 180.;
    XYZ ret;
    if (useFieldAngles) {
        const Plane p =
            Plane::Through(Eigen::Vector3d(0, 0, -depth), Eigen::Vector3d(1, 0, -depth), Eigen::Vector3d(0, 1, -depth));
        // These are expressed as angles with respect to a screen that is
        // straight ahead, independent in X and Y.  The -Z axis is straight
        // ahead.  Positive rotation in longitude points towards +X,
        // positive rotation in latitude points towards +Y.
        ret = XYZ{depth * tan(longLat.longitude()), // X
                  depth * tan(longLat.latitude()),  // Y
                  -depth};                          // Z
    } else {
        // Compute the 3D coordinate of the point w.r.t. the eye at the origin.
        // longitude = 0, latitude = 0 points along the -Z axis in eye space.
        // Positive rotation in longitude is towards -X and positive rotation in
        // latitude points towards +Y.
        double theta = longLat.longitude();
        double phi = MY_PI / 2. - longLat.latitude();
        ret = XYZ{-depth * (-sin(theta)) * sin(phi), // X
                  depth * std::cos(phi),             // Y
                  -depth * cos(theta) * sin(phi)};   // Z
    }
    return ret;
}

XYZList convertAdditionalAngles(std::vector<LongLat> const& additionalAngles, double depth, bool useFieldAngles) {
    XYZList ret;

    for (auto& additionalLongLat : additionalAngles) {
        ret.push_back(longLatToWorldSpace(additionalLongLat, useFieldAngles, depth));
    }
    return ret;
}

bool convert_to_normalized_and_meters(std::vector<Mapping>& mapping, double toMeters, double depth, double left,
                                      double bottom, double right, double top, bool useFieldAngles) {
    for (auto& thisMapping : mapping) {
        //  Convert the input coordinates from its input space into meters
        // and then convert (using the screen dimensions) into normalized screen
        // units.
        thisMapping.xyLatLong.x *= toMeters;
        thisMapping.xyLatLong.x = (thisMapping.xyLatLong.x - left) / (right - left);
        thisMapping.xyLatLong.y *= toMeters;
        thisMapping.xyLatLong.y = (thisMapping.xyLatLong.y - bottom) / (top - bottom);

        LongLat tempLongLat = {thisMapping.xyLatLong.longitude, thisMapping.xyLatLong.latitude};
        // Convert the input latitude and longitude from degrees to radians.
        // Then, compute 3d coordinates of a point.
        thisMapping.xyz = longLatToWorldSpace(tempLongLat, useFieldAngles, depth);
    }

    // Make sure that the normalized screen coordinates are all within the range 0 to 1.
    for (size_t i = 0; i < mapping.size(); i++) {
        if ((mapping[i].xyLatLong.x < 0) || (mapping[i].xyLatLong.x > 1)) {
            std::cerr << "Warning: Point " << i << " (line " << i + 1 << " in the file if no filtration took place):"
                      << " x out of range [0,1]: " << mapping[i].xyLatLong.x
                      << " (increase bounds on -screen or don't specify it)" << std::endl;
        }
        if ((mapping[i].xyLatLong.y < 0) || (mapping[i].xyLatLong.y > 1)) {
            std::cerr << "Warning: Point " << i << " (line " << i + 1 << " in the file if no filtration took place):"
                      << " y out of range [0,1]: " << mapping[i].xyLatLong.y
                      << " (increase bounds on -screen or don't specify it)" << std::endl;
        }
    }

    return true;
}

NormalizedMeasurements convert_to_normalized_and_meters(InputMeasurements const& input, double toMeters, double depth,
                                                        RectBoundsd screenDims, bool useFieldAngles) {
    const auto screenWidth = (screenDims.right - screenDims.left);
    const auto screenHeight = (screenDims.top - screenDims.bottom);
    NormalizedMeasurements ret;
    ret.inputSource = input.inputSource;
    for (auto& inputMeas : input.measurements) {
        //  Convert the input coordinates from its input space into meters
        // and then convert (using the screen dimensions) into normalized screen
        // units.
        auto x = (inputMeas.screen[0] * toMeters - screenDims.left) / screenWidth;
        auto y = (inputMeas.screen[1] * toMeters - screenDims.bottom) / screenHeight;
        auto screen = Point2d{x, y};

        // Convert the input latitude and longitude from degrees to radians.
        // Then, compute 3d coordinates of a point.
        auto pt = longLatToWorldSpace(inputMeas.viewAnglesDegrees, useFieldAngles, depth);
        auto pointFromView = Point3d{pt.x, pt.y, pt.z};
        ret.measurements.push_back(NormalizedMeasurement{screen, pointFromView, inputMeas.lineNumber});
    }

    // Make sure that the normalized screen coordinates are all within the range 0 to 1.
    for (auto const& outputMeas : ret.measurements) {
        auto& screen = outputMeas.screen;
        if ((screen[0] < 0) || (screen[0] > 1)) {
            std::cerr << "Warning: Point from " << ret.inputSource << ":" << outputMeas.lineNumber << ":"
                      << " normalized x out of range [0,1]: " << screen[0]
                      << " (increase bounds on screen or don't specify it)" << std::endl;
        }
        if ((screen[1] < 0) || (screen[1] > 1)) {
            std::cerr << "Warning: Point from " << ret.inputSource << ":" << outputMeas.lineNumber << ":"
                      << " normalized y out of range [0,1]: " << screen[1]
                      << " (increase bounds on screen or don't specify it)" << std::endl;
        }
    }
    return ret;
}

bool findScreen(ScreenDescription& outScreen, const std::vector<Mapping>& mapping,
                XYZList const& additionalPointsFromAngles, bool verbose) {
    if (mapping.empty()) {
        std::cerr << "findScreen(): Error: No points in mapping" << std::endl;
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
    XYZ& screenLeft = outScreen.screenLeft;
    XYZ& screenRight = outScreen.screenRight;

    screenLeft = screenRight = mapping[0].xyz;
    if (verbose) {
        std::cerr << "First point rotation about Y (degrees): " << radToDegree(screenLeft.rotationAboutY())
                  << std::endl;
    }
    /// Debug print functor.
    auto horizontalAngleDebugPrint = [&](const char* when) {
        if (verbose) {
            std::cerr << "[" << when << "] Horizontal angular range: "
                      << radToDegree(screenLeft.rotationAboutY() - screenRight.rotationAboutY()) << std::endl;
            std::cerr << "Screen left: ";
            screenLeft.debugPrint(std::cerr);
            std::cerr << std::endl;
            std::cerr << "Screen right: ";
            screenRight.debugPrint(std::cerr);
            std::cerr << std::endl;
        }
    };
    /// Little utility functor so we can go through both the mappings as well as the points that came from bare angles
    /// with the same code.
    auto considerXYZForScreenBound = [&](XYZ const& xyz) {
        if (xyz.rotationAboutY() > screenLeft.rotationAboutY()) {
            screenLeft = xyz;
        }
        if (xyz.rotationAboutY() < screenRight.rotationAboutY()) {
            screenRight = xyz;
        }
    };
    for (const auto& i : mapping) {
        considerXYZForScreenBound(i.xyz);
    }
    horizontalAngleDebugPrint("mappings only");
    for (auto& extraPoint : additionalPointsFromAngles) {
        considerXYZForScreenBound(extraPoint);
    }
    horizontalAngleDebugPrint("full");
    if (screenLeft.rotationAboutY() - screenRight.rotationAboutY() >= MY_PI) {
        std::cerr << "findScreen(): Error: Field of view > 180 degrees: found "
                  << radToDegree(screenLeft.rotationAboutY() - screenRight.rotationAboutY()) << std::endl;
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

    /// @todo Note this will always give us the same plane if we feed in field angles, since field angles are
    /// transformed into xyz by projecting onto a straight-ahead screen
    double dx = screenRight.x - screenLeft.x;
    double dz = screenRight.z - screenLeft.z;
    double& A = outScreen.A;
    double& B = outScreen.B;
    double& C = outScreen.C;
    double& D = outScreen.D;
    A = -dz;
    B = 0;
    C = dx;
    double len = std::sqrt(A * A + B * B + C * C);
    A /= len;
    B /= len;
    C /= len;
    D = -(A * screenRight.x + B * screenRight.y + C * screenRight.z);
    if (verbose) {
        std::cerr << "Plane of the screen A,B,C, D: " << A << "," << B << "," << C << ", " << D << std::endl;
    }
    const Plane screenPlane(Eigen::Vector3d(A, B, C), D);

    //====================================================================
    // Figure out the Y screen-space extents.
    // The Y screen-space extents are symmetric and correspond to the lines parallel
    //  to the screen X axis that are within the plane of the X line specifying the
    //  axis extents at the largest magnitude angle up or down from the horizontal.
    // Find the highest-magnitude Y value of all points when they are
    // projected into the plane of the screen.
    double& maxY = outScreen.maxY;
#if 0
    auto computeYMagnitude = [&](XYZ const& xyz) { return std::abs(xyz.projectOntoPlane(A, B, C, D).y); };
#else
    auto computeYMagnitude = [&](XYZ const& xyz) { return std::abs(screenPlane.projection(toEigen(xyz)).y()); };
#endif
    auto considerXYZForMaxYMagnitude = [&](XYZ const& xyz) {
        double Y = computeYMagnitude(xyz);
        if (Y > maxY) {
            maxY = Y;
        }
    };
    /// initial value
    maxY = computeYMagnitude(mapping[0].xyz);

    /// Go through all mappings
    for (const auto& i : mapping) {
        considerXYZForMaxYMagnitude(i.xyz);
    }
    if (verbose) {
        std::cerr << "Maximum-magnitude Y projection after just mappings: " << maxY << std::endl;
    }

    /// Go through all points from extra angles
    for (auto& extraPoint : additionalPointsFromAngles) {
        considerXYZForMaxYMagnitude(extraPoint);
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
    const Eigen::Vector3d leftProjVec = toEigen(leftProj);
    const Eigen::Vector3d rightProjVec = toEigen(rightProj);
#ifdef ASSUME_CENTER_OF_SCREEN
    double screenWidth = leftProj.distanceFrom(rightProj);
    double hFOVRadians = 2 * std::atan((screenWidth / 2) / std::abs(D));
#else
    // get these at unit distance (if not there already)
    // which makes them rather like unit-distance left and right clipping planes
    double leftEdge = leftProj.x / (-leftProj.z);
    double rightEdge = rightProj.x / (-rightProj.z);
    // this isn't actually terribly useful here...
    double screenWidth = -leftEdge + rightEdge;
    auto leftHalfFOV = std::atan(-leftEdge);
    auto rightHalfFOV = std::atan(rightEdge);
    if (verbose) {
        std::cerr << "Left half-FOV: " << radToDegree(leftHalfFOV) << std::endl;
        std::cerr << "Right half-FOV: " << radToDegree(rightHalfFOV) << std::endl;
    }
    auto hFOVRadians = leftHalfFOV + rightHalfFOV;
#endif
    if (verbose) {
        std::cerr << "Screen width: " << screenWidth << std::endl;
    }
    double hFOVDegrees = radToDegree(hFOVRadians);
    if (verbose) {
        std::cerr << "Horizontal field of view (degrees): " << hFOVDegrees << std::endl;
    }

    //====================================================================
    // Figure out the monocular vertical field of view for the screen.
    // The FOV is twice the arctangent of half of the Y
    // distance divided by the distance to the screen.
    double vFOVRadians = 2 * std::atan(maxY / std::abs(D));
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
    double angleRadians = std::abs(atan2(A, C));
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
#if 0
    double yCOP = 0.5;
    XYZ projection;
    projection.x = -D * A;
    projection.y = -D * B;
    projection.z = -D * C;
    double xCOP = leftProj.distanceFrom(projection) / leftProj.distanceFrom(rightProj);
#else
    Eigen::Vector3d copOnPlane = screenPlane.projection(Eigen::Vector3d::Zero());
    double xCOP = (copOnPlane.x() - leftProjVec.x()) / screenWidth;
    double yCOP = 0.5;
#endif
    if (verbose) {
        std::cerr << "Center of projection x,y: " << xCOP << ", " << yCOP << std::endl;
    }

    //====================================================================
    // Fill in the screen parameters.
    outScreen.hFOVDegrees = hFOVDegrees;
    outScreen.vFOVDegrees = vFOVDegrees;
    outScreen.overlapPercent = overlapPercent;
    outScreen.xCOP = xCOP;
    outScreen.yCOP = yCOP;

    return true;
}

template <typename Derived> inline double tanAboutY(Eigen::MatrixBase<Derived> const& vec) {
    return vec.x() / -vec.z();
}
template <typename Derived> inline double rotationAboutY(Eigen::MatrixBase<Derived> const& vec) {
    return std::atan2(vec.x(), -vec.z());
}
inline std::string to_string(Eigen::Vector3d const& vec) {

    static const auto PRECISION = 4;
    static const auto WIDTH = PRECISION + 3;
    std::ostringstream ss;
    ss << std::setprecision(PRECISION);
    ss << "(" << std::setw(WIDTH) << vec.x();
    ss << ", " << std::setw(WIDTH) << vec.y();
    ss << ", " << std::setw(WIDTH) << vec.z() << ")";
    return ss.str();
}
struct Vec3dAndOrigin {
    Eigen::Vector3d vec;
    DataOrigin origin;
};
template <typename Comparator> class ViewExtrema : public GenericExtremaFinder<Vec3dAndOrigin, Comparator> {
  public:
    using base = GenericExtremaFinder<Vec3dAndOrigin, Comparator>;
    ViewExtrema(Comparator&& comp) : base(std::move(comp)) {}

    /// Process a NormalizedMeasurements struct's contents.
    void process(NormalizedMeasurements const& data) {
        for (const auto& meas : data.measurements) {
            base::process(Vec3dAndOrigin{ei::map(meas.pointFromView), meas.getOrigin(data)});
        }
    }
    /// Process a vector of NormalizedMeasurements structs
    void process(std::vector<NormalizedMeasurements> const& channels) {
        for (const auto& data : channels) {
            process(data);
        }
    }

    /// Process the XYZList of extra angles loaded from a single file.
    ///
    /// @note Synthesizes DataOrigin objects for the elements.
    void processExtraAngles(XYZList const& list) {
        const auto n = list.size();
        const auto src = std::string("ExtraAnglesFile");
        for (std::size_t i = 0; i < n; ++i) {
            base::process(Vec3dAndOrigin{toEigen(list[i]), DataOrigin{src, i + 1}});
        }
    }
};
struct HorizontalExtremaComparator {
#if 1
    /// This ought to be equivalent, but faster (and marginally more accurate) than the alternative.
    bool operator()(Vec3dAndOrigin const& lhs, Vec3dAndOrigin const& rhs) const {
        return tanAboutY(lhs.vec) < tanAboutY(rhs.vec);
    }
#else
    bool operator()(Vec3dAndOrigin const& lhs, Vec3dAndOrigin const& rhs) const {
        return rotationAboutY(lhs.vec) < rotationAboutY(rhs.vec);
    }

#endif
};
class ScreenHorizontalExtrema : public ViewExtrema<HorizontalExtremaComparator> {
  public:
    ScreenHorizontalExtrema() : ViewExtrema<HorizontalExtremaComparator>(HorizontalExtremaComparator()) {}

    Eigen::Vector3d const& getLeft() const { return getMin().vec; }
    DataOrigin const& getLeftOrigin() const { return getMin().origin; }
    Eigen::Vector3d const& getRight() const { return getMax().vec; }
    DataOrigin const& getRightOrigin() const { return getMax().origin; }

    double angularRangeRadians() const { return rotationAboutY(getRight()) - rotationAboutY(getLeft()); }

    double dx() const { return (getRight() - getLeft()).x(); }
    double dz() const { return (getRight() - getLeft()).z(); }

    void debugAnglePrint(bool verbose, const char* when) const {
        if (!verbose) {
            return;
        }
        std::cerr << "[" << when << "] Horizontal angular range: " << radToDegree(angularRangeRadians()) << std::endl;
        std::cerr << "Screen left: " << to_string(getLeft()) << " from " << getLeftOrigin() << std::endl;
        std::cerr << "Screen right: " << to_string(getRight()) << " from " << getRightOrigin() << std::endl;
    }
};

static double computeYMagnitude(Plane const& screenPlane, Eigen::Vector3d const& xyz) {
    return std::abs(screenPlane.projection(xyz).y());
}

static double computeYMagnitude(Plane const& screenPlane, XYZ const& xyz) {
    return computeYMagnitude(screenPlane, toEigen(xyz));
}
struct VerticalProjectionComparator {
    VerticalProjectionComparator(Plane const& screenPlane) : screenPlane_(screenPlane) {}
    bool operator()(Vec3dAndOrigin const& lhs, Vec3dAndOrigin const& rhs) const {
        return getProjectedY(lhs.vec) < getProjectedY(rhs.vec);
    }

  private:
    double getProjectedY(Eigen::Vector3d const& vec) const { return screenPlane_.projection(vec).y(); }
    Plane screenPlane_;
};
class ScreenVerticalExtrema : public ViewExtrema<VerticalProjectionComparator> {
  public:
    using base = ViewExtrema<VerticalProjectionComparator>;

    ScreenVerticalExtrema(Plane const& screenPlane)
        : base(VerticalProjectionComparator(screenPlane)), screenPlane_(screenPlane) {}

    Vec3dAndOrigin const& getMaxMagnitudeObj() const {
        return (computeMagnitude(getMin().vec) > computeMagnitude(getMax().vec)) ? getMin() : getMax();
    }
    double getMaxMagnitude() const { return std::max(computeMagnitude(getMin().vec), computeMagnitude(getMax().vec)); }
    DataOrigin const& getMaxMagnitudeOrigin() const { return getMaxMagnitudeObj().origin; }
    void debugAnglePrint(bool verbose, const char* when) const {
        if (!verbose) {
            return;
        }
        auto& maxMag = getMaxMagnitudeObj();
        std::cerr << "[" << when << "] Maximum-magnitude Y projection: " << getMaxMagnitude() << std::endl;
        std::cerr << "Point: " << to_string(maxMag.vec) << " from " << maxMag.origin << std::endl;
    }

  private:
    double computeMagnitude(Eigen::Vector3d const& vec) const { return computeYMagnitude(screenPlane_, vec); }
    Plane const screenPlane_;
};

static ScreenHorizontalExtrema computeXScreenSpaceExtents(std::vector<NormalizedMeasurements> const& dataSets,
                                                          XYZList const& additionalPointsFromAngles, bool verbose) {

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

    ScreenHorizontalExtrema horizExtrema;
    horizExtrema.process(dataSets);
    horizExtrema.debugAnglePrint(verbose, "mappings only");
    horizExtrema.processExtraAngles(additionalPointsFromAngles);
    horizExtrema.debugAnglePrint(verbose, "full");
    return horizExtrema;
}

static Plane computeScreenPlane(ScreenHorizontalExtrema const& horizExtrema, bool verbose) {

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

    /// @todo Note this will always give us the same plane if we feed in field angles, since field angles are
    /// transformed into xyz by projecting onto a straight-ahead screen
    Eigen::Vector3d planeNormal = Eigen::Vector3d(-horizExtrema.dz(), 0, horizExtrema.dx()).normalized();
    auto screenPlane = Plane(planeNormal, horizExtrema.getRight());
    if (verbose) {
        std::cerr << "Plane of the screen A,B,C, D: " << PlaneA::get(screenPlane) << "," << PlaneB::get(screenPlane)
                  << "," << PlaneC::get(screenPlane) << ", " << PlaneC::get(screenPlane) << std::endl;
    }
    return screenPlane;
}

#if 0
static double computeMaximumYMagnitude(std::vector<NormalizedMeasurements> const& dataSets,
                                       XYZList const& additionalPointsFromAngles, bool verbose,
                                       Plane const& screenPlane) {
    //====================================================================
    // Figure out the Y screen-space extents.
    // The Y screen-space extents are symmetric and correspond to the lines parallel
    //  to the screen X axis that are within the plane of the X line specifying the
    //  axis extents at the largest magnitude angle up or down from the horizontal.
    // Find the highest-magnitude Y value of all points when they are
    // projected into the plane of the screen.
    double maxY = computeYMagnitude(screenPlane, ei::map(data.measurements.front().pointFromView));
    for (auto& meas : data.measurements) {
        auto yMag = computeYMagnitude(screenPlane, ei::map(meas.pointFromView));
        if (yMag > maxY) {
            maxY = yMag;
        }
    }
    if (verbose) {
        std::cerr << "Maximum-magnitude Y projection after just mappings: " << maxY << std::endl;
    }
    for (auto& extraPoint : additionalPointsFromAngles) {
        auto yMag = computeYMagnitude(screenPlane, extraPoint);
        if (yMag > maxY) {
            maxY = yMag;
        }
    }
    if (verbose) {
        std::cerr << "Maximum-magnitude Y projection: " << maxY << std::endl;
    }
    return maxY;
}
#else
static double computeMaximumYMagnitude(std::vector<NormalizedMeasurements> const& dataSets,
                                       XYZList const& additionalPointsFromAngles, bool verbose,
                                       Plane const& screenPlane) {
    //====================================================================
    // Figure out the Y screen-space extents.
    // The Y screen-space extents are symmetric and correspond to the lines parallel
    //  to the screen X axis that are within the plane of the X line specifying the
    //  axis extents at the largest magnitude angle up or down from the horizontal.
    // Find the highest-magnitude Y value of all points when they are
    // projected into the plane of the screen.
    ScreenVerticalExtrema extrema(screenPlane);
    extrema.process(dataSets);

    extrema.debugAnglePrint(verbose, "mappings only");
    extrema.processExtraAngles(additionalPointsFromAngles);
    extrema.debugAnglePrint(verbose, "full");
    return extrema.getMaxMagnitude();
}
#endif

struct HorizontalOutputs {
    /// left clipping plane at unit distance
    double left;
    /// right clipping plane at unit distance
    double right;

    /// @name Parameters used by the v1 schema model
    /// @{
    /// horizontal FOV in radians
    double hFOVradians;
    /// normalized (to [0, 1]) horizontal center of projection
    double xCOP;
    /// @}
};

static inline double projectToYPlaneAndZDivide(Eigen::Vector3d const& vec) {
    // get these at unit distance (if not there already)
    // which makes them rather like unit-distance left and right clipping planes
    return vec.x() / (-vec.z());
}

static HorizontalOutputs computeHFOV(ScreenHorizontalExtrema const& xExtrema, Plane const& screenPlane, bool verbose) {
    HorizontalOutputs ret;
    ret.left = projectToYPlaneAndZDivide(xExtrema.getLeft());
    ret.right = projectToYPlaneAndZDivide(xExtrema.getRight());
    if (verbose) {
        std::cerr << "Left clip: " << ret.left << std::endl;
        std::cerr << "Right clip: " << ret.right << std::endl;
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
    auto leftHalfFOV = std::atan(-ret.left);
    auto rightHalfFOV = std::atan(ret.right);
    if (verbose) {
        std::cerr << "Left half-FOV: " << radToDegree(leftHalfFOV) << std::endl;
        std::cerr << "Right half-FOV: " << radToDegree(rightHalfFOV) << std::endl;
    }
    ret.hFOVradians = leftHalfFOV + rightHalfFOV;
    auto screenWidth = -ret.left + ret.right;
    if (verbose) {
        std::cerr << "Screen width: " << screenWidth << std::endl;
        double hFOVDegrees = radToDegree(ret.hFOVradians);
        std::cerr << "Horizontal field of view (degrees): " << hFOVDegrees << std::endl;
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
    double xCOPOnPlane = screenPlane.projection(Eigen::Vector3d::Zero()).x();
    ret.xCOP = (xCOPOnPlane - ret.left) / screenWidth;

    return ret;
}

/// @returns FOV in radians
double computeVFOV(Plane const& screenPlane, double maxY, bool verbose) {
    //====================================================================
    // Figure out the monocular vertical field of view for the screen.
    // The FOV is twice the arctangent of half of the Y
    // distance divided by the distance to the screen.
    const auto vFOVRadians = 2 * std::atan(maxY / std::abs(PlaneD::get(screenPlane)));
    if (verbose) {
        auto vFOVDegrees = radToDegree(vFOVRadians);
        std::cerr << "Vertical field of view (degrees): " << vFOVDegrees << std::endl;
    }
    return vFOVRadians;
}

/// @returns so-called "overlap percent" (in [0, 100])
static double findOverlapPercent(Plane const& screenPlane, double hFOVradians, bool verbose) {

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
    double angleRadians = std::abs(std::atan2(PlaneA::get(screenPlane), PlaneC::get(screenPlane)));
    if (verbose) {
        std::cerr << "Angle degrees: " << radToDegree(angleRadians) << std::endl;
    }
    auto overlapFrac = 1. - 2. * angleRadians / hFOVradians;
    auto overlapPercent = overlapFrac * 100;
    if (verbose) {
        std::cerr << "Overlap percent: " << overlapPercent << std::endl;
    }
    return overlapPercent;
}

bool findScreen(ProjectionDescription& outProjection, ScreenDetails& outScreen,
                std::vector<NormalizedMeasurements> const& dataSets, XYZList const& additionalPointsFromAngles,
                bool verbose) {
    if (dataSets.empty()) {
        std::cerr << "findScreen(): Error: No channels of data!" << std::endl;
        return false;
    }
    for (auto& data : dataSets) {
        if (data.empty()) {
            std::cerr << "findScreen(): Error: No points in mapping" << std::endl;
            return false;
        }
    }
    const auto horizExtrema = computeXScreenSpaceExtents(dataSets, additionalPointsFromAngles, verbose);

    if (horizExtrema.angularRangeRadians() >= MY_PI) {
        std::cerr << "findScreen(): Error: Field of view > 180 degrees: found "
                  << radToDegree(horizExtrema.angularRangeRadians()) << std::endl;
        return false;
    }

    const auto screenPlane = computeScreenPlane(horizExtrema, verbose);

    const auto horizData = computeHFOV(horizExtrema, screenPlane, verbose);

    const auto maxY = computeMaximumYMagnitude(dataSets, additionalPointsFromAngles, verbose, screenPlane);

    const auto vFOVRadians = computeVFOV(screenPlane, maxY, verbose);
    // enforced by method of computing VFOV using max magnitude.
    auto yCOP = 0.5;

    if (verbose) {
        std::cerr << "Center of projection x,y: " << horizData.xCOP << ", " << yCOP << std::endl;
    }

    //====================================================================
    // Fill in the screen parameters.
    ei::map(outScreen.screenLeft) = horizExtrema.getLeft();
    ei::map(outScreen.screenRight) = horizExtrema.getRight();
    outScreen.screenPlane = screenPlane;
    outProjection.hFOVDegrees = radToDegree(horizData.hFOVradians);
    outProjection.vFOVDegrees = radToDegree(vFOVRadians);
    outProjection.overlapPercent = findOverlapPercent(screenPlane, horizData.hFOVradians, verbose);
    outProjection.cop = Point2d{horizData.xCOP, yCOP};
    return true;
}

bool findMesh(const std::vector<Mapping>& mapping, ScreenDescription const& screen, MeshDescription& mesh,
              bool /*verbose*/) {
    if (mapping.empty()) {
        std::cerr << "findMesh(): Error: No points in mapping" << std::endl;
        return false;
    }

    const double& A = screen.A;
    const double& B = screen.B;
    const double& C = screen.C;
    const double& D = screen.D;

    //====================================================================
    // Map each incoming mesh coordinate into the corresponding output
    // coordinate, storing them into the output mesh.
    mesh.clear();

    // Scale and offset to apply to the points projected onto the plane
    // to convert their values into normalized screen coordinates.  This
    // checks the assumption that the screen has some width in the X
    // coordinate (not rotated 90 degrees) and uses only the X value to
    // scale X.
    const XYZ& leftProj = screen.screenLeft;
    const XYZ& rightProj = screen.screenRight;
    if (leftProj.x == rightProj.x) {
        std::cerr << "Error computing mesh: screen has no X extent" << std::endl;
        return false;
    }
    double xOutOffset = -leftProj.x;
    double xOutScale = 1 / (rightProj.x - leftProj.x);
    double yOutOffset = screen.maxY; // Negative of negative maxY is maxY
    double yOutScale = 1 / (2 * screen.maxY);

    for (const auto& i : mapping) {

        // Input point coordinates are already normalized.
        double xNormIn = i.xyLatLong.x;
        double yNormIn = i.xyLatLong.y;
        std::array<double, 2> in;
        in[0] = xNormIn;
        in[1] = yNormIn;

        // Project the 3D points back into the plane of the screen and determine
        // the normalized coordinates in the coordinate system with the lower left
        // corner at (0,0) and the upper right at (1,1).  Because we oversized the
        // screen, these will all be in this range.  Otherwise, they might not be.
        XYZ onScreen = i.xyz.projectOntoPlane(A, B, C, D);
        double xNormOut = (onScreen.x + xOutOffset) * xOutScale;
        double yNormOut = (onScreen.y + yOutOffset) * yOutScale;
        std::array<double, 2> out;
        out[0] = xNormOut;
        out[1] = yNormOut;

        // Add a new entry onto the mesh
        std::array<std::array<double, 2>, 2> element;
        element[0] = in;
        element[1] = out;

        mesh.push_back(element);
    }

    return true;
}

MeshDescription findMesh(const NormalizedMeasurements& data, ScreenDetails const& screen, bool /*verbose*/) {
    MeshDescription ret;
    if (data.empty()) {
        std::cerr << "findMesh(): Error: No points in mapping!" << std::endl;
        return ret;
    }

    //====================================================================
    // Map each incoming mesh coordinate into the corresponding output
    // coordinate, storing them into the output mesh.

    // Scale and offset to apply to the points projected onto the plane
    // to convert their values into normalized screen coordinates.  This
    // checks the assumption that the screen has some width in the X
    // coordinate (not rotated 90 degrees) and uses only the X value to
    // scale X.
    Eigen::Vector3d leftProj(screen.screenLeft[0], 0, screen.screenLeft[2]);
    Eigen::Vector3d rightProj(screen.screenRight[0], 0, screen.screenRight[2]);
    if (leftProj.x() == rightProj.x()) {
        std::cerr << "Error computing mesh: screen has no X extent" << std::endl;
        return ret;
    }

    // Negative of negative maxY is maxY
    Eigen::Array2d outOffset = Eigen::Array2d(-leftProj.x(), screen.maxY);
    double xOutScale = 1 / (rightProj.x() - leftProj.x());
    double yOutScale = 1 / (2 * screen.maxY);
    Eigen::Array2d outScale = Eigen::Array2d(xOutScale, yOutScale);

    for (const auto& meas : data.measurements) {

        // Input (screen) point coordinates are already normalized.

        // Project the 3D points back into the plane of the screen and determine
        // the normalized coordinates in the coordinate system with the lower left
        // corner at (0,0) and the upper right at (1,1).  Because we oversized the
        // screen, these will all be in this range.  Otherwise, they might not be.
        Eigen::Vector3d onScreen = screen.screenPlane.projection(ei::map(meas.pointFromView));

        Point2d out;
        ei::mapArray(out) = (onScreen.head<2>().array() + outOffset) * outScale;

        // Add a new entry onto the mesh
        ret.push_back(MeshDescriptionRow{meas.screen, out});
    }

    return ret;
}

static void normalize(std::array<double, 2>& v) {
    double len = std::sqrt(v[0] * v[0] + v[1] * v[1]);
    if (len > 0) {
        v[0] /= len;
        v[1] /= len;
    }
}

/// Finds out whether the neighbor of the specified index
/// in the mapping violates the strictures of the
/// remove_invalid_points_based_on_angle function.
/// @return True if the neighbor angle difference is too large.
static bool neighbor_error(const std::vector<Mapping>& mapping, size_t index0, size_t index1, double xx, double xy,
                           double yx, double yy, double minDotProduct) {
    // Map the difference in angle space between the point and
    // its neighbor into screen space.
    std::array<double, 2> angleVec;
    angleVec[0] = mapping[index1].xyLatLong.longitude - mapping[index0].xyLatLong.longitude;
    angleVec[1] = mapping[index1].xyLatLong.latitude - mapping[index0].xyLatLong.latitude;
    std::array<double, 2> screenMappedVec;
    screenMappedVec[0] = angleVec[0] * xx + angleVec[1] * yx;
    screenMappedVec[1] = angleVec[0] * xy + angleVec[1] * yy;

    // Find the screen-space difference between the point and its
    // neighbor.
    std::array<double, 2> screenVec;
    screenVec[0] = mapping[index1].xyLatLong.x - mapping[index0].xyLatLong.x;
    screenVec[1] = mapping[index1].xyLatLong.y - mapping[index0].xyLatLong.y;

    // Normalize the two vectors
    normalize(screenMappedVec);
    normalize(screenVec);

    // Find the dot product between the two vectors.
    double dotProduct = screenMappedVec[0] * screenVec[0] + screenMappedVec[1] * screenVec[1];

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
static size_t neighbor_errors(const std::vector<Mapping>& mapping, size_t index, double xx, double xy, double yx,
                              double yy, double minDotProduct) {
    // Sort all points other than the one we're looking for
    // based on distance from the one we are looking for.
    typedef std::multimap<double, size_t> PointDistanceIndexMap;
    PointDistanceIndexMap map;
    for (size_t i = 0; i < mapping.size(); i++) {
        if (i != index) {
            map.insert(
                std::make_pair(point_distance(mapping[index].xyLatLong.longitude, mapping[index].xyLatLong.latitude,
                                              mapping[i].xyLatLong.longitude, mapping[i].xyLatLong.latitude),
                               i));
        }
    }

    // Check the first 8 points on the list (if there are that many)
    // and count up how many violate the angle condition.
    size_t ret = 0;
    PointDistanceIndexMap::const_iterator it = map.begin();
    for (size_t i = 0; (i < map.size()) && (i < 8); i++) {
        if (neighbor_error(mapping, index, it->second, xx, xy, yx, yy, minDotProduct)) {
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
static size_t find_index_of_angle_worst_offender(std::vector<Mapping>& mapping, double xx, double xy, double yx,
                                                 double yy, double minDotProduct) {
    size_t worstIndex = 0;
    size_t worstCount = neighbor_errors(mapping, 0, xx, xy, yx, yy, minDotProduct);
    for (size_t i = 1; i < mapping.size(); i++) {
        size_t count = neighbor_errors(mapping, i, xx, xy, yx, yy, minDotProduct);
        if (count > worstCount) {
            worstCount = count;
            worstIndex = i;
        }
    }
    if (worstCount == 0) {
        return mapping.size();
    }
    return worstIndex;
}

int remove_invalid_points_based_on_angle(std::vector<Mapping>& mapping, double xx, double xy, double yx, double yy,
                                         double maxAngleDegrees) {
    int ret = 0; // No points yet removed from the mesh.

    // Find the dot product associated with two unit vectors
    // separated by the angle specified.  This is the cosine
    // of the angle.
    double minDotProduct = std::cos(maxAngleDegrees / 180.0 * MY_PI);

    // We remove the worst offender from the list each time,
    // then re-start.  Assuming that we get the actual outlier,
    // as opposed to one of its neighbors, this avoids trimming
    // too many points from the vector.
    bool foundOutlier;
    do {
        foundOutlier = false;
        size_t off = find_index_of_angle_worst_offender(mapping, xx, xy, yx, yy, minDotProduct);
        if (off < mapping.size()) {
            mapping.erase(mapping.begin() + off);
            foundOutlier = true;
            ret++;
        }
    } while (foundOutlier);

    return ret;
}

XYZ reflect(XYZ input) {
    input.x *= -1;
    return input;
}

XYLatLong reflect(XYLatLong input) {
    input.longitude *= -1;
    input.x *= -1;
    return input;
}

XYLatLong reflect_normalized(XYLatLong input) {
    input.longitude *= -1;
    input.x = 1 - input.x;
    return input;
}

Mapping reflect(Mapping const& entry) {
    Mapping ret;
    ret.xyLatLong = reflect(entry.xyLatLong);
    ret.xyz = reflect(entry.xyz);
    return ret;
}

Mapping reflect_normalized(Mapping const& entry) {
    Mapping ret;
    ret.xyLatLong = reflect_normalized(entry.xyLatLong);
    ret.xyz = reflect(entry.xyz);
    return ret;
}

std::vector<Mapping> reflect_mapping(std::vector<Mapping> const& mapping) {
    std::vector<Mapping> ret;
    ret.reserve(mapping.size());
    for (auto& thisMapping : mapping) {
        ret.push_back(reflect(thisMapping));
    }

    return ret;
}

std::vector<Mapping> reflect_normalized_mapping(std::vector<Mapping> const& mapping) {
    std::vector<Mapping> ret;
    ret.reserve(mapping.size());
    for (auto& thisMapping : mapping) {
        ret.push_back(reflect_normalized(thisMapping));
    }

    return ret;
}

XYZList reflectPoints(XYZList const& input) {
    XYZList ret;
    ret.reserve(input.size());
    for (auto& pt : input) {
        ret.push_back(reflect(pt));
    }
    return ret;
}
