/** @file
    @brief Implementation

    @date 2017

    @author
    Sensics, Inc.
    <http://sensics.com/osvr>
*/

// Copyright 2017 Sensics, Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//        http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Internal Includes
#include "Subproblems.h"
#include "EigenStdArrayInterop.h"
#include "InvalidPointRemoval.h"
#include "helper.h"

// Library/third-party includes
#include <Eigen/Core>
#include <Eigen/Geometry>

// Standard includes
#include <iostream>
#include <tuple>

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

static inline double getNeighborError(InputMeasurement const& a, InputMeasurement const& b, Point2d const& xxxy,
                                      Point2d const& yxyy, bool verbose = false) {
    // Map the difference in angle space between the point and
    // its neighbor into screen space.
    Eigen::Vector2d angleVec = ei::map(a.viewAnglesDegrees.longLat) - ei::map(b.viewAnglesDegrees.longLat);
    Eigen::Vector2d screenMappedVec = angleVec[0] * ei::map(xxxy) + angleVec[1] * ei::map(yxyy);

    // Find the screen-space difference between the point and its
    // neighbor.
    Eigen::Vector2d screenVec = ei::map(a.screen) - ei::map(b.screen);

    if (verbose) {
        std::cerr << "\nangleVec:        " << angleVec.transpose()
                  << " normalized: " << angleVec.normalized().transpose() << std::endl;
        std::cerr << "screenMappedVec: " << screenMappedVec.transpose()
                  << " normalized: " << screenMappedVec.normalized().transpose() << std::endl;
        std::cerr << "screenVec:       " << screenVec.transpose()
                  << " normalized: " << screenVec.normalized().transpose() << std::endl;
        std::cerr << "dot product: " << screenMappedVec.normalized().dot(screenVec.normalized()) << std::endl;
    }
    // Normalize the two vectors
    // Find the dot product between the two vectors.
    // If the dot product is too small, the angle is too large
    // and so we fail.
    return screenMappedVec.normalized().dot(screenVec.normalized());
}
/// Finds out how many neighbors of the specified index
/// in the mapping violate the strictures of the
/// remove_invalid_points_based_on_angle function.
/// Checks up to 8 nearest neighbors in angle space.
/// @return How many neighbor angle differences are too large.
static std::pair<std::size_t, double> neighbor_errors(InputMeasurements const& input, size_t index, Point2d const& xxxy,
                                                      Point2d yxyy, double minDotProduct, bool verbose) {
    // Sort all points other than the one we're looking for
    // based on distance from the one we are looking for.
    // Start by just making an unsorted vector with them.
    using InputMeasurementDistanceIndex = std::pair<double, std::size_t>;
    using InputMeasurementDistances = std::vector<InputMeasurementDistanceIndex>;
    auto refMeas = input.measurements[index];
    InputMeasurementDistances dists;
    const auto n = input.size();
    dists.reserve(n - 1);
    for (size_t i = 0; i < n; ++i) {
        if (i == index) {
            continue;
        }
        auto dist = InputMeasurementDistanceIndex{
            angleSquaredDistance(input.measurements[i].viewAnglesDegrees, refMeas.viewAnglesDegrees), i};
        dists.push_back(dist);
    }

    // Now, we use nth_element to find exactly the 8th element and as a side effect sort those smaller elements ahead of
    // it.
    static const auto desired_element = 8;
    auto endIter = dists.end();
    if (dists.size() > desired_element) {
        auto midElt = dists.begin();
        std::advance(midElt, desired_element - 1);
        std::nth_element(dists.begin(), midElt, dists.end());

        /// Adjust the range we'll iterate through below so we only iterate thru the ones we've split.
        endIter = midElt;
        endIter++;
    }

    // Check the first 8 points on the list (if there are that many)
    // and count up how many violate the angle condition.
    double product = 1;
    std::size_t count = 0;
    std::for_each(dists.begin(), endIter, [&](InputMeasurementDistanceIndex const& elt) {
        auto err = getNeighborError(input.measurements[elt.second], refMeas, xxxy, yxyy, verbose);
        if (err < minDotProduct) {
            count++;
            product *= err;
        }
    });
    return std::make_pair(count, product);
}

/// Looks through the mapping to see if there is a point whose
/// neighbor angles violate the strictures of the
/// remove_invalid_points_based_on_angle function.
/// @return index of the point that has the largest number
/// of invalid neighbor angles if one is found, size of mapping
/// if not.
static size_t find_index_of_angle_worst_offender(InputMeasurements const& input, Point2d const& xxxy,
                                                 Point2d const& yxyy, double minDotProduct, bool verbose = false) {
    size_t worstIndex = 0;
    size_t worstCount = 0;
    double worstProduct = minDotProduct * 2; // arbitrarily too large to matter.
    const auto n = input.size();
    for (size_t i = 0; i < n; ++i) {
#if 1
        const bool childVerbose = false;
#else
        const bool childVerbose = verbose && (i == 0);
#endif
        size_t count;
        double product;
        std::tie(count, product) = neighbor_errors(input, i, xxxy, yxyy, minDotProduct, childVerbose);
        if (count > worstCount || (count != 0 && count == worstCount && product < worstProduct)) {
            worstCount = count;
            worstProduct = product;
            worstIndex = i;
        }
    }
    if (worstCount == 0) {
        return n;
    }
    if (verbose) {
        std::cerr << "Worst remaining: " << input.measurements[worstIndex].getOrigin(input) << " with " << worstCount
                  << " invalid neighbor angles." << std::endl;
        std::cerr << "\tproduct: " << worstProduct << std::endl;
    }
    return worstIndex;
}

int remove_invalid_points_based_on_angle(InputMeasurements& input, double maxAngleDegrees, Point2d const& xxxy,
                                         Point2d const& yxyy, bool verbose) {
//#define USE_OLD_REMOVE_INVALID
#ifdef USE_OLD_REMOVE_INVALID
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
        size_t off = find_index_of_angle_worst_offender(input, xxxy, yxyy, minDotProduct, verbose);
        if (off < input.size()) {
            input.measurements.erase(input.measurements.begin() + off);
            foundOutlier = true;
            ret++;
        }
    } while (foundOutlier);
#else
    detail::InvalidPointRemover remover(input, maxAngleDegrees, xxxy, yxyy, verbose);
    auto ret = remover();
#endif
    return ret;
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

XYZList convertAdditionalAngles(std::vector<LongLat> const& additionalAngles, double depth, bool useFieldAngles) {
    XYZList ret;

    for (auto& additionalLongLat : additionalAngles) {
        ret.push_back(longLatToWorldSpace(additionalLongLat, useFieldAngles, depth));
    }
    return ret;
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

ScreenDetails::ScreenDetails(Plane const& scrPlane, Eigen::Vector3d const& left, Eigen::Vector3d const& right,
                             double maxYMagnitude)
    : valid(true), screenPlane(scrPlane), screenYBasis(Eigen::Vector3d::UnitY()) {

    // Scale and offset to apply to the points projected onto the plane
    // to convert their values into normalized screen coordinates.  This
    // checks the assumption that the screen has some width in the X
    // coordinate (not rotated 90 degrees) and uses only the X value to
    // scale X.
    auto leftProjX = left[0];
    auto rightProjX = right[0];
    if (leftProjX == rightProjX) {
        std::cerr << "Error computing mesh: screen has no X extent" << std::endl;
        valid = false;
        return;
    }
    screenOrigin = left - Eigen::Vector3d(0, maxYMagnitude, 0);
    Eigen::Vector3d screenHorizExtent = right - left;
    screenXBasis = screenHorizExtent.stableNormalized();
    scale = Eigen::Array2d(1. / screenHorizExtent.norm(), 1. / (2. * maxYMagnitude));
}

Eigen::Array2d ScreenDetails::projectAndNormalize(Point3d const& angleViewPoint, bool verbose) const {
    Eigen::Vector3d onScreenPlane = screenPlane.projection(ei::map(angleViewPoint));
    Eigen::Vector3d onScreenPlaneAtOrigin = onScreenPlane - screenOrigin;
    Eigen::Array2d unscaled =
        Eigen::Array2d(onScreenPlaneAtOrigin.dot(screenXBasis), onScreenPlaneAtOrigin.dot(screenYBasis));

    if (verbose) {
        std::cerr << "Origin: " << screenOrigin.transpose() << std::endl;
        std::cerr << "Point: " << ei::map(angleViewPoint).transpose() << std::endl;
        std::cerr << "onScreenPlane: " << onScreenPlane.transpose() << std::endl;
        std::cerr << "onScreenPlaneAtOrigin: " << onScreenPlaneAtOrigin.transpose() << std::endl;
        std::cerr << "unscaled: " << unscaled.transpose() << std::endl;
    }
    Eigen::Array2d ret = unscaled * scale;
    return ret;
}

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

struct HorizontalOutputs {
    /// left clipping plane at unit distance
    double left;
    /// right clipping plane at unit distance
    double right;

    /// minimum x value point, projected to screen plane and y plane.
    Eigen::Vector3d leftPoint;
    /// maximum x value point, projected to screen plane and y plane.
    Eigen::Vector3d rightPoint;

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
#define USE_REAL_COP
#ifdef USE_REAL_COP
    ret.left = projectToYPlaneAndZDivide(xExtrema.getLeft());
    ret.right = projectToYPlaneAndZDivide(xExtrema.getRight());
    {
        Eigen::Vector3d pt = screenPlane.projection(xExtrema.getLeft());
        pt.y() = 0;
        ret.leftPoint = pt;
    }
    {
        Eigen::Vector3d pt = screenPlane.projection(xExtrema.getRight());
        pt.y() = 0;
        ret.rightPoint = pt;
    }

#else
    /// @todo temporary override to force COP=0.5
    auto left = projectToYPlaneAndZDivide(xExtrema.getLeft());
    auto right = projectToYPlaneAndZDivide(xExtrema.getRight());
    auto maxX = std::max(std::abs(left), std::abs(right));

    ret.left = -maxX;
    ret.right = maxX;
    assert(xExtrema.getLeft().z() == -1 && "Code below assumes that screen is at -1");
    ret.leftPoint = Eigen::Vector3d(-maxX, 0, -1);
    ret.rightPoint = Eigen::Vector3d(maxX, 0, -1);
#endif
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
#if 0
    double xCOPOnPlane = projectToYPlaneAndZDivide(screenPlane.projection(Eigen::Vector3d::Zero()));
#else
    double xCOPOnPlane = screenPlane.projection(Eigen::Vector3d::Zero()).x();
#endif
    ret.xCOP = (xCOPOnPlane - ret.left) / screenWidth;

    return ret;
}

/// @returns FOV in radians
double computeVFOV(Plane const& screenPlane, double maxY, bool verbose) {
    //====================================================================
    // Figure out the monocular vertical field of view for the screen.
    // The FOV is twice the arctangent of half of the Y
    // distance divided by the distance to the screen.
    const auto vFOVRadians = 2 * std::atan(maxY / std::abs(screenPlane.offset()));
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
    outScreen = ScreenDetails(screenPlane, horizData.leftPoint, horizData.rightPoint, maxY);
    if (!outScreen.valid) {
        std::cerr << "Screen details are not valid!" << std::endl;
        return false;
    }
    if (verbose) {
        std::cerr << "Screen origin: " << outScreen.screenOrigin.transpose() << std::endl;
        std::cerr << "Screen x basis: " << outScreen.screenXBasis.transpose() << std::endl;
        Point3d zero = {0., 0., 0.};
        Eigen::Array2d testCOP = outScreen.projectAndNormalize(zero, true);
        std::cerr << "Screen normalized CoP:" << testCOP.transpose() << std::endl;
    }

    outProjection.hFOVDegrees = radToDegree(horizData.hFOVradians);
    outProjection.vFOVDegrees = radToDegree(vFOVRadians);
    outProjection.overlapPercent = findOverlapPercent(screenPlane, horizData.hFOVradians, verbose);
    outProjection.cop = Point2d{horizData.xCOP, yCOP};
    return true;
}

MeshDescription findMesh(const NormalizedMeasurements& data, ScreenDetails const& screen, bool verbose) {
    MeshDescription ret;
    if (data.empty()) {
        std::cerr << "findMesh(): Error: No points in mapping!" << std::endl;
        return ret;
    }

    //====================================================================
    // Map each incoming mesh coordinate into the corresponding output
    // coordinate, storing them into the output mesh.

    for (const auto& meas : data.measurements) {

        // Input (screen) point coordinates are already normalized.

        // Project the 3D points back into the plane of the screen and determine
        // the normalized coordinates in the coordinate system with the lower left
        // corner at (0,0) and the upper right at (1,1).  Because we oversized the
        // screen, these will all be in this range.  Otherwise, they might not be.
        /// @todo How/where did we oversize the screen? (Y axis?)
        Point2d out;
        ei::mapArray(out) = screen.projectAndNormalize(meas.pointFromView);

        // Add a new entry onto the mesh
        ret.push_back(MeshDescriptionRow{meas.screen, out});
    }

    return ret;
}
