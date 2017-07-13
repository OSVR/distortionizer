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

// Library includes
#include <Eigen/Core>
#include <Eigen/Geometry>

// Standard includes
#include <array>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

// Global constants and variables
#define MY_PI (4.0 * std::atan(1.0))

template <typename T> struct RectBounds {
    using value_type = T;
    value_type left;
    value_type right;
    value_type top;
    value_type bottom;
    RectBounds<T> reflectedHorizontally() const { return RectBounds<T>{-right, -left, top, bottom}; }
};
using RectBoundsd = RectBounds<double>;

using Point2d = std::array<double, 2>;
using Point3d = std::array<double, 3>;
namespace ei {
/// Function to map std::array as an Eigen vector type
template <typename Scalar, size_t Dim>
inline Eigen::Map<Eigen::Matrix<Scalar, Dim, 1> > map(std::array<Scalar, Dim>& v) {
    return Eigen::Matrix<Scalar, Dim, 1>::Map(v.data());
}
/// Function to map constant std::array as an Eigen vector type
template <typename Scalar, size_t Dim>
inline Eigen::Map<Eigen::Matrix<Scalar, Dim, 1> const> map(std::array<Scalar, Dim> const& v) {
    return Eigen::Matrix<Scalar, Dim, 1>::Map(v.data());
}

/// Function to map std::array as an Eigen array-vector type
template <typename Scalar, size_t Dim>
inline Eigen::Map<Eigen::Array<Scalar, Dim, 1> > mapArray(std::array<Scalar, Dim>& v) {
    return Eigen::Array<Scalar, Dim, 1>::Map(v.data());
}
/// Function to map constant std::array as an Eigen array-vector type
template <typename Scalar, size_t Dim>
inline Eigen::Map<Eigen::Array<Scalar, Dim, 1> const> mapArray(std::array<Scalar, Dim> const& v) {
    return Eigen::Array<Scalar, Dim, 1>::Map(v.data());
}
} // namespace ei

/// Gentle wrapper around Point2d assigning longitude and latitude meaning (respectively) to the elements.
struct LongLat {
    Point2d longLat;
    /// angle in x
    double& longitude() { return longLat[0]; }
    /// angle in x (read-only)
    double longitude() const { return longLat[0]; }

    /// angle in y
    double& latitude() { return longLat[1]; }
    /// angle in y (read-only)
    double latitude() const { return longLat[1]; }

    /// Access underlying array for things like Eigen::Map
    double* data() { return longLat.data(); }
    /// Access underlying array for things like Eigen::Map (read-only)
    const double* data() const { return longLat.data(); }
#if 0
    /// Access as Eigen data structure.
    Eigen::Map<Eigen::Vector2d> ei() { return Eigen::Vector2d::Map(data()); }
    /// Access as Eigen data structure (constant)
    Eigen::Map<Eigen::Vector2d const> ei() const { return Eigen::Vector2d::Map(data()); }
#endif
};

struct Config {
    bool useRightEye = false;
    bool computeScreenBounds = true;
    RectBoundsd suppliedScreenBounds;
    bool useFieldAngles = true;
    double toMeters = 1.0;
    double depth = 2.0;
    bool verifyAngles = false;
    /// parameters to verify_angles
    double xx, xy, yx, yy, maxAngleDiffDegrees;

    bool verbose = false;
};

// Screen-space to/from angle-space map entry
class XYLatLong {
  public:
    double x = 0;
    double y = 0;
    double latitude = 0;
    double longitude = 0;

    XYLatLong(double px, double py, double plat, double plong) {
        x = px;
        y = py;
        latitude = plat;
        longitude = plong;
    }
    XYLatLong() { x = y = latitude = longitude = 0; }
};

struct DataOrigin {
    std::string inputSource;
    std::size_t lineNumber;
    bool known() const { return !inputSource.empty(); }
};
inline std::ostream& operator<<(std::ostream& os, DataOrigin const& orig) {
    if (orig.known()) {
        std::ostringstream oss;
        oss << orig.inputSource << ":" << orig.lineNumber;
        os << oss.str();
    } else {
        os << "(unknown)";
    }
    return os;
}

struct InputMeasurements;
struct InputMeasurement {
    /// In arbitrary units
    Point2d screen;

#if 0
    using ScreenEigenType = Eigen::Vector2d;
    Eigen::Map<ScreenEigenType> screenEi() { return ScreenEigenType::Map(screen.data()); }
    Eigen::Map<ScreenEigenType const> screenEi() const { return ScreenEigenType::Map(screen.data()); }
#endif

    /// in degrees
    LongLat viewAnglesDegrees;
#if 0
    using LongLatEigenType = Eigen::Vector2d;
    Eigen::Map<LongLatEigenType> viewAnglesEi() { return LongLatEigenType::Map(viewAnglesDegrees.data()); }
    Eigen::Map<LongLatEigenType const> screenEi() const { return LongLatEigenType::Map(viewAnglesDegrees.data()); }
#endif

    std::size_t lineNumber;

    DataOrigin getOrigin(InputMeasurements const& parent) const;
};

struct InputMeasurements {
    std::string inputSource;
    std::vector<InputMeasurement> measurements;
    bool empty() const { return measurements.empty(); }
    std::size_t size() const { return measurements.size(); }
};

inline DataOrigin InputMeasurement::getOrigin(InputMeasurements const& parent) const {
    return DataOrigin{parent.inputSource, lineNumber};
}

struct NormalizedMeasurements;
struct NormalizedMeasurement {
    /// Normalized screen units, in [0, 1]
    Point2d screen;

    /// In arbitrary units in 3d space, based on the view angles from the corresponding input measurement.
    Point3d pointFromView;
    std::size_t lineNumber;
    DataOrigin getOrigin(NormalizedMeasurements const& parent) const;
};

struct NormalizedMeasurements {
    std::string inputSource;
    std::vector<NormalizedMeasurement> measurements;
    bool empty() const { return measurements.empty(); }
    std::size_t size() const { return measurements.size(); }
};

inline DataOrigin NormalizedMeasurement::getOrigin(NormalizedMeasurements const& parent) const {
    return DataOrigin{parent.inputSource, lineNumber};
}

// 3D coordinate
class XYZ {
  public:
    double x;
    double y;
    double z;

    XYZ(double px, double py, double pz) {
        x = px;
        y = py;
        z = pz;
    }
    XYZ() { x = y = z = 0; }

    /// Return the rotation about the Y axis, where 0 rotation points along
    /// the -Z axis and positive rotation heads towards the -X axis.
    /// The X axis in atan space corresponds to the -z axis in head space,
    /// and the Y axis in atan space corresponds to the -x axis in head space.
    double rotationAboutY() const { return std::atan2(-x, -z); }

    /// Project from the origin through our point onto a plane whose
    /// equation is specified.
    XYZ projectOntoPlane(double A, double B, double C, double D) const {
        XYZ ret;

        // Solve for the value of S that satisfies:
        //    Asx + Bsy + Csz + D = 0,
        //    s = -D / (Ax + By + Cz)
        // Then for the location sx, sy, sz.

        double s = -D / (A * x + B * y + C * z);
        ret.x = s * x;
        ret.y = s * y;
        ret.z = s * z;

        return ret;
    }

    /// Return the rotation distance from another point.
    double distanceFrom(const XYZ& p) const {
        return std::sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y) + (z - p.z) * (z - p.z));
    }

    void debugPrint(std::ostream& os) const {
        static const auto PRECISION = 4;
        static const auto WIDTH = PRECISION + 3;
        std::ostringstream ss;
        ss << std::setprecision(PRECISION);
        ss << "(" << std::setw(WIDTH) << x;
        ss << ", " << std::setw(WIDTH) << y;
        ss << ", " << std::setw(WIDTH) << z << ")";
        os << ss.str();
    }
};

using XYZList = std::vector<XYZ>;

/// Mapping entry, along with its associated 3D coordinate
class Mapping {
  public:
    /// eye/camera space
    XYLatLong xyLatLong;
    /// Screen/world space
    XYZ xyz;

    Mapping(XYLatLong const& ll, XYZ const& x) : xyLatLong(ll), xyz(x) {}
    Mapping() = default;
};

// Description of a screen
struct ScreenDescription {
    double hFOVDegrees;
    double vFOVDegrees;
    double overlapPercent;
    double xCOP;
    double yCOP;

    // These are quantities computed along the way to getting the
    // screen that are needed by the mesh calculations, so they
    // are stored in the screen to pass from the findScreen to
    // the findMesh functions.
    double A, B, C, D;           //!< Ax + By + Cz + D = 0 screen plane
    XYZ screenLeft, screenRight; //!< Left-most and right-most points on screen
    double maxY;                 //!< Maximum absolute value of Y for points on screen
};

using Plane = Eigen::Hyperplane<double, 3>;

struct PlaneA;
struct PlaneB;
struct PlaneC;
struct PlaneD;
namespace detail {
template <typename T> struct PlaneCoefficientIndexTrait;
template <> struct PlaneCoefficientIndexTrait<PlaneA> : std::integral_constant<std::size_t, 0> {};

template <> struct PlaneCoefficientIndexTrait<PlaneB> : std::integral_constant<std::size_t, 1> {};

template <> struct PlaneCoefficientIndexTrait<PlaneC> : std::integral_constant<std::size_t, 2> {};

template <> struct PlaneCoefficientIndexTrait<PlaneD> : std::integral_constant<std::size_t, 3> {};

template <typename Derived> struct PlaneAccessProxy {
  public:
    /// Const accessor for this coefficient of the plane equation.
    static double get(Plane const& p) { return p.coeffs()[PlaneCoefficientIndexTrait<Derived>()]; }
#if 0
        /// Reference accessor for this coefficient of the plane equation.
        static double& get(Plane& p) { return p.coeffs()[PlaneCoefficientIndexTrait<Derived>()]; }
#endif
};
} // namespace detail
#if 0
  /// Get a coefficient of a plane by name (using the PlaneA, PlaneB, PlaneC, PlaneD tag types)
template <typename PlaneComponent> inline double get(Plane const& p) {
    return p.coeffs()[detail::PlaneCoefficientIndexTrait<PlaneComponent>()];
}
template <typename PlaneComponent> inline double& get(Plane& p) {
    return p.coeffs()[detail::PlaneCoefficientIndexTrait<PlaneComponent>()];
}
#endif

struct PlaneA : detail::PlaneAccessProxy<PlaneA> {};
struct PlaneB : detail::PlaneAccessProxy<PlaneB> {};
struct PlaneC : detail::PlaneAccessProxy<PlaneC> {};
struct PlaneD : detail::PlaneAccessProxy<PlaneD> {};

inline Eigen::Vector3d toEigen(XYZ const& p) { return Eigen::Vector3d(p.x, p.y, p.z); }
inline XYZ toXYZ(Eigen::Vector3d const& p) { return XYZ{p.x(), p.y(), p.z()}; }
#if 0
inline Eigen::Vector3d projectOntoPlane(Plane const& plane, XYZ const& p) { return plane.projection(toEigen(p)); }
#endif

/// Output from find_screen that is used to generate the configuration
struct ProjectionDescription {
    double hFOVDegrees;
    double vFOVDegrees;
    double overlapPercent = 100.;
    /// Center of Projection
    Point2d cop = {0.5, 0.5};
};
/// Output from find_screen that is only needed by the mesh computation.
struct ScreenDetails {
    Plane screenPlane;               //!< Ax + By + Cz + D = 0 screen plane
    Point3d screenLeft, screenRight; //!< Left-most and right-most points on screen
    double maxY;                     //!< Maximum absolute value of Y for points on screen
};

using MeshDescriptionRow = std::array< //!< 2-vector of from, to coordinates
    std::array<double, 2>,             //!< 2-vector of unit coordinates (x,y)
    2>;
/// Holds a list of mappings from physical-display normalized
/// coordinates to canonical-display normalized coordinates.
typedef std::vector< //!< Vector of mappings
    MeshDescriptionRow>
    MeshDescription;

template <typename T> class InclusiveBounds {
  public:
    using value_type = T;
    InclusiveBounds() = default;
    InclusiveBounds(value_type minVal, value_type maxVal) : valid_(true), minVal_(minVal), maxVal_(maxVal) {
        if (maxVal_ < minVal_) {
            using std::swap;
            swap(minVal_, maxVal_);
        }
    }
    explicit operator bool() const { return valid_; }
    bool contains(value_type val) const { return (!valid_) || (val >= minVal_ && val <= maxVal_); }
    bool outside(value_type val) const { return valid_ && (val < minVal_ || val > maxVal_); }

    value_type getMin() const { return minVal_; }
    value_type getMax() const { return maxVal_; }

  private:
    bool valid_ = false;
    value_type minVal_;
    value_type maxVal_;
};
template <typename T> inline std::ostream& operator<<(std::ostream& os, InclusiveBounds<T> const& bounds) {
    std::ostringstream oss;
    oss << "[";
    if (bounds) {
        const auto oldFlags = oss.flags();
        const auto streamFlags = os.flags();
        oss.flags(streamFlags);
        oss.precision(os.precision());
        oss << bounds.getMin();
        oss.flags(oldFlags);
        oss << ", ";
        oss.precision(os.precision());
        oss.flags(streamFlags);
        oss << bounds.getMax();
        oss.flags(oldFlags);
    } else {
        oss << "unbounded";
    }
    oss << "]";
    os << oss.str();
    return os;
}

using InclusiveBoundsd = InclusiveBounds<double>;
using InclusiveBoundsf = InclusiveBounds<float>;

template <typename T> struct XYInclusiveBounds {
    // Do we have any bounds?
    explicit operator bool() const { return static_cast<bool>(x) || static_cast<bool>(y); }
    InclusiveBounds<T> x;
    InclusiveBounds<T> y;
};

template <typename T> inline std::ostream& operator<<(std::ostream& os, XYInclusiveBounds<T> const& xyBounds) {
    if (!xyBounds) {
        os << "unbounded";
        return os;
    }
    if (xyBounds.x) {
        os << "x: " << xyBounds.x;
    }
    if (xyBounds.x && xyBounds.y) {
        os << ", ";
    }
    if (xyBounds.y) {
        os << "y: " << xyBounds.y;
    }
    return os;
}

using XYInclusiveBoundsd = XYInclusiveBounds<double>;
using XYInclusiveBoundsf = XYInclusiveBounds<float>;
