/** @file
    @brief Header

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

#ifndef INCLUDED_Subproblems_h_GUID_07251D99_D976_4FBA_E5CF_5512B0598BE8
#define INCLUDED_Subproblems_h_GUID_07251D99_D976_4FBA_E5CF_5512B0598BE8

// Internal Includes
#include "EigenStdArrayInterop.h"
#include "types.h"

// Library/third-party includes
// - none

// Standard includes
#include <cassert>

/// Reads the four-column whitespace-delimited mapping file.
/// Returns empty vector if it fails to read anything.
InputMeasurements readInputMeasurements(std::string const& inputSource, std::istream& in);

/// Reads a two-column whitespace-delimited file of additional angles
/// Returns empty if it fails to read anything.
std::vector<LongLat> readAdditionalAngles(std::istream& in);

/// This removes invalid points from the mesh if the angle
/// between the vector from a point to its neighbor in lat/long
/// space (when transformed by the specified mapping into screen
/// space) is more than the maximum specified number of degrees
/// different than the vector to that same neighbor in screen
/// space.  The mapping: +longitude (left) points in (xx, xy) in screen
/// space and +latitude (up) points in (yx, yy) in screen space.
///
/// @param xxxy is a unit vector in screen space in the direction of +longitude (to the right)
/// @param yxyy is a unit vector in screen space in the direction of +latitude (up)
///
/// @return -1 on error, the number of points that were removed from the mesh otherwise.
int remove_invalid_points_based_on_angle(InputMeasurements& input, double maxAngleDegrees,
                                         Point2d const& xxxy = {1., 0.}, Point2d const& yxyy = {0., 1},
                                         bool verbose = false);

NormalizedMeasurements convert_to_normalized_and_meters(InputMeasurements const& input, double toMeters, double depth,
                                                        RectBoundsd rect, bool useFieldAngles = false);

/// Converts angles (in degrees) to eye-space 3d points.
XYZList convertAdditionalAngles(std::vector<LongLat> const& additionalAngles, double depth, bool useFieldAngles);

/// Output from find_screen that is only needed by the mesh computation.
struct ScreenDetails {
    ScreenDetails() = default;
    ScreenDetails(Plane const& scrPlane, Eigen::Vector3d const& left, Eigen::Vector3d const& right,
                  double maxYMagnitude);
    bool valid = false;

    Plane screenPlane; //!< Ax + By + Cz + D = 0 screen plane
    Eigen::Vector3d screenOrigin;
    Eigen::Vector3d screenYBasis;
    Eigen::Vector3d screenXBasis;
    Eigen::Array2d scale;
    Eigen::Array2d projectAndNormalize(Point3d const& angleViewPoint, bool verbose = false) const;
};

bool findScreen(ProjectionDescription& outProjection, ScreenDetails& outScreen,
                std::vector<NormalizedMeasurements> const& dataSets, XYZList const& additionalPointsFromAngles,
                bool verbose = false);

MeshDescription findMesh(const NormalizedMeasurements& data, ScreenDetails const& screen, bool verbose = false);

template <typename ValueType, typename Comparator> class GenericExtremaFinder {
  public:
    using comparator_type = Comparator;
    using value_type = ValueType;
    GenericExtremaFinder() : compare_() {}
    GenericExtremaFinder(comparator_type&& compare) : compare_(std::move(compare)) {}
    GenericExtremaFinder(comparator_type const& compare) : compare_(compare) {}

    value_type const& getMin() const {
        assert(valid_ && "Only makes sense to get min value if you've actually processed any elements!");
        return minVal_;
    }
    value_type const& getMax() const {
        assert(valid_ && "Only makes sense to get max value if you've actually processed any elements!");
        return maxVal_;
    }
    void process(value_type const& val) {
        if (!valid_) {
            minVal_ = val;
            maxVal_ = val;
            valid_ = true;
            return;
        }

        if (compare_(val, minVal_)) {
            /// new is less than min
            minVal_ = val;
        }

        if (compare_(maxVal_, val)) {
            /// max is less than new
            maxVal_ = val;
        }
    }

  private:
    comparator_type compare_;
    bool valid_ = false;

    value_type minVal_;
    value_type maxVal_;
};

static inline double angleSquaredDistance(LongLat const& a, LongLat const& b) {
    return (ei::map(a.longLat) - ei::map(b.longLat)).squaredNorm();
}

#endif // INCLUDED_Subproblems_h_GUID_07251D99_D976_4FBA_E5CF_5512B0598BE8
