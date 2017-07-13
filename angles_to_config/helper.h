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

#include "types.h"
#include <iostream>
#include <vector>

/// Reads the four-column whitespace-delimited mapping file.
/// Returns empty mapping if it fails to read anything.
std::vector<Mapping> read_from_infile(std::istream& in);

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
///   @return -1 on error, the number of points that were
/// removed from the mesh otherwise.
int remove_invalid_points_based_on_angle(std::vector<Mapping>& mapping, double xx, double xy, double yx, double yy,
                                         double maxAngleDegrees);
/// Converts angles (in degrees) to eye-space 3d points.
XYZList convertAdditionalAngles(std::vector<LongLat> const& additionalAngles, double depth, bool useFieldAngles);

bool convert_to_normalized_and_meters(std::vector<Mapping>& mapping, double toMeters, double depth, double left,
                                      double bottom, double right, double top, bool useFieldAngles = false);

NormalizedMeasurements convert_to_normalized_and_meters(InputMeasurements const& input, double toMeters, double depth,
                                                        RectBoundsd rect, bool useFieldAngles = false);

bool findScreen(ScreenDescription& outScreen, const std::vector<Mapping>& mapping,
                XYZList const& additionalPointsFromAngles, bool verbose = false);

/// Compatibility wrapping.
inline bool findScreen(const std::vector<Mapping>& mapping, ScreenDescription& outScreen, bool verbose = false) {
    return findScreen(outScreen, mapping, XYZList(), verbose);
}

/// Output from find_screen that is only needed by the mesh computation.
struct ScreenDetails {
    ScreenDetails() = default;
    ScreenDetails(Plane const& scrPlane, Eigen::Vector3d const& left, Eigen::Vector3d const& right, double maxYMagnitude);
    bool valid = false;

    Plane screenPlane;               //!< Ax + By + Cz + D = 0 screen plane
    Eigen::Vector3d screenLeft, screenRight; //!< Left-most and right-most points on screen
    double maxY;                     //!< Maximum absolute value of Y for points on screen

    Eigen::Array2d offset;
    Eigen::Array2d scale;
    Eigen::Array2d projectAndNormalize(Point3d const& angleViewPoint) const;
};

bool findScreen(ProjectionDescription& outProjection, ScreenDetails& outScreen,
                std::vector<NormalizedMeasurements> const& dataSets, XYZList const& additionalPointsFromAngles,
                bool verbose = false);

bool findMesh(const std::vector<Mapping>& mapping, ScreenDescription const& screen, MeshDescription& mesh,
              bool verbose = false);

MeshDescription findMesh(const NormalizedMeasurements& data, ScreenDetails const& screen, bool verbose = false);

/// Reflect a point around X=0
XYZ reflect(XYZ input);

/// Reflect a (non-normalized) screen location and angles around X=0
XYLatLong reflect(XYLatLong input);

/// Reflect a (normalized) screen location and angles around X=0
XYLatLong reflect_normalized(XYLatLong input);

/// Produce a mapping that is reflected around X=0 in both angles and
/// screen coordinates.
Mapping reflect(Mapping const& entry);

/// Same as above except for mappings that have already gone through convert_to_normalized_and_meters()
Mapping reflect_normalized(Mapping const& entry);

/// Produce a mapping that is reflected around X=0 in both angles and
/// screen coordinates.
std::vector<Mapping> reflect_mapping(std::vector<Mapping> const& mapping);

/// Same as reflect_mapping except for mappings that have already gone through convert_to_normalized_and_meters()
std::vector<Mapping> reflect_normalized_mapping(std::vector<Mapping> const& mapping);

/// Reflect a list of points around X=0
XYZList reflectPoints(XYZList const& input);

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
