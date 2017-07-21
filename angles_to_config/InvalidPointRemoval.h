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

#ifndef INCLUDED_InvalidPointRemoval_h_GUID_1FF303B7_FC7F_448F_E416_A278BCAD6B0A
#define INCLUDED_InvalidPointRemoval_h_GUID_1FF303B7_FC7F_448F_E416_A278BCAD6B0A

// Internal Includes
#include "Subproblems.h"
#include "types.h"

// Library/third-party includes
// - none

// Standard includes
#include <algorithm>
#include <cstddef>
#include <tuple>
#include <vector>

namespace detail {
class InvalidPointRemover;
enum class MeasurementStatus { Active = 0, Removed = 1 };
/// Represents the distance between an (implied) element and the explicitly specified (here by index) measurement.
struct StatusDistanceIndex {
    StatusDistanceIndex(double dist, std::size_t i);
    MeasurementStatus status = MeasurementStatus::Active;

    double distance;
    std::size_t index;

    /// Memoizes the computation.
    double getError(InputMeasurements const& input, InputMeasurement const& refMeas, Eigen::Vector2d const& xxxy,
                    Eigen::Vector2d const& yxyy, bool verbose = false);

  private:
    bool populatedError = false;
    double error;
};
/// Comparison operator for StatusDistanceIndex: compares status first (with active < removed) so we tend to ignore
/// removed points.
static inline bool operator<(StatusDistanceIndex const& lhs, StatusDistanceIndex const& rhs) {
    return std::tie(lhs.status, lhs.distance, lhs.index) < std::tie(rhs.status, rhs.distance, rhs.index);
}

using StatusDistanceIndexList = std::vector<StatusDistanceIndex>;

class MeasurementNeighborhood {
  public:
    MeasurementNeighborhood(InvalidPointRemover const& parent, std::size_t i);
    /// own index.
    const std::size_t index;

    /// Is this index/measurement valid/active?
    bool isValid() const { return MeasurementStatus::Active == status; }

    void markIndexRemoved(std::size_t removedIndex);

    std::pair<std::size_t, double> getNeighborErrorCount(bool verbose);

  private:
    InvalidPointRemover const* const parent_;
    /// own status.
    MeasurementStatus status = MeasurementStatus::Active;

    /// called once during constructor
    void populateNeighborData(InputMeasurements const& input);

    using iterator = StatusDistanceIndexList::iterator;
    /// @param n How many neighbors do you want to consider for "bad ones" judgement?
    iterator findEndOfUpToNClosest(const std::size_t n = 8);

    std::pair<std::size_t, double> countNeighborErrors(bool verbose);

    bool dirty_ = true;

    /// This list is in no particular order...
    /// EXCEPT immediately following findEndOfUpToNWorst, in which case the first
    /// up to n elements are closest to this point.
    /// Note that not all of them may be unremoved, though in practice (since we compare first by status, and we
    /// shouldn't remove nearly every measurement) the removed ones should be after the split point, not before.
    StatusDistanceIndexList neighborData_;

    /// Memoized result of countNeighborErrors
    std::pair<std::size_t, double> neighborErrors_;
};

class InvalidPointRemover {
  public:
    InvalidPointRemover(InputMeasurements& input, double maxAngleDegrees, Point2d const& xxxyVec,
                        Point2d const& yxyyVec, bool verbose);
    /// @name For use by related (contained) objects
    /// @{
    const Eigen::Vector2d xxxy;
    const Eigen::Vector2d yxyy;
    const double minDotProduct;
    InputMeasurements const& getInput() const { return input_; }
    /// @}

    /// Performs the actual processing of the input for invalid angle detection.
    /// Returns the number of elements removed.
    std::size_t operator()();

  private:
    void populateMeasurementNeighborhoods();
    void markIndexAsRemoved(std::size_t i);
    std::size_t getWorstOffenderIndex();
    void iterativelyFindAndMarkForRemoval();

    static double computeMinDotProduct(double maxAngleDegrees) {
        // Find the dot product associated with two unit vectors
        // separated by the angle specified.  This is the cosine
        // of the angle.
        return std::cos(maxAngleDegrees / 180.0 * MY_PI);
    }
    InputMeasurements& input_;
    bool verbose_;
    bool done_ = false;
    std::vector<MeasurementNeighborhood> measNeighborhoods_;
    std::vector<std::size_t> removedLineNumbers_;
};

} // namespace detail
#endif // INCLUDED_InvalidPointRemoval_h_GUID_1FF303B7_FC7F_448F_E416_A278BCAD6B0A
