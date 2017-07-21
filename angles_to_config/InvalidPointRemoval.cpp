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
#include "InvalidPointRemoval.h"

// Library/third-party includes
// - none

// Standard includes
#include <iostream>

static inline double getNeighborError(InputMeasurement const& a, InputMeasurement const& b, Eigen::Vector2d const& xxxy,
                                      Eigen::Vector2d const& yxyy, bool verbose = false) {
    // Map the difference in angle space between the point and
    // its neighbor into screen space.
    Eigen::Vector2d angleVec = ei::map(a.viewAnglesDegrees.longLat) - ei::map(b.viewAnglesDegrees.longLat);
    Eigen::Vector2d screenMappedVec = angleVec[0] * xxxy + angleVec[1] * yxyy;

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

namespace detail {

MeasurementNeighborhood::MeasurementNeighborhood(InvalidPointRemover const& parent, std::size_t i)
    : index(i), parent_(&parent) {
    populateNeighborData(parent_->getInput());
}

void MeasurementNeighborhood::populateNeighborData(InputMeasurements const& input) {
    auto refMeas = input.measurements[index];
    const auto n = input.size();
    neighborData_.reserve(n - 1);
    for (std::size_t i = 0; i < n; ++i) {
        if (i == index) {
            continue;
        }
        auto dist = StatusDistanceIndex{
            angleSquaredDistance(input.measurements[i].viewAnglesDegrees, refMeas.viewAnglesDegrees), i};
        neighborData_.push_back(dist);
    }
}
void MeasurementNeighborhood::markIndexRemoved(const std::size_t removedIndex) {
    if (MeasurementStatus::Removed == status) {
        // we can ignore this, we're already a non-entity.
        return;
    }
    if (index == removedIndex) {
        status = MeasurementStatus::Removed;
        /// won't be doing anything else with this object, so we can exit now.
        return;
    }
    for (auto& neighbor : neighborData_) {
        if (removedIndex == neighbor.index) {
            neighbor.status = MeasurementStatus::Removed;
            /// Only one entry for each neighbor.
            return;
        }
    }
    assert(0 && "Index not found in list or self: self and list should contain all indices!");
}
std::pair<std::size_t, double> MeasurementNeighborhood::getNeighborErrorCount(bool verbose) {
    if (MeasurementStatus::Removed == status) {
        /// we've already been removed, can't be better than that.
        return std::make_pair(0, 100.);
    }
#ifdef SENSICS_JUST_TESTING_TODO
    dirty_ = true;
#endif

    if (dirty_) {
        neighborErrors_ = countNeighborErrors(verbose);
        dirty_ = false;
    }
    return neighborErrors_;
}
std::pair<std::size_t, double> MeasurementNeighborhood::countNeighborErrors(bool verbose) {
    if (MeasurementStatus::Removed == status) {
        /// we've already been removed, can't be better than that.
        return std::make_pair(0, 100.);
    }
    auto& input = parent_->getInput();
    auto& refMeas = input.measurements[index];
    auto endIter = findEndOfUpToNClosest();
    auto minDotProduct = parent_->minDotProduct;
    auto xxxy = parent_->xxxy;
    auto yxyy = parent_->yxyy;
#if 0
    return std::count_if(neighborData_.begin(), endIter, [&](StatusDistanceIndex& elt) {
        /// count if element is active
        return (MeasurementStatus::Active == elt.status) &&
               /// and if error is less than minDotProduct.
               (elt.getError(input, refMeas, xxxy, yxyy, verbose) < minDotProduct);
    });
#else

    double product = 1;
    std::size_t count = 0;
    std::for_each(neighborData_.begin(), endIter, [&](StatusDistanceIndex& elt) {
        if (MeasurementStatus::Active == elt.status) {
            auto err = elt.getError(input, refMeas, xxxy, yxyy, verbose);
            if (err < minDotProduct) {
                count++;
                product *= err;
            }
        }
    });

    return std::make_pair(count, product);
#endif
}

MeasurementNeighborhood::iterator MeasurementNeighborhood::findEndOfUpToNClosest(const std::size_t n) {
    const auto beginIter = neighborData_.cbegin();
    auto endIter = neighborData_.end();
    if (neighborData_.size() > n) {
        auto midElt = neighborData_.begin();
        std::advance(midElt, n - 1);
        std::nth_element(neighborData_.begin(), midElt, neighborData_.end());

        /// Adjust the range we'll iterate through below so we only iterate thru the ones we've split.
        if (midElt == endIter) {
            // not sure how this happened, but makes it unsafe to dereference midElt.
            return endIter;
        }
#if 0
            // Back up to exclude any removed elements
            while (MeasurementStatus::Removed == midElt->status && midElt != beginIter) {
                midElt--;
            }
#endif
        // OK, so then we want the "one past the end"
        endIter = midElt;
        endIter++;
    }
    return endIter;
}

StatusDistanceIndex::StatusDistanceIndex(double dist, std::size_t i) : distance(dist), index(i) {}

double StatusDistanceIndex::getError(InputMeasurements const& input, InputMeasurement const& refMeas,
                                     Eigen::Vector2d const& xxxy, Eigen::Vector2d const& yxyy, bool verbose) {
    if (!populatedError) {
        populatedError = true;
        error = getNeighborError(input.measurements[index], refMeas, xxxy, yxyy, verbose);
    }
    return error;
}

InvalidPointRemover::InvalidPointRemover(InputMeasurements& input, double maxAngleDegrees, Point2d const& xxxyVec,
                                         Point2d const& yxyyVec, bool verbose)
    : xxxy(ei::map(xxxyVec)),
      yxyy(ei::map(yxyyVec)),
      minDotProduct(computeMinDotProduct(maxAngleDegrees)),
      input_(input),
      verbose_(verbose) {}

std::size_t InvalidPointRemover::operator()() {
    assert(!done_ && "Can only process once!");
    done_ = true;

    populateMeasurementNeighborhoods();

    iterativelyFindAndMarkForRemoval();

    // Now that we've found and recorded which ones we want to remove, actually remove them before returning.
    /// Easiest way to do this is with an .erase(remove_if( idiom, but we don't get indices that way, so this object has
    /// actually been storing source line numbers.

    /// Sort removed line numbers so we can do a quick check for membership.
    std::sort(removedLineNumbers_.begin(), removedLineNumbers_.end());

    input_.measurements.erase(
        std::remove_if(input_.measurements.begin(), input_.measurements.end(), [&](InputMeasurement const& meas) {
            /// if we find our line number in the list, we're a goner.
            return std::binary_search(removedLineNumbers_.begin(), removedLineNumbers_.end(), meas.lineNumber);
        }));
    return removedLineNumbers_.size();
}

void InvalidPointRemover::populateMeasurementNeighborhoods() {
    const auto n = input_.size();
    for (std::size_t i = 0; i < n; ++i) {
        measNeighborhoods_.emplace_back(*this, i);
    }
}

void InvalidPointRemover::markIndexAsRemoved(std::size_t i) {
    removedLineNumbers_.push_back(input_.measurements[i].lineNumber);
    for (auto& measNeigh : measNeighborhoods_) {
        measNeigh.markIndexRemoved(i);
    }
}

std::size_t InvalidPointRemover::getWorstOffenderIndex() {
    size_t worstIndex = 0;
    size_t worstCount = 0;
    double worstProduct = minDotProduct * 2; // arbitrarily too large to matter.

    for (auto& neighborhood : measNeighborhoods_) {
        if (!neighborhood.isValid()) {
            continue;
        }
#if 1
        const bool childVerbose = false;
#else
        const bool childVerbose = verbose_ && (neighborhood.index == 0);
#endif
        std::size_t count;
        double product;
        std::tie(count, product) = neighborhood.getNeighborErrorCount(childVerbose);
        if (count > worstCount || (count != 0 && count == worstCount && product < worstProduct)) {
            worstCount = count;
            worstProduct = product;
            worstIndex = neighborhood.index;
        }
    }
    if (worstCount == 0) {
        return input_.size();
    }
    if (verbose_) {
        std::cerr << "Worst remaining: " << input_.measurements[worstIndex].getOrigin(input_) << " with " << worstCount
                  << " invalid neighbor angles." << std::endl;
        std::cerr << "\tproduct: " << worstProduct << std::endl;
    }
    return worstIndex;
}
void InvalidPointRemover::iterativelyFindAndMarkForRemoval() {
    // We remove the worst offender from the list each time,
    // then re-start. (Well, actually, we just tag it as removed, but...)
    // Assuming that we get the actual outlier,
    // as opposed to one of its neighbors, this avoids trimming
    // too many points from the vector.
    bool foundOutlier;
    do {
        foundOutlier = false;
        std::size_t worstIdx = getWorstOffenderIndex();
        if (worstIdx < input_.size()) {
            // OK, we found one. "Remove" it.
            markIndexAsRemoved(worstIdx);
            foundOutlier = true;
        }
    } while (foundOutlier);
}
} // namespace detail
