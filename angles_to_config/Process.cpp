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
#include "Process.h"
#include "helper.h"

// Library/third-party includes
// - none

// Standard includes
#include <cassert>
#include <iostream>

AnglesToConfigSingleEyeProcess::AnglesToConfigSingleEyeProcess(Config const& c) : config_(c) {}

int AnglesToConfigSingleEyeProcess::supplyInputMapping(std::vector<Mapping>&& mapping) {

    assert(
        (status_ == Status::Empty || status_ == Status::HasSomeMapping) &&
        "Should only call this function on an empty object or one with only calls to supplyInputMapping() completed");
    mappings_.emplace_back(std::move(mapping));
    status_ = Status::HasSomeMapping;
    if (!config_.verifyAngles) {
        return 0;
    }
    auto& latestMapping = mappings_.back();
    auto m = mappings_.size();

    //====================================================================
    // If we've been asked to verify the angles on the meshes, do so now.
    // This makes sure that the direction between neighbors in angle space
    // is consistent with their direction in screen space, removing points
    // that don't satisfy the criterion.  This removes inconsistent points
    // from the simulation (caused by multiple ray bounces or other
    // singularities in the simulation).
    int ret = remove_invalid_points_based_on_angle(latestMapping, config_.xx, config_.xy, config_.yx, config_.yy,
                                                   config_.maxAngleDiffDegrees);
    if (ret < 0) {
        std::cerr << "Error verifying angles for mesh " << m << std::endl;
        return 60;
    }
    if (config_.verbose) {
        std::cerr << "Removed " << ret << " points from mesh " << m << std::endl;
    }
    return 0;
}
void AnglesToConfigSingleEyeProcess::computeBounds() {
    assert(status_ == Status::HasSomeMapping &&
           "Should only call this function after one or more calls to supplyInputMapping()");
    status_ = Status::HasBoundsComputed;
    if (config_.computeScreenBounds) {
        //====================================================================
        // If we've been asked to auto-range the screen coordinates, compute
        // them here.  Look at all of the points from all of the colors and
        // make a bound on all of them.
        screenBounds_.left = screenBounds_.right = mappings_[0][0].xyLatLong.x;
        screenBounds_.top = screenBounds_.bottom = mappings_[0][0].xyLatLong.y;
    } else {
        screenBounds_ = config_.suppliedScreenBounds;
    }

    if (config_.verbose) {
        std::cerr << "Left, bottom, right, top = " << screenBounds_.left << ", " << screenBounds_.bottom << ", "
                  << screenBounds_.right << ", " << screenBounds_.top << std::endl;
    }
}

void AnglesToConfigSingleEyeProcess::normalizeMappings() {
    assert(status_ == Status::HasBoundsComputed && "Should only call this function after calling computeBounds()");
    status_ = Status::HasMappingsNormalized;
    for (auto& mapping : mappings_) {
        //====================================================================
        // Convert the input values into normalized coordinates and into 3D
        // locations.
        convert_to_normalized_and_meters(mapping, config_.toMeters, config_.depth, screenBounds_.left,
                                         screenBounds_.bottom, screenBounds_.right, screenBounds_.top,
                                         config_.useFieldAngles);
    }

    /// Prepare additional data structure for "full mapping"
    prepareFullMapping();
}

int AnglesToConfigSingleEyeProcess::computeScreenAndMeshes(SingleEyeOutput& outResults) const {
    if (!::findScreen(getFullMapping(), screenBounds_.left, screenBounds_.bottom, screenBounds_.right,
                      screenBounds_.top, outResults.screen, config_.verbose)) {
        std::cerr << "Error: Could not find screen" << std::endl;
        return 3;
    }
    for (auto& mapping : mappings_) {
        MeshDescription mesh;
        if (!findMesh(mapping, screenBounds_.left, screenBounds_.bottom, screenBounds_.right, screenBounds_.top,
                      outResults.screen, mesh, config_.verbose)) {

            std::cerr << "Error: Could not find mesh" << std::endl;
            return 30;
        }
        if (mesh.size() != mapping.size()) {
            std::cerr << "Error: Mesh size " << mesh.size() << " does not match mapping size" << mapping.size()
                      << std::endl;
            return 4;
        }
        outResults.meshes.push_back(std::move(mesh));
    }

    return 0;
}

std::vector<Mapping> reflect_normalized_mapping(std::vector<Mapping> const& mapping) {
    std::vector<Mapping> ret;
    for (size_t i = 0; i < mapping.size(); i++) {
        ret.push_back(mapping[i]);
        ret[i].xyLatLong.longitude = 1 - ret[i].xyLatLong.longitude;
        ret[i].xyLatLong.x = 1 - ret[i].xyLatLong.x;
    }

    return ret;
}

AnglesToConfigSingleEyeProcess AnglesToConfigSingleEyeProcess::reflectedHorizontally() const {
    AnglesToConfigSingleEyeProcess ret(config_);
    ret.status_ = status_;
    auto reflect = &reflect_mapping;
    if (status_ == Status::HasMappingsNormalized) {
        reflect = &reflect_normalized_mapping;
    }
    for (auto& mapping : mappings_) {
        ret.mappings_.push_back(reflect(mapping));
    }
    ret.screenBounds_ = screenBounds_.reflectedHorizontally();
    ret.fullMapping_ = reflect(fullMapping_);

    return ret;
}
void AnglesToConfigSingleEyeProcess::prepareFullMapping() {
    /// if we only have one mapping, then it's the full mapping
    if (mappings_.size() < 2) {
        return;
    }

    /// Otherwise, we need to combine the mappings.
    if (fullMapping_.empty()) {
        for (auto const& mapping : mappings_) {
            fullMapping_.insert(fullMapping_.end(), mapping.begin(), mapping.end());
        }
    }
}

std::vector<Mapping> const& AnglesToConfigSingleEyeProcess::getFullMapping() const {
    /// if we only have one mapping, then it's the full mapping
    if (mappings_.size() == 1) {
        return mappings_[0];
    }

    /// Otherwise, we need to return the combined mappings.
    return fullMapping_;
}
