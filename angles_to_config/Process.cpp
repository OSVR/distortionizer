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
#include "Subproblems.h"
#include "helper.h" // for a reflection function

// Library/third-party includes
// - none

// Standard includes
#include <cassert>
#include <iostream>

AnglesToConfigSingleEyeProcess::AnglesToConfigSingleEyeProcess(Config const& c) : config_(c) {}

void AnglesToConfigSingleEyeProcess::setScreenSpaceTrimBounds(XYInclusiveBoundsd const& bounds) {
    assert(status_ == Status::Empty && "Should only call this function on an empty object");
    screenTrim_ = bounds;
}

void AnglesToConfigSingleEyeProcess::setInputAngleBounds(XYInclusiveBoundsd const& bounds) {
    assert(status_ == Status::Empty && "Should only call this function on an empty object");
    angleBounds_ = bounds;
}

template <typename ValueType, typename F>
inline bool trimVectorToBounds(std::vector<ValueType>& vec, InclusiveBoundsd bounds, F&& memberGetter) {
    if (!bounds) {
        return false;
    }
    vec.erase(
        std::remove_if(vec.begin(), vec.end(), [&](ValueType const& m) { return bounds.outside(memberGetter(m)); }));
    return true;
}
template <typename ValueType, typename F1, typename F2>
inline bool trimVectorToBounds(std::vector<ValueType>& vec, XYInclusiveBoundsd bounds, F1&& memberGetterX,
                               F2&& memberGetterY) {

    bool ret = false;
    if (!bounds) {
        return ret;
    }
    if (trimVectorToBounds(vec, bounds.x, std::forward<F1>(memberGetterX))) {
        ret = true;
    }
    if (trimVectorToBounds(vec, bounds.y, std::forward<F2>(memberGetterY))) {
        ret = true;
    }
    return ret;
}

int AnglesToConfigSingleEyeProcess::supplyInputMeasurements(InputMeasurements&& meas) {
    assert((status_ == Status::Empty || status_ == Status::HasSomeMapping) && "Should only call this function on an "
                                                                              "empty object or one with only calls to "
                                                                              "supplyInputMeasurements() completed");
    if (config_.verbose) {
        std::cerr << "supplyInputMeasurements provided with an input measurement collection of size " << meas.size()
                  << std::endl;
    }
    bool trimmed = false;
    if (trimVectorToBounds(meas.measurements, screenTrim_,
                           /// get x coordinate of screen
                           [](InputMeasurement const& m) { return m.screen[0]; },
                           /// get y coordinate of screen
                           [](InputMeasurement const& m) { return m.screen[1]; })) {
        trimmed = true;
    }

    /// At this point, we're still in degrees (conversion happens with convert_to_normalized_and_meters)
    if (trimVectorToBounds(meas.measurements, angleBounds_,
                           [](InputMeasurement const& m) { return m.viewAnglesDegrees.longitude(); },
                           [](InputMeasurement const& m) { return m.viewAnglesDegrees.latitude(); })) {
        trimmed = true;
    }

    if (config_.verbose && trimmed) {
        std::cerr << "Size after input trimming: " << meas.size() << std::endl;
    }
    if (config_.verifyAngles) {

        //====================================================================
        // If we've been asked to verify the angles on the meshes, do so now.
        // This makes sure that the direction between neighbors in angle space
        // is consistent with their direction in screen space, removing points
        // that don't satisfy the criterion.  This removes inconsistent points
        // from the simulation (caused by multiple ray bounces or other
        // singularities in the simulation).
        int ret = remove_invalid_points_based_on_angle(meas, config_.maxAngleDiffDegrees, {config_.xx, config_.xy},
                                                       {config_.yx, config_.yy}, config_.verbose);
        if (ret < 0) {
            std::cerr << "Error verifying angles for mesh " << (inputMeasurementChannels_.size() + 1) << std::endl;
            return 60;
        }
        if (config_.verbose) {
            std::cerr << "Removed " << ret << " points from mesh " << (inputMeasurementChannels_.size() + 1)
                      << std::endl;
        }
    }

    inputMeasurementChannels_.emplace_back(std::move(meas));
    status_ = Status::HasSomeMapping;
    return 0;
}

void AnglesToConfigSingleEyeProcess::supplyAdditionalAngles(std::vector<LongLat> const& additionalAngles) {
    additionalAnglePoints_ = convertAdditionalAngles(additionalAngles, config_.depth, config_.useFieldAngles);
}
// for x positive to the right, y positive up.
inline void extendBounds(RectBoundsd& bounds, double x, double y) {
    if (x > bounds.right) {
        bounds.right = x;
    }
    if (x < bounds.left) {
        bounds.left = x;
    }
    if (y > bounds.top) {
        bounds.top = y;
    }
    if (y < bounds.bottom) {
        bounds.bottom = y;
    }
}
// for x positive to the right, y positive up.
inline void extendBounds(RectBoundsd& bounds, Point2d const& pt) { extendBounds(bounds, pt[0], pt[1]); }
void AnglesToConfigSingleEyeProcess::computeBounds() {
    assert(status_ == Status::HasSomeMapping &&
           "Should only call this function after one or more calls to supplyInputMapping()");
    status_ = Status::HasBoundsComputed;
    if (config_.computeScreenBounds) {
        //====================================================================
        // If we've been asked to auto-range the screen coordinates, compute
        // them here.  Look at all of the points from all of the colors and
        // make a bound on all of them.
        screenBounds_.left = screenBounds_.right = inputMeasurementChannels_.front().measurements.front().screen[0];
        screenBounds_.top = screenBounds_.bottom = inputMeasurementChannels_.front().measurements.front().screen[1];
        for (auto const& chan : inputMeasurementChannels_) {
            for (auto const& meas : chan.measurements) {
                extendBounds(screenBounds_, meas.screen);
            }
        }
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
    for (auto const& chan : inputMeasurementChannels_) {
        //====================================================================
        // Convert the input values into normalized coordinates and into 3D
        // locations.
        normalizedMeasurementChannels_.push_back(convert_to_normalized_and_meters(
            chan, config_.toMeters, config_.depth, screenBounds_, config_.useFieldAngles));
    }
}

int AnglesToConfigSingleEyeProcess::computeScreenAndMeshes(SingleEyeOutput& outResults,
                                                           OutputOptions const& outOpts) const {
    ScreenDetails screenDetails;
    if (!::findScreen(outResults.projection, screenDetails, normalizedMeasurementChannels_, additionalAnglePoints_,
                      config_.verbose)) {
        std::cerr << "Error: Could not find screen" << std::endl;
        return 3;
    }
    for (auto& normalizedMeasChan : normalizedMeasurementChannels_) {
        MeshDescription mesh = findMesh(normalizedMeasChan, screenDetails, config_.verbose);
        if (mesh.empty()) {
            std::cerr << "Error: Could not find mesh" << std::endl;
            return 30;
        }
        if (mesh.size() != normalizedMeasChan.size()) {
            std::cerr << "Error: Mesh size " << mesh.size() << " does not match mapping size"
                      << normalizedMeasChan.size() << std::endl;
            return 4;
        }
        if (outOpts.u1) {
            trimVectorToBounds(mesh, outOpts.u1, [](MeshDescriptionRow const& row) {
                // 0 is for the first uv pair
                // then 0 is for u.
                return row[0][0];
            });
        }
        outResults.meshes.push_back(std::move(mesh));
    }

    return 0;
}
inline static InputMeasurements reflect(InputMeasurements const& in) {
    InputMeasurements ret;
    ret.inputSource = in.inputSource + "[reflected]";
    for (auto& meas : in.measurements) {
        InputMeasurement reflected = meas;

        /// negate x
        /// @todo subtract x from right?
        reflected.screen[0] *= 1;

        /// negate x/longitude
        reflected.viewAnglesDegrees.longitude() *= -1;
        ret.measurements.push_back(reflected);
    }
    return ret;
}

inline static NormalizedMeasurements reflect(NormalizedMeasurements const& in) {
    NormalizedMeasurements ret;
    ret.inputSource = in.inputSource + "[reflected]";
    for (auto& meas : in.measurements) {
        NormalizedMeasurement reflected = meas;

        /// subtract x from 1
        reflected.screen[0] = 1. - reflected.screen[0];
        /// negate x
        reflected.pointFromView[0] *= -1;

        ret.measurements.push_back(reflected);
    }
    return ret;
}
AnglesToConfigSingleEyeProcess AnglesToConfigSingleEyeProcess::reflectedHorizontally() const {
    AnglesToConfigSingleEyeProcess ret(config_);
    ret.status_ = status_;
    for (auto& chan : inputMeasurementChannels_) {
        ret.inputMeasurementChannels_.push_back(reflect(chan));
    }

    for (auto& chan : normalizedMeasurementChannels_) {
        ret.normalizedMeasurementChannels_.push_back(reflect(chan));
    }

    ret.additionalAnglePoints_ = reflectPoints(additionalAnglePoints_);
    ret.screenBounds_ = screenBounds_.reflectedHorizontally();

    return ret;
}
