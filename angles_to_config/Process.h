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

#ifndef INCLUDED_Process_h_GUID_B979F1F0_7103_458A_77CB_BF5B10B5C3A0
#define INCLUDED_Process_h_GUID_B979F1F0_7103_458A_77CB_BF5B10B5C3A0

// Internal Includes
#include "types.h"

// Library/third-party includes
// - none

// Standard includes
#include <vector>

struct OutputOptions {
    InclusiveBoundsd u1;
};

struct SingleEyeOutput {
    ProjectionDescription projection;
    std::vector<MeshDescription> meshes;
};

class AnglesToConfigSingleEyeProcess {
  public:
    AnglesToConfigSingleEyeProcess(Config const& c);

    /// Set an offset/translation in angle space (degrees) to apply to input data.
    void setAngleOffset(Point2d const& v);

    /// Set a twist/rotation (in degrees) about the normal to 0, 0 in angle space to apply to input data.
    void setAngleTwist(double v);

    /// @name Set input trimming bounds
    /// @{
    void setScreenSpaceTrimBounds(XYInclusiveBoundsd const& bounds);
    void setInputAngleBounds(XYInclusiveBoundsd const& bounds);
    /// @}

    /// Configure the screen that is used when projecting input angles.
    void setScreenYRotation(double degrees);
    /// Set configuration option that only has an effect if screen y rotation is nonzero.
    void setAnglesAreScreenRelative(bool val);

    /// Call this as many times as you have input mappings (1 for mono, 3 for RGB)
    int supplyInputMeasurements(InputMeasurements&& meas);

    /// Optional: call at most once with additional angles (in degrees, and same type [field angles/lat-long] as the
    /// input mapping angles) that should be considered visible even though there isn't a mapping to screen space known
    /// for them.
    void supplyAdditionalAngles(std::vector<LongLat> const& additionalAngles);
    /// Call this after finishing all calls to supplyInputMapping() - even if the config contains a supplied bounds.
    void computeBounds();
    /// Call after computeBounds()
    void normalizeMappings();
    /// Call after normalizeMappings() to compute and retrieve final output.
    int computeScreenAndMeshes(SingleEyeOutput& outResults, OutputOptions const& outOpts = OutputOptions{}) const;

    /// Creates a copy of this object but with the mappings and bounds reflected horizontally.
    AnglesToConfigSingleEyeProcess reflectedHorizontally() const;

  private:
    Config config_;
    enum class Status { Empty, HasSomeMapping, HasBoundsComputed, HasMappingsNormalized };

    Status status_ = Status::Empty;
    Point2d angleOffset_;
    double angleTwist_ = 0.;
    double screenRotateYRadians_ = 0.;
    bool anglesScreenRelative_ = true;
    XYInclusiveBoundsd screenTrim_;
    XYInclusiveBoundsd angleBounds_;
    std::vector<InputMeasurements> inputMeasurementChannels_;
    std::vector<NormalizedMeasurements> normalizedMeasurementChannels_;
    XYZList additionalAnglePoints_;
    RectBoundsd screenBounds_;
};

#endif // INCLUDED_Process_h_GUID_B979F1F0_7103_458A_77CB_BF5B10B5C3A0
