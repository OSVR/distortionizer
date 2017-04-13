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

struct SingleEyeOutput {

    ScreenDescription screen;
    std::vector<MeshDescription> meshes;
};
class AnglesToConfigSingleEyeProcess {
  public:
    AnglesToConfigSingleEyeProcess(Config const& c);
    /// Call this as many times as you have input mappings (1 for mono, 3 for RGB)
    int supplyInputMapping(std::vector<Mapping>&& mapping);
    /// Call this after finishing all calls to supplyInputMapping() - even if the config contains a supplied bounds.
    void computeBounds();
    /// Call after computeBounds()
    void normalizeMappings();
    /// Call after normalizeMappings() to compute and retrieve final output.
    int computeScreenAndMeshes(SingleEyeOutput& outResults) const;

    /// Creates a copy of this object but with the mappings and bounds reflected horizontally.
    AnglesToConfigSingleEyeProcess reflectedHorizontally() const;

  private:
    std::vector<Mapping> const& getFullMapping() const;
    void prepareFullMapping();
    Config config_;
    enum class Status { Empty, HasSomeMapping, HasBoundsComputed, HasMappingsNormalized };

    Status status_ = Status::Empty;
    std::vector<std::vector<Mapping>> mappings_;
    RectBoundsd screenBounds_;
    /// only populated when more than one element in mappings_
    std::vector<Mapping> fullMapping_;
};

#endif // INCLUDED_Process_h_GUID_B979F1F0_7103_458A_77CB_BF5B10B5C3A0
