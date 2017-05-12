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

bool convert_to_normalized_and_meters(std::vector<Mapping>& mapping, double toMeters, double depth, double left,
                                      double bottom, double right, double top, bool useFieldAngles = false);

bool findScreen(const std::vector<Mapping>& mapping, double left, double bottom, double right, double top,
                ScreenDescription& screen, bool verbose = false);

bool findMesh(const std::vector<Mapping>& mapping, double left, double bottom, double right, double top,
              ScreenDescription const& screen, MeshDescription& mesh, bool verbose = false);

// Produce a mapping that is reflected around X=0 in both angles and
// screen coordinates.
std::vector<Mapping> reflect_mapping(std::vector<Mapping> const& mapping);
