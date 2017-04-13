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

#ifndef INCLUDED_GenerateOutput_h_GUID_463EC8A0_A36C_4598_123C_7239DE78853E
#define INCLUDED_GenerateOutput_h_GUID_463EC8A0_A36C_4598_123C_7239DE78853E

// Internal Includes
#include "types.h"

// Library/third-party includes
// - none

// Standard includes
#include <iosfwd>

void writeMesh(std::ostream& s, MeshDescription const& mesh);

#endif // INCLUDED_GenerateOutput_h_GUID_463EC8A0_A36C_4598_123C_7239DE78853E
