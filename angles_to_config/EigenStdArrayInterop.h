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

#ifndef INCLUDED_EigenStdArrayInterop_h_GUID_86D6CCCF_AE53_4DBA_1E69_400EF6110109
#define INCLUDED_EigenStdArrayInterop_h_GUID_86D6CCCF_AE53_4DBA_1E69_400EF6110109


// Internal Includes
// - none

// Library/third-party includes
#include <Eigen/Core>

// Standard includes
// - none

namespace ei {
/// Function to map std::array as an Eigen vector type
template <typename Scalar, size_t Dim>
inline Eigen::Map<Eigen::Matrix<Scalar, Dim, 1> > map(std::array<Scalar, Dim>& v) {
    return Eigen::Matrix<Scalar, Dim, 1>::Map(v.data());
}
/// Function to map constant std::array as an Eigen vector type
template <typename Scalar, size_t Dim>
inline Eigen::Map<Eigen::Matrix<Scalar, Dim, 1> const> map(std::array<Scalar, Dim> const& v) {
    return Eigen::Matrix<Scalar, Dim, 1>::Map(v.data());
}

/// Function to map std::array as an Eigen array-vector type
template <typename Scalar, size_t Dim>
inline Eigen::Map<Eigen::Array<Scalar, Dim, 1> > mapArray(std::array<Scalar, Dim>& v) {
    return Eigen::Array<Scalar, Dim, 1>::Map(v.data());
}
/// Function to map constant std::array as an Eigen array-vector type
template <typename Scalar, size_t Dim>
inline Eigen::Map<Eigen::Array<Scalar, Dim, 1> const> mapArray(std::array<Scalar, Dim> const& v) {
    return Eigen::Array<Scalar, Dim, 1>::Map(v.data());
}
} // namespace ei

#endif // INCLUDED_EigenStdArrayInterop_h_GUID_86D6CCCF_AE53_4DBA_1E69_400EF6110109

