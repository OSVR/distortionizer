/** @file
    @brief Header

    Based on code from GetOptionalParameter in OSVR.

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

#ifndef INCLUDED_JsonUtils_h_GUID_F3BC6DA9_BF2B_49A5_A324_8E3A999758DC
#define INCLUDED_JsonUtils_h_GUID_F3BC6DA9_BF2B_49A5_A324_8E3A999758DC

// Internal Includes
// - none

// Library/third-party includes
#include <json/value.h>

// Standard includes
#include <array>
#include <string>

namespace osvr {
namespace util {
    namespace detail {
        template <typename T> struct JsonTypeGetter;
        template <typename T> struct JsonTypeContainment;
#define OSVR_DECLARE_JSON_TYPE_TRAITS(TYPENAME, GET_METHOD, IS_METHOD)                                                 \
    template <> struct JsonTypeGetter<TYPENAME> {                                                                      \
        static TYPENAME apply(Json::Value const& val) { return val.GET_METHOD(); }                                     \
    };                                                                                                                 \
    template <> struct JsonTypeContainment<TYPENAME> {                                                                 \
        static bool apply(Json::Value const& val) { return val.IS_METHOD(); }                                          \
    };
        OSVR_DECLARE_JSON_TYPE_TRAITS(bool, asBool, isBool)
        OSVR_DECLARE_JSON_TYPE_TRAITS(float, asFloat, isDouble)
        OSVR_DECLARE_JSON_TYPE_TRAITS(double, asDouble, isDouble)
        OSVR_DECLARE_JSON_TYPE_TRAITS(int, asInt, isInt)
        OSVR_DECLARE_JSON_TYPE_TRAITS(unsigned int, asUInt, isUInt)
        OSVR_DECLARE_JSON_TYPE_TRAITS(std::string, asString, isString)

#undef OSVR_DECLARE_JSON_TYPE_TRAITS
    } // namespace detail
    template <typename T> inline T json_cast(Json::Value const& val) { return detail::JsonTypeGetter<T>::apply(val); }
    /* End of code from GetOptionalParameter */
    template <typename T> inline bool json_is(Json::Value const& val) {
        return detail::JsonTypeContainment<T>::apply(val);
    }

    template <typename T> inline T getWithDefault(Json::Value const& root, const char* memName, T def) {
        const bool wasMember = root.isMember(memName);
        auto& elt = root[memName];
        if (wasMember) {
            if (json_is<T>(elt)) {
                return json_cast<T>(elt);
            }
        }
        return def;
    }

    template <typename T> inline T getIntrusiveDefault(Json::Value& root, const char* memName, T def) {
        const bool wasMember = root.isMember(memName);
        auto& elt = root[memName];
        if (wasMember) {
            if (json_is<T>(elt)) {
                return json_cast<T>(elt);
            }
        }
        elt = def;
        return def;
    }

    template <typename T> inline T getIntrusiveNonzeroDefault(Json::Value& root, const char* memName, T def) {
        const bool wasMember = root.isMember(memName);
        auto& elt = root[memName];
        if (wasMember) {
            if (json_is<T>(elt)) {
                auto ret = json_cast<T>(elt);
                if (ret == 0) {
                    elt = def;
                    ret = def;
                }
                return ret;
            }
        }
        elt = def;
        return def;
    }
    template <typename T, std::size_t N>
    inline std::array<T, N> getIntrusiveDefaultArray(Json::Value& root, const char* memName,
                                                     std::array<T, N> const& def) {
        std::array<T, N> ret;
        auto& elt = root[memName];
        bool success = false;
        if (elt.isArray() && elt.size() == N) {
            success = true;
            for (Json::Value::ArrayIndex i = 0; i < N; ++i) {
                auto& member = elt[i];
                if (!json_is<T>(member)) {
                    success = false;
                    break;
                }
                ret[i] = json_cast<T>(elt[i]);
            }
        }
        if (success) {
            return ret;
        }
        /// Wasn't successful in loading, will set instead.
        elt.resize(N);
        for (Json::Value::ArrayIndex i = 0; i < N; ++i) {
            elt[i] = def[i];
        }
        return def;
    }
    template <typename T, std::size_t N>
    inline std::array<T, N> getDefaultArray(Json::Value const& root, const char* memName, std::array<T, N> const& def) {
        std::array<T, N> ret;
        auto& elt = root[memName];
        bool success = false;
        if (elt.isArray() && elt.size() == N) {
            success = true;
            for (Json::Value::ArrayIndex i = 0; i < N; ++i) {
                auto& member = elt[i];
                if (!json_is<T>(member)) {
                    success = false;
                    break;
                }
                ret[i] = json_cast<T>(elt[i]);
            }
        }
        if (success) {
            return ret;
        }
        /// Wasn't successful in loading
        return def;
    }

    template <typename T, std::size_t N, typename F>
    inline std::array<T, N> getIntrusiveDefaultArrayWithValidity(Json::Value& root, const char* elemName,
                                                                 std::array<T, N> const& defaultValue,
                                                                 F elementValidityPredicate) {
        std::array<T, N> ret;
        if (root.isMember(elemName)) {
            auto& elt = root[elemName];
            if (elt.isArray() && elt.size() == N) {
                bool didFail = false;
                for (Json::ArrayIndex i = 0; i < N; ++i) {
                    auto& entry = elt[i];
                    if (!json_is<T>(entry)) {
                        // error in parsing, break before completion
                        didFail = true;
                        break;
                    } else if (!elementValidityPredicate(entry)) {
                        // user supplied predicate said no
                        didFail = true;
                        break;
                    } else {
                        // no failure so far.
                        ret[i] = json_cast<T>(entry);
                    }
                }
                // at end of for loop: did we make it through without failure?
                if (!didFail) {
                    return ret;
                }
            }
        }
        // Somehow not valid - we'll put in the default value instead.
        ret = defaultValue;
        auto& elt = root[elemName];
        elt = Json::Value(Json::arrayValue);
        for (Json::ArrayIndex i = 0; i < N; ++i) {
            elt[i] = ret[i];
        }
        return ret;
    }
} // namespace util
} // namespace osvr
#endif // INCLUDED_JsonUtils_h_GUID_F3BC6DA9_BF2B_49A5_A324_8E3A999758DC
