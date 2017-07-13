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
#include "GenerateOutput.h"
#include "JsonUtils.h"
#include "Process.h"
#include "helper.h"

// Library/third-party includes
#include <json/reader.h>
#include <json/value.h>

// Standard includes
#include <cassert>
#include <fstream>
#include <iostream>

using namespace osvr::util;
/// Returns true and changes outval only if a full set of bounds is found.
bool getBounds(Json::Value const& obj, RectBoundsd& outVal) {
    if (!obj.isObject()) {
        return false;
    }
    RectBoundsd ret;
    {
        auto& left = obj["left"];
        if (!json_is<double>(left)) {
            return false;
        }
        ret.left = json_cast<double>(left);
    }

    {
        auto& right = obj["right"];
        if (!json_is<double>(right)) {
            return false;
        }
        ret.right = json_cast<double>(right);
    }
    {
        auto& top = obj["top"];
        if (!json_is<double>(top)) {
            return false;
        }
        ret.top = json_cast<double>(top);
    }
    {
        auto& bottom = obj["bottom"];
        if (!json_is<double>(bottom)) {
            return false;
        }
        ret.bottom = json_cast<double>(bottom);
    }
    outVal = ret;
    return true;
}

void basicConfigParsing(Json::Value const& root, Config& conf) {
    conf.depth = getWithDefault(root, "depth", conf.depth);
    {
        auto& screen = root["screen"];

        auto gotBounds = getBounds(screen, conf.suppliedScreenBounds);
        conf.computeScreenBounds = !gotBounds;
    }
    {
        auto& inputType = root["inputType"];
        if (json_is<std::string>(inputType)) {
            auto inputTypeStr = json_cast<std::string>(inputType);
            if (inputTypeStr == "fieldAngles") {
                conf.useFieldAngles = true;
            } else if (inputTypeStr == "latitudeLongitude") {
                conf.useFieldAngles = false;
            } else {
                std::cerr << "Warning: Unrecognized value for 'inputType'" << std::endl;
            }
        }
    }
    /// @todo verify angles
    conf.verbose = getWithDefault(root, "verbose", conf.verbose);
}

bool supplySingleFileData(std::string const& fn, AnglesToConfigSingleEyeProcess& process) {
    std::ifstream inFile(fn);
    if (!inFile.good()) {
        std::cerr << "Error: Could not open file " << fn << std::endl;
        return false;
    }
    auto measurements = readInputMeasurements(fn, inFile);
    if (measurements.empty()) {
        std::cerr << "Error: Failure to read records from " << fn << std::endl;
        return false;
    }
    auto success = process.supplyInputMeasurements(std::move(measurements));
    if (0 != success) {
        std::cerr << "Error: supplyInputMeasurements failed, return value: " << success << std::endl;
        return false;
    }
    return true;
}
bool withDataJsonSchemaError() {
    std::cerr << "Warning: Elements of 'input' must include a 'mapping' property containing either be a string or an "
                 "object with red, green, blue as keys, and may include an 'additionalVisibleAngles' property"
              << std::endl;
    return false;
}

static const auto RED_STR = "red";
static const auto GREEN_STR = "green";
static const auto BLUE_STR = "blue";
static const auto CHANNELS_STRS = {RED_STR, GREEN_STR, BLUE_STR};

bool processMappingElement(Json::Value const& mapping, AnglesToConfigSingleEyeProcess& process) {
    if (json_is<std::string>(mapping)) {
        // OK this is a mono channel
        auto fn = json_cast<std::string>(mapping);
        if (!supplySingleFileData(fn, process)) {
            return false;
        }
        return true;
    }
    if (mapping.isObject()) {
        for (auto& channel : CHANNELS_STRS) {
            auto chanFnVal = mapping[channel];
            if (!json_is<std::string>(chanFnVal)) {
                return withDataJsonSchemaError();
            }
            if (!supplySingleFileData(json_cast<std::string>(chanFnVal), process)) {
                return false;
            }
        }
        return true;
    }
    return withDataJsonSchemaError();
}

bool attemptSingleEyeProcessing(Json::Value const& inputData, AnglesToConfigSingleEyeProcess& process) {
    if (inputData.isNull()) {
        return false;
    }
    if (!inputData.isObject() || inputData["mapping"].isNull()) {
        /// OK, so this is the old schema, where the "mapping" level didn't exist.
        /// Deal with it, as a fall back.
        if (!processMappingElement(inputData["mapping"], process)) {
            std::cerr << "Failed loading mapping data in fallback-schema path." << std::endl;
            return false;
        }
    } else {
        if (!processMappingElement(inputData["mapping"], process)) {
            std::cerr << "Failed loading mapping data." << std::endl;
            return false;
        }
    }

    auto& additionalFile = inputData["additionalVisibleAngles"];
    if (json_is<std::string>(additionalFile)) {
        auto newFn = json_cast<std::string>(additionalFile);
        auto stream = std::ifstream(newFn);
        auto addAngles = readAdditionalAngles(stream);
        if (addAngles.empty()) {
            std::cerr << "Warning: additionalVisibleAngles filename " << newFn
                      << " supplied but no pairs of angles could be read." << std::endl;
        } else {
            process.supplyAdditionalAngles(addAngles);
        }
    }

    // OK, if we made it here we have at least one file loaded in properly.
    process.computeBounds();
    process.normalizeMappings();

    return true;
}

int outputClientMeshData(std::ostream& os, SingleEyeOutput const& left, SingleEyeOutput const& right) {
    //====================================================================
    // Construct Json screen description.
    // We do this by hand rather than using JsonCPP because we need
    // to control the printed precision of the numbers to avoid making
    // a huge file.
    os << "{" << std::endl;
    os << R"( "display": {)" << std::endl;
    os << R"(  "hmd": {)" << std::endl;

    os << R"(   "field_of_view": {)" << std::endl;
    os << R"(    "monocular_horizontal": )" << right.projection.hFOVDegrees << "," << std::endl;
    os << R"(    "monocular_vertical": )" << right.projection.vFOVDegrees << "," << std::endl;
    os << R"(    "overlap_percent": )" << right.projection.overlapPercent << "," << std::endl;
    os << R"(    "pitch_tilt": 0)" << std::endl;
    os << "   }," << std::endl; // field_of_view

    os << R"(   "eyes": [)" << std::endl;
    os << "    {" << std::endl;
    os << R"(     "center_proj_x": )" << left.projection.cop[0] << "," << std::endl;
    os << R"(     "center_proj_y": )" << left.projection.cop[1] << "," << std::endl;
    os << R"(     "rotate_180": 0)" << std::endl;
    os << "    }," << std::endl;
    os << "    {" << std::endl;
    os << R"(     "center_proj_x": )" << right.projection.cop[0] << "," << std::endl;
    os << R"(     "center_proj_y": )" << right.projection.cop[1] << "," << std::endl;
    os << R"(     "rotate_180": 0)" << std::endl;
    os << "    }" << std::endl;
    os << "   ]," << std::endl; // eyes

    os << R"(   "distortion": {)" << std::endl;
    switch (left.meshes.size()) {
    case 1:
        os << R"(    "type": "mono_point_samples",)" << std::endl;
        os << R"(    "mono_point_samples": [)" << std::endl;
        writeMesh(os, left.meshes[0]);
        os << "," << std::endl;
        writeMesh(os, right.meshes[0]);
        os << "    ]" << std::endl; // mono_point_samples
        os << "   }" << std::endl;  // distortion
        break;
    case 3:
        os << R"(    "type": "rgb_point_samples",)" << std::endl;
        os << R"(    "red_point_samples": [)" << std::endl;
        writeMesh(os, left.meshes[0]);
        os << "," << std::endl;
        writeMesh(os, right.meshes[0]);
        os << "    ]," << std::endl; // red_point_samples
        os << R"(    "green_point_samples": [)" << std::endl;
        writeMesh(os, left.meshes[1]);
        os << "," << std::endl;
        writeMesh(os, right.meshes[1]);
        os << "    ]," << std::endl; // green_point_samples
        os << R"(    "blue_point_samples": [)" << std::endl;
        writeMesh(os, left.meshes[2]);
        os << "," << std::endl;
        writeMesh(os, right.meshes[2]);
        os << "    ]" << std::endl; // blue_point_samples
        os << "   }" << std::endl;  // distortion
        break;
    default:
        std::cerr << "Error: Unexpected number of meshes: " << left.meshes.size() << std::endl;
        return 3;
    }

    os << "  }" << std::endl; // hmd
    os << " }" << std::endl;  // display
    os << "}" << std::endl;   // Closes outer object
    return 0;
}
inline InclusiveBoundsd parseBounds(Json::Value const& trimRange) {
    InclusiveBoundsd ret;
    if (trimRange.isArray() && (trimRange.size() == 2) && trimRange[0].isDouble() && trimRange[1].isDouble()) {
        ret = InclusiveBoundsd{json_cast<double>(trimRange[0]), json_cast<double>(trimRange[1])};
    }
    return ret;
}

inline XYInclusiveBoundsd parseXYBounds(Json::Value const& trimRangeXY) {
    XYInclusiveBoundsd ret;
    ret.x = parseBounds(trimRangeXY["x"]);
    ret.y = parseBounds(trimRangeXY["y"]);
    return ret;
}

inline void parseTrimmingForEye(Json::Value const& trimming, const char eyeName[],
                                AnglesToConfigSingleEyeProcess& process) {
    {
        auto bounds = parseXYBounds(trimming["screenPointRange"][eyeName]);
        if (bounds) {
            std::cerr << "Trimming " << eyeName << " eye screen point components: " << bounds << std::endl;
            process.setScreenSpaceTrimBounds(bounds);
        }
    }

    {
        auto bounds = parseXYBounds(trimming["fieldAngleRange"][eyeName]);
        if (bounds) {
            std::cerr << "Trimming " << eyeName << " eye input degrees X longitude/Y latitude components: " << bounds
                      << std::endl;
            process.setInputAngleBounds(bounds);
        }
    }
}
int main(int argc, char* argv[]) {
    std::cerr << "Using config file " << argv[1] << std::endl;
    Json::Value root;
    {
        Json::Reader reader;
        std::ifstream infile(argv[1]);
        if (!reader.parse(infile, root)) {
            std::cerr << "Error parsing config file " << argv[1] << "!" << std::endl;
            return -1;
        }
    }
    Config conf;
    basicConfigParsing(root, conf);

    auto& input = root["input"];
    AnglesToConfigSingleEyeProcess leftEyeProcess(conf);
    AnglesToConfigSingleEyeProcess rightEyeProcess(conf);

    auto& trimming = root["trimming"];
    parseTrimmingForEye(trimming, "left", leftEyeProcess);
    parseTrimmingForEye(trimming, "right", rightEyeProcess);

    auto haveLeft = attemptSingleEyeProcessing(input["left"], leftEyeProcess);

    auto haveRight = attemptSingleEyeProcessing(input["right"], rightEyeProcess);
    if (!haveLeft && !haveRight) {
        std::cerr << "Error: must have valid input data for at least one eye, and none provided." << std::endl;
        return -1;
    } else if (haveLeft && !haveRight) {
        std::cerr << "Note: Provided only left eye, assuming right is symmetrical." << std::endl;
        rightEyeProcess = leftEyeProcess.reflectedHorizontally();
    } else if (haveRight && !haveLeft) {
        std::cerr << "Note: Provided only right eye, assuming left is symmetrical." << std::endl;
        leftEyeProcess = rightEyeProcess.reflectedHorizontally();
    } else {
        assert(haveLeft && haveRight && "only remaining option");
        std::cerr << "Note: Data provided for both eyes." << std::endl;
    }

    OutputOptions leftOutOpts{};
    OutputOptions rightOutOpts{};
    {
        auto& leftEyeU1 = trimming["outputRange"]["u1"]["left"];
        auto bounds = parseBounds(trimming["outputRange"]["u1"]["left"]);
        if (bounds) {
            leftOutOpts.u1 = bounds;
            std::cerr << "Trimming left eye u1 to range " << leftOutOpts.u1 << std::endl;
            /// mirror to right eye
            rightOutOpts.u1 = InclusiveBoundsd(1. - bounds.getMax(), 1. - bounds.getMin());
        }
    }

    SingleEyeOutput leftOutput;
    {
        auto success = leftEyeProcess.computeScreenAndMeshes(leftOutput, leftOutOpts);
        if (0 != success) {
            std::cerr << "Error: Failure in computeScreenAndMeshes for left eye, code " << success << std::endl;
            return success;
        }
    }

    SingleEyeOutput rightOutput;
    {
        auto success = rightEyeProcess.computeScreenAndMeshes(rightOutput, rightOutOpts);
        if (0 != success) {
            std::cerr << "Error: Failure in computeScreenAndMeshes for right eye, code " << success << std::endl;
            return success;
        }
    }

    if (rightOutput.meshes.size() != leftOutput.meshes.size()) {
        std::cerr << "Error: Mismatch in number of meshes between eyes - either both are monochromatic or both are RGB"
                  << std::endl;
        return -1;
    }
    {
        auto success = outputClientMeshData(std::cout, leftOutput, rightOutput);
        if (0 != success) {
            std::cerr << "Error: failure in outputClientMeshData: " << success << std::endl;
            return success;
        }
    }
    std::cerr << "Processing complete!" << std::endl;
    return 0;
}
