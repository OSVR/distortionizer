cmake_minimum_required(VERSION 3.1.0)
project(ImageBasedDistortion)

#-----------------------------------------------------------------------------
# Local CMake Modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set (CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH};F:/Packages/nVidia/VR_SDK_under_NDA/2015-11-19-gameworks-vr-for-headset-developers-1.0/direct_mode/DirectModeSample;F:/Packages/SDL/SDL2-2.0.3;F:/Packages/OpenCV/opencv-2.4.10/build;"C:/Program Files/jsoncpp";"C:/Program Files/jsoncpp/lib/cmake/jsoncpp";C:/Program Files/Microsoft SDKs/Windows/v6.0;F:/taylorr/Personal/Work/consulting/Sensics/OSVR/src/OSVR-Core/vendor/eigen-3.2.7;C:/usr/local CACHE PATH "Path that CMake uses to find packages")

#-----------------------------------------------------------------------------
find_package(jsoncpp REQUIRED)
if(TARGET jsoncpp_lib_static AND NOT TARGET jsoncpp_lib)
	add_library(jsoncpp_lib ALIAS jsoncpp_lib_static)
endif()

#-----------------------------------------------------------------------------
add_executable(AnglesToConfig AnglesToConfig.cpp)
target_link_libraries(AnglesToConfig PRIVATE jsoncpp_lib)

