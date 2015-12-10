cmake_minimum_required(VERSION 3.1.0)
project(ImageBasedDistortion)

#-----------------------------------------------------------------------------
# Local CMake Modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set (CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH};F:/Packages/nVidia/VR_SDK_under_NDA/2015-11-19-gameworks-vr-for-headset-developers-1.0/direct_mode/DirectModeSample;F:/Packages/SDL/SDL2-2.0.3;F:/Packages/OpenCV/opencv-2.4.10/build;C:/Program Files/jsoncpp;C:/Program Files/Microsoft SDKs/Windows/v6.0;F:/taylorr/Personal/Work/consulting/Sensics/OSVR/src/OSVR-Core/vendor/eigen-3.2.7;C:/usr/local CACHE PATH "Path that CMake uses to find packages")

#-----------------------------------------------------------------------------
# This looks for an osvrConfig.cmake file - most of the time it can be
# autodetected but you might need to specify osvr_DIR to be something like
# C:/Users/Ryan/Desktop/build/OSVR-Core-vc12 or
# C:/Users/Ryan/Downloads/OSVR-Core-Snapshot-v0.1-406-gaa55515-build54-vs12-32bit
# in the CMake GUI or command line
find_package(osvr REQUIRED)
find_package(jsoncpp REQUIRED)
if(TARGET jsoncpp_lib_static AND NOT TARGET jsoncpp_lib)
	add_library(jsoncpp_lib ALIAS jsoncpp_lib_static)
endif()

find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(SDL2 REQUIRED)
find_package(osvrRenderManager REQUIRED)
include_directories(${OSVRRENDERMANAGER_INCLUDE_DIRS})

#-----------------------------------------------------------------------------
add_executable(PresentPatternRenderManager PresentPatternRenderManager.cpp)
target_link_libraries(PresentPatternRenderManager PRIVATE ${OSVRRENDERMANAGER_LIBRARIES} osvr::osvrClientKitCpp osvr::osvrClient osvr::osvrCommon ${OPENGL_LIBRARY} GLEW::GLEW SDL2::SDL2 jsoncpp_lib)

