cmake_minimum_required(VERSION 3.1.0)
project(distortUsingRenderManager)

#-----------------------------------------------------------------------------
# Local CMake Modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set (CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH};C:/Program Files/jsoncpp;C:/Program Files/Microsoft SDKs/Windows/v6.0;C:/usr/local CACHE PATH "Path that CMake uses to find packages")

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

find_package(VRPN REQUIRED)	# For time-manipulation routines
find_package(quatlib REQUIRED)	# For transformation construction
include_directories(${VRPN_INCLUDE_DIRS})

find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(SDL2 REQUIRED)
find_package(osvrRenderManager REQUIRED)
include_directories(${OSVRRENDERMANAGER_INCLUDE_DIRS})

#-----------------------------------------------------------------------------
# OpenGL Example program, which should eventually be open source
add_executable(DistortionCorrectRenderManager DistortionCorrectRenderManager.cpp font.c)
# Surprisingly, this also lets it know where to find the header files.
target_link_libraries(DistortionCorrectRenderManager PRIVATE ${OSVRRENDERMANAGER_LIBRARIES} osvr::osvrClientKitCpp osvr::osvrClient osvr::osvrCommon ${OPENGL_LIBRARY} GLEW::GLEW SDL2::SDL2 jsoncpp_lib ${VRPN_LIBRARIES})
