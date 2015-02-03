cmake_minimum_required(VERSION 2.8.12)
project(distortionizer)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
###
# Configuration Options
###

###
# CMake Modules
###
include(EnableExtraCompilerWarnings)

include(SetDefaultBuildType)
set_default_build_type(RelWithDebInfo)

include(MSVCMultipleProcessCompile)

include(InstallRequiredSystemLibraries)

###
# Dependencies
###
# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
find_package(Qt5Widgets REQUIRED)
find_package(Qt5OpenGL REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)

###
# Building
###
add_subdirectory(calibration)
