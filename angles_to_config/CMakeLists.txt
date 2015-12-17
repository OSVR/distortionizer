cmake_minimum_required(VERSION 3.1.0)
project(AnglesToConfig)

#-----------------------------------------------------------------------------
# Local CMake Modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#-----------------------------------------------------------------------------
add_executable(AnglesToConfig AnglesToConfig.cpp)

