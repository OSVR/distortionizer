cmake_minimum_required(VERSION 3.1.0)
project(PolynomialFit)

#-----------------------------------------------------------------------------
# Local CMake Modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(jsoncpp)
if(TARGET jsoncpp_lib_static AND NOT TARGET jsoncpp_lib)
    add_library(jsoncpp_lib INTERFACE)
    target_link_libraries(jsoncpp_lib INTERFACE jsoncpp_lib_static)
endif()

#-----------------------------------------------------------------------------
add_executable(PolynomialFit PolynomialFit.cpp)
