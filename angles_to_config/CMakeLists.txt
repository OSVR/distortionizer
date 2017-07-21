cmake_minimum_required(VERSION 3.1.0)
project(AnglesToConfig)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#-----------------------------------------------------------------------------
# Local CMake Modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake" "${CMAKE_CURRENT_SOURCE_DIR}/../cmake-local")
###
# Configuration Options
###

option(BUILD_WITH_RENDERMANAGER "Build apps that use RM" OFF)

find_package(OpenGL)
find_package(GLEW)
#find_package(SDL2)
find_package(osvrRenderManager)
find_package(osvr)
find_package(JsonCpp)

include(GNUInstallDirs)

add_subdirectory(../vendor vendorlibs)

#-----------------------------------------------------------------------------
add_library(AnglesToConfig-core STATIC
    EigenStdArrayInterop.h
    GenerateOutput.cpp
    GenerateOutput.h
    helper.cpp
    helper.h
    InvalidPointRemoval.cpp
    InvalidPointRemoval.h
    Process.cpp
    Process.h
    Subproblems.cpp
    Subproblems.h
    types.h)

target_link_libraries(AnglesToConfig-core PUBLIC distort-eigen-vendored)

add_executable(AnglesToConfig
    AnglesToConfig.cpp)

target_link_libraries(AnglesToConfig PRIVATE AnglesToConfig-core)

add_executable(MakeExampleMesh
    MakeExampleMesh.cpp)

install(TARGETS AnglesToConfig MakeExampleMesh
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(JSONCPP_FOUND)
    add_executable(AnglesToConfigModular
        AnglesToConfigModular.cpp)
    target_link_libraries(AnglesToConfigModular PRIVATE AnglesToConfig-core JsonCpp::JsonCpp)
    install(TARGETS AnglesToConfigModular
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

if (OPENGL_FOUND AND GLEW_FOUND AND SDL2_FOUND AND BUILD_WITH_RENDERMANAGER AND OSVRRENDERMANAGER_FOUND)
    add_executable(DebugAnglesToConfig DebugAnglesToConfig.cpp)
    target_include_directories(DebugAnglesToConfig PRIVATE ${OSVRRENDERMANAGER_INCLUDE_DIRS})
    target_link_libraries(DebugAnglesToConfig PRIVATE AnglesToConfig-core osvrRenderManager::osvrRenderManagerCpp osvr::osvrClientKitCpp osvr::osvrClient osvr::osvrCommon ${OPENGL_LIBRARY} GLEW::GLEW SDL2::SDL2 JsonCpp::JsonCpp)
endif()
