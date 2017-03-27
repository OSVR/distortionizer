cmake_minimum_required(VERSION 3.1.0)
project(ImageBasedDistortion)

#-----------------------------------------------------------------------------
# Local CMake Modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake" "${CMAKE_CURRENT_SOURCE_DIR}/../cmake-local")

#-----------------------------------------------------------------------------
# This looks for an osvrConfig.cmake file - most of the time it can be
# autodetected but you might need to specify osvr_DIR to be something like
# C:/Users/Ryan/Desktop/build/OSVR-Core-vc12 or
# C:/Users/Ryan/Downloads/OSVR-Core-Snapshot-v0.1-406-gaa55515-build54-vs12-32bit
# in the CMake GUI or command line
find_package(osvr REQUIRED)
find_package(JsonCpp REQUIRED)

find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
#find_package(SDL2 REQUIRED)
find_package(osvrRenderManager REQUIRED)

#-----------------------------------------------------------------------------

add_executable(PresentLinesRenderManager PresentLinesRenderManager.cpp)
target_link_libraries(PresentLinesRenderManager PRIVATE osvrRenderManager::osvrRenderManagerCpp osvr::osvrClientKitCpp osvr::osvrClient osvr::osvrCommon ${OPENGL_LIBRARY} GLEW::GLEW SDL2::SDL2 JsonCpp::JsonCpp)

add_executable(PresentPatternRenderManager PresentPatternRenderManager.cpp)
target_link_libraries(PresentPatternRenderManager PRIVATE osvrRenderManager::osvrRenderManagerCpp osvr::osvrClientKitCpp osvr::osvrClient osvr::osvrCommon ${OPENGL_LIBRARY} GLEW::GLEW SDL2::SDL2 JsonCpp::JsonCpp)
