cmake_minimum_required(VERSION 3.1.0)
project(AnglesToConfig)

#-----------------------------------------------------------------------------
# Local CMake Modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(OpenGL)
find_package(GLEW)
find_package(SDL2)
find_package(osvrRenderManager)
find_package(osvr)
find_package(jsoncpp)
if(TARGET jsoncpp_lib_static AND NOT TARGET jsoncpp_lib)
    add_library(jsoncpp_lib INTERFACE)
    target_link_libraries(jsoncpp_lib INTERFACE jsoncpp_lib_static)
endif()

#-----------------------------------------------------------------------------
add_executable(AnglesToConfig AnglesToConfig.cpp helper.cpp)

if (OPENGL_FOUND AND GLEW_FOUND AND SDL2_FOUND)
  include_directories(${OSVRRENDERMANAGER_INCLUDE_DIRS})

  add_executable(DebugAnglesToConfig DebugAnglesToConfig.cpp helper.cpp)
  target_link_libraries(DebugAnglesToConfig PRIVATE ${OSVRRENDERMANAGER_LIBRARIES} osvr::osvrClientKitCpp osvr::osvrClient osvr::osvrCommon ${OPENGL_LIBRARY} GLEW::GLEW SDL2::SDL2 jsoncpp_lib)
endif()
