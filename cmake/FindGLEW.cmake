#.rst:
# FindGLEW
# --------
#
# Find the OpenGL Extension Wrangler Library (GLEW)
#
# IMPORTED Targets
# ^^^^^^^^^^^^^^^^
#
# This module defines the :prop_tgt:`IMPORTED` target ``GLEW::GLEW``,
# if GLEW has been found.
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module defines the following variables:
#
# ::
#
#   GLEW_INCLUDE_DIRS - include directories for GLEW
#   GLEW_LIBRARIES - libraries to link against GLEW
#   GLEW_FOUND - true if GLEW has been found and can be used

#=============================================================================
# Copyright 2012 Benjamin Eikel, 2015 Ryan Pavlik
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

find_path(GLEW_INCLUDE_DIR GL/glew.h)
if(WIN32 AND GLEW_INCLUDE_DIR)
  set(GLEW_ROOT "${GLEW_INCLUDE_DIR}/..")
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(GLEW_PLATFORM x64)
  else()
    set(GLEW_PLATFORM Win32)
  endif()

  find_library(GLEW_SHARED_LIBRARY
    NAMES glew32
    HINTS "${GLEW_ROOT}"
    PATH_SUFFIXES "lib/Release/${GLEW_PLATFORM}")
  find_file(GLEW_DLL glew32.dll
    HINTS "${GLEW_ROOT}"
    PATH_SUFFIXES "bin/Release/${GLEW_PLATFORM}")

  find_library(GLEW_STATIC_LIBRARY
    NAMES glew32s
    HINTS "${GLEW_ROOT}"
    PATH_SUFFIXES "lib/Release/${GLEW_PLATFORM}")

  find_library(GLEW_MX_SHARED_LIBRARY
    NAMES glew32mx
    HINTS "${GLEW_ROOT}"
    PATH_SUFFIXES "lib/Release MX/${GLEW_PLATFORM}")
  find_file(GLEW_MX_DLL glew32mx.dll
    HINTS "${GLEW_ROOT}"
    PATH_SUFFIXES "bin/Release MX/${GLEW_PLATFORM}")

  find_library(GLEW_MX_STATIC_LIBRARY
    NAMES glew32mxs
    HINTS "${GLEW_ROOT}"
    PATH_SUFFIXES "lib/Release MX/${GLEW_PLATFORM}")
  if(GLEW_SHARED_LIBRARY)
    set(GLEW_LIBRARY ${GLEW_SHARED_LIBRARY})
  elseif(GLEW_LIBRARY_STATIC)
    set(GLEW_LIBRARY ${GLEW_STATIC_LIBRARY})
  endif()
else()
  find_library(GLEW_LIBRARY NAMES GLEW glew32 glew glew32s PATH_SUFFIXES lib64)
  set(GLEW_LIBRARIES ${GLEW_LIBRARY})
endif()

set(GLEW_INCLUDE_DIRS ${GLEW_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GLEW
                                  REQUIRED_VARS GLEW_INCLUDE_DIR GLEW_LIBRARY)

if(GLEW_FOUND AND NOT TARGET GLEW::GLEW)
  if(WIN32)
    if(GLEW_SHARED_LIBRARY AND GLEW_DLL AND NOT TARGET GLEW::GLEW)
      add_library(GLEW::GLEW SHARED IMPORTED)
      set_target_properties(GLEW::GLEW PROPERTIES
        IMPORTED_IMPLIB "${GLEW_SHARED_LIBRARY}"
        IMPORTED_LOCATION "${GLEW_DLL}"
        INTERFACE_INCLUDE_DIRECTORIES "${GLEW_INCLUDE_DIRS}")
    endif()
    if(GLEW_MX_SHARED_LIBRARY AND GLEW_MX_DLL AND NOT TARGET GLEW::GLEW_MX)
      add_library(GLEW::GLEW_MX SHARED IMPORTED)
      set_target_properties(GLEW::GLEW_MX PROPERTIES
        IMPORTED_IMPLIB "${GLEW_MX_SHARED_LIBRARY}"
        IMPORTED_LOCATION "${GLEW_MX_DLL}"
        INTERFACE_INCLUDE_DIRECTORIES "${GLEW_INCLUDE_DIRS}")
    endif()

    if(GLEW_STATIC_LIBRARY AND NOT TARGET GLEW::GLEW_static)
      add_library(GLEW::GLEW_static STATIC IMPORTED)
      set_target_properties(GLEW::GLEW_static PROPERTIES
        IMPORTED_LOCATION "${GLEW_STATIC_LIBRARY}"
        INTERFACE_COMPILE_DEFINITIONS "GLEW_STATIC"
        INTERFACE_INCLUDE_DIRECTORIES "${GLEW_INCLUDE_DIRS}")
    endif()

    if(GLEW_MX_STATIC_LIBRARY AND NOT TARGET GLEW::GLEW_MX_static)
      add_library(GLEW::GLEW_MX_static STATIC IMPORTED)
      set_target_properties(GLEW::GLEW_MX_static PROPERTIES
        IMPORTED_LOCATION "${GLEW_MX_STATIC_LIBRARY}"
        INTERFACE_COMPILE_DEFINITIONS "GLEW_STATIC"
        INTERFACE_INCLUDE_DIRECTORIES "${GLEW_INCLUDE_DIRS}")
    endif()
  endif()
  if(NOT TARGET GLEW::GLEW)
    add_library(GLEW::GLEW UNKNOWN IMPORTED)
    set_target_properties(GLEW::GLEW PROPERTIES
      IMPORTED_LOCATION "${GLEW_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${GLEW_INCLUDE_DIRS}")
  endif()
endif()

mark_as_advanced(GLEW_INCLUDE_DIR GLEW_LIBRARY)
