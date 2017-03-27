cmake_minimum_required(VERSION 3.1.0)
project(distortionizer)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" "${CMAKE_CURRENT_SOURCE_DIR}/cmake-local")
###
# Configuration Options
###

option(BUILD_WITH_RENDERMANAGER "Build apps that use RM" OFF)

###
# CMake Modules
###
include(EnableExtraCompilerWarnings)

include(SetDefaultBuildType)
set_default_build_type(RelWithDebInfo)

include(MSVCMultipleProcessCompile)
include(CopyImportedTarget)

include(InstallRequiredSystemLibraries)
if(WIN32)
    # Flat install for docs on Windows.
    set(CMAKE_INSTALL_DOCDIR .)
endif()
include(GNUInstallDirs)
if(NOT CMAKE_INSTALL_DOCDIR)
    set(CMAKE_INSTALL_DOCDIR ${CMAKE_INSTALL_DATAROOTDIR}/doc/${PROJECT_NAME})
endif()

###
# Dependencies
###
# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
find_package(Qt5 COMPONENTS Widgets OpenGL)
find_package(OpenGL REQUIRED)
find_package(GLEW)
if(NOT CMAKE_VERSION VERSION_LESS 3.6.0)
    find_program(IWYU_PROGRAM NAMES include-what-you-use)
    if(IWYU_PROGRAM)
        set(IWYU_DO TRUE)
    endif()
    find_program(
        CLANG_TIDY_PROGRAM
        NAMES "clang-tidy"
        DOC "Path to clang-tidy executable")
    if(CLANG_TIDY_PROGRAM)
        set(CLANG_TIDY_DO TRUE)
    endif()
endif()
find_package(osvrRenderManager)

if(IWYU_DO)
	set(CMAKE_CXX_INCLUDE_WHAT_YOU_USE "${IWYU_PROGRAM}")
	if(MSVC)
		list(APPEND CMAKE_CXX_INCLUDE_WHAT_YOU_USE "-fdiagnostics-format=msvc")
	endif()
endif()

# For generating documentation in HTML
if(WIN32)
    list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/third-party/discount-windows-bins")
endif()
find_package(Markdown)

###
# Docs
###
# Convert the README to HTML if we can
set(MARKDOWN_INPUT README.md CONTRIBUTING.md)
if(MARKDOWN_FOUND)
    include(UseMarkdown)
    add_markdown_target(markdown_docs "${CMAKE_CURRENT_BINARY_DIR}" ${MARKDOWN_INPUT})

    install_markdown_target(markdown_docs DESTINATION ${CMAKE_INSTALL_DOCDIR})
else()
    install(FILES ${MARKDOWN_INPUT} DESTINATION ${CMAKE_INSTALL_DOCDIR})
endif()

install(FILES LICENSE NOTICE DESTINATION ${CMAKE_INSTALL_DOCDIR})

# Install the fragment shaders as documentation.
install(FILES vizard/ShaderTest.frag vizard/ShaderTest.vert DESTINATION ${CMAKE_INSTALL_DOCDIR})

###
# Building
###
add_subdirectory(calibration)
add_subdirectory(angles_to_config)
