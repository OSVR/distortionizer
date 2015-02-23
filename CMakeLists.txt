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
list(APPEND CMAKE_PREFIX_PATH "C:\\Qt\\5.4\\msvc2013_64_opengl\\")
set(CMAKE_AUTOMOC ON)
find_package(Qt5Widgets REQUIRED)
find_package(Qt5OpenGL REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)


# For generating documentation in HTML
if(WIN32)
    list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/third-party/discount-windows-bins")
endif()
find_package(Markdown)

###
# Docs
###
# Convert the README to HTML if we can
if(MARKDOWN_FOUND)
    include(UseMarkdown)
    add_markdown_target(markdown_docs "${CMAKE_CURRENT_BINARY_DIR}" README.md)

    install_markdown_target(markdown_docs DESTINATION ${CMAKE_INSTALL_DOCDIR})
else()
    install(FILES README.md DESTINATION ${CMAKE_INSTALL_DOCDIR})
endif()

# Install the fragment shaders as documentation.
install(FILES vizard/ShaderTest.frag vizard/ShaderTest.vert DESTINATION ${CMAKE_INSTALL_DOCDIR})

###
# Building
###
add_subdirectory(calibration)
