include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../shaders")

if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

set(SOURCES
    main.cpp
    mainwindow.cpp
    mainwindow.h
    opengl_widget.cpp
    opengl_widget.h)
set(SHADERS_SOURCES
    ../shaders/undistort_shader.cpp
    ../shaders/undistort_shader.h)
source_group(shaders FILES ${SHADERS_SOURCES})
qt5_wrap_ui(UI_HEADERS mainwindow.ui)

add_executable(distortionizer-calibration ${SOURCES} ${SHADERS_SOURCES} ${UI_HEADERS})

target_link_libraries(distortionizer-calibration Qt5::Widgets Qt5::OpenGL ${OPENGL_LIBRARIES})
install(TARGETS distortionizer-calibration
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

set(deps Qt5::Core Qt5::Gui Qt5::Widgets Qt5::OpenGL)

if(TARGET GLEW::GLEW_static)
    target_link_libraries(distortionizer-calibration GLEW::GLEW_static)
elseif(TARGET GLEW::GLEW)
    target_link_libraries(distortionizer-calibration GLEW::GLEW)
    list(APPEND deps GLEW::GLEW)
endif()

foreach(dep ${deps})
    copy_imported_targets(distortionizer-calibration ${dep})
    install_imported_target(${dep} DESTINATION ${CMAKE_INSTALL_BINDIR})
endforeach()

if(WIN32)
    install_imported_target(Qt5::QWindowsIntegrationPlugin DESTINATION ${CMAKE_INSTALL_BINDIR}/platforms)
    list(APPEND deps Qt5::QWindowsIntegrationPlugin)
endif()

include(${CMAKE_CURRENT_SOURCE_DIR}/GetSearchDirs.cmake)

get_search_dirs(libSearchDirs ${deps})

set(APPS
    "\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/distortionizer-calibration${CMAKE_EXECUTABLE_SUFFIX}")
install(CODE "include(BundleUtilities)
    fixup_bundle(\"${APPS}\" \"\" \"${libSearchDirs}\")")
