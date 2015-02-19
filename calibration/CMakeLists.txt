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
if(TARGET GLEW::GLEW_static)
    target_link_libraries(distortionizer-calibration GLEW::GLEW_static)
elseif(TARGET GLEW::GLEW)
    target_link_libraries(distortionizer-calibration GLEW::GLEW)
endif()

install(TARGETS distortionizer-calibration
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
