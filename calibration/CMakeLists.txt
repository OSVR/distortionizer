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
add_executable(distortionizer-calibration ${SOURCES} ${SHADERS_SOURCES})
target_link_libraries(distortionizer-calibration Qt5::Widgets Qt5::OpenGL)
