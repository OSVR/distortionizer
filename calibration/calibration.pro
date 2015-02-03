#-------------------------------------------------
#
# Project created by QtCreator 2014-12-22T17:42:44
#
#-------------------------------------------------

# Modified to work with Qt5.

QT       += core gui opengl widgets
CONFIG += console

TARGET = calibration
TEMPLATE = app

INCLUDEPATH += /usr/local/include
INCLUDEPATH += C:/usr/local/include

INCLUDEPATH += ../shaders

# Avoid some warnings on Windows
DEFINES += _CRT_SECURE_NO_WARNINGS=1

LIBS += -LC:/usr/local/lib
LIBS += -lglew32s
DEFINES += GLEW_STATIC

SOURCES += main.cpp\
        mainwindow.cpp \
    opengl_widget.cpp \
    ../shaders/undistort_shader.cpp

HEADERS  += mainwindow.h \
    opengl_widget.h \
    ../shaders/undistort_shader.h

FORMS    += mainwindow.ui
