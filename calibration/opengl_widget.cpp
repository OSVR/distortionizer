/** @file
    @brief Implementation

    @date 2014

    @author
    Russell Taylor
    <russ@reliasolve.com>
    <http://sensics.com>
*/

// Copyright 2014 Sensics, Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "opengl_widget.h"

#include <QColor>
#include <QFileDialog>
#include <QtGui>
#include <QtOpenGL>
#include <iostream>
#include <math.h>
#include <stdio.h>

#ifndef GL_MULTISAMPLE
#define GL_MULTISAMPLE 0x809D
#endif

#define CONFIG_FILE "HMD_Config.json"

//----------------------------------------------------------------------
// Helper functions

OpenGL_Widget::OpenGL_Widget(QWidget* parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent),
      d_cop_l(QPoint(0, 0)),
      d_cop_r(QPoint(0, 0)),
      d_cop(QPoint(0, 0)),
      d_k1_red(0),
      d_k1_green(0),
      d_k1_blue(0),
      fullscreen(false) {
    using namespace std;
    cout << "Distortion estimation for HMD using K1 (quadratic) term" << endl
         << "The program always runs on the last screen, full screen" << endl
         << "Click with the mouse in the left eye to locate center" << endl
         << "Keyboard controls:" << endl
         << "  r/R: Increase/Decrease distortion in R+G+B" << endl
         << "  g/G: Increase/Decrease distortion in G+B" << endl
         << "  b/B: Increase/Decrease distortion in B only" << endl
         << "  S/L: Save/Load state from JSON config file (" << CONFIG_FILE << " by default)" << endl
         << "  Left,Right,Up,Down: Move the center of projection by one pixel" << endl
         << "  f/F:    Toggle fullscreen on/off" << endl
         << "  c/C:    Reset center of projection" << endl
         << "  v/V:    Reset distortion values to 0" << endl
         << "  ESC/Q: Quit the application" << endl
         << endl;
}

OpenGL_Widget::~OpenGL_Widget() {}

void OpenGL_Widget::initializeGL() {
    qglClearColor(Qt::black);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_MULTISAMPLE);
    static GLfloat lightPosition[4] = {0.5, 5.0, 7.0, 1.0};
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);

    // Makes the colors for the primitives be what we want.
    glDisable(GL_LIGHTING);
}

// The distortion is with respect to a center of projection, which
// is common for all three colors.  The first-order correction that
// can be applied assumes that there is an additional radial shift
// that is proportionaly to the square of the distance of a pixel
// from the center of projection, with the coefficient called K1:
//    RcorrR = Rinit + K1R * (Rinit * Rinit)
//    RcorrG = Rinit + K1G * (Rinit * Rinit)
//    RcorrB = Rinit + K1B * (Rinit * Rinit)
//    0 <= K1R <= K1G <= K1B
//
// We can solve for a pre-distorted position that would produce
// a desired end location after passing through the optical system
// by solving for Rinit in the above equations (here done for a
// single color):
//    K1 * Rinit^2 + Rinit - Rcorr = 0
//    Solving for the positive root of the quadratic equation
//       A = K1, B = 1, C = -Rcorr
//    Rinit = (-B + sqrt(B^2 - 4AC)) / 2A
//          = (-1 + sqrt(1 + 4*K1*Rcorr)) / (2*K1)
//    K1 > 0

QPointF OpenGL_Widget::transformPoint(QPointF p, QPoint cop, unsigned color) {
    QPointF ret = p;
    QPointF offset = p - cop;
    float r2 = offset.x() * offset.x() + offset.y() * offset.y();
    float r = sqrt(r2);
    float k1;
    switch (color) {
    case 0:
        k1 = d_k1_red;
        break;
    case 1:
        k1 = d_k1_green;
        break;
    case 2:
        k1 = d_k1_blue;
        break;
    }
    k1 = k1 / ((d_width / 4.0) * (d_width / 4.0) * 16);

    // We will calculate the transformed point
    // by calculating the new location using the
    // following formula
    // x_d = distorted x; y_d = distorted y
    // x_u = undistorted x; y_u = undistorted y
    // x_c = center x; y_c = center y
    // r = sqrt( (x_u - x_c)^2 +  (y_u - y_c)^2 )
    // x_d = x_u [(1 + k1*r^2) * (x_u - x_c) / r]
    // y_d = y_u [(1 + k1*r^2) * (y_u - y_c) / r]
    ret = cop + (1 - k1 * r * r) * offset;
    return ret;
}

void OpenGL_Widget::drawCorrectedLine(QPoint begin, QPoint end, QPoint cop, unsigned color) {
    QPointF offset = end - begin;
    float len = sqrt(offset.x() * offset.x() + offset.y() * offset.y());
    QPointF offset_dir = offset / len;
    glBegin(GL_LINE_STRIP);
    for (float s = 0; s <= len; s++) {
        QPointF p = begin + s * offset_dir;
        QPointF tp = transformPoint(p, cop, color);
        glVertex2f(tp.x(), tp.y());
    }
    glEnd();
}

void OpenGL_Widget::drawCorrectedCircle(QPoint center, float radius, QPoint cop, unsigned color) {
    glBegin(GL_LINE_STRIP);
    float step = 1 / radius;
    for (float r = 0; r <= 2 * M_PI; r += step) {
        QPointF p(center.x() + radius * cos(r), center.y() + radius * sin(r));
        QPointF tp = transformPoint(p, cop, color);
        glVertex2f(tp.x(), tp.y());
    }
    glEnd();
}

void OpenGL_Widget::drawCorrectedLines(QPoint begin, QPoint end, QPoint cop) {
    float bright = 0.5f;
    glColor3f(bright, 0.0, 0.0);
    drawCorrectedLine(begin, end, cop, 0);

    glColor3f(0.0, bright, 0.0);
    drawCorrectedLine(begin, end, cop, 1);

    glColor3f(0.0, 0.0, bright);
    drawCorrectedLine(begin, end, cop, 2);
}

void OpenGL_Widget::drawCorrectedCircles(QPoint center, float radius, QPoint cop) {
    float bright = 0.5f;
    glColor3f(bright, 0.0, 0.0);
    drawCorrectedCircle(center, radius, cop, 0);

    glColor3f(0.0, bright, 0.0);
    drawCorrectedCircle(center, radius, cop, 1);

    glColor3f(0.0, 0.0, bright);
    drawCorrectedCircle(center, radius, cop, 2);
}

void OpenGL_Widget::drawCrossHairs() {
    // Draw two perpendicular lines through the center of
    // projection on the left eye, and the right eye.
    glColor3f(1.0, 1.0, 1.0);

    if (fullscreen) {
        glBegin(GL_LINES);
        glVertex2f(0, d_cop.y());
        glVertex2f(d_width, d_cop.y());
        glVertex2f(d_cop.x(), 0);
        glVertex2f(d_cop.x(), d_height);
        glEnd();
    } else {
        glBegin(GL_LINES);
        glVertex2f(0, d_cop_l.y());
        glVertex2f(d_width / 2, d_cop_l.y());
        glVertex2f(d_cop_l.x(), 0);
        glVertex2f(d_cop_l.x(), d_height);

        glVertex2f(d_width / 2, d_cop_r.y());
        glVertex2f(d_width, d_cop_r.y());
        glVertex2f(d_cop_r.x(), 0);
        glVertex2f(d_cop_r.x(), d_height);
        glEnd();
    }
}

void OpenGL_Widget::drawGrid() {
    // Draw a set of vertical grid lines to the right and left
    // of the center of projection for each eye.  Draw a red,
    // green, and blue line at each location with less than
    // full brightness.  Draw from the top of the screen to
    // the bottom.
    int spacing = 40;

    if (fullscreen) {
        // Vertical lines
        for (int r = spacing; r < d_width; r += spacing) {
            // Vertical lines, left eye
            if (d_cop.x() + r < d_width) {
                QPoint begin(d_cop.x() + r, 0);
                QPoint end(d_cop.x() + r, d_height);
                drawCorrectedLines(begin, end, d_cop);
            }
            if (d_cop.x() - r >= 0) {
                QPoint begin(d_cop.x() - r, 0);
                QPoint end(d_cop.x() - r, d_height);
                drawCorrectedLines(begin, end, d_cop);
            }
        }

        // Horizontal lines
        for (int r = spacing; r < d_height; r += spacing) {
            // Horizontal lines, left eye
            if (d_cop.y() + r < d_height) {
                QPoint begin(0, d_cop.y() + r);
                QPoint end(d_width, d_cop.y() + r);
                drawCorrectedLines(begin, end, d_cop);
            }
            if (d_cop.y() - r >= 0) {
                QPoint begin(0, d_cop.y() - r);
                QPoint end(d_width, d_cop.y() - r);
                drawCorrectedLines(begin, end, d_cop);
            }
        }
    }

    else {
        // Vertical lines
        for (int r = spacing; r < d_width / 2; r += spacing) {
            // Vertical lines, left eye
            if (d_cop_l.x() + r < d_width / 2) {
                QPoint begin(d_cop_l.x() + r, 0);
                QPoint end(d_cop_l.x() + r, d_height);
                drawCorrectedLines(begin, end, d_cop_l);
            }
            if (d_cop_l.x() - r >= 0) {
                QPoint begin(d_cop_l.x() - r, 0);
                QPoint end(d_cop_l.x() - r, d_height);
                drawCorrectedLines(begin, end, d_cop_l);
            }

            // Vertical lines, right eye
            if (d_cop_r.x() + r < d_width) {
                QPoint begin(d_cop_r.x() + r, 0);
                QPoint end(d_cop_r.x() + r, d_height - 1);
                drawCorrectedLines(begin, end, d_cop_r);
            }
            if (d_cop_r.x() - r >= d_width / 2) {
                QPoint begin(d_cop_r.x() - r, 0);
                QPoint end(d_cop_r.x() - r, d_height - 1);
                drawCorrectedLines(begin, end, d_cop_r);
            }
        }

        // Horizontal lines
        for (int r = spacing; r < d_height; r += spacing) {
            // Horizontal lines, left eye
            if (d_cop_l.y() + r < d_height) {
                QPoint begin(0, d_cop_l.y() + r);
                QPoint end(d_width / 2 - 1, d_cop_l.y() + r);
                drawCorrectedLines(begin, end, d_cop_l);
            }
            if (d_cop_l.y() - r >= 0) {
                QPoint begin(0, d_cop_l.y() - r);
                QPoint end(d_width / 2 - 1, d_cop_l.y() - r);
                drawCorrectedLines(begin, end, d_cop_l);
            }

            // Horizontal lines, right eye
            if (d_cop_r.y() + r < d_height) {
                QPoint begin(d_width / 2, d_cop_r.y() + r);
                QPoint end(d_width - 1, d_cop_r.y() + r);
                drawCorrectedLines(begin, end, d_cop_r);
            }
            if (d_cop_r.y() - r >= 0) {
                QPoint begin(d_width / 2, d_cop_r.y() - r);
                QPoint end(d_width - 1, d_cop_r.y() - r);
                drawCorrectedLines(begin, end, d_cop_r);
            }
        }
    }
}

void OpenGL_Widget::drawCircles() {
    if (fullscreen) {
        drawCorrectedCircles(d_cop, 0.1 * d_width / 4, d_cop);
        drawCorrectedCircles(d_cop, 0.7 * d_width / 4, d_cop);
    } else {
        drawCorrectedCircles(d_cop_l, 0.1 * d_width / 4, d_cop_l);
        drawCorrectedCircles(d_cop_r, 0.1 * d_width / 4, d_cop_r);
        drawCorrectedCircles(d_cop_l, 0.7 * d_width / 4, d_cop_l);
        drawCorrectedCircles(d_cop_r, 0.7 * d_width / 4, d_cop_r);
    }
}

void OpenGL_Widget::paintGL() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(0.0, 0.0, -10.0);

    // Set up rendering state.
    // Turn on blending, so that we'll get white
    // lines when we draw three different-colored lines
    // in the same location.  Also turn off Z-buffer
    // test so we get all of the lines drawn.  Also,
    // turn off texture-mapping.
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_TEXTURE_2D);

    drawCrossHairs();
    drawGrid();
    drawCircles();
}

void OpenGL_Widget::resizeGL(int width, int height) {
    d_width = width;
    d_height = height;
    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    // Make the window one unit high (-0.5 to 0.5) and have an aspect ratio that matches
    // the aspect ratio of the window.  We also make the left side of the window be at
    // the origin.
    float aspect;
    if ((height <= 0) || (width < 0)) {
        aspect = 1.0;
    } else {
        aspect = static_cast<float>(width) / height;
    }
    glOrtho(0, d_width - 1, 0, d_height - 1, 5.0, 15.0);
    glMatrixMode(GL_MODELVIEW);

    setDeftCOPVals();
}

void OpenGL_Widget::setDeftCOPVals() {

    // Default center of projection is the center of the left half
    // of the screen.
    d_cop_l.setX(d_width / 4);
    d_cop_l.setY(d_height / 2);

    // Find the mirror of the left-eye's center of projection
    // around the screen center to find the right eye's COP.
    d_cop_r = QPoint(d_width - d_cop_l.x(), d_cop_l.y());

    // Default center of projection for fullscreen mode is center of the screen
    d_cop.setX(d_width / 2);
    d_cop.setY(d_height / 2);
}

void OpenGL_Widget::keyPressEvent(QKeyEvent* event) {
    // Color shift equivalent to one pixel at the edge of
    // the screen if the center of projection is in the
    // middle of the screen.
    // float color_shift = 0.001 / ((d_width/4.0)*(d_width/4.0));
    float color_shift = 0.001f;
    // std::string fileName;
    // QStringList arguments = qApp->arguments();
    // if (arguments.size() > 1){
    //    fileName = qPrintable(arguments.at(1));
    //    //printf("Supplied file name as %s", fileName);
    //    std::cout << "file name is " << fileName << endl;
    //}
    switch (event->key()) {
    case Qt::Key_Escape:
    case Qt::Key_Q:
        QApplication::quit();
        break;
    case Qt::Key_S: // Save the state to an output file.
        // XXX Would like to throw a dialog box, but it shows in HMD
        // and cannot be moved.
        saveConfigToJson(CONFIG_FILE);
        break;
    case Qt::Key_L: // Load the state from an output file.
        // XXX Would like to throw a dialog box, but it shows in HMD
        // and cannot be moved.
        loadConfigFromJson(CONFIG_FILE);
        break;
    case Qt::Key_Left:

        if (fullscreen) {
            d_cop.setX(d_cop.x() - 1);
        } else {
            d_cop_l.setX(d_cop_l.x() - 1);

            // Find the mirror of the left-eye's center of projection
            // around the screen center to find the right eye's COP.
            d_cop_r = QPoint(d_width - d_cop_l.x(), d_cop_l.y());
        }
        break;
    case Qt::Key_Right:

        if (fullscreen) {
            d_cop.setX(d_cop.x() + 1);
        } else {
            d_cop_l.setX(d_cop_l.x() + 1);

            // Find the mirror of the left-eye's center of projection
            // around the screen center to find the right eye's COP.
            d_cop_r = QPoint(d_width - d_cop_l.x(), d_cop_l.y());
        }
        break;
    case Qt::Key_Down:

        if (fullscreen) {
            d_cop.setY(d_cop.y() - 1);
        } else {
            d_cop_l.setY(d_cop_l.y() - 1);

            // Find the mirror of the left-eye's center of projection
            // around the screen center to find the right eye's COP.
            d_cop_r = QPoint(d_width - d_cop_l.x(), d_cop_l.y());
        }
        break;
    case Qt::Key_Up:

        if (fullscreen) {
            d_cop.setY(d_cop.y() + 1);
        } else {
            d_cop_l.setY(d_cop_l.y() + 1);

            // Find the mirror of the left-eye's center of projection
            // around the screen center to find the right eye's COP.
            d_cop_r = QPoint(d_width - d_cop_l.x(), d_cop_l.y());
        }
        break;
    case Qt::Key_R:
        if (event->modifiers() & Qt::ShiftModifier) {
            d_k1_red -= color_shift;
        } else {
            d_k1_red += color_shift;
        }
    // Do not break here: Red shift also affects green and blue
    // break;
    case Qt::Key_G:
        if (event->modifiers() & Qt::ShiftModifier) {
            d_k1_green -= color_shift;
        } else {
            d_k1_green += color_shift;
        }
    // Do not break here: Green shift also affects blue
    // break;
    case Qt::Key_B:
        if (event->modifiers() & Qt::ShiftModifier) {
            d_k1_blue -= color_shift;
        } else {
            d_k1_blue += color_shift;
        }
        break;
    case Qt::Key_F:
        if (fullscreen) {
            fullscreen = false;
        } else {
            fullscreen = true;
        }
        break;
    case Qt::Key_C:
        // reset center to default values
        setDeftCOPVals();
        break;
    case Qt::Key_V:
        // reset distortion values to original
        d_k1_red = 0.0;
        d_k1_green = 0.0;
        d_k1_blue = 0.0;

        break;
    }

    printf("R = %g, G = %g, B = %g;\n", d_k1_red, d_k1_green, d_k1_blue);
    if (fullscreen) {
        printf("Center coords: x: %d, y: %d\n", d_cop.x(), d_cop.y());
    } else {
        printf("Left eye coords: x: %d, y: %d; ", d_cop_l.x(), d_cop_l.y());
        printf("Right eye coords: x: %d, y: %d;\n", d_cop_r.x(), d_cop_r.y());
    }
    updateGL();
}

void OpenGL_Widget::mousePressEvent(QMouseEvent* event) {
    if (fullscreen) {
        d_cop = event->pos();
        d_cop.setY(d_height - d_cop.y());
    } else {
        if (event->pos().x() < d_width / 2) {
            d_cop_l = event->pos();
            d_cop_l.setY(d_height - d_cop_l.y());

            // Find the mirror of the left-eye's center of projection
            // around the screen center to find the right eye's COP.
            d_cop_r = QPoint(d_width - d_cop_l.x(), d_cop_l.y());
        }
    }
    updateGL();
    //    d_last_pos = event->pos();
}

void OpenGL_Widget::mouseMoveEvent(QMouseEvent* event) {

    if (event->buttons() & Qt::LeftButton) {
        // XXX
    } else if (event->buttons() & Qt::RightButton) {
        // XXX
    }
}

QPointF OpenGL_Widget::pixelToRelative(QPointF cop) {
    QPointF relative_cop;
    float cop_x, cop_y;
    cop_x = (float)cop.x() / (float)d_width;
    cop_y = (float)cop.y() / (float)d_height;
    relative_cop.setX(cop_x);
    relative_cop.setY(cop_y);

    return relative_cop;
}

QPoint OpenGL_Widget::relativeToPixel(QPointF cop) {

    QPoint pixel_cop;
    int cop_x, cop_y;
    cop_x = cop.x() * d_width;
    cop_y = cop.y() * d_height;
    pixel_cop.setX(cop_x);
    pixel_cop.setY(cop_y);

    return pixel_cop;
}

bool OpenGL_Widget::saveConfigToJson(QString filename) {
    FILE* f = fopen(filename.toStdString().c_str(), "w");
    if (f == NULL) {
        fprintf(stderr, "OpenGL_Widget::saveConfigToJson(): Can't save to %s", filename.toStdString().c_str());
        return false;
    }
    // convert from pixels to Relative screen size to use by shader
    QPointF relative_cop;
    if (fullscreen) {
        relative_cop = pixelToRelative(d_cop);
    } else {
        relative_cop = pixelToRelative(d_cop_l);
    }

    fprintf(f, "{\n");
    fprintf(f, "    \"hmd\": {\n");
    fprintf(f, "        \"distortion\": {\n");
    fprintf(f, "            \"k1_red\": %g,\n", d_k1_red);
    fprintf(f, "            \"k1_green\": %g,\n", d_k1_green);
    fprintf(f, "            \"k1_blue\": %g\n", d_k1_blue);
    fprintf(f, "        },\n");
    fprintf(f, "        \"eyes\": [\n");
    fprintf(f, "            {\n");
    fprintf(f, "                \"center_proj_x\": %f,\n", relative_cop.x());
    fprintf(f, "                \"center_proj_y\": %f\n", relative_cop.y());
    fprintf(f, "            }\n");
    fprintf(f, "        ],\n");
    fprintf(f, "        \"fullscreen\": %d\n", (int)fullscreen);
    /*
    else{
        fprintf(f, "    {\n");
        fprintf(f, "      \"Center_of_projection_pixels_x\": %d,\n", d_cop_l.x() / d_width);
        fprintf(f, "      \"Center_of_projection_pixels_y\": %d\n",  d_cop_l.y() / d_height);
        fprintf(f, "    },\n");
        fprintf(f, "    {\n");
        // Center of projection with respect to the lower-left corner of the frame.
        fprintf(f, "      \"Center_of_projection_pixels_x\": %d,\n", d_cop_r.x() - d_width / 2);
        fprintf(f, "      \"Center_of_projection_pixels_y\": %d\n", d_cop_r.y());
        fprintf(f, "    }\n");

    }
    */

    fprintf(f, "    }\n");
    fprintf(f, "}\n");

    fclose(f);
    return true;
}

bool OpenGL_Widget::loadConfigFromJson(QString filename) {
    FILE* f = fopen(filename.toStdString().c_str(), "r");
    if (f == NULL) {
        fprintf(stderr, "OpenGL_Widget::loadConfigFromJson(): Can't load from to %s", filename.toStdString().c_str());
        return false;
    }

    // Lines and parameters to read
    char line[1024];
    char param[1024];
    float val;
    QPointF cop;

    // Skip the first three lines.
    char* unused;
    unused = fgets(line, sizeof(line), f);
    unused = fgets(line, sizeof(line), f);
    unused = fgets(line, sizeof(line), f);

    // Try and read the red term
    if (fgets(line, sizeof(line), f) == NULL) {
        fprintf(stderr, "OpenGL_Widget::loadConfigFromJson(): Can't read red line");
        return false;
    }
    if (sscanf(line, "%s %g", param, &val) != 2) {
        fprintf(stderr, "OpenGL_Widget::loadConfigFromJson(): Bad red line: %s", line);
        return false;
    }
    d_k1_red = val;

    // Try and read the green term
    if (fgets(line, sizeof(line), f) == NULL) {
        fprintf(stderr, "OpenGL_Widget::loadConfigFromJson(): Can't read green line");
        return false;
    }
    if (sscanf(line, "%s %g", param, &val) != 2) {
        fprintf(stderr, "OpenGL_Widget::loadConfigFromJson(): Bad green line: %s", line);
        return false;
    }
    d_k1_green = val;

    // Try and read the blue term
    if (fgets(line, sizeof(line), f) == NULL) {
        fprintf(stderr, "OpenGL_Widget::loadConfigFromJson(): Can't read blue line");
        return false;
    }
    if (sscanf(line, "%s %g", param, &val) != 2) {
        fprintf(stderr, "OpenGL_Widget::loadConfigFromJson(): Bad blue line: %s", line);
        return false;
    }
    d_k1_blue = val;

    // Skip the next three lines
    unused = fgets(line, sizeof(line), f);
    unused = fgets(line, sizeof(line), f);
    unused = fgets(line, sizeof(line), f);

    // Try and read the eye X COP term
    if (fgets(line, sizeof(line), f) == NULL) {
        fprintf(stderr, "OpenGL_Widget::loadConfigFromJson(): Can't read COP x line");
        return false;
    }
    if (sscanf(line, "%s %f", param, &val) != 2) {
        fprintf(stderr, "OpenGL_Widget::loadConfigFromJson(): Bad COP x line: %s", line);
        return false;
    }
    cop.setX(val);

    // Try and read the eye Y COP term
    if (fgets(line, sizeof(line), f) == NULL) {
        fprintf(stderr, "OpenGL_Widget::loadConfigFromJson(): Can't read COP y line");
        return false;
    }
    if (sscanf(line, "%s %f", param, &val) != 2) {
        fprintf(stderr, "OpenGL_Widget::loadConfigFromJson(): Bad COP y line: %s", line);
        return false;
    }
    cop.setY(val);

    if (fullscreen) {
        d_cop = relativeToPixel(cop);
    } else {
        d_cop_l = relativeToPixel(cop);
        // Find the mirror of the left-eye's center of projection
        // around the screen center to find the right eye's COP.
        d_cop_r = QPoint(d_width - d_cop_l.x(), d_cop_l.y());
    }

    fclose(f);
    return true;
}
