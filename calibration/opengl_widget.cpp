#include <QtGui>
#include <QtOpenGL>
#include <QColor>
#include <math.h>

#include "opengl_widget.h"

#ifndef GL_MULTISAMPLE
#define GL_MULTISAMPLE  0x809D
#endif

//----------------------------------------------------------------------
// Helper functions

OpenGL_Widget::OpenGL_Widget(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
    , d_cop(QPoint(0,0))
    , k1_red(1)
    , dk_green(1)
    , dk_blue(1)
{
}

OpenGL_Widget::~OpenGL_Widget()
{
}

void OpenGL_Widget::initializeGL()
{
    qglClearColor(Qt::black);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_MULTISAMPLE);
    static GLfloat lightPosition[4] = { 0.5, 5.0, 7.0, 1.0 };
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);

    // Makes the colors for the primitives be what we want.
    glDisable(GL_LIGHTING);

    // Tell the version number.  If we do this in the constructor, it doesn't
    // get displayed.
//    emit newVersionLabel(tr("Version 1.0.0"));
}

void OpenGL_Widget::drawCircle(QPoint center, float radius)
{
    glBegin(GL_LINE_STRIP);
    float step = 1 / radius;
    for (float r = 0; r <= 2*M_PI; r += step) {
        glVertex2f(center.x() + radius * cos(r),
                   center.y() + radius * sin(r));
    }
    glEnd();
}

void OpenGL_Widget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(0.0, 0.0, -10.0);

    // Set up rendering state.
    glDisable(GL_TEXTURE_2D);

    // Turn on blending, so that we'll get white
    // lines when we draw three different-colored lines
    // in the same location.  Also turn off Z-buffer
    // test so we get all of the lines drawn.
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);
    glDisable(GL_DEPTH_TEST);

    // Find the mirror of the left-eye's center of projection
    // around the screen center to find the right eye's COP.
    QPoint cop_r = QPoint(d_width - d_cop.x(), d_cop.y());

    // Draw two perpendicular lines through the center of
    // projection on the left eye, and the right eye.
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINES);
        glVertex2f(0, d_cop.y());
        glVertex2f(d_width/2, d_cop.y());
        glVertex2f(d_cop.x(), 0);
        glVertex2f(d_cop.x(), d_height);

        glVertex2f(d_width/2, cop_r.y());
        glVertex2f(d_width, cop_r.y());
        glVertex2f(cop_r.x(), 0);
        glVertex2f(cop_r.x(), d_height);
    glEnd();

    // Draw three colored circles around the center of projection
    // that go most of the way across each half of the screen.  The
    // spacing is determined by the dk_green and dk_blue factors.
    // Because we're using blending, this will result in a white
    // circle when all three have the same radius and there is no
    // lens color aberration.
    float r_red = 0.75 * d_width/4;
    float r_green = r_red * dk_green;
    float r_blue = r_red * dk_blue;
    glColor3f(1.0, 0.0, 0.0);
    drawCircle(d_cop, r_red);
    drawCircle(cop_r, r_red);

    glColor3f(0.0, 1.0, 0.0);
    drawCircle(d_cop, r_green);
    drawCircle(cop_r, r_green);

    glColor3f(0.0, 0.0, 1.0);
    drawCircle(d_cop, r_blue);
    drawCircle(cop_r, r_blue);
}

void OpenGL_Widget::resizeGL(int width, int height)
{
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
        aspect = static_cast<float>(width)/height;
    }
    glOrtho(0, d_width-1, 0, d_height-1, 5.0, 15.0);
    glMatrixMode(GL_MODELVIEW);

    // Default center of projection is the center of the left half
    // of the screen.
    d_cop.setX(d_width / 4);
    d_cop.setY(d_height / 2);
}


void OpenGL_Widget::keyPressEvent(QKeyEvent *event)
{
    float color_shift_scale = 1 + 1.0 / (d_width/2);
    switch (event->key()) {
    case Qt::Key_Left:
        d_cop.setX(d_cop.x()-1);
        break;
    case Qt::Key_Right:
        d_cop.setX(d_cop.x()+1);
        break;
    case Qt::Key_Down:
        d_cop.setY(d_cop.y()-1);
        break;
    case Qt::Key_Up:
        d_cop.setY(d_cop.y()+1);
        break;
    case Qt::Key_G:
        if(event->modifiers() & Qt::ShiftModifier) {
            dk_green *= color_shift_scale;
        } else {
            dk_green /= color_shift_scale;
        }
        break;
    case Qt::Key_B:
        if(event->modifiers() & Qt::ShiftModifier) {
            dk_blue *= color_shift_scale;
        } else {
            dk_blue /= color_shift_scale;
        }
        break;
    }
    updateGL();
}

void OpenGL_Widget::mousePressEvent(QMouseEvent *event)
{
    if (event->pos().x() < d_width/2) {
        d_cop = event->pos();
        d_cop.setY(d_height - d_cop.y());
    }
    updateGL();
//    lastPos = event->pos();
}

void OpenGL_Widget::mouseMoveEvent(QMouseEvent *event)
{
    //int dx = event->x() - lastPos.x();
    //int dy = event->y() - lastPos.y();

    if (event->buttons() & Qt::LeftButton) {
        // XXX
    } else if (event->buttons() & Qt::RightButton) {
        // XXX
    }
//    lastPos = event->pos();
}
