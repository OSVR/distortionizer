#pragma once
#include <QGLWidget>

class OpenGL_Widget : public QGLWidget
{
    Q_OBJECT

public:
    OpenGL_Widget(QWidget *parent = 0);
    ~OpenGL_Widget();

public slots:

signals:

protected:
    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);

private:
    int d_width, d_height;  //< Size of the window we're rendering into
    QPoint  d_cop;          //< Center of projection
};

