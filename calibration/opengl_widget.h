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
};

