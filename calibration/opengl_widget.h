/** @file
    @brief Header

    @date 2014

    @author
    Russell Taylor
    <russ@reliasolve.com>
    <http://sensics.com>
*/

// Copyright 2014 Sensics, Inc.
//
// All rights reserved.
//
// (Final version intended to be licensed under
// the Apache License, Version 2.0)

#pragma once
#include "opengl_widget.h"
#include <QGLWidget>

// There are three different indices of refraction for the three
// different wavelengths in the head-mounted display (R, G, B).
// This is equivalent to having lenses with three different
// powers, one per color.  This produces effectively three different
// radial distortion patterns, one per color.  The distortion of
// red is the least, green next, and blue the most.
//
// The relative values (and ratios) of these distortions depends
// on the particular material that the lenses are made of, and is
// based on an emperical fit.  This means that we must determine
// each of them independently, subject to the above inequality
// constraint.
//
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
//
// When there are nonzero coefficients for higher-order terms
// (K2 and above), the result is a fourth-order polynomial that
// is challenging to invert analytically.

class OpenGL_Widget : public QGLWidget
{
    Q_OBJECT

public:
    OpenGL_Widget(QWidget *parent = 0);
    ~OpenGL_Widget();

public slots:

signals:

protected:
    //------------------------------------------------------
    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void keyPressEvent(QKeyEvent *event);

    //------------------------------------------------------
    bool saveConfigToJson(QString filename);
    bool loadConfigFromJson(QString filename);

    //------------------------------------------------------
    // Used as options in the rendering, depending on our
    // mode.
    void drawCrossHairs();
    void drawGrid();
    void drawCircles();

    //------------------------------------------------------
    // Helper functions for the draw routines.

    /// Draw a line from the specified begin point to the
    // specified end, doing distortion correcton.  The line
    // is drawn in short segments, with the correction applied
    // to each segment endpoint.  The color index tells whether
    // we use red (0), green (1), or blue (2) correction factors.
    // It uses the specified center of projection.
    void drawCorrectedLine(QPoint begin, QPoint end,
                           QPoint cop, unsigned color);
    void drawCorrectedCircle(QPoint center, float radius,
                           QPoint cop, unsigned color);

    /// Draw a set of 3 colored lines from the specified begin
    // point to the specified end, doing distortion correcton.
    void drawCorrectedLines(QPoint begin, QPoint end, QPoint cop);
    void drawCorrectedCircles(QPoint center, float radius, QPoint cop);

    /// Transform the specified pixel coordinate by the
    // color-correction distortion matrix using the appropriate
    // distortion correction.  The color index tells whether
    // we use red (0), green (1), or blue (2) correction factors.
    // It uses the specified center of projection.
    QPointF transformPoint(QPointF p, QPoint cop, unsigned color);

private:
    int d_width, d_height;  //< Size of the window we're rendering into
    QPoint d_cop_l;         //< Center of projection for the left eye
    QPoint d_cop_r;         //< Center of projection for the right eye
    float  d_k1_red;        //< Quadratic term for distortion of red
    float  d_k1_green;      //< Quadratic term for distortion of green
    float  d_k1_blue;       //< Quadratic term for distortion of blue
};
