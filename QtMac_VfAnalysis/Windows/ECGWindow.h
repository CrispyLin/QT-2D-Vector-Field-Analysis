#pragma once
#ifndef ECGWINDOW_H
#define ECGWINDOW_H

#include "GL_LIB/glew.h"

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QSurfaceFormat>
#include <QDebug>
#include <QColor>
#include <QDir>
#include <stdlib.h>
#include <OpenGL/GL.h>
#include <OpenGl/GLU.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>

#include "Predefined.h"
#include "BuildGeom/Geometry.h"
#include "Analysis/MorseDecomp.h"
#include "VField.h"
#include "StreamlineCalculate/EvenStreamlines.h"

class ECGWindow : public QOpenGLWidget, public QOpenGLFunctions
{
    Q_OBJECT
public:
    // constructor and desctuctor
    ECGWindow(QWidget *parent = nullptr);
    ~ECGWindow();

    // public member variables
    bool ShowECGOn = false;

    // public member functions
    void draw_nodes(GLenum mode);
    void draw_solid_circle(double cx, double cy, double R);
    void set_color(int nodetype);
    void draw_edges();
    void draw_wings(double head[2], icVector2 direct);
    void display_label(int x, int y, char* string);
    void draw_labels();
    void DrawGLScene(GLenum mode);
    void draw_circle(double cx, double cy, double radius, int type);
    void update_scene();
protected:
    void    initializeGL() override;
    void	paintGL() override;
    void	resizeGL(int w, int h) override;
};

#endif // ECGWINDOW_H
