#ifndef MCGWINDOW_H
#define MCGWINDOW_H

#define SELECTBUFFERSIZE 128

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
#include <QMouseEvent>

#include "Predefined.h"
#include "BuildGeom/Geometry.h"
#include "Analysis/MorseDecomp.h"
#include "VField.h"
#include "StreamlineCalculate/EvenStreamlines.h"
#include "utility_functions.h"

class MCGWindow : public QOpenGLWidget, public QOpenGLFunctions
{
    Q_OBJECT
public:
    // constructor and desctuctor
    MCGWindow(QWidget *parent = nullptr);
    ~MCGWindow();


    // public member variables
    double s2_old = -1, t2_old = -1;
    bool ShowMCGOn = false;
    bool ShowConleyCircle = false;

    // public member functions
    void draw_nodes(GLenum mode);
    void draw_solid_circle(double cx, double cy, double R);
    void set_color(int nodetype);
    void draw_edges();
    void draw_wings(double head[2], icVector2 direct);
    void display_label(int x, int y, char* string);
    void draw_labels();
    void DrawGLScene(GLenum mode);
    void draw_highlights();
    void draw_circle(double cx, double cy, double radius, int type);
    void update_scene();

    void HitProcessforGraph(double ss, double tt);

protected:
    void    initializeGL() override;
    void	paintGL() override;
    void	resizeGL(int w, int h) override;

private:
    void mouseMoveEvent(QMouseEvent * event) override;
};

#endif // MCGWINDOW_H
