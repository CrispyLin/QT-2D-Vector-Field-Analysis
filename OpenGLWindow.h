#pragma once
#ifndef OPENGLWINDOW_H
#define OPENGLWINDOW_H
#include "openGL_lib/glew.h" // include glew.h first is necessary

#include <mainwindow.h>
#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QSurfaceFormat>
#include <QDebug>
#include <QColor>

#include "BuildGeom/Geometry.h"
#include "Analysis/MorseDecomp.h"
#include "VField.h"

class OpenGLWindow : public QOpenGLWidget, public QOpenGLFunctions
{
    Q_OBJECT
public:
    // public member variables
    int MoveOrStop ;

    icVector3 rot_center;
    double zoom_factor;
    float rotmat[4][4]; //storing rotation parameters
    float ObjXmat[16]; //Storing current modelview transformation

    // public member functions
    MainWindow* mainWindow = nullptr;
    Polyhedron* polyhedron = nullptr;

    // constructor and desctuctor
    OpenGLWindow(QWidget *parent = nullptr);
    ~OpenGLWindow();

    // public member functions
    void    set_up_MainWindow_ptr(MainWindow* MW_ptr);
    void    build_Polyhedron();
    void    draw();
    int     InitGL();
    void    init_flags();
protected:
    void    initializeGL() override;
    void	paintGL() override;
    void	resizeGL(int w, int h) override;
};

#endif // OPENGLWINDOW_H
