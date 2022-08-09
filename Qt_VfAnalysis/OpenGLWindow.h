#pragma once
#ifndef OPENGLWINDOW_H
#define OPENGLWINDOW_H

#include "GL/glew.h"

#include <mainwindow.h>
#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QSurfaceFormat>
#include <QDebug>
#include <QColor>

#include "Predefined.h"
#include "BuildGeom/Geometry.h"
#include "Analysis/MorseDecomp.h"
#include "VField.h"

class OpenGLWindow : public QOpenGLWidget, public QOpenGLFunctions
{
    Q_OBJECT
public:
    // public member variables
    int MoveOrStop;
    int ShowFixedPtOn;
    int ShowSeparatricesOn;
    int ShowSCCsOn;
    int ShowPeriodicOrbitsOn;
    int EvenStreamlinePlacement;
    int ShowColorVFMagOn;
    int ShowMorseConn;
    int ShowTriangleConn;
    int morseC;
    int triC;
    int sccC;
    bool ShowEdgeSamplesOn;
    bool ShowTriMappingOn;
    int selected_triangle;
    int ShowBackward;
    int sampling_edge;
    int ShowConnectionRegion;
    int ShowBoundary;
    bool FlipNormalOn;
    bool IBFVOff;

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
