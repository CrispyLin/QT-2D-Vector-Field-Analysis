#pragma once

#ifndef OPENGLWINDOW_H
#define OPENGLWINDOW_H


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

#include "MainWindow.h"
#include "Predefined.h"
#include "BuildGeom/Geometry.h"
#include "Analysis/MorseDecomp.h"
#include "VField.h"
#include "StreamlineCalculate/EvenStreamlines.h"
#include "RegionTauMap.h"
#include "utility_functions.h"
#include "Others/TraceBall.h"

#define	NPN 64
#define NMESH  100
#define DM  ((double) (1.0/(NMESH-1.0)))
#define NPIX  800 // number of pixels


typedef float Matrix[4][4];

class VectorFieldWindow : public QOpenGLWidget, public QOpenGLFunctions
{
    Q_OBJECT
public:
    // constructor and desctuctor
    VectorFieldWindow(QWidget *parent = nullptr);
    ~VectorFieldWindow();

protected:
    void    initializeGL() override;

    void	paintGL() override;
    void	resizeGL(int w, int h) override;

private:
    // for mouse events
    QPoint pos;
    double s_old = -1, t_old = -1;
    double last_x = -1, last_y = -1;

    void mousePressEvent(QMouseEvent * event) override ;
    void mouseReleaseEvent(QMouseEvent * event) override ;
    void mouseMoveEvent(QMouseEvent * event) override ;


    void rightButtonUp(QMouseEvent * event);
    void leftButtonMoved(QMouseEvent *event);
    void middleButtonDown(QMouseEvent *event);
    void middleButtonMoved(QMouseEvent *event);

    void wheelEvent(QWheelEvent * event) override;

public:
    // public member variables
    CTraceBall traceball;
    Quaternion rvec;

    bool    MoveOrStop;
    int     ShowFixedPtOn;
    int     ShowSeparatricesOn;
    int     ShowSCCsOn;
    int     ShowPeriodicOrbitsOn;
    bool    EvenStreamlinePlacement;
    bool    ShowColorVFMagOn;
    int     ShowMorseConn;
    int     ShowTriangleConn;
    int     morseC;
    int     triC;
    int     sccC;
    bool    ShowEdgeSamplesOn;
    bool    ShowTriMappingOn;
    int     selected_triangle;
    bool    ShowBackward;
    int     sampling_edge;
    int     ShowConnectionRegion;
    int     ShowBoundary;
    bool    FlipNormalOn;
    bool    IBFVOff;
    double  min_pri = 0;
    int    iter_max = 0;
    double  init_tau = 0;
    double  tau_max = 0;


    icVector3   rot_center;
    double  zoom_factor;
    float   rotmat[4][4]; //storing rotation parameters
    float   ObjXmat[16]; //Storing current modelview transformation

    MainWindow* mainWindowPtr = nullptr;
    Polyhedron* polyhedron = nullptr;

    // public member functions
    void    initializeGL2(QString file);

    void    set_up_MainWindow_ptr(MainWindow* MW_ptr);
    void    draw();
    int     InitGL();
    void    init_flags();
    void    init_tex();
    void    makePatterns();
    int     DrawGLScene(GLenum mode); // Here's Where We Do All The Drawing
    void    vis_rot_sum();
    void    ReCalTexcoord();
    void    IBFVSEffect(GLenum mode);
    void    draw_shadedObj(GLenum mode);
    void    without_antialiasing(GLenum mode);
    void    set_view(GLenum mode);
    void    set_scene(GLenum mode);
    void    set_ColorByType(unsigned char type);
    void    draw_ASolidSphere(double x, double y, double z, double r);

    void    multmatrix(const Matrix m);
    void    mat_ident(Matrix m);

    void    detect_FixedPts();

    void    display_MCG_connections(void);
    void    display_SCCs(GLenum mode);
    void    display_color_VFMag();
    void    display_separatrices(GLenum mode);
    void    display_periodicorbits(GLenum mode);
    void    display_even_streamlines(GLenum mode);
    void    display_FixedPtsIcon(GLenum mode);
    void    display_single_traj(Trajectory *traj, int sep_id);
    void    display_tau_distribution();
    void    display_diff_ECG_MCG();
    void    display_diff_MCGs();

    int    get_num_fixedPts();
    void   HitProcess(double ss, double st);
    int    SelecteObj(int hits, GLuint buffer[]);
};



#endif // OPENGLWINDOW_H
