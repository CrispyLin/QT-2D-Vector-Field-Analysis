/*
CGraphView.h
This header file contains the declaration of the routines used for displaying the
ECG/MCG

Created and Modified by Guoning Chen
        copyright @2007
    */

#pragma once
#include "GL_LIB/glew.h"
#include <stdlib.h>
#include "ExternalDependencies/icVector.h"

class CGraphView
{
public:
    CGraphView(void);

    virtual ~CGraphView(void);

    ////attributes
public:
    ////For singularities pair cancellation
    int PairCancelOn;
    int SelectRepellerOrAttractor;

    int cancel_node1, cancel_node2;
    int *cancel_nodes;
    int num_cancel_nodes;

    int PairCounter;

    ////methods
public:
    //GLvoid ReSizeGLScene(GLsizei width, GLsizei height);
    int InitGL(GLvoid);
    void DrawGLScene(GLenum mode);
    int OnCreate();

    ////Draw the components of the graph
    void draw_nodes(GLenum mode);
    void draw_edges();
    void draw_wings(double head[2], icVector2 direct);

    void draw_solid_circle(double cx, double cy, double R);
    void draw_solid_square(double cx, double cy);
    void draw_circle(double cx, double cy, double radius,int type);
    void draw_node_for_periodicOrbit(double cx, double cy, int type);
    void draw_highlights();
    void set_color(int nodetype);

    //void DrawLimitCycleLegend(double cx, double cy);

    void display_label(int x, int y, char *string); ////display the labels
    void draw_labels();

    void HitProcessforGraph(double ss, double tt);
    void HitProcessforCancel(double ss, double tt);

    void get_highlighted_node_and_edges();
    void get_highlighted_separatrices();
    int find_connected_separatrices(int saddle, int othersing, int type);

    void init_flags();
    int ShowMCGOn;
    int ShowConleyCircle;


    ////For multiple repellers and attractors selection  11/20/05
    ////Note that these routines are different from the routines in class "CGlView"
    void add_to_repellerList(int);
    void add_to_attractorList(int);
    void clear_repell_attractLists();

    //DECLARE_MESSAGE_MAP()
    //afx_msg void OnMouseMove(UINT nFlags, CPoint point);

    void get_rotation(double p[3], icVector3 rot_axis, double rotang, double rotmatrix[16]);
    void get_rotation(double rotang, icVector3 rot_axis, double rotmatrix[16]);
    void rightMultiply16(double old_rotmat[16], double rot_mat[16]);

};
