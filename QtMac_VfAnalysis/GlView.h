#pragma once
#ifndef GLVIEW_H
#define GLVIEW_H

#include "GL_LIB/glew.h"
#include <stdlib.h>

#include "ExternalDependencies/icVector.h"
#include "Analysis/MorseDecomp.h"

typedef float Matrix[4][4];

// CGlView

class CGlView
{
public:
    CGlView();
    virtual ~CGlView();

    // Attributes
public:
    int MoveOrStop ;

    icVector3 rot_center;
    double zoom_factor;
    float rotmat[4][4]; //storing rotation parameters
    float ObjXmat[16]; //Storing current modelview transformation

    // Operations for OpenGl window setting and displaying
public:
    int InitGL(GLvoid);	           //initialize the OpenGl envrionment
    int DrawGLScene(/*GLvoid*/GLenum);       //The same function as display routine to show the visual effect in the opengl window
    int OnCreate();

    /*--------------------------------------------------------------------------*/
    //Operations for ibfv tool visual effect displaying
public:
    void getDP(double x, double y, double *px, double *py); ////just for testing 06/23/05
    void makePatterns(void);  ////To initialize the texture pattern for flow effect visualization

    void IBFVSEffect(GLenum mode);
    void IBFVSEffect_noshiny(GLenum mode);   // without shiny effect
    void getInitTex(void);
    void set_view(GLenum mode);
    void set_scene(GLenum mode);

    void set_light();

    void multmatrix(const Matrix m); ////This routine should be put in the matrix operation file not here!!
    void mat_ident( Matrix m);       ////Set back to identity matrix

    /*get the rotation matrix using OpenGL*/
    void get_rotation(double p[3], icVector3 rot_axis, double rotang, double rotmatrix[16]);
    void get_rotation(double rotang, icVector3 rot_axis, double rotmatrix[16]);

    void ReCalTexcoord(void);      //recalculate the texture coordinates after user changing the view point

    //Display the features of the vector field

    /**Display the fixed points**/
    int ShowFixedPtOn;
    void display_FixedPtsIcon(GLenum mode);
    void set_ColorByType(unsigned char type);
    void draw_ASolidSphere(double, double, double, double);

    /**Display separatrices**/
    int ShowSeparatricesOn;
    void display_separatrices(GLenum mode);
    void display_single_traj(Trajectory *traj, int sep_id);

    /**Display the SCCs**/
    int ShowSCCsOn;
    void display_SCCs(GLenum mode);
    void get_Color(int num, float rgb[3]);
    void get_Color_frac(int num, int frac, float rgb[3]);

    /**Display the calculated attachement and separation points **/
    void display_sepandatt_pts(GLenum mode);

    /**Display the extracted periodic orbits**/
    int ShowPeriodicOrbitsOn;
    void display_periodicorbits(GLenum mode);

    /**Display the evenly placed streamlines**/
    int EvenStreamlinePlacement;
    void display_even_streamlines(GLenum mode);
    void draw_shadedObj(GLenum mode);

    /**Display the color map of the vector field magnitude**/
    int ShowColorVFMagOn;
    void display_color_VFMag();

    void init_flags();

    //Display Connection
    int ShowMorseConn;
    int ShowTriangleConn;
    int morseC;
    int triC;
    int sccC;
    void display_morse_connection(GLenum mode);
    void display_triangle_connection(GLenum mode);

    bool ShowEdgeSamplesOn;
    bool ShowTriMappingOn;

    /*  02/10/2010  */
    void HitProcess(double ss, double st);
    int SelecteObj(int hits, GLuint buffer[]);
    int selected_triangle;
    int ShowBackward;
    int sampling_edge;

    //display connection regions
    void display_MCG_connections(void);
    int ShowConnectionRegion;

    //display boundary
    int ShowBoundary;

    // display tau value distribution
    void display_tau_distribution();
    void display_diff_ECG_MCG();
    void display_diff_MCGs();

    void without_antialiasing(GLenum mode);
    void with_antialiasing(GLenum mode);

    // visualize the total rotation
    void vis_rot_sum();

    bool FlipNormalOn;

    bool IBFVOff;
};

#endif
