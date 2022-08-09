/*
Created and modified by Guoning Chen
copyright @2007
*/
#pragma once
#ifndef __PREDEFINED_H__
#define __PREDEFINED_H__

#define M_PI 3.14159265358979324      ////Move to the datastructure.h file 06/23/05

#define DistanceThreshold 0.00001
///////////////////////////////////////////////////
enum Elementtype{
    NOTHING,
    SOURCE,
    SINK,
    SADDLE,
    CWCENTER,
    CCWCENTER,
    AFOCUS,
    RFOCUS,
    REGULAR
};

/*current color scheme  08/25/05
source: pure green (0, 1, 0)
repeller: light green (0, 1, 0.7)
sink:   pure red (1, 0, 0)
attractor:  orange (1, 0.5, 0)
saddle: pure blue (0, 0, 1)
center: light red (1, 0, 1)
*/

///////////////////////////////////////////////////
enum which_point{
    NON,
    LOWLEFT,
    UPPERLEFT,
    UPPERRIGHT,
    LOWRIGHT,
    LEFT,
    UPPER,
    RIGHT,
    BUTTOM,
    UPPERROTATE,
    ARROWBASE,
    ARROWHEAD
};


////////////////////////////////////////////////////////
////Load name for object selection
#define NAMEOFSINGELEM      1           ////name of singular element for mouse selection
#define NAMEOFREGELEM       2001        ////name of regular element for mouse selection
#define NAMEOFSINGCONTROL   3000        ////name of singular element control points for mouse selection
#define NAMEOFREGCONTROL    4000        ////name of regular element control points for mouse selection
#define NAMEOFSINGULARITY   5000        ////name of singularities for mouse selection of topology editing
#define NAMEOFLIMITCYCLES   7000        ////name of limit cycles for limit cycle editing
#define NAMEOFSHAPECONTROL  8000        ////name of shape control point for limit cycle shape controlling
#define NAMEOFTRIANGLE      10001       ////name of triangles for mouse selection



////some constants
const int NUMTRACINGTRIANGLE = 1000;
const double RegularStrength = 0.5;  //the scaler for the strength of regular element
const double ARROWSCALE = 0.07;
const double EDITBOXSIZE = 0.04;
const double INOROUTJUDGE = 0.015;     //for single separatrices
/* modify the step on 01/23/07 */
const double SEPARATRIXSTEP = 0.09;     //for separatrices, maybe a better setting for separatrices drawing
const double INTERSECTIONERROR = 1e-5; //the error threshold between beginning point and next intersection point
                                       //for closed streamline tracing
const double SCALESTEP = 15.;   //The scaler step control for scale
const double ROTATESTEP = 10.;   //The scaler step control for rotation

const int InitSaddleRegionLength = 100;  ////tracing n's triangles for each direction


#define	NPN 64
#define NMESH  100
#define DM  ((double) (1.0/(NMESH-1.0)))
#define NPIX  /*512*/800 // number of pixels

#endif /* __PREDEFINED_H__ */

