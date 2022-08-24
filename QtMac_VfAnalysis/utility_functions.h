#ifndef UTILITY_FUNCTIONS_H
#define UTILITY_FUNCTIONS_H
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

void    display_sel_tri(int tri);
void    matrix_ident( float m[4][4]);

void ScreenToSecondWin(
        int px, int py,
        int screen_leftx, int screen_bottomy,
        int win_screen_sizex, int win_screen_sizey,
        double world_leftx, double world_bottomy,
        double win_world_sizex, double win_world_sizey,
        double &s, double &t
);

void display_Obj(GLenum mode);

#endif // UTILITY_FUNCTIONS_H
