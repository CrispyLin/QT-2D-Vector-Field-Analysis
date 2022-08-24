#include "VectorFieldWindow.h"

extern Polyhedron* object;

void    display_sel_tri(int tri)
{
    Triangle *t = object->tlist.tris[tri];
    int i;
    glDepthFunc(GL_LEQUAL);
    glBegin(GL_LINE_LOOP);
    glColor3f(1, 1, 1);
    for (i=0; i<t->nverts; i++)
    {
        glVertex3f(t->verts[i]->x, t->verts[i]->y, t->verts[i]->z);
    }
    glEnd();
    glDepthFunc(GL_LESS);
}


void    matrix_ident( float m[4][4])
{
    int i;

    for (i = 0; i <= 3; i++) {
        m[i][0] = 0.0;
        m[i][1] = 0.0;
        m[i][2] = 0.0;
        m[i][3] = 0.0;
        m[i][i] = 1.0;
    }
}


void ScreenToSecondWin(
        int px, int py,
        int screen_leftx, int screen_bottomy,
        int win_screen_sizex, int win_screen_sizey,
        double world_leftx, double world_bottomy,
        double win_world_sizex, double win_world_sizey,
        double &s, double &t
){
    double ratiox = (double)(px-screen_leftx)/ (double)win_screen_sizex;
    double ratioy = (double)(screen_bottomy - py)/(double)win_screen_sizey;

    s = (double) world_leftx + ratiox * win_world_sizex;
    t = (double) world_bottomy + ratioy * win_world_sizey;
}

void display_Obj(GLenum mode)
{
    int i, j;
    Triangle *face;
    int *verts;

    for (i=0; i<object->tlist.ntris; i++) {
        if (mode == GL_SELECT)
            glLoadName(i+1);

        face = object->tlist.tris[i];

        glBegin(GL_POLYGON);
        for (j=0; j<face->nverts; j++) {
            glVertex3f(face->verts[j]->x, face->verts[j]->y, face->verts[j]->z);
        }
        glEnd();
    }

}
