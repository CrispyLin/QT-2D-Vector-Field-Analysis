/*
 *
    EdgeSamplePts.cpp
*/
#include "GL_LIB/glew.h"
#include "BuildGeom/Geometry.h"
#include "Others/TraceBall.h"



extern Polyhedron *object;
extern double g_zoom_factor;

extern void  HsvRgb( float hsv[3], float rgb[3] );
GLuint startList;


void
cal_rot_mat(icVector3 vec1, icVector3 vec2, float rotmat[4][4])
{
    int i, j;

    for (i=0; i<4; i++)
    {
        for (j=0; j<4; j++)
            rotmat[i][j]=0.;
    }
    rotmat[0][0] = rotmat[1][1] = rotmat[2][2] = rotmat[3][3] = 1.;

    normalize(vec1);
    normalize(vec2);

    icVector3 vs = cross(vec1, vec2);

    icVector3 v(vs);
    normalize(v);

    float ca = dot(vec1, vec2);

    icVector3 vt((1.-ca)*v);

    //rotM.M11 = vt.x * v.x + ca;
    //rotM.M22 = vt.y * v.y + ca;
    //rotM.M33 = vt.z * v.z + ca;

    rotmat[0][0] = vt.entry[0]*v.entry[0]+ca;
    rotmat[1][1] = vt.entry[1]*v.entry[1]+ca;
    rotmat[2][2] = vt.entry[2]*v.entry[2]+ca;

    vt.entry[0] *= v.entry[1];
    vt.entry[2] *= v.entry[0];
    vt.entry[1] *= v.entry[2];

    //rotM.M12 = vt.x - vs.z;
    //rotM.M13 = vt.z + vs.y;
    //rotM.M21 = vt.x + vs.z;
    //rotM.M23 = vt.y - vs.x;
    //rotM.M31 = vt.z - vs.y;
    //rotM.M32 = vt.y + vs.x;

    rotmat[0][1] = vt.entry[0] - vs.entry[2];
    rotmat[0][2] = vt.entry[2] + vs.entry[1];
    rotmat[1][0] = vt.entry[0] + vs.entry[2];
    rotmat[1][2] = vt.entry[1] - vs.entry[0];
    rotmat[2][0] = vt.entry[2] - vs.entry[1];
    rotmat[2][1] = vt.entry[1] + vs.entry[0];
}


void
get_rot_between_two_vecs(double vec1[3], double vec2[3], float rot_mat[4][4])
{
    /*  compute the cross product of these two vectors   */
    double x, y, z, w;
    double cross_vec[3];
    cross_vec[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
    cross_vec[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
    cross_vec[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];

    double dot_rs = vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];

    double sq_len1 = vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2];
    double sq_len2 = vec2[0]*vec2[0]+vec2[1]*vec2[1]+vec2[2]*vec2[2];

    x = cross_vec[0];
    y = cross_vec[1];
    z = cross_vec[2];
    w = sqrt(sq_len1)*sqrt(sq_len2) + dot_rs;

    float convert_quat[4];
    convert_quat[0] = x;
    convert_quat[1] = y;
    convert_quat[2] = z;
    convert_quat[3] = w;

    icVector3 pre_vec(vec1);
    icVector3 cur_vec(vec2);

    cal_rot_mat(pre_vec, cur_vec, rot_mat);

}
void
transpose_mat(float mat[4][4])
{
    float temp[4][4];

    int i, j;

    for (i=0; i<4; i++)
    {
        for (j=0; j<4; j++)
            temp[i][j] = mat[j][i];
    }

    for (i=0; i<4; i++)
        for (j=0; j<4; j++)
            mat[i][j]=temp[i][j];
}

void multmatrix(float m[4][4])
{
    int i,j, index = 0;

    GLfloat mat[16];

    for ( i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            mat[index++] = m[i][j];

    glMultMatrixf (mat);
}

void
EdgeSamplePt::get_sample_coords(double x[3])
{
    Edge *e = object->elist.edges[which_edge];

    x[0] = e->verts[0]->x + alpha_*(e->verts[1]->x - e->verts[0]->x);
    x[1] = e->verts[0]->y + alpha_*(e->verts[1]->y - e->verts[0]->y);
    x[2] = e->verts[0]->z + alpha_*(e->verts[1]->z - e->verts[0]->z);
}


void
EdgeSamplePt_List::display()
{
    if (samples == NULL) return;

    int i;

    float hsv[3], rgb[3];

    glDisable(GL_LIGHTING);

    glDepthFunc(GL_LEQUAL);
    glBegin(GL_LINES);
    glLineWidth(0.5);
    glColor3f(0.3, 0.3, 0.3);
    for (i=0; i<object->elist.nedges; i++)
    {
        Edge *cur_e = object->elist.edges[i];
        glVertex3f(cur_e->verts[0]->x, cur_e->verts[0]->y, cur_e->verts[0]->z);
        glVertex3f(cur_e->verts[1]->x, cur_e->verts[1]->y, cur_e->verts[1]->z);
    }
    glEnd();
    glLineWidth(1.);

    GLUquadricObj *q;
    q = gluNewQuadric();

    icVector3 z(0, 0, 1);
    float rot_mat[4][4];

    glBegin(GL_POINTS);
    for (i=0; i<num; i++)
    {
        EdgeSamplePt *cur_p = samples[i];
        Edge *e = object->elist.edges[cur_p->which_edge];
        icVector3 ave_norm = 0.5*(e->verts[0]->normal+e->verts[1]->normal);
        double x[3]={0.};
        cur_p->get_sample_coords(x);

        hsv[1]=0.8;
        hsv[2]=1;
        hsv[0]=360.-360.*cur_p->alpha_;

        HsvRgb(hsv, rgb);
        glColor3fv(rgb);

        cal_rot_mat(z, ave_norm, rot_mat);
        transpose_mat(rot_mat);

        glPushMatrix();
        glTranslatef(x[0]+0.002*ave_norm.entry[0], x[1]+0.002*ave_norm.entry[1], x[2]+0.002*ave_norm.entry[2]);
        // we need to define a rotation here
        multmatrix(rot_mat);
        gluDisk (q, 0.0, 0.008/g_zoom_factor, 20, 1);
        glPopMatrix();

        //gluDisk(
    }
    glEnd();

    //gluDeleteQuadric (q);

    //glCallList(startList+1);

    glDepthFunc(GL_LESS);
}


void
EdgeSamplePt_List::display_sel_edges(int tri, bool backward,int sampling_edge)
{
    if (samples == NULL) return;

    if (tri < 0 || tri >= object->tlist.ntris) return;

    Triangle *t = object->tlist.tris[tri];

    int i, j;

    float hsv[3], rgb[3];
    glDisable(GL_LIGHTING);

    glDepthFunc(GL_LEQUAL);
    glBegin(GL_LINES);
    glLineWidth(0.5);
    glColor3f(0.3, 0.3, 0.3);
    for (i=0; i<t->nverts; i++)
    {
        Edge *cur_e = t->edges[i];
        glVertex3f(cur_e->verts[0]->x, cur_e->verts[0]->y, cur_e->verts[0]->z);
        glVertex3f(cur_e->verts[1]->x, cur_e->verts[1]->y, cur_e->verts[1]->z);
    }
    glEnd();
    glLineWidth(1.);

    GLUquadricObj *q;
    q = gluNewQuadric();

    icVector3 z(0, 0, 1);
    float rot_mat[4][4];


    glBegin(GL_POINTS);
    for (i=0; i<num; i++)
    {
        EdgeSamplePt *cur_p = samples[i];
        Edge *e = object->elist.edges[cur_p->which_edge];
        //Edge *e = object->elist.edges[ sampling_edge];

        if (e!=t->edges[0] && e!=t->edges[1] && e!=t->edges[2])
            continue;
        if(e!=t->edges[sampling_edge])continue;

        icVector3 ave_norm = 0.5*(e->verts[0]->normal+e->verts[1]->normal);
        normalize (ave_norm);
        double x[3]={0.};
        cur_p->get_sample_coords(x);

        hsv[1]=0.8;
        hsv[2]=1;
        hsv[0]=360.-360.*cur_p->alpha_;

        HsvRgb(hsv, rgb);
        //glColor3fv(rgb);
        qInfo() << "displaying edges";
        glColor4f(rgb[0],rgb[1],rgb[2],1.);
        /*if(cur_p->backward)glColor3f(1,0,1);
        else glColor3f(1,1,0);*/

        if(backward)
        {
            if(!cur_p->backward)continue;
        }
        else
        {
            if(cur_p->backward)continue;

        }


        cal_rot_mat(z, ave_norm, rot_mat);

        glPushMatrix();
        glTranslatef(x[0]+0.002*ave_norm.entry[0], x[1]+0.002*ave_norm.entry[1], x[2]+0.002*ave_norm.entry[2]);
        // we need to define a rotation here
        multmatrix(rot_mat);
        gluDisk (q, 0.0, 0.008/g_zoom_factor, 20, 1);
        glPopMatrix();

        //gluDisk(

        /*
           Draw the mapping triangle now
        */

        Triangle *map_t = object->tlist.tris[cur_p->end_tri];
        if (cur_p->end_tri < 0 || cur_p->end_tri >= object->tlist.ntris) continue;

        float center[3];
        center[0]=(map_t->verts[0]->x+map_t->verts[1]->x+map_t->verts[2]->x)/3;
        center[1]=(map_t->verts[0]->y+map_t->verts[1]->y+map_t->verts[2]->y)/3;
        center[2]=(map_t->verts[0]->z+map_t->verts[1]->z+map_t->verts[2]->z)/3;
        float h[3];
        h[0]=map_t->verts[0]->normal.entry[0];
        h[1]=map_t->verts[0]->normal.entry[1];
        h[2]=map_t->verts[0]->normal.entry[2];


        glBegin(GL_LINES);

        glVertex3f(x[0]+h[0]*0.003,x[1]+h[1]*0.003,x[2]+h[2]*0.003);
        glVertex3f(center[0]+h[0]*0.003,center[1]+h[1]*0.003,center[2]+h[2]*0.003);
        glEnd();


        glBegin(GL_TRIANGLES);
        glVertex3f(map_t->verts[0]->x, map_t->verts[0]->y, map_t->verts[0]->z);
        glVertex3f(map_t->verts[1]->x, map_t->verts[1]->y, map_t->verts[1]->z);
        glVertex3f(map_t->verts[2]->x, map_t->verts[2]->y, map_t->verts[2]->z);
        glEnd();
    }
    glEnd();

    gluDeleteQuadric (q);

    glDepthFunc(GL_LESS);
}



/*
   create a display list for the edge sample points
*/
void
EdgeSamplePt_List::create_display_list()
{
    startList = glGenLists(2);

    glNewList(startList+1, GL_COMPILE);
    GLUquadricObj *q;
    q = gluNewQuadric();
    float hsv[3], rgb[3];

    glBegin(GL_POINTS);
    for (int i=0; i<num; i++)
    {
        EdgeSamplePt *cur_p = samples[i];
        Edge *e = object->elist.edges[cur_p->which_edge];
        icVector3 ave_norm = 0.5*(e->verts[0]->normal+e->verts[1]->normal);
        double x[3]={0.};
        cur_p->get_sample_coords(x);

        hsv[1]=0.8;
        hsv[2]=1;
        hsv[0]=360.-360.*cur_p->alpha_;

        HsvRgb(hsv, rgb);
        glColor3fv(rgb);

        /*if(cur_p->backward)glColor3f(1,0,1);
        else glColor3f(1,1,0);*/

        glPushMatrix();
        glTranslatef(x[0]+0.002*ave_norm.entry[0], x[1]+0.002*ave_norm.entry[1], x[2]+0.002*ave_norm.entry[2]);
        // we need to define a rotation here
        gluDisk (q, 0.0, 0.008/g_zoom_factor, 20, 1);
        glPopMatrix();

        //gluDisk(
    }
    glEnd();

    gluDeleteQuadric (q);
    glEndList();

}
