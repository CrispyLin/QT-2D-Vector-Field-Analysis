/*
This file contains routines of class Edge, Triangle that we use to compute the
Separation and attachment points for periodic orbit extraction

Created and Modified by Guoning Chen
        copyright @2007
    */
#include "BuildGeom/Geometry.h"

extern Polyhedron *object;

/*
Implementation of the member functions of class Triangle
*/

void Triangle::cal_Jacobian()
{
    Triangle *face = this;
    icVector2 vec[3];
    icVector3 gPos;
    Vertex *v;
    double alpha[3];

    int i;
    for(i = 0; i < 3; i++)
    {
        v = face->verts[i];

        if(v->angle_deficit == 0)
        {
            gPos.entry[0] = v->x;
            gPos.entry[1] = v->y;
            gPos.entry[2] = v->z;
            alpha[i] = 1;
            alpha[(i+1)%3] = alpha[(i+2)%3] = 0;

            vec[i] = object->get_Vector_At_Point(index, gPos.entry, alpha, 0, 0);
        }

        else
        {
            /* move away from the vertex a little bit */
            alpha[i]=1-0.001;
            alpha[(i+1)%3]=0.0005;
            alpha[(i+2)%3]=0.0005;

            double cur_point[2];

            /* Get the new cur_point */
            cur_point[0] = alpha[1]*face->x1+alpha[2]*face->x2;
            cur_point[1] = alpha[2]*face->y2;

            icVector3 globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

            gPos.entry[0] = face->verts[0]->x + globalv.entry[0];
            gPos.entry[1] = face->verts[0]->y + globalv.entry[1];
            gPos.entry[2] = face->verts[0]->z + globalv.entry[2];

            vec[i] = object->get_Vector_At_Point(index, gPos.entry, alpha, 0, 0);
        }
    }

    double x[3], y[3], vx[3], vy[3];
    for (int j=0; j<face->nverts; j++) {

        if(j == 0)
        {
            x[j] = y[j] = 0;
        }

        else if(j == 1)
        {
            x[j] = face->x1; y[j] = 0;
        }
        else
        {
            x[j] = face->x2; y[j] = face->y2;
        }

        vx[j] = vec[j].entry[0];
        vy[j] = vec[j].entry[1];
    }

    //new method to get the coefficients 07/19/06
    double ta, tb, tc, td, in_ta, in_tb, in_tc, in_td, det;
    double la, lb,  ld, le;

    //solve la, lb, lc
    ta = x[1]-x[0];
    tb = y[1]-y[0];
    tc = x[2]-x[0];
    td = y[2]-y[0];
    det = ta*td-tb*tc;

    in_ta = td/det;
    in_tb = -tb/det;
    in_tc = -tc/det;
    in_td = ta/det;

    la = in_ta*(vx[1]-vx[0])+in_tb*(vx[2]-vx[0]);
    lb = in_tc*(vx[1]-vx[0])+in_td*(vx[2]-vx[0]);

    //solve ld, le, lf

    ld = in_ta*(vy[1]-vy[0])+in_tb*(vy[2]-vy[0]);
    le = in_tc*(vy[1]-vy[0])+in_td*(vy[2]-vy[0]);

    double temp1[][2]= {{la, lb}, {ld, le}};

    //set the general tensor for the triangle
    face->Jacobian.set(temp1);
}


void Triangle::cal_Jacobian2x2_neighborhood(icVector2 v[3], icVector2 vec[3], icMatrix2x2 &mat)
{
    double x[3] = {0.};
    double y[3] = {0.};
    double vx[3] = {0.};
    double vy[3] = {0.};

    int i;

    for(i = 0; i < 3; i++)
    {
        x[i] = v[i].entry[0];
        y[i] = v[i].entry[1];

        vx[i] = vec[i].entry[0];
        vy[i] = vec[i].entry[1];
    }

    //new method to get the coefficients 07/19/06
    double ta, tb, tc, td, in_ta, in_tb, in_tc, in_td, det, la, lb, le, ld;

    //solve la, lb
    ta = x[1]-x[0];
    tb = y[1]-y[0];
    tc = x[2]-x[0];
    td = y[2]-y[0];
    det = ta*td-tb*tc;

    in_ta = td/det;
    in_tb = -tb/det;
    in_tc = -tc/det;
    in_td = ta/det;

    la = in_ta*(vx[1]-vx[0])+in_tb*(vx[2]-vx[0]);
    lb = in_tc*(vx[1]-vx[0])+in_td*(vx[2]-vx[0]);

    //solve ld, le

    ld = in_ta*(vy[1]-vy[0])+in_tb*(vy[2]-vy[0]);
    le = in_tc*(vy[1]-vy[0])+in_td*(vy[2]-vy[0]);

    double temp1[][2]= {{la, lb}, {ld, le}};

    mat.set(temp1);
}


void Triangle::cal_localvec_of_vertex_in_tri(int vertID1, int vertID2, icVector2 &local_vec)
{
    int i;
    Triangle *face = this;

    int vert_angledeficit = -1;
    int v1, v2;

    for(i = 0; i < face->nverts; i++)
    {
        if(face->verts[i]->angle_deficit == 1)
        {
            vert_angledeficit = i;
        }

        if(face->verts[i]->index == vertID1)
        {
            v1 = i;
        }

        if(face->verts[i]->index == vertID2)
            v2 = i;
    }

    ////There is no angle deficit for this triangle
    if(vert_angledeficit == -1)
    {
        for(i = 0; i < face->nverts; i++)
        {
            if(face->verts[i]->index == vertID1)
            {
                local_vec = face->dir_vec[i];
            }
        }
    }

    else
        cal_localvec_with_angledeficit(v1, v2, vert_angledeficit, local_vec);
}


void Triangle::cal_localvec_with_angledeficit(int v1, int v2, int deficitvert, icVector2 &local_vec)
{
    Triangle *face = this;

    if(deficitvert == 0)         //v0 has angle deficit
    {
        if(v1 == 0)
        {
            if(v2 == 1)
                local_vec = face->dir_vec[1];
            else
                local_vec = face->dir_vec[0];
        }

        local_vec = face->dir_vec[v1+1];
    }

    else if(deficitvert == 1)    //v1 has angle deficit
    {
        if(v1 == 1)
        {
            if(v2 == 0)
                local_vec = face->dir_vec[1];
            else
                local_vec = face->dir_vec[2];
        }

        else if(v1 == 0)
            local_vec = face->dir_vec[0];

        else
            local_vec = face->dir_vec[3];
    }

    else                                      //v2 has angle deficit
    {
        if(v1 == 2)
        {
            if(v2 == 0)
                local_vec = face->dir_vec[3];
            else
                local_vec = face->dir_vec[2];
        }
        else if(v1 == 0)
            local_vec = face->dir_vec[0];

        else
            local_vec = face->dir_vec[1];
    }
}

/*
Implementation of member function of the class Polyhedron for computing Jacobians of all triangles
*/
void Polyhedron::cal_all_tri_Jacobians()
{
    int i;

    for(i = 0; i < tlist.ntris; i++)
    {
        tlist.tris[i]->cal_Jacobian();
    }
}

/*
Implementation of member functions of the class Edge for computing the separation and attachement points
*/

void Edge::cal_Jacobian_edge(int edge_tri)
{
    //get the middle point of the edge
    double middle[3], lm[2];
    Triangle *face;
    Vertex *v1, *v2;
    icVector3 VP;

    if(edge_tri == tris[0]->index)
        face = tris[0];
    else
        face = tris[1];

    v1 = verts[0];
    v2 = verts[1];

    middle[0] = (v1->x + v2->x)/2.;
    middle[1] = (v1->y + v2->y)/2.;
    middle[2] = (v1->z + v2->z)/2.;

    VP.entry[0] = middle[0] - face->verts[0]->x;
    VP.entry[1] = middle[1] - face->verts[0]->y;
    VP.entry[2] = middle[2] - face->verts[0]->z;

    lm[0] = dot(VP, face->LX);
    lm[1] = dot(VP, face->LY);

    double alpha[3] = {0.};

    //get the position of the middle inside the triangle
    object->get_2D_Barycentric_Facters(face->index, lm[0], lm[1], alpha);


    //get the small neighborhood of the middle point
    icVector2 v[3], vec[3];
    cal_neighborhood_with_vectors(face->index, alpha, lm[0], lm[1], v, vec);

    //get the approximate Jacobian for the edge
    face->cal_Jacobian2x2_neighborhood(v, vec, Jacobian);
    triangle = face->index;  //save which triangle we calculate the Jacobian inside
}


void Edge::cal_neighborhood_with_vectors(int faceid, double alpha[3],
                                         double sincx, double sincy,
                                         icVector2 Coord[3], icVector2 vec[3])
{
    icVector3 glPos;

    icVector2 vec2D;

    double t_alpha[3];
    Triangle *face = object->tlist.tris[faceid];

    Coord[0].entry[0] = sincx;
    Coord[0].entry[1] = sincy;
    object->local_To_Global(faceid, Coord[0], glPos);

    vec[0] = object->get_Vector_At_Point(faceid, glPos.entry, alpha, 0, 0);

    ////On the edge v1v2
    if(fabs(alpha[0]-0) < 1e-8)
    {
        ///calculate the other two points along the two edges incident to v0 respectively
        vec2D.entry[0] = face->x2 - sincx;
        vec2D.entry[1] = face->y2 - sincy;

        Coord[1] = Coord[0] + 0.5 * vec2D;

        ////Get the vector associated with the point
        object->get_2D_Barycentric_Facters(faceid, Coord[1].entry[0], Coord[1].entry[1], t_alpha);
        object->local_To_Global(faceid, Coord[1], glPos);
        vec[1] = object->get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);

        ////
        vec2D.entry[0] = - face->x1;
        vec2D.entry[1] = 0;

        Coord[2] = Coord[0] + 0.2 * vec2D;

        ////Get the vector associated with the point
        object->get_2D_Barycentric_Facters(faceid, Coord[2].entry[0], Coord[2].entry[1], t_alpha);
        object->local_To_Global(faceid, Coord[2], glPos);
        vec[2] = object->get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);
    }

    ////On the edge v2v0
    else if(fabs(alpha[1]-0) < 1e-8)
    {
        ///calculate the other two points along the two edges incident to v0 respectively
        vec2D.entry[0] = - sincx;
        vec2D.entry[1] = - sincy;

        Coord[1] = Coord[0] + 0.5 * vec2D;

        ////Get the vector associated with the point
        object->get_2D_Barycentric_Facters(faceid, Coord[1].entry[0], Coord[1].entry[1], t_alpha);
        object->local_To_Global(faceid, Coord[1], glPos);
        vec[1] = object->get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);

        ////
        vec2D.entry[0] = face->x1;
        vec2D.entry[1] = 0;

        Coord[2] = Coord[0] + 0.2 * vec2D;

        ////Get the vector associated with the point
        object->get_2D_Barycentric_Facters(faceid, Coord[2].entry[0], Coord[2].entry[1], t_alpha);
        object->local_To_Global(faceid, Coord[2], glPos);
        vec[2] = object->get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);
    }

    ////On the edge v0v1
    else
    {
        ///calculate the other two points along the two edges incident to v0 respectively
        vec2D.entry[0] = face->x1 - sincx;
        vec2D.entry[1] =  - sincy;

        Coord[1] = Coord[0] + 0.5 * vec2D;

        ////Get the vector associated with the point
        object->get_2D_Barycentric_Facters(faceid, Coord[1].entry[0], Coord[1].entry[1], t_alpha);
        object->local_To_Global(faceid, Coord[1], glPos);
        vec[1] = object->get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);

        ////
        vec2D.entry[0] = face->x2;
        vec2D.entry[1] = face->y2;

        Coord[2] = Coord[0] + 0.2 * vec2D;

        ////Get the vector associated with the point
        object->get_2D_Barycentric_Facters(faceid, Coord[2].entry[0], Coord[2].entry[1], t_alpha);
        object->local_To_Global(faceid, Coord[2], glPos);
        vec[2] = object->get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);
    }
}


void Edge::cal_eigenvals_edge()
{
    icMatrix2x2 mat = Jacobian;

    double evalues[2] = {0.};

    double trace = (mat.entry[0][0] + mat.entry[1][1])/2.;
    double antidiag = (mat.entry[0][1] + mat.entry[1][0])/2.;
    double la = mat.entry[0][0] - trace;
    double lb = antidiag;
    double lc = antidiag;
    double ld = mat.entry[1][1] - trace;

    mat.entry[0][0] = la;
    mat.entry[0][1] = lb;
    mat.entry[1][0] = lc;
    mat.entry[1][1] = ld;

    cal_eigenval_for_Matrix2x2(mat, evec);
}


bool Edge::cal_attachment_pt(int triangle, double attp[3])
{
    Vertex *v1, *v2;
    v1 = verts[0];
    v2 = verts[1];

    icMatrix2x2 mat;
    mat.set(Jacobian);

    ////if the R3 component of the Jacobian of the edge is smaller than some threshold, ignore that
    double r1, r2;

    icVector3 r_coef;
    r_coef.entry[0] = (mat.entry[0][0]+mat.entry[1][1])/2.;
    r_coef.entry[1] = (mat.entry[0][1]-mat.entry[1][0])/2.;

    r1 = (mat.entry[0][0]-mat.entry[1][1])/2.;
    r2 = (mat.entry[0][1]+mat.entry[1][0])/2.;

    r_coef.entry[2] = sqrt(r1*r1+r2*r2);

    ////Remove those edges that have small R1, R2 and R3 magnitude
    if(::length(r_coef) < 5)  //modified from 10 at 08/19/06
    {
        return false;
    }

    normalize(r_coef);

    ////Remove those edges that have weak R3 component
    //if(r_coef.entry[2] < R3Filter)
    if(r_coef.entry[2] < 0.5)//modified from 0.7 at 08/19/06
    {
        return false;
    }

    ////Another filter to remove some edges 07/21/06
    double cr = evec[0].entry[0]*evec[1].entry[1]\
                -evec[1].entry[0]*evec[0].entry[1];

    if(abs(cr) < 1e-8)
    {
        return false;
    }

    ////we need to use the local vector here
    icVector2 lvec1, lvec2;
    Triangle *face;
    if(triangle == this->tris[0]->index)
        face = tris[0];
    else
        face = tris[1];

    face->cal_localvec_of_vertex_in_tri(verts[0]->index, verts[1]->index, lvec1);
    face->cal_localvec_of_vertex_in_tri(verts[1]->index, verts[0]->index, lvec2);

    double re1 = evec[1].entry[0]*lvec1.entry[1] - evec[1].entry[1]*lvec1.entry[0];
    double re2 = evec[1].entry[0]*lvec2.entry[1] - evec[1].entry[1]*lvec2.entry[0];

    if(re1 * re2 >= 0)
        //if(re1 * re2 >= -1e-8)
        return false;

    if(abs(re1) <= 1e-40)
    {
        attp[0] = v1->x;
        attp[1] = v1->y;
        return true;
        //return false;
    }

    if(abs(re2) <= 1e-40)
    {
        attp[0] = v2->x;
        attp[1] = v2->y;
        return true;
        //return false;
    }

    /* The following is not an accurate estimation */
    double ratio = abs(re1)/(abs(re2)+abs(re1));

    attp[0] = (1-ratio)*v1->x + ratio*v2->x;
    attp[1] = (1-ratio)*v1->y + ratio*v2->y;
    attp[2] = (1-ratio)*v1->z + ratio*v2->z;
    return true;
}


bool Edge::cal_separation_pt(int triangle, double sep[3])
{
    Vertex *v1, *v2;
    v1 = verts[0];
    v2 = verts[1];

    icMatrix2x2 mat;
    mat.set(Jacobian);

    ////if the R3 component of the Jacobian of the edge is smaller than some threshold, ignore that
    double r1, r2;

    icVector3 r_coef;
    r_coef.entry[0] = (mat.entry[0][0]+mat.entry[1][1])/2.;
    r_coef.entry[1] = (mat.entry[0][1]-mat.entry[1][0])/2.;

    r1 = (mat.entry[0][0]-mat.entry[1][1])/2.;
    r2 = (mat.entry[0][1]+mat.entry[1][0])/2.;

    r_coef.entry[2] = sqrt(r1*r1+r2*r2);

    ////Remove those edges that have small R1, R2 and R3 magnitude
    if(::length(r_coef) < 5)//modified from 10 at 08/19/06
    {
        return false;
    }

    normalize(r_coef);

    ////Remove those edges that have weak R3 component
    //if(r_coef.entry[2] < R3Filter)
    if(r_coef.entry[2] < 0.5)//modified from 0.7 at 08/19/06
    {
        return false;
    }

    ////Calculate the eigen vector of the symmetric matrix of the Jacobian of the edge
    //double trace = (mat.entry[0][0] + mat.entry[1][1])/2.;
    //double antidiag = (mat.entry[0][1] + mat.entry[1][0])/2.;
    //mat.entry[0][0] = mat.entry[0][0] - trace;
    //mat.entry[1][1] = mat.entry[1][1] - trace;
    //mat.entry[0][1] = antidiag;
    //mat.entry[1][0] = antidiag;
    //CalEigenForMatrix2x2(mat, cur_e->evec);

    ////Another filter to remove some edges 07/21/06
    double cr = evec[0].entry[0]*evec[1].entry[1]\
                -evec[1].entry[0]*evec[0].entry[1];

    if(abs(cr) < 1e-8)
    {
        return false;
    }


    ////we need to use the local vector here
    icVector2 lvec1, lvec2;

    Triangle *face;
    if(triangle == this->tris[0]->index)
        face = tris[0];
    else
        face = tris[1];

    face->cal_localvec_of_vertex_in_tri(verts[0]->index, verts[1]->index, lvec1);
    face->cal_localvec_of_vertex_in_tri(verts[1]->index, verts[0]->index, lvec2);

    double re1 = evec[0].entry[0]*lvec1.entry[1] - evec[0].entry[1]*lvec1.entry[0];
    double re2 = evec[0].entry[0]*lvec2.entry[1] - evec[0].entry[1]*lvec2.entry[0];

    if(re1*re2 >= 0)
    {
        return false;
    }

    if(abs(re1) <= 1e-40)
    {
        sep[0] = v1->x;
        sep[1] = v1->y;
        return false;
    }

    if(abs(re2) <= 1e-40)
    {
        sep[0] = v2->x;
        sep[1] = v2->y;
        return false;
    }

    /* The following is not an accurate estimation */
    double ratio = abs(re1)/(abs(re2)+abs(re1));

    sep[0] = (1-ratio)*v1->x + ratio*v2->x;
    sep[1] = (1-ratio)*v1->y + ratio*v2->y;
    sep[2] = (1-ratio)*v1->z + ratio*v2->z;
    return true;
}


/*
Implementation of member function of class SCComponent for computing the separation and attachment points
inside a SCC
*/

void SCComponent::cal_sep_attp_pts()
{
    int j, k;
    Triangle *face;
    Edge *cur_e;
    double p[3] = {0.};

    //num_seppts = 0;
    //num_attpts = 0;

    /*initialize the edges*/
    for(j = 0; j < nnodes; j++)
    {
        ////First, you need to mark all the boundary edges, and don't calculate the special points on them!
        ////Now, just let me ignore that step and experiment a little bit further

        face = object->tlist.tris[nodes[j]];

        for(k = 0; k < 3; k++)
        {
            face->edges[k]->visited = false;
        }
    }


    for(j = 0; j < nnodes; j++)
    {
        ////First, you need to mark all the boundary edges, and don't calculate the special points on them!
        ////Now, just let me ignore that step and experiment a little bit further

        face = object->tlist.tris[nodes[j]];

        for(k = 0; k < 3; k++)
        {
            cur_e = face->edges[k];

            if(cur_e->visited)
                continue;

            cur_e->cal_Jacobian_edge(face->index);
            cur_e->cal_eigenvals_edge();

            cur_e->valid = 0;  //the default setting is invalide

            ////calculate the separation point
            if(cur_e->cal_separation_pt(cur_e->triangle, p))
            {
                cur_e->sep = new Point3Dd(p[0],p[1],p[2]);
                cur_e->find_sep = true;
                cur_e->valid = true;
                cur_e->sep_visit = 0;

                nseppts ++;
                //Num_SpecialPts++;  ////08/11/06
            }

            ////calculate the attachment point
            if(cur_e->cal_attachment_pt(cur_e->triangle, p))
            {
                cur_e->attp = new Point3Dd(p[0],p[1],p[2]);
                cur_e->find_attp = true;
                cur_e->valid = true;
                cur_e->att_visit = 0;

                nattpts ++;
                //Num_SpecialPts++;  ////08/11/06
            }

            cur_e->visited = true;
        }
    }
}
