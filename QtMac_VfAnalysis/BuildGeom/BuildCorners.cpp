/*
This file contains routines of class Polyhedron that we use to build the corner list for the
being loaded geometrical domain

Created and Modified by Guoning Chen
        copyright @2007
    */

#include "VField.h"

extern Polyhedron *object;

void Polyhedron::build_Corners()
{
    int i, j, m;

    ////////////////Define variables for corner information building
    Vertex **vert = vlist.verts;

    clist.corners = new Corner *[tlist.ntris * 3]; //allocate the memory for corners
    clist.ncorners = tlist.ntris * 3;

    for( i = 0; i < tlist.ntris * 3; i++)
    {
        clist.corners[i] = (Corner *) malloc( sizeof(Corner));
        clist.corners[i]->index = 0;
        clist.corners[i]->Edge_count = 0;
    }

    for( i = 0; i < vlist.nverts; i++)
        vlist.verts[i]->ncorners = 0;

    if(elist.nedges == 0)
        build_Edges();

    for( i = 0; i < tlist.ntris; i++)
    {
        for( j = 0; j < tlist.tris[i]->nverts; j++)
        {
            //////////////////////////////////////////////////////////////////////////////
            ////////Add some information to corner structure 1/21/05
            ///Get the current id of triangle to the corner
            clist.corners[i*3 + j]->index = i*3 + j;
            clist.corners[i*3 + j]->t = i;
            clist.corners[i*3 + j]->v = tlist.tris[i]->verts[j]->index; //Get the ID of current vertex


            /////At the same time, we need to add this corner to current vertex
            add_Corner_to_Vertex(clist.corners[i*3 + j], tlist.tris[i]->verts[j]);

            //Get its previous corner ID
            //if( (j - 1) < 0)
            //	clist.corners[i*3 + j]->p = i*3 + 2;
            //else
            clist.corners[i*3 + j]->p = i*3 + (3+(j-1))%3;

            //Get its next corner ID
            clist.corners[i*3 + j]->n = i*3 + (j+1)%3;


            ///Get the angle of the corner
            //clist.corners[i*3 + j]->angle = get_Angle(tlist.tris[i]->verts[j],
            //	                                    tlist.tris[i]->verts[(3+(j-1))%3],
            //			                            tlist.tris[i]->verts[(j+1)%3]);


            ///Get the edges constitute the corner
            for( m = 0; m < 3; m++)
            {

                if( tlist.tris[i]->edges[m]->verts[0] == tlist.tris[i]->verts[j]
                    ||tlist.tris[i]->edges[m]->verts[1] == tlist.tris[i]->verts[j])
                {
                    ////If one of the end points of the edge is the current vertex
                    ////It means that the edge is incident to the vertex
                    ///Add this edge to the corner associate with the vertex
                    if( clist.corners[i*3 + j]->Edge_count == 0)
                    {
                        clist.corners[i*3 + j]->edge[0] = tlist.tris[i]->edges[m];
                        clist.corners[i*3 + j]->Edge_count++;
                    }
                    else
                    {
                        clist.corners[i*3 + j]->edge[1] = tlist.tris[i]->edges[m];
                        clist.corners[i*3 + j]->Edge_count++;
                    }
                }
                else //it means that the edge is not incident to the vertex,
                    //so it must be the opposite edge of that vertex
                    clist.corners[i*3 + j]->e = tlist.tris[i]->edges[m];
            }

            /////////////////////////////////////////////////////////////////////////
            ////Initialize the  opposite corner o to -1
            clist.corners[i*3 + j]->o = -1;
            clist.corners[i*3 + j]->ot = -1;
        }
    }

    find_Opposite();     //finding the opposite corner of each corner

    ////Sort the edges for each corner, make edge[0] store the edge that is the opposite edge of corner c.n
    ////edge[1]  store the edge that is the opposite edge of corner c.p;

    Corner *c;
    for( i = 0 ; i < clist.ncorners; i++)
    {
        c = clist.corners[clist.corners[i]->n];
        clist.corners[i]->edge[0] = c->e;

        c = clist.corners[clist.corners[i]->p];
        clist.corners[i]->edge[1] = c->e;
    }
}



/************************************************************************
This routine implements add an corner into the vertex
************************************************************************/
void Polyhedron::add_Corner_to_Vertex(Corner *c, Vertex *v)
{
    Corner **temp = v->corners;
    v->corners = new Corner *[v->ncorners + 1];
    if(v->ncorners > 0)
    {
        for(int i = 0; i < v->ncorners; i++)
            v->corners[i] = temp[i];
        delete []temp;
    }

    v->corners[v->ncorners] = c;

    v->ncorners++;
}


/*********************************************************************************************
This routine is designed to find the opposite corner
It must be called after all the corners have been traversed and basic information (max, min) have been set up
There may be some corner withour opposite corner, because their opposite edge is the boundary of the object
But for a closed 3D object, all the corners should have an opposite corner
***************************************************************************/
void Polyhedron::find_Opposite()
{
    int i, j, k;

    ////Get the opposite corner without sorting
    Vertex **vert = vlist.verts;
    int face_index;

    for( i = 0; i < tlist.ntris; i++)
    {
        for( j = 0; j < tlist.tris[i]->nverts; j++)
        {
            if(clist.corners[i*3 + j]->o < 0)
            {
                /////Getting the opposite face ID through 'c.e' information
                //if(face_index < 0) ////the oposite triangle does not exist!
                //	continue;

                if(clist.corners[i*3 + j]->e->tris[0] != NULL &&
                    clist.corners[i*3 + j]->e->tris[0]->index != i)
                    face_index = clist.corners[i*3 + j]->e->tris[0]->index;
                else if(clist.corners[i*3 + j]->e->tris[1] != NULL &&
                         clist.corners[i*3 + j]->e->tris[1]->index != i)
                    face_index = clist.corners[i*3 + j]->e->tris[1]->index;
                else
                    continue;


                for( k = 0; k < tlist.tris[face_index]->nverts; k++)
                {
                    ///Search all the corners belong to the opposite face
                    ///Finding the corner that c'.e has
                    if( clist.corners[face_index*3 + k]->e == clist.corners[i*3 + j]->e)
                    { // This corner is the opposite corner we want to find
                        //Because the opposite corners are pairs
                        //So we found one pair of opposite corner here
                        clist.corners[i*3 + j]->o = face_index*3 + k;
                        clist.corners[i*3 + j]->ot = face_index; //get the opposite triangle

                        clist.corners[face_index*3 + k]->o = i*3 + j;
                        clist.corners[face_index*3 + k]->ot = i;
                    }
                }
            }
        }
    }
}


double Polyhedron::get_Angle(Corner *c)
{
    icVector3 v1, v2;

    Vertex **verts = vlist.verts;

    v1.entry[0] = verts[clist.corners[c->p]->v]->x - verts[c->v]->x;
    v1.entry[1] = verts[clist.corners[c->p]->v]->y - verts[c->v]->y;
    v1.entry[2] = verts[clist.corners[c->p]->v]->z - verts[c->v]->z;

    v2.entry[0] = verts[clist.corners[c->n]->v]->x - verts[c->v]->x;
    v2.entry[1] = verts[clist.corners[c->n]->v]->y - verts[c->v]->y;
    v2.entry[2] = verts[clist.corners[c->n]->v]->z - verts[c->v]->z;

    ///Getting the angle of this two vectors
    normalize(v1);
    normalize(v2);

    double ang = acos(dot(v1, v2));

    return ang;
}


/*
An implementation of quick sort algorithm
*/
//void QuickSort(int array[], int left, int right)
//{
//        if(left>=right)
//        return;
//
//        int i = left;
//        int j = right+1;
//        int t = array[left];//left most element is the pivot element
//
//        while(1)
//        {
//                //elements to the left of t should be less than t and elements to the right of t should be greater than t
//                do
//                {
//                        i++;
//                }while(i<=right && array[i]<t);
//                do
//                {
//                        j--;
//                }while(array[j]>t);
//                if(i>j)//if the indices cross
//                        break;
//                int temp = array[i];
//                array[i] = array[j];
//                array[j] = temp;
//        }
//        //swap the pivot element and array[j]
//        int temp = array[left];
//        array[left] = array[j];
//        array[j] = temp;
//
//        QuickSort(array,left,j-1);
//        QuickSort(array,j+1,right);
//}


void QuickSortforCorner(int corner[], int min_v[], int max_v[], int left, int right)
{
    if(left>=right)
        return;

    int i = left;
    int j = right+1;
    int t = max_v[left];//left most element is the pivot element

    while(1)
    {
        //elements to the left of t should be less than t and elements to the right of t should be greater than t
        do
        {
            i++;
        }while(i<=right && max_v[i]<t);
        do
        {
            j--;
        }while(max_v[j]>t);
        if(i>j)//if the indices cross
            break;
        int temp = max_v[i];
        max_v[i] = max_v[j];
        max_v[j] = temp;

        //swap the corner and min_v array as well
        temp = corner[i];
        corner[i] = corner[j];
        corner[j] = temp;

        temp = min_v[i];
        min_v[i] = min_v[j];
        min_v[j] = temp;
    }

    //swap the pivot element and array[j]
    int temp = max_v[left];
    max_v[left] = max_v[j];
    max_v[j] = temp;

    //swap the corner and min_v array as well
    temp = corner[left];
    corner[left] = corner[j];
    corner[j] = temp;

    temp = min_v[left];
    min_v[left] = min_v[j];
    min_v[j] = temp;

    QuickSortforCorner(corner,min_v, max_v, left,j-1);
    QuickSortforCorner(corner,min_v, max_v, j+1,right);
}



/**********************************************************
Select the first corner that most aligh with tangent plane
**********************************************************/
int Polyhedron::get_First_Corner(Vertex *p)
{
    int i, firstID = 0;
    double dot1, maxdot = -1;
    Triangle *face;
    for(i = 0; i < p->ncorners; i++)
    {
        face = tlist.tris[p->corners[i]->t];
        dot1 = dot(p->normal, face->normal);

        if(dot1 > maxdot)
        {
            firstID = i;
            maxdot = dot1;
        }
    }

    return firstID;
}


void Polyhedron::sort_Corner_on_Vertex()
{
    Corner **Newlist;
    int i, j;
    Vertex *vert;
    Corner *c;

    for(i = 0; i < vlist.nverts; i++)
    {
        vert = vlist.verts[i];

        /* consider the boundary of the planar domains */
        //if(vert->ncorners <= 3)
        //	continue;

        /*if one of its edge is boundary edge, this vertex is a boundary vertex, continue*/
        bool flag = false;
        for(j=0; j<vert->nedges; j++)
        {
            //if(vert->edges[j]->tris[0]->index == -1
            //	||vert->edges[j]->tris[1]->index == -1)
            if(vert->edges[j]->tris[0] == NULL
                ||vert->edges[j]->tris[1] == NULL)
            {
                flag = true;
                break;
            }
        }

        Newlist = new Corner *[vert->ncorners];

        //// Modified by Guoning at 09/28/2010 to handle the boundary vertices sorting
        if(flag)
        {
            flag = false;
            //continue;

            // we still need to sort them !!! However, we need to start with a corner having a boundary edge

            for (j=0; j<vert->ncorners; j++)
            {
                Corner *c=vert->corners[j];

                if (c->edge[0]->tris[1]==NULL || c->edge[1]->tris[1]==NULL)
                {
                    Newlist[0] = vert->corners[j];

                    if (clist.corners[Newlist[0]->p]->o>=0)
                        break;
                }
            }
        }

        else
            Newlist[0] = vert->corners[get_First_Corner(vert)];

        bool keep_original = false;
        for(j = 1; j < vert->ncorners; j++)
        {
            //Need to deal with boundary for open manifold
            //if(Newlist[j-1]->index < 0) //07/18/06
            if(Newlist[j-1]== NULL) //07/18/06
                continue;

            // The following is not correct if the vertex is non manifold (i.e., the touch point of two disks)
            if(clist.corners[Newlist[j-1]->p]->o < 0) //07/18/06
            {
                keep_original = true;
                Newlist[j] = NULL;
                continue;
            }

            c = clist.corners[clist.corners[Newlist[j-1]->p]->o];
            Newlist[j] = clist.corners[c->p];
        }

        if (keep_original)
            continue;

        //store the sorted corners
        for(j = 0; j < vert->ncorners; j ++)
        {
            vert->corners[j] = Newlist[j];
        }

        delete[] Newlist;
    }
}


int Polyhedron::get_orientation(Vertex *p, double cur_ang, Edge *e1)
{
    double a, b;
    icVector3 edge1, normal;
    Vertex *adjv;
    double ang1, ang_diff;

    normal = p->normal;

    if( e1->verts[0] != p)
        adjv = e1->verts[0];
    else
        adjv = e1->verts[1];

    edge1.entry[0] = adjv->x - p->x;
    edge1.entry[1] = adjv->y - p->y;
    edge1.entry[2] = adjv->z - p->z;


    ////project to the tangent plane
    edge1 = edge1 - dot(edge1, normal)*normal;
    normalize(edge1);

    ////calculate the angle of the projected vector
    a = dot(edge1, p->T);
    b = dot(edge1, p->B);

    ang1 = atan2(b, a);
    if(ang1 < 0)
        ang1 += 2*M_PI;

    ang_diff = ang1 - cur_ang;

    if(ang_diff < -M_PI)
        ang_diff += 2 * M_PI;
    if(ang_diff > M_PI)
        ang_diff -= 2 * M_PI;

    if(ang_diff >= 0)
        return 1;
    else
        return -1;
}


/*
    Allocate the angles for the corners associated with each vertex
*/
void Polyhedron::alloc_Corner_Angle()
{
    int i, j;
    Vertex *vert, *adjv;
    icVector3 normals, edgevec, projv;
    Corner *c;
    double sum_ang, cur_ang, r;
    double a, b;
    int orient;

    orient = 0;
    r = sum_ang = cur_ang = 0;
    a = b = 0;

    double *Anglist;



    for(i = 0; i < vlist.nverts; i++)
    {
        //_cprintf ("Processing the corners at vertex %d\n", i);

        //if (i==3729)
        //{
        //	int stop = 1;
        //}

        vert = vlist.verts[i];
        sum_ang = 0;
        Anglist = new double[vert->ncorners];

        ////calculate the total angle around current vertex
        bool defect_corner = false;
        for(j = 0; j < vert->ncorners; j++)
        {
            if(vert->corners[j] == NULL) //07/18/06
            {
                defect_corner = true;
                //_cprintf ("Defect corner at Vertex %d\n", i);
                continue;
            }

            c = vert->corners[j];
            //sum_ang += c->angle;
            sum_ang += get_Angle(c);
        }

        if (defect_corner) continue;

        r = 2 * M_PI / sum_ang;

        //////////////////////////////////////////////////////////////////////


        ////choose first edge of first corner to do the projection to get the initial angle
        c = vert->corners[0];
        if( c->edge[0]->verts[0]->index != i)
            adjv = c->edge[0]->verts[0];
        else
            adjv = c->edge[0]->verts[1];

        edgevec.entry[0] = adjv->x - vert->x;
        edgevec.entry[1] = adjv->y - vert->y;
        edgevec.entry[2] = adjv->z - vert->z;

        normalize(edgevec);

        normals = vert->normal;

        ////project to the tangent plane
        double dot1 = dot(edgevec, normals);
        projv = edgevec - dot1*normals;
        normalize(projv);

        ////calculate the angle of the projected vector
        a = dot(projv, vert->T);
        b = dot(projv, vert->B);

        ////adjust the angle to make sure that it is between 0 and 2Pi
        Anglist[0] = atan2(b,a);
        if(Anglist[0] < 0)
            Anglist[0] += 2*M_PI;

        cur_ang = Anglist[0];  ////set it as current angle position

        //////Store it as the beginning angle for the first corner in the corners' list
        vert->corners[0]->r = r;
        vert->corners[0]->BeginAng = Anglist[0];
        vert->corners[vert->ncorners-1]->EndAng = cur_ang;


        ////we may need to first decide the orientation of the angle
        orient = get_orientation(vert, cur_ang, c->edge[1]);

        /*   modified at 07/01/2009   */
        if(orient>0)
            c->orient = true;
        else
            c->orient = false;

        ////calculate the angles for other corners/triangles
        for(j = 1 ; j < vert->ncorners; j++)
        {
            c = vert->corners[j-1];

            if(orient > 0)
            {
                Anglist[j] = r * get_Angle(c) + cur_ang;

                if(Anglist[j] >= 2*M_PI)
                    Anglist[j] -= 2*M_PI;
            }
            else{
                Anglist[j] = cur_ang - r * get_Angle(c);

                if(Anglist[j] < 0)  ////modified on 04/23/05
                    Anglist[j] += 2*M_PI;
            }

            cur_ang = Anglist[j];

            //////Store it as the beginning angle for the first corner in the corners' list
            vert->corners[j-1]->EndAng = cur_ang;
            vert->corners[j]->BeginAng = cur_ang;
            vert->corners[j]->r = r;
            //vert->corners[j]]->orient = orient;
        }
    }
}



void Polyhedron::find_all_corners_around_vertex(int vertid)
{
    Vertex *temp_v = vlist.verts[vertid];
    int i;
    Triangle *temp_t;
    Corner *temp_c, *start_c, *current_c;

    //temp_v->incident_corners = new CornerList(10);
    if (temp_v->corners!=NULL)
        free (temp_v->corners);
    temp_v->corners = (Corner **) malloc( sizeof(Corner*)*temp_v->nedges);  // assume that edge list has been added
    temp_v->ncorners = 0;

    temp_t = temp_v->edges[0/*temp_v->tangent_x_edgeid*/]->tris[0];
    for (i=0; i<3; i++) {
        if (clist.corners[3*temp_t->index+i]->v == temp_v->index) {
            temp_c = clist.corners[3*temp_t->index+i];
            break;
        }
    }
    if (orient == 0) {
        if (clist.corners[temp_c->p]->e == temp_v->edges[0/*temp_v->tangent_x_edgeid*/]) {
            temp_t = temp_v->edges[0/*temp_v->tangent_x_edgeid*/]->tris[1];
            for (i=0; i<3; i++) {
                if (clist.corners[3*temp_t->index+i]->v == temp_v->index) {
                    temp_c = clist.corners[3*temp_t->index+i];
                    break;
                }
            }
        }
    }
    current_c = start_c = temp_c;

    //if (orient == 1) {
    //	if (temp_v->boundary_vertex == 0) {
    //		temp_v->Corners[temp_v->Num_corners]=current_c->index;
    //		temp_v->Num_corners ++;
    //		while (Object.clist[Object.clist[Object.clist[current_c->n]->o]->n] != start_c) {
    //			current_c = Object.clist[Object.clist[Object.clist[current_c->n]->o]->n];
    //			temp_v->Corners[temp_v->Num_corners]=current_c->index;
    //			temp_v->Num_corners ++;
    //		}
    //	} else {
    //		while (Object.clist[current_c->p]->o >=0 )
    //			current_c = Object.clist[Object.clist[Object.clist[current_c->p]->o]->p];
    //		temp_v->Corners[temp_v->Num_corners]=current_c->index;
    //		temp_v->Num_corners ++;
    //		while (Object.clist[current_c->n]->o >= 0) {
    //			current_c = Object.clist[Object.clist[Object.clist[current_c->n]->o]->n];
    //			temp_v->Corners[temp_v->Num_corners]=current_c->index;
    //			temp_v->Num_corners ++;
    //		}
    //	}
    //} else {  // plantoic and bunny
    //	if (temp_v->boundary_vertex == 0) {
    //		temp_v->Corners[temp_v->Num_corners]=current_c->index;
    //		temp_v->Num_corners ++;
    //		while (Object.clist[Object.clist[Object.clist[current_c->p]->o]->p] != start_c) {
    //			current_c = Object.clist[Object.clist[Object.clist[current_c->p]->o]->p];
    //			temp_v->Corners[temp_v->Num_corners]=current_c->index;
    //			temp_v->Num_corners ++;
    //		}
    //	} else {
    //		while (Object.clist[current_c->n]->o >= 0)
    //			current_c = Object.clist[Object.clist[Object.clist[current_c->n]->o]->n];
    //		temp_v->Corners[temp_v->Num_corners]=current_c->index;
    //		temp_v->Num_corners ++;
    //		while (Object.clist[current_c->p]->o >= 0) {
    //			current_c = Object.clist[Object.clist[Object.clist[current_c->p]->o]->p];
    //			temp_v->Corners[temp_v->Num_corners]=current_c->index;
    //			temp_v->Num_corners ++;
    //		}
    //	}
    //}
}
