/*
This file contains routines of class Polyhedron that we use to build the edge list for the
being loaded geometrical domain

Created and Modified by Guoning Chen
        */

#include "VField.h"

extern Polyhedron *object;


void Polyhedron::add_Edge_to_Vertex(Edge *e, Vertex *v)
{

    Edge **temp = v->edges;
    v->edges = new Edge *[v->nedges + 1];

    if(v->nedges > 0)
    {
        for(int i = 0; i < v->nedges; i++)
            v->edges[i] = temp[i];
        delete [] temp;
        //free(temp);
    }

    v->edges[v->nedges] = e;

    v->nedges++;
}



void Polyhedron::calc_Edge_Length()
{
    int i;
    icVector3 v1, v2;
    shortest_edgelength = 1e49;

    for (i=0; i<elist.nedges; i++) {
        v1.set(elist.edges[i]->verts[0]->x, elist.edges[i]->verts[0]->y, elist.edges[i]->verts[0]->z);
        v2.set(elist.edges[i]->verts[1]->x, elist.edges[i]->verts[1]->y, elist.edges[i]->verts[1]->z);
        //v2.set(elist.edges[i]->verts[1]->x, elist.edges[i]->verts[1]->y, elist.edges[i]->verts[1]->z);
        elist.edges[i]->length = length(v1-v2);
        //v1.set(elist.edges[i]->verts[1]->x-elist.edges[i]->verts[0]->x,
        //	elist.edges[i]->verts[1]->y-elist.edges[i]->verts[0]->y,
        //	elist.edges[i]->verts[1]->z-elist.edges[i]->verts[0]->z);
        //elist.edges[i]->length=length(v1);

        if(elist.edges[i]->length<shortest_edgelength)
            shortest_edgelength = elist.edges[i]->length;
    }
}


void Polyhedron::build_Edges()
{
    int i;
    unsigned char j;
    Triangle *f;
    Vertex *v1,*v2;
    double count = 0;

    /* count up how many edges we may require */

    for (i = 0; i < tlist.ntris; i++) {
        f = tlist.tris[i];
        for (j = 0; j < f->nverts; j++) {
            v1 = f->verts[j];
            v2 = f->verts[(j+1) % f->nverts];
            Triangle *result = find_common_Edge (f, v1, v2);
            if (result)
                count += 0.5;
            else
                count += 1;
        }
    }

    /*
    printf ("counted %f edges\n", count);
    */

    /* create space for edge list */

    elist.curMaxNumEdges = (int) (count + 10);  /* leave some room for expansion */
    elist.edges = new Edge *[elist.curMaxNumEdges];
    elist.nedges = 0;

    /* zero out all the pointers from faces to edges */

    for (i = 0; i < tlist.ntris; i++)
        for (j = 0; j < 3; j++)
            tlist.tris[i]->edges[j] = NULL;

    /* create all the edges by examining all the triangles */

    for (i = 0; i < tlist.ntris; i++) {
        f = tlist.tris[i];
        for (j = 0; j < 3; j++) {
            /* skip over edges that we've already created */
            if (f->edges[j])
                continue;
            v1 = f->verts[j];
            v2 = f->verts[(j+1) % f->nverts];
            create_Edge (v1, v2);
        }
    }

    /* add the edges to the corresponding vertices */
    for(i = 0; i < elist.nedges; i++)
    {
        add_Edge_to_Vertex(elist.edges[i], elist.edges[i]->verts[0]);
        add_Edge_to_Vertex(elist.edges[i], elist.edges[i]->verts[1]);
    }
}


void Polyhedron::create_Edge(Vertex *v1, Vertex *v2)
{
    int i,j;
    Triangle *f;

    /* make sure there is enough room for a new edge */

    if (elist.nedges >= elist.curMaxNumEdges) {

        elist.curMaxNumEdges += 100;
        Edge **list = new Edge *[elist.curMaxNumEdges];

        /* copy the old list to the new one */
        for (i = 0; i < elist.nedges; i++)
            list[i] = elist.edges[i];

        /* replace list */
        free (elist.edges);
        elist.edges = list;
    }

    /* create the edge */

    elist.edges[elist.nedges] = new Edge;
    Edge *e = elist.edges[elist.nedges];
    e->index = elist.nedges;
    e->verts[0] = v1;
    e->verts[1] = v2;
    elist.nedges++;

    /* count all triangles that will share the edge, and do this */
    /* by looking through all faces of the first vertex */

    e->ntris = 0;

    for (i = 0; i < v1->ntris; i++) {
        f = v1->tris[i];
        /* examine the vertices of the face for a match with the second vertex */
        for (j = 0; j < 3; j++) {
            /* look for a match */
            if (f->verts[j] == v2) {
                e->ntris++;
                break;
            }
        }
    }

    /* make room for the face pointers (at least two) */
    if (e->ntris < 2)
    {
        e->tris = new Triangle *[2];
        e->tris[0] = e->tris[1] = NULL;
    }

    else
    {
        e->tris = new Triangle *[e->ntris];
        e->tris[0] = e->tris[1] = NULL;
    }

    /* create pointers from edges to faces and vice-versa */

    e->ntris = 0; /* start this out at zero again for creating ptrs to tris */

    for (i = 0; i < v1->ntris; i++) {

        f = v1->tris[i];

        /* examine the vertices of the face for a match with the second vertex */
        for (j = 0; j < 3; j++)
            if (f->verts[j] == v2) {

                e->tris[e->ntris] = f;
                e->ntris++;

                if (f->verts[(j+1)%3] == v1)
                    f->edges[j] = e;
                else if (f->verts[(j+2)%3] == v1)
                    f->edges[(j+2)%3] = e;
                else {
                    fprintf (stderr, "Non-recoverable inconsistancy in create_edge()\n");
                    exit (-1);
                }

                break;  /* we'll only find one instance of v2 */
            }

    }

    /*initialize the edge content*/

    e->attp = NULL;
    e->sep = NULL;
    e->att_visit = 0;
    e->sep_visit = 0;
    e->valid = false;
    e->find_attp = false;
    e->find_sep = false;
}





