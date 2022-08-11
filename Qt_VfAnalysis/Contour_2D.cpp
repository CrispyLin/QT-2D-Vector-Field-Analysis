
#include "VField.h"

#include "FileLoader/PlyLoader.h"
#include "Graph2D.h"

#include "GL_LIB/glew.h"


extern Polyhedron *object;


/**************************************************************************
    Implementation of contour (level set) extraction
*/

void
Contour_Graph::extract_iso_contour_from_distance_field(double iso_val)
{
    /* First, we need to find out all the edges that may have this iso_val
     We use MidPointID of each edge to record the corresponding node id in the contour,
     use find_attp to marked the edges with the isoline crossing it
     */

    this->iso_val = iso_val;

    int i, j;

    for (i=0; i<object->elist.nedges; i++)
    {
        object->elist.edges[i]->find_attp = false;
        object->elist.edges[i]->MidPointID = -1;
        object->elist.edges[i]->visited = false;
    }

    for (i=0; i<object->elist.nedges; i++)
    {
        Edge *e=object->elist.edges[i];
        Vertex *v1=e->verts[0];
        Vertex *v2=e->verts[1];

        if ((v1->vf_mag <= iso_val && v2->vf_mag >= iso_val)
            || (v1->vf_mag >= iso_val && v2->vf_mag <= iso_val))
        {
            e->find_attp = true;

            // Assume v1 has smaller value
            if (v2->vf_mag < v1->vf_mag)
            {
                Vertex *tmp = v1;
                v1 = v2;
                v2 = tmp;
            }

            double pos[3];
            double t = (iso_val-v1->vf_mag)/(v2->vf_mag-v1->vf_mag);
            pos[0] = v1->x + t*(v2->x-v1->x);
            pos[1] = v1->y + t*(v2->y-v1->y);
            pos[2] = v1->z + t*(v2->z-v1->z);

            // Add a new node to the vlist
            ContourVertex one_v;
            one_v.pos.set(pos);
            one_v.edges.clear();
            one_v.at_edge = i;
            one_v.visited = false;
            one_v.stop = false;
            one_v.index = vlist.size();
            add_one_node(one_v);

            //
            e->MidPointID = one_v.index;
        }
    }


    // next, we need to connect these discrete points to get some isolines
    for (i=0; i<object->tlist.ntris; i++)
    {
        Triangle *f = object->tlist.tris[i];

        /*
            There could have two, three or non isolines lying within a triangle
        */

        std::vector<Edge*> with_isolines;
        for (j=0; j<f->nverts; j++)
        {
            if (f->edges[j]->find_attp)
                with_isolines.push_back(f->edges[j]);
        }

        if (with_isolines.size()<2)
            continue;

        else if (with_isolines.size()==2)
        {
            // add one edge to the contour
            ContourEdge one_e;
            one_e.at_triangle = i;
            one_e.index = elist.size();
            one_e.verts[0] = with_isolines[0]->MidPointID;
            one_e.verts[1] = with_isolines[1]->MidPointID;
            one_e.visited = false;
            icVector3 dist_vec=vlist[one_e.verts[0]].pos - vlist[one_e.verts[1]].pos;
            one_e.length = length(dist_vec);
            add_one_edge(one_e);
        }

        else
        {
            // They should form a loop
            for (j=0; j<with_isolines.size(); j++)
            {
                // add one edge to the contour
                ContourEdge one_e;
                one_e.at_triangle = i;
                one_e.index = elist.size();
                one_e.verts[0] = with_isolines[j]->MidPointID;
                one_e.verts[1] = with_isolines[(j+1)%3]->MidPointID;
                one_e.visited = false;
                icVector3 dist_vec=vlist[one_e.verts[0]].pos - vlist[one_e.verts[1]].pos;
                one_e.length = length(dist_vec);
                add_one_edge(one_e);
            }
        }
    }

}


void
Contour_Graph::display()
{
    int i;

    glColor3f (0, 0, 0);
    icVector3 norm;
    glBegin(GL_LINES);
    for (i=0; i<elist.size(); i++)
    {
        ContourEdge &e=elist[i];
        ContourVertex &v1=vlist[e.verts[0]];
        ContourVertex &v2=vlist[e.verts[1]];
        norm = object->tlist.tris[e.at_triangle]->normal;

        icVector3 offset = v1.pos+0.0012*norm;
        glVertex3dv(offset.entry);
        offset = v2.pos+0.0012*norm;
        glVertex3dv(offset.entry);
    }
    glEnd();
}
