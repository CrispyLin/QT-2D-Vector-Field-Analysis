/*
This file contains routines for Morse Decomposition

Created and Modified by Guoning Chen
        copyright @2007
    */

#include "Analysis/MorseDecomp.h"

extern Polyhedron *object;

MorseDecomp *morse_decomp = nullptr;
MorseDecomp *local_decomp = nullptr;

SCCList *tscclist=nullptr;

extern MCG_Graph *mcg;
extern PeriodicOrbitList *periodic_orbits;
extern bool RemoveDisconnMSOn;

/*
   A variable used to record the tau value that is used to compute a FG locally or globally
*/
double used_tau = 0;

double edge_sample_error = 1.e-15;

extern  bool EnConleyIndexComp;
extern int ndisplay_POs;

extern bool enFastRefinement;

extern int Integrator_opt;

/*
   07/21/2010 For debugging purpose
*/
MorseDecomp *L1_morse=nullptr;
MorseDecomp *L2_morse=nullptr;

/*
    To store the sample point list  02/09/2010
*/
EdgeSamplePt_List *edge_samples = nullptr;

/*
Entry: Edge *cur_edge
       int triangle
*/
bool Edge::is_exit_edge(int tri)
{
    double dot1, dot2, dotresult = 0;
    icVector2 vec[2];

    get_vecs_at_edge_endingpts(this, tri, vec);

    /////
    dot1 = dot(this->normal_2d, vec[0]);
    dot2 = dot(this->normal_2d, vec[1]);

    dotresult = max(dot1, dot2); /*this already consider the mixed edges*/

    if(dotresult >= 0) return true;

    return false;
}


bool
Edge::is_pure_exit_edge(int tri)
{
    double dot1, dot2, dotresult = 0;
    icVector2 vec[2];

    //get_vecs_at_edge_endingpts(this, tri, vec);
    Triangle *t = this->tris[0];

    if (t->index != tri)
        t = this->tris[1];

    vec[0].entry[0] = dot(this->verts[0]->t_vec, t->LX);
    vec[0].entry[1] = dot(this->verts[0]->t_vec, t->LY);

    vec[1].entry[0] = dot(this->verts[1]->t_vec, t->LX);
    vec[1].entry[1] = dot(this->verts[1]->t_vec, t->LY);

    /////
    dot1 = dot(this->normal_2d, vec[0]);
    dot2 = dot(this->normal_2d, vec[1]);

    dotresult = min(dot1, dot2); /*this already consider the mixed edges*/

    if(dotresult > 0) return true;

    return false;
}


bool
Edge::is_pure_entrance_edge(int tri)
{
    double dot1, dot2, dotresult = 0;
    icVector2 vec[2];

    //get_vecs_at_edge_endingpts(this, tri, vec);
    Triangle *t = this->tris[0];

    if (t->index != tri)
        t = this->tris[1];

    vec[0].entry[0] = dot(this->verts[0]->t_vec, t->LX);
    vec[0].entry[1] = dot(this->verts[0]->t_vec, t->LY);

    vec[1].entry[0] = dot(this->verts[1]->t_vec, t->LX);
    vec[1].entry[1] = dot(this->verts[1]->t_vec, t->LY);

    dot1 = dot(this->normal_2d, vec[0]);
    dot2 = dot(this->normal_2d, vec[1]);

    dotresult = max(dot1, dot2); //pure entrance edge growing

    if(dotresult < 0) return true;

    return false;
}


//bool IsEntranceEdge(Edge *cur_edge, int triangle)
//{
//	double dot1, dot2, dotresult = 0;
//	icVector2 vec[2];
//
//	GetVecsAtEdgeEndingPts(cur_edge, triangle, vec);
//
//	dot1 = dot(cur_edge->attract_normal, vec[0]);
//	dot2 = dot(cur_edge->attract_normal, vec[1]);
//
//	dotresult = max(dot1, dot2); //pure entrance edge growing
//
//	if(dotresult < 0) return true;
//
//	return false;
//}


bool Edge::is_mixed_edge(int triangle)
{
    double dot1, dot2, dotresult = 0;
    icVector2 vec[2];

    get_vecs_at_edge_endingpts(this, triangle, vec);

    dot1 = dot(this->normal_2d, vec[0]);
    dot2 = dot(this->normal_2d, vec[1]);

    if(dot1*dot2 < 0)  //do we need to consider tangent case here "=" 07/10/06
        return true;

    //double ratio;
    //ratio = (this->normal_2d.entry[0]*vec[0].entry[0]+vec[0].entry[1]*this->normal_2d.entry[1])/
    //	((vec[0].entry[0]-vec[1].entry[0])*this->normal_2d.entry[0]+(vec[0].entry[1]-vec[1].entry[1])*this->normal_2d.entry[1]);

    //if(ratio < 0.5 || ratio > 0.5)
    //if(ratio < 0.01 || ratio > 0.99)
    ////if(ratio < MixedEdgeRatio || ratio > 1-MixedEdgeRatio)
    //	return false;

    return false;
}


////Get the corresponding vectors on the ending vertices of the specific edge under the
////local frame of the given triangle 1/21/06
void Edge::get_vecs_at_edge_endingpts(Edge *cur_edge, int triangle, icVector2 vec[2])
{
    Triangle *face = object->tlist.tris[triangle];
    int i, direct_id = 0;         // index to read the directional vectors
    int deficit_index = -1;
    icVector2 direct_vec[6];

    Vertex *v1, *v2;
    v1 = cur_edge->verts[0];
    v2 = cur_edge->verts[1];

    ////Find the corresponding vectors for the two vertices
    for(i = 0; i < face->nverts; i++)
    {
        direct_vec[direct_id] = face->dir_vec[direct_id];
        direct_id ++;

        v1 = face->verts[i];
        if(v1->angle_deficit == 1) //if angle difict happens, we need to store the second vector
        {
            direct_vec[direct_id] = face->dir_vec[direct_id];
            direct_id ++;

            deficit_index = i;    //record the index of the vertex that has angle deficit
        }
    }

    if(deficit_index == -1)  //there is no angle deficit for this triangle
    {
        ////Get the vector at the first vertex
        for(i = 0; i < face->nverts; i++)
        {
            if(v1 == face->verts[i])
            {
                vec[0] = direct_vec[i];
                break;
            }
        }

        ////Get the vector at the second vertex
        for(i = 0; i < face->nverts; i++)
        {
            if(v2 == face->verts[i])
            {
                vec[1] = direct_vec[i];
                return;
            }
        }
    }

    else   ////if there is angle deficit existing, things will be somewhat complexed
    {
        if(!face->orient) //counter clockwise orientation
            get_vecs_withangledeficit(v1, v2, face, deficit_index, vec);
        else
            get_vecs_withangledeficit_2(v1, v2, face, deficit_index, vec);
    }
}

////Here we suppose there is only one vertex for each triangle having angle deficit
////This routine may have some problems for assign the correct vectors at the vertex with angle deficit
////1/22/06
void Edge::get_vecs_withangledeficit(Vertex *v1, Vertex *v2, Triangle *face,
                                     int deficit_index, icVector2 vec[2])
{
    //The current codes suppose that the orientation is counter-clockwise
    switch(deficit_index)
    {
    case 0:  //if the first vertex has angle deficit
        if(v1 != face->verts[0] && v2 != face->verts[0])
        {
            if(v1 == face->verts[1]) //v1 = vert[1], v2 = vert[2]
            {
                vec[0] = face->dir_vec[2];
                vec[1] = face->dir_vec[3];
            }
            else                             //v1 = vert[2], v2 = vert[1]
            {
                vec[0] = face->dir_vec[3];
                vec[1] = face->dir_vec[2];
            }
        }

        else if(v1 == face->verts[0])
        {
            if(v2 == face->verts[2]) //v1 = vert[0], v2 = vert[2]
            {
                vec[0] = face->dir_vec[0];
                vec[1] = face->dir_vec[3];
            }
            else                             //v1 = vert[0], v2 = vert[1]
            {
                vec[0] = face->dir_vec[1];
                vec[1] = face->dir_vec[2];
            }
        }

        else if(v2 == face->verts[0])
        {
            if(v1 == face->verts[2]) //v1 = vert[2], v2 = vert[0]
            {
                vec[1] = face->dir_vec[0];
                vec[0] = face->dir_vec[3];
            }
            else                             //v1 = vert[1], v2 = vert[0]
            {
                vec[1] = face->dir_vec[1];
                vec[0] = face->dir_vec[2];
            }
        }

        break;

    case 1:
        if(v1 != face->verts[1] && v2 != face->verts[1])
        {
            if(v1 == face->verts[0])
            {
                vec[0] = face->dir_vec[0];
                vec[1] = face->dir_vec[3];
            }
            else
            {
                vec[0] = face->dir_vec[3];
                vec[1] = face->dir_vec[0];
            }
        }

        else if(v1 == face->verts[1])
        {
            if(v2 == face->verts[0])
            {
                vec[0] = face->dir_vec[1];
                vec[1] = face->dir_vec[0];
            }
            else
            {
                vec[0] = face->dir_vec[2];
                vec[1] = face->dir_vec[3];
            }
        }

        else if(v2 == face->verts[1])
        {
            if(v1 == face->verts[0])
            {
                vec[1] = face->dir_vec[1];
                vec[0] = face->dir_vec[0];
            }
            else
            {
                vec[1] = face->dir_vec[2];
                vec[0] = face->dir_vec[3];
            }
        }
        break;

    case 2: //verts[2] has angle deficit
        if(v1 != face->verts[2] && v2 != face->verts[2])
        {
            if(v1 == face->verts[0])  //v1 = verts[0], v2 = verts[1]
            {
                vec[0] = face->dir_vec[0];
                vec[1] = face->dir_vec[1];
            }
            else                              //v1 = verts[1], v2 = verts[0]
            {
                vec[0] = face->dir_vec[1];
                vec[1] = face->dir_vec[0];
            }
        }

        else if(v1 == face->verts[2])
        {
            if(v2 == face->verts[0])  //v1 = verts[2], v2 = verts[0]
            {
                vec[0] = face->dir_vec[3];
                vec[1] = face->dir_vec[0];
            }
            else                              //v1 = verts[2], v2 = verts[1]
            {
                vec[0] = face->dir_vec[2];
                vec[1] = face->dir_vec[1];
            }
        }

        else if(v2 == face->verts[2])
        {
            if(v1 == face->verts[0])  //v1 = verts[0], v2 = verts[2]
            {
                vec[1] = face->dir_vec[3];
                vec[0] = face->dir_vec[0];
            }
            else                              //v1 = verts[1], v2 = verts[2]
            {
                vec[1] = face->dir_vec[2];
                vec[0] = face->dir_vec[1];
            }
        }

        break;
    }
}


/*
Modified at 06/10/09
Consider clockwise orientation triangle case
*/
void Edge::get_vecs_withangledeficit_2(Vertex *v1, Vertex *v2, Triangle *face,
                                       int deficit_index, icVector2 vec[2])
{
    //The current codes suppose that the orientation is counter-clockwise
    switch(deficit_index)
    {
    case 0:  //if the first vertex has angle deficit
        if(v1!= face->verts[0] && v2 != face->verts[0])
        {
            vec[0] = face->dir_vec[2];
            vec[1] = face->dir_vec[3];
        }

        else if(v1 == face->verts[0])
        {
            if(v2 == face->verts[2]) //v1 = vert[0], v2 = vert[2]
            {
                vec[0] = face->dir_vec[1];
                vec[1] = face->dir_vec[3];
            }
            else                             //v1 = vert[0], v2 = vert[1]
            {
                vec[0] = face->dir_vec[0];
                vec[1] = face->dir_vec[2];
            }
        }

        else if(v2 == face->verts[0])
        {
            if(v1 == face->verts[2]) //v1 = vert[2], v2 = vert[0]
            {
                vec[1] = face->dir_vec[1];
                vec[0] = face->dir_vec[3];
            }
            else                             //v1 = vert[1], v2 = vert[0]
            {
                vec[1] = face->dir_vec[0];
                vec[0] = face->dir_vec[2];
            }
        }

        break;

    case 1: //verts[1] has angle deficit
        if(v1 != face->verts[1] && v2 != face->verts[1])
        {
            vec[0] = face->dir_vec[0];
            vec[1] = face->dir_vec[3];
        }

        else if(v1 == face->verts[1])
        {
            if(v2 == face->verts[0])
            {
                vec[0] = face->dir_vec[2];
                vec[1] = face->dir_vec[0];
            }
            else
            {
                vec[0] = face->dir_vec[1];
                vec[1] = face->dir_vec[3];
            }
        }

        else if(v2 == face->verts[1])
        {
            if(v1 == face->verts[0])
            {
                vec[1] = face->dir_vec[2];
                vec[0] = face->dir_vec[0];
            }
            else
            {
                vec[1] = face->dir_vec[1];
                vec[0] = face->dir_vec[3];
            }
        }
        break;

    case 2: //verts[2] has angle deficit
        if(v1 != face->verts[2] && v2 != face->verts[2])//v1 = verts[1], v2 = verts[0]
        {
            vec[0] = face->dir_vec[0];
            vec[1] = face->dir_vec[1];
        }

        else if(v1 == face->verts[2])
        {
            if(v2 == face->verts[0])  //v1 = verts[2], v2 = verts[0]
            {
                vec[0] = face->dir_vec[2];
                vec[1] = face->dir_vec[0];
            }
            else                              //v1 = verts[2], v2 = verts[1]
            {
                vec[0] = face->dir_vec[3];
                vec[1] = face->dir_vec[1];
            }
        }

        else if(v2 == face->verts[2])
        {
            if(v1 == face->verts[0])  //v1 = verts[0], v2 = verts[2]
            {
                vec[1] = face->dir_vec[2];
                vec[0] = face->dir_vec[0];
            }
            else                              //v1 = verts[1], v2 = verts[2]
            {
                vec[1] = face->dir_vec[3];
                vec[0] = face->dir_vec[1];
            }
        }

        break;
    }
}



/*------------------------------------------------------------------*/

/*The implementation of the Morse Decomposition*/


/*
~Geometry-based Morse Decomposition
Build the directed graph based on the vector field without using multi-valued map
Note: we assume the directed graph has been initialized before calling this routine
*/
void MorseDecomp::build_direted_graph()
{
    int i;

    /*reset all the edges*/
    for(i = 0; i < object->elist.nedges; i++)
    {
        object->elist.edges[i]->visited = false;
    }

    for(i = 0; i < object->tlist.ntris; i++)
    {
        //Add the face to the node list

        dg->nlist->dirnodes[i]->node_index = i;

        build_DirGraph_Tri_no_Tau(i,0);
    }

    /*write the graph (edge list) into a file for comparison */
    //FILE *fp = fopen("dirgraph1.txt", "w");
    //for(i = 0; i < dg->elist->nedges; i++)
    //{
    //	fprintf(fp, "%d: %d, %d\n", i, dg->elist->edges[i]->node_index1, dg->elist->edges[i]->node_index2);
    //}
    //fclose(fp);
}


/*
Compute the directed edges for one particular triangle according to
the flow behaviour cross the edges of this triangle
*/
void MorseDecomp::build_DirGraph_Tri_no_Tau(int tri, int backward)
{
    Triangle *face = object->tlist.tris[tri];
    Edge *cur_e;
    int j, node1, node2;

    int cur_end_index;

    for(j = 0; j < 3; j++)
    {
        cur_e = face->edges[j];

        if(cur_e->visited)
            continue;

        ////Calculate the outward normal of the edge
        if(cur_e->tris[0]==nullptr || cur_e->tris[1]==nullptr) /*don't consider boundary edges*/
            continue;

        if(face == cur_e->tris[0])
            cur_e->cal_normal_at_edge(cur_e, face, 0);
        else
            cur_e->cal_normal_at_edge(cur_e, face, 1);

        //build the directed connection according to the type of the edge
        if(cur_e->is_mixed_edge(tri))
        {
            node1 = tri;   //go from currrent triangle
            node2 = cur_e->tris[0]->index;
            if(node2 == node1)
                node2 = cur_e->tris[1]->index;

            if(node2 < 0)  // we reach boundary
                continue;

            //if(dg->is_repeated_edge(node1, node2))
            //	continue;

            //if(dg->is_repeated_edge_2(node1, node2))  /*  Modified by Guoning 09/24/09 */
            //	continue;

            //add to the edge list
            dg->add_to_edgelist(node1, node2, cur_end_index);
            ////add the edge to the nodes
            dg->add_edge_to_node(node1, cur_end_index-1);
            dg->add_edge_to_node(node2, cur_end_index-1);

            /*add the second directed edge*/
            int temp = node1;
            node1 = node2;
            node2 = temp;
            //add to the edge list
            dg->add_to_edgelist(node1, node2, cur_end_index);
            ////add the edge to the nodes
            dg->add_edge_to_node(node1, cur_end_index-1);
            dg->add_edge_to_node(node2, cur_end_index-1);

            cur_e->visited = true;
        }

        else{
            if(backward == 0)
            {
                if(cur_e->is_exit_edge(tri))
                {
                    node1 = tri;   //go from currrent triangle
                    node2 = cur_e->tris[0]->index;
                    if(node2 == node1)
                        node2 = cur_e->tris[1]->index;

                    if(node2 < 0)  // we reach boundary
                        continue;

                    //if(dg->is_repeated_edge(node1, node2))
                    //	continue;

                    //add to the edge list
                    dg->add_to_edgelist(node1, node2, cur_end_index);
                    ////add the edge to the nodes
                    dg->add_edge_to_node(node1, cur_end_index-1);
                    dg->add_edge_to_node(node2, cur_end_index-1);

                    //cur_e->visited = true;
                }
            }
        }

    }
}



void
MorseDecomp::build_direted_graph_2()
{
    int i;

    /*reset all the edges*/
    for(i = 0; i < object->elist.nedges; i++)
    {
        object->elist.edges[i]->visited = false;
    }

    for (i=0; i<object->tlist.ntris; i++)
    {
        dg->nlist->dirnodes[i]->node_index = i;
    }

    for(i = 0; i < object->elist.nedges; i++)
    {
        //Add the face to the node list

        //build_DirGraph_Tri_no_Tau(i,0);
        build_DirEdges_edge (i);
    }

    /*write the graph (edge list) into a file for comparison */
    //FILE *fp = fopen("dirgraph1.txt", "w");
    //for(i = 0; i < dg->elist->nedges; i++)
    //{
    //	fprintf(fp, "%d: %d, %d\n", i, dg->elist->edges[i]->node_index1, dg->elist->edges[i]->node_index2);
    //}
    //fclose(fp);
}


void
MorseDecomp::build_DirEdges_edge(int which_edge)
{
    Edge *cur_e = object->elist.edges[which_edge];
    int cur_end_index;

    int node1, node2;
    if (cur_e->tris[1] == nullptr) return;

    node1 = cur_e->tris[0]->index;
    node2 = cur_e->tris[1]->index;

    //if (node1 < 0 || node2 < 0) return;  // boundary edge

    cur_e->cal_normal_through_edge(cur_e, cur_e->tris[0], 0);  // compute only outer normal

    if (cur_e->is_pure_exit_edge(node1))
    {
        //if(dg->is_repeated_edge(node1, node2))
        //	return;

        //add to the edge list
        dg->add_to_edgelist(node1, node2, cur_end_index);
        ////add the edge to the nodes
        dg->add_edge_to_node(node1, cur_end_index-1);
        dg->add_edge_to_node(node2, cur_end_index-1);
    }

    else if (cur_e->is_pure_entrance_edge(node1))
    {
        //add to the edge list
        dg->add_to_edgelist(node2, node1, cur_end_index);
        ////add the edge to the nodes
        dg->add_edge_to_node(node1, cur_end_index-1);
        dg->add_edge_to_node(node2, cur_end_index-1);
    }

    else
    {
        //add to the edge list
        dg->add_to_edgelist(node1, node2, cur_end_index);
        ////add the edge to the nodes
        dg->add_edge_to_node(node1, cur_end_index-1);
        dg->add_edge_to_node(node2, cur_end_index-1);

        //add to the edge list
        dg->add_to_edgelist(node2, node1, cur_end_index);
        ////add the edge to the nodes
        dg->add_edge_to_node(node1, cur_end_index-1);
        dg->add_edge_to_node(node2, cur_end_index-1);
    }

}

/*
*/
void MorseDecomp::init_graph()
{
    int i;

    FILE *fp;
    //fp = fopen("detect_porbit_cooling.txt", "a");
    //fprintf(fp, "Start initializing Directed Graph...\n");
    //fclose(fp);

    if(dg != nullptr)
        delete dg;

    //fp = fopen("detect_porbit_cooling.txt", "a");
    //fprintf(fp, "Start initializing dg...\n");
    //fprintf(fp, "nodes %d, edges %d...\n", object->tlist.ntris, object->elist.nedges);
    //fclose(fp);

    /*this initialization is suitable for Morse decomposition without using \tau */
    //dg = new DirGraph(object->tlist.ntris, object->elist.nedges*2);
    dg = new DirGraph(object->tlist.ntris, object->tlist.ntris*10); // For HCCI data set

    if(dg->elist == nullptr || dg->nlist->dirnodes == nullptr || dg->elist->edges == nullptr)
    {
        char rout[256], var[256];
        sprintf(rout, "%s", "MorseDecomp::init_graph");
        if(dg->nlist == nullptr)
            sprintf(var, "%s", "dg->nlist");
        else
            sprintf(var, "%s", "dg->elist");

        //MessageBox(nullptr, "failed to reallocate memory for the directed graph!", "Error", MB_OK);
        //write_mem_error(rout, var, 0);
        exit(-1);
    }

    //fp = fopen("detect_porbit_cooling.txt", "a");
    //fprintf(fp, "The address for nlist is %f...\n", dg->nlist);
    //fprintf(fp, "The address for elist is %f...\n", dg->elist);
    //fprintf(fp, "Start initializing nlist...\n");
    //fclose(fp);

    /*allocate the real space for the node list and edge list*/
    for(i = 0; i < dg->nlist->curMaxNumENodes; i++)
    {
        dg->nlist->dirnodes[i] = new DirGraph_Node();
        //dg->nlist->dirnodes[i]->edges = nullptr;
        //dg->nlist->dirnodes[i]->nedges = 0;
    }

    //fp = fopen("detect_porbit_cooling.txt", "a");
    //fprintf(fp, "Start initializing elist...\n");
    //fclose(fp);

    for(i = 0; i < dg->elist->curMaxNumGedges; i++)
    {
        dg->elist->edges[i] = new Graph_Edge();
    }

    /*initialize the SCC list as well here*/

    //if(scclist != nullptr)
    //	delete scclist;

    //scclist = new SCCList();  /*initial size is 1000 now*/

}

/*
build the SCC list according to the obtained SCC labels of the nodes
*/

void MorseDecomp::build_SCCElemList()
{
    int i;
    ////first, we need to search the sub-trees to figure how many scc in the graph and
    ////what are they
    for(i = 0; i < dg->nlist->ndirnodes; i++)
    {
        dg->nlist->dirnodes[i]->visited = 0;
    }

    /* To get the elements in each component */
    int cur_sccomp = 0;
    int num_morsesets = 0;

    //allocate space for the SCC components

    if(scclist != nullptr)
    {
        delete scclist;
    }

    scclist = new SCCList(dg->num_sccomps);
    scclist->nsccs = dg->num_sccomps;

    for(i = 0; i < dg->num_sccomps; i++)
        scclist->scccomponents[i] = new SCComponent();

    int *numtriangles_eachscc = (int*)malloc(sizeof(int)*(dg->num_sccomps));

    for(i = 0; i < dg->num_sccomps; i++)
        numtriangles_eachscc[i] = 0;

    //first, we search the whole mesh once, and get the number of triangles in each scc components
    for(i = 0; i < dg->nlist->ndirnodes; i++)
    {
        numtriangles_eachscc[dg->nlist->dirnodes[i]->sscomp_index] += 1;
    }

    /*write the number of triangles in each SCC into a file 04/22/07*/
    //FILE *fp = fopen("ntris_SCCs1.txt", "w");
    //for(i = 0; i < dg->num_sccomps; i++)
    //{
    //	if(numtriangles_eachscc[i]>100)
    //	fprintf(fp, "%d: %d\n", i, numtriangles_eachscc[i]);
    //}
    //fclose(fp);

    //allocate space for each scc
    for(i = 0; i < dg->num_sccomps; i++)
    {
        scclist->scccomponents[i]->nodes = (int*) malloc(sizeof(int)*numtriangles_eachscc[i]);
        scclist->scccomponents[i]->nnodes = 0;

        ////initialize
        scclist->scccomponents[i]->singular_tri = nullptr;
        scclist->scccomponents[i]->nfixedpoints = 0;
        scclist->scccomponents[i]->periodicorbits = nullptr;
        scclist->scccomponents[i]->nperiodicorbits = 0;
        scclist->scccomponents[i]->nattpts = 0;
        scclist->scccomponents[i]->nseppts = 0;

        scclist->scccomponents[i]->node_index = -1;
    }

    //record the scc components
    if(dg->nlist->ndirnodes!=object->tlist.ntris)//local
    {
        for(i = 0; i < dg->nlist->ndirnodes /*the number of total triangles*/; i++)
        {

            scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nodes[
                scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nnodes] = i;
            scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nnodes++;

            int gi;//global i
            gi=dg->nlist->dirnodes[i]->global_index;

            /* count the number of fixed points in the SCC 02/21/07 */
            if(object->tlist.tris[gi]->singularityID >= 0)
            {
                scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->singular_tri
                    =extend_link(scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->singular_tri,
                                  scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nfixedpoints);

                scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->singular_tri[
                    scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nfixedpoints] = gi;
                scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nfixedpoints++;
            }
        }

    }
    else//global
    {

        for(i = 0; i < dg->nlist->ndirnodes /*the number of total triangles*/; i++)
        {

            scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nodes[
                scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nnodes] = i;
            scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nnodes++;

            /* count the number of fixed points in the SCC 02/21/07 */
            if(object->tlist.tris[i]->singularityID >= 0)
            {
                scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->singular_tri
                    =extend_link(scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->singular_tri,
                                  scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nfixedpoints);

                scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->singular_tri[
                    scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nfixedpoints] = i;
                scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nfixedpoints++;
            }
        }
    }

    //Mark the boundary edges of each SCC component 07/18/06
    //Edge *cur_e;
    //Triangle *face;
    //int j, k;
    //Triangle *other_t;
    //for(i = 0; i < scclist->nsccs; i++)
    //{
    //	for(j = 0; j < scclist->scccomponents[i]->nnodes; j++)
    //	{
    //		face = object->tlist.tris[scclist->scccomponents[i]->nodes[j]];

    //		for(k = 0; k < 3; k++)
    //		{
    //			cur_e = face->edges[k];

    //			other_t = cur_e->tris[0];

    //			if(other_t == face)
    //				other_t = cur_e->tris[1];

    //			//Test whether it is a boundary edge
    //			if(IsRepeated(scclist->scccomponents[i]->nodes, other_t->index,
    //				scclist->scccomponents[i]->nnodes))
    //			{
    //				cur_e->OnBoundary = false;
    //				continue;  //not a boundary edge
    //			}

    //			cur_e->OnBoundary = true;
    //		}
    //	}
    //}

    if(numtriangles_eachscc != nullptr)
        free(numtriangles_eachscc);
}



void MorseDecomp::build_SCCElemList_local()
{
    int i;
    ////first, we need to search the sub-trees to figure how many scc in the graph and
    ////what are they
    for(i = 0; i < dg->nlist->ndirnodes; i++)
    {
        dg->nlist->dirnodes[i]->visited = 0;
    }

    /* To get the elements in each component */
    int cur_sccomp = 0;
    int num_morsesets = 0;

    //allocate space for the SCC components

    if(scclist != nullptr)
    {
        delete scclist;
    }

    scclist = new SCCList(dg->num_sccomps);
    scclist->nsccs = dg->num_sccomps;

    for(i = 0; i < dg->num_sccomps; i++)
        scclist->scccomponents[i] = new SCComponent();

    int *numtriangles_eachscc = (int*)malloc(sizeof(int)*(dg->num_sccomps));

    for(i = 0; i < dg->num_sccomps; i++)
        numtriangles_eachscc[i] = 0;

    //first, we search the whole mesh once, and get the number of triangles in each scc components
    for(i = 0; i < dg->nlist->ndirnodes; i++)
    {
        numtriangles_eachscc[dg->nlist->dirnodes[i]->sscomp_index] += 1;
    }

    /*write the number of triangles in each SCC into a file 04/22/07*/
    //FILE *fp = fopen("ntris_SCCs1.txt", "w");
    //for(i = 0; i < dg->num_sccomps; i++)
    //{
    //	if(numtriangles_eachscc[i]>100)
    //	fprintf(fp, "%d: %d\n", i, numtriangles_eachscc[i]);
    //}
    //fclose(fp);

    //allocate space for each scc
    for(i = 0; i < dg->num_sccomps; i++)
    {
        scclist->scccomponents[i]->nodes = (int*) malloc(sizeof(int)*numtriangles_eachscc[i]);
        scclist->scccomponents[i]->nnodes = 0;

        ////initialize
        scclist->scccomponents[i]->singular_tri = nullptr;
        scclist->scccomponents[i]->nfixedpoints = 0;
        scclist->scccomponents[i]->periodicorbits = nullptr;
        scclist->scccomponents[i]->nperiodicorbits = 0;
        scclist->scccomponents[i]->nattpts = 0;
        scclist->scccomponents[i]->nseppts = 0;

        scclist->scccomponents[i]->node_index = -1;
    }

    //record the scc components
    if(dg->nlist->ndirnodes!=object->tlist.ntris)//local
    {
        for(i = 0; i < dg->nlist->ndirnodes /*the number of total triangles*/; i++)
        {

            scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nodes[
                scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nnodes] = i;
            scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nnodes++;

            int gi;//global i
            gi=dg->nlist->dirnodes[i]->global_index;

            /* count the number of fixed points in the SCC 02/21/07 */
            if(object->tlist.tris[gi]->singularityID >= 0)
            {
                scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->singular_tri
                    =extend_link(scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->singular_tri,
                                  scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nfixedpoints);

                scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->singular_tri[
                    scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nfixedpoints] = gi;
                scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nfixedpoints++;
            }
        }

    }
    else//global
    {

        for(i = 0; i < dg->nlist->ndirnodes /*the number of total triangles*/; i++)
        {

            scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nodes[
                scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nnodes] = i;
            scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nnodes++;

            /* count the number of fixed points in the SCC 02/21/07 */
            if(object->tlist.tris[i]->singularityID >= 0)
            {
                scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->singular_tri
                    =extend_link(scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->singular_tri,
                                  scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nfixedpoints);

                scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->singular_tri[
                    scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nfixedpoints] = i;
                scclist->scccomponents[dg->nlist->dirnodes[i]->sscomp_index]->nfixedpoints++;
            }
        }
    }

    //Mark the boundary edges of each SCC component 07/18/06
    //Edge *cur_e;
    //Triangle *face;
    //int j, k;
    //Triangle *other_t;
    //for(i = 0; i < scclist->nsccs; i++)
    //{
    //	for(j = 0; j < scclist->scccomponents[i]->nnodes; j++)
    //	{
    //		face = object->tlist.tris[scclist->scccomponents[i]->nodes[j]];

    //		for(k = 0; k < 3; k++)
    //		{
    //			cur_e = face->edges[k];

    //			other_t = cur_e->tris[0];

    //			if(other_t == face)
    //				other_t = cur_e->tris[1];

    //			//Test whether it is a boundary edge
    //			if(IsRepeated(scclist->scccomponents[i]->nodes, other_t->index,
    //				scclist->scccomponents[i]->nnodes))
    //			{
    //				cur_e->OnBoundary = false;
    //				continue;  //not a boundary edge
    //			}

    //			cur_e->OnBoundary = true;
    //		}
    //	}
    //}

    if(numtriangles_eachscc != nullptr)
        free(numtriangles_eachscc);
}

/*
Judge whether a strongly connected component is valid or not
*/
bool MorseDecomp::is_valid_SCC(int scc_index)
{
    if(scclist->scccomponents[scc_index]->nnodes < 2
        || scclist->scccomponents[scc_index]->nnodes >= object->tlist.ntris)
        return false;  //the minimum number of triangles that can contain a limit cycle is 2

    /* this calculation may have problem 07/22/06 */

    Region_Tri_Index *temp_region = new Region_Tri_Index(scclist->scccomponents[scc_index]->nodes,
                                                         scclist->scccomponents[scc_index]->nnodes, 0);

    /*   check the connectivity of the Morse set 03/08/2010   */
    //if (!is_connected_reg (temp_region->tris, temp_region->ntris))
    //	return false;

    //if (!is_connected_reg(temp_region->tris, temp_region->ntris))
    //{
    //	delete temp_region;
    //	return false;
    //}
    int Euler_num = temp_region->cal_euler_value();


    //scclist->scccomponents[scc_index]->num_boundaries = 1 + (2-Euler_num)/2;


    // check the connectivity of the Morse set (it should be connected!)

    ////if the scc forms a topological disk
    if(Euler_num == 1)
    {
        if(!temp_region->contain_fixedpts())
        {
            delete temp_region;
            return false;  //do not find any singularity inside this SCC
        }
    }

    delete temp_region;
    return true;
}


void MorseDecomp::mark_all_valid_SCCS()
{

    for(int i = 0; i < scclist->nsccs; i++)
        scclist->scccomponents[i]->valid = true;


    for(int i = 0; i < scclist->nsccs; i++)
    {

        //if(scclist->scccomponents[i]->nnodes < 0
        //	|| scclist->scccomponents[i]->nnodes >= object->tlist.ntris)
        //{
        //	scclist->scccomponents[i]->valid = false;
        //	continue;
        //}

        if(scclist->scccomponents[i]->nnodes < 2)
        {
            if(scclist->scccomponents[i]->nfixedpoints == 0)
            {
                scclist->scccomponents[i]->valid = false;
                continue;
            }
            scclist->scccomponents[i]->valid = true;
            continue;

        }
        else if(!is_valid_SCC(i))
        {
            scclist->scccomponents[i]->valid = false;
            continue;
        }

        int *temp_tris = new int [scclist->scccomponents[i]->nnodes];
        for (int j=0; j<scclist->scccomponents[i]->nnodes; j++)
        {
            temp_tris[j]=scclist->scccomponents[i]->nodes[j];
        }

        /*
            Added by Guoning 07/01/2010
        */

        if (RemoveDisconnMSOn&&!is_connected_reg(temp_tris, scclist->scccomponents[i]->nnodes))
        {
            scclist->scccomponents[i]->valid = false;
            continue;
        }

        delete [] temp_tris;

        /*
           How to remove the boundary Morse sets?? 03/01/2010
        */
        //else
        //scclist->scccomponents[i]->valid = true;   //this is a valid SCC

    }
}


/*
  Modified by Guoning at 02/16/2010
*/
bool
MorseDecomp::mark_all_valid_SCCS(int *tris, int ntri)  // for local region refinement
{
    for(int i = 0; i < scclist->nsccs; i++)
        scclist->scccomponents[i]->valid = true;

    for(int i = 0; i < object->tlist.ntris; i++)
        object->tlist.tris[i]->MS_type = 0;    // set default label value for the triangle

    int *temp_tris = new int [ntri];


    for(int i = 0; i < scclist->nsccs; i++)
    {
        //if(i>=587)
        //{
        //FILE *fp=fopen("detect_porbit.txt", "a");
        //fprintf(fp, "start %d SCC with %d triangles. \n", i,
        //	scclist->scccomponents[i]->nnodes);
        //fprintf(fp, "They are %d and %d. \n", scclist->scccomponents[i]->nodes[0],
        //	scclist->scccomponents[i]->nodes[1]);
        //fclose(fp);
        //}


        //if(scclist->scccomponents[i]->nnodes <= 0
        //	|| scclist->scccomponents[i]->nnodes >= object->tlist.ntris)
        //{
        //	scclist->scccomponents[i]->valid = false;
        //	continue;
        //}

        for (int j=0; j<scclist->scccomponents[i]->nnodes; j++)
        {
            temp_tris[j]=tris[scclist->scccomponents[i]->nodes[j]];
        }

        if(scclist->scccomponents[i]->nnodes < 2)
        {
            if(scclist->scccomponents[i]->nfixedpoints == 0)
            {
                scclist->scccomponents[i]->valid = false;
                continue;
            }
            scclist->scccomponents[i]->valid = true;
            continue;
        }

        else if(!is_valid_scc(tris, ntri, i))
        {
            scclist->scccomponents[i]->valid = false;

            continue;
        }

        //scclist->scccomponents[i]->valid = is_valid_scc(tris, ntri, i);


        //// Added by Guoning: we disregard the refined Morse sets that become disconnected 06/29/2010
        //// but if we do it here, we will probably missing some features :(
        else if (RemoveDisconnMSOn&&!is_connected_reg(temp_tris, scclist->scccomponents[i]->nnodes))
        {
            scclist->scccomponents[i]->valid = false;
            //delete temp_tris;
            return false;  // this decomposition is not valid since we get disconnected strongly connected component
            //continue;
        }

        //scclist->scccomponents[i]->valid = true;   //this is a valid SCC
    }

    delete [] temp_tris;

    return true;
}



bool
MorseDecomp::is_valid_scc(int *tris, int ntri, int scc_index)
{
    if(scclist->scccomponents[scc_index]->nnodes < 2
        || scclist->scccomponents[scc_index]->nnodes >= object->tlist.ntris)
        return false;  //the minimum number of triangles that can contain a limit cycle is 2

    /* this calculation may have problem 07/22/06 */

    Region_Tri_Index *temp_region = new Region_Tri_Index(tris,
                                                         ntri, 0);

    /*   check the connectivity of the Morse set 03/08/2010   */
    //if (!is_connected_reg (temp_region->tris, temp_region->ntris))
    //	return false;

    int Euler_num = temp_region->cal_euler_value();


    //scclist->scccomponents[scc_index]->num_boundaries = 1 + (2-Euler_num)/2;

    ////if the scc forms a topological disk
    if(Euler_num == 1)
    {
        if(!temp_region->contain_fixedpts())
        {
            delete temp_region;
            //if(scc_index==587)
            //{
            //	FILE *fp=fopen("detect_porbit.txt", "a");
            //	fprintf(fp, "finish calculating the Euler characteristics\n");
            //	fprintf(fp, "the E number is %d.\n",Euler_num);
            //	fclose(fp);
            //}
            return false;  //do not find any singularity inside this SCC
        }
    }

    delete temp_region;
    return true;
}


void MorseDecomp::cal_sepsandatts_valid_SCCs()
{
    int i;

    for(i = 0; i < scclist->nsccs; i++)
    {
        if(scclist->scccomponents[i]->valid && scclist->scccomponents[i]->nnodes >= 2)
            scclist->scccomponents[i]->cal_sep_attp_pts();
    }
}



void MorseDecomp::morse_decomp()
{
    //FILE *fp;
    //fp = fopen("detect_porbit.txt", "w");
    //fprintf(fp, "Start initializing Directed Graph...\n");
    //fclose(fp);

    init_graph();

    //fp = fopen("detect_porbit.txt", "a");
    ////fp = fopen("detect_porbit_swirl.txt", "a");
    //fprintf(fp, "Start constructing Directed Graph...\n");
    //fclose(fp);

    //build_direted_graph();
    build_direted_graph_2();
    used_tau = 0;

    //fp = fopen("detect_porbit_cooling.txt", "a");
    //fp = fopen("detect_porbit.txt", "a");
    //fprintf(fp, "Finish constructing Directed Graph...\n");
    //fprintf(fp, "Start finding SCC...\n");
    ////fprintf(fp, "current date and time are :  %s. \n", ctime (&g_rawtime) );
    //fclose(fp);

    dg->find_SCCS();

    //fp = fopen("detect_porbit_cooling.txt", "a");
    //fp = fopen("detect_porbit.txt", "a");
    //fprintf(fp, "Finish finding SCC...\n");
    ////fprintf(fp, "current date and time are :  %s. \n", ctime (&g_rawtime) );
    //fprintf(fp, "Start building SCC List...\n");
    //fclose(fp);

    build_SCCElemList();

    //fp = fopen("detect_porbit.txt", "a");
    //fprintf(fp, "Finish finding SCC...\n");
    ////fprintf(fp, "current date and time are :  %s. \n", ctime (&g_rawtime) );
    //fprintf(fp, "finish building SCC List...\n");
    //fprintf(fp, "marking the valid SCCs...\n");
    //fprintf(fp, "the obtained SCCs %d\n", scclist->nsccs);
    //fclose(fp);

    mark_all_valid_SCCS();

    //fp = fopen("detect_porbit.txt", "a");
    //fprintf(fp, "finish building SCC List...\n");
    //fprintf(fp, "finish marking the valid SCCs...\n");
    //fclose(fp);
}



/*use the idea of \tau maps*/
void MorseDecomp::morse_decomp_tau(double tau)
{
    FILE *fp;
    clock_t start, finish;

    start = clock();

    init_graph();

    used_tau = tau;

    if(fabs(tau)<1.e-8)
        //build_direted_graph();
        build_direted_graph_2();   //geometric based

    else
        build_multivalued_graph(tau);   // tau-map based

    finish = clock();

    fp = fopen("tvcg_performance.txt", "w");
    fprintf(fp, "time for constructing FG is %f. \n",(double)(finish - start)/CLOCKS_PER_SEC);
    fclose(fp);

    start = clock();
    dg->find_SCCS();


    finish = clock();

    fp = fopen("tvcg_performance.txt", "a");
    fprintf(fp, "time for extracting SCC is %f. \n",(double)(finish - start)/CLOCKS_PER_SEC);
    fclose(fp);

    //FILE *fp1=fopen("cooling_test.txt", "a");
    //fprintf(fp1, "finish finding SCCs.\n");
    //fclose(fp1);

    build_SCCElemList();

    //fp1=fopen("cooling_test.txt", "a");
    //fprintf(fp1, "finish building the SCC list.\n");
    //fclose(fp1);

    mark_all_valid_SCCS();

    //fp1=fopen("cooling_test.txt", "a");
    //fprintf(fp1, "finish marking the valid SCCs.\n");
    //fclose(fp1);
}


void MorseDecomp::build_multivalued_graph(double tau)
{
    /*to visualize the relationship between spatial tau and temporal tau
    08/01/07*/
    //init_samplepts_tautracing();

    /***1. first trace backward*/
    trace_all_Verts(tau, 1);

    //FILE *fp=fopen("cooling_test.txt", "w");
    //fprintf(fp, "finish vertex tracing (backward).\n");
    //fclose(fp);

    _cprintf("finish vertex tracing (backward).\n");

    /*1.2. build the edges according to the result*/
    build_edges_all_Vers(1);

    //fp=fopen("cooling_test.txt", "a");
    //fprintf(fp, "finish building edges for vertex tracing results.\n");
    //fclose(fp);

    _cprintf("finish building edges for vertex tracing results (backward).\n");

    /************************************************************/
    /*   initialize the sample point list  02/09/2010  */
    if (edge_samples != nullptr)
        delete edge_samples;
    edge_samples = new EdgeSamplePt_List(object->elist.nedges/**100*/);
    /************************************************************/

    trace_all_centers_tris_build_edges(tau, 1);

    _cprintf("start adaptive edge sampling (backward).\n");

    trace_all_edges_build_di_edges_adp(tau, 1);

    _cprintf("finish adaptive edge sampling (backward).\n");

    //fp=fopen("cooling_test.txt", "a");
    //fprintf(fp, "finish edge sampling.\n");
    //fclose(fp);

    /***2. perform forward tracing*/
    trace_all_Verts(tau, 0);

    //fopen("cooling_test.txt", "a");
    //fprintf(fp, "finish vertex tracing (forward).\n");
    //fclose(fp);
    _cprintf("finish vertex tracing (forward).\n");

    /*2.2. build the edges according to the result*/
    build_edges_all_Vers(0);

    //fp=fopen("cooling_test.txt", "a");
    //fprintf(fp, "finish building edges for vertex tracing results.\n");
    //fclose(fp);

    trace_all_centers_tris_build_edges(tau, 0);
    trace_all_edges_build_di_edges_adp(tau, 0);

    //fp=fopen("cooling_test.txt", "a");
    //fprintf(fp, "finish edge sampling.\n");
    //fclose(fp);

    /*to visualize the relationship between spatial tau and temporal tau
    08/01/07*/
    ////assign_color();

    ////compute_density();

    ////assign_density_colors();

    ////test_pro(tau, 0.1); /*test Konstantin's ideas*/

    //edge_samples->create_display_list();
}


void MorseDecomp::trace_all_Verts(double tau, int backward)
{
    int i, k;
    int tri_id, end_tris;
    Vertex *v;
    icVector3 stP, newP;

    double lp[2];
    double alpha[3];
    Triangle *face;

    Trajectory *temp= new Trajectory(-1, 1);

    for(i = 0; i < object->vlist.nverts; i++)
    {
        v = object->vlist.verts[i];

        /*   removed by Guoning 07/22/2010 */
        //if(length(v->t_vec) < 1e-10)
        //{
        //	v->imgtri = -1;
        //	continue;   //probably no vector value on it
        //}

        // we current ignore the non-magnifold vertex
        //if (v->ncorners + 1 != v->nedges) continue;

        //_cprintf ("Processing vertex %d\n", i);

        //if (i==222)
        //{
        //	int stop = 0;
        //}

        /* we can choose the triangle that the flow will lead the vertex go into */
        temp->pass_vertex(v->index, tri_id, backward);

        //
        //FILE *fp=fopen("tracing_tris.txt", "a");
        //fprintf(fp, "%d\n", tri_id);
        //fclose(fp);

        if(tri_id < 0)
        {
            tri_id = v->corners[0]->t;
        }

        stP.entry[0] = v->x;
        stP.entry[1] = v->y;
        stP.entry[2] = v->z;
        //face = object->tlist.tris[tri_id];
        //for(k = 0; k < 3; k++)
        //{
        //	if(i == face->verts[k]->index)
        //		break;
        //}

        //alpha[k] = 1-0.00001;
        //alpha[(k+1)%3]=alpha[(k+2)%3]=0.000005;

        //lp[0] = alpha[1]*face->x1+alpha[2]*face->x2;
        //lp[1] = alpha[2]*face->y2;

        //icVector3 glv = lp[0]*face->LX + lp[1]*face->LY;
        //stP.entry[0] = face->verts[0]->x+glv.entry[0];
        //stP.entry[1] = face->verts[0]->y+glv.entry[1];
        //stP.entry[2] = face->verts[0]->z+glv.entry[2];

        //trace_Ver(tri_id, stP.entry, newP.entry, end_tris, tau, backward);

        if (backward == 0)
            trace_Ver_f(tri_id, stP.entry, newP.entry, end_tris, tau);
        else
            trace_Ver_b(tri_id, stP.entry, newP.entry, end_tris, tau);

        if (end_tris > -1)
        {
            int stp = 0;
        }

        v->img_tau[0] = newP.entry[0];
        v->img_tau[1] = newP.entry[1];
        v->img_tau[2] = newP.entry[2];

        v->imgtri = end_tris;

        //if(backward == 0)
        //	v->end_tri[0] = end_tris;
        //else
        //	v->end_tri[1] = end_tris;

        /*save the information for this vertex*/
        //if(backward == 0)
        //	v->tau[0] = trace_time;
        //else
        //	v->tau[1] = trace_time;

    }

    delete temp;
}


void MorseDecomp::trace_Ver(int tri, double st[3], double en[3], int &en_tri, double tau, int backward)
{
    int i;
    int flag = 0;
    double globalp[3];
    int pre_face, cur_face;

    int steps = (int)(tau*20);

    pre_face = cur_face = tri;


    ////if the beginning point is quite close to a vertex of the triangle, move it a little bit away from the vertex

    globalp[0] = st[0];
    globalp[1] = st[1];
    globalp[2] = st[2];

    gt_tau = 0;
    Trajectory *temp = new Trajectory(-2, 1);

    for(i = 0; i < NUMTRACINGTRIANGLE /*&& gt_tau < tau*/; i++)
    {

        if(cur_face < 0)
        {
            en[0] = globalp[0];
            en[1] = globalp[1];
            en[2] = globalp[2];
            en_tri = cur_face;

            delete temp;
            return;
        }

        if (object->tlist.tris[cur_face]->has_zero_vec)
        {
            en[0] = globalp[0];
            en[1] = globalp[1];
            en[2] = globalp[2];
            en_tri = cur_face;

            delete temp;
            return;
        }

        pre_face = cur_face;
        //cur_face = temp->trace_in_triangle_tau(cur_face, globalp, backward, tau, gt_tau, flag);

        if (backward == 0)
            cur_face = temp->trace_in_triangle_tau_f(cur_face, globalp, tau, gt_tau, flag);
        else
            cur_face = temp->trace_in_triangle_tau_b(cur_face, globalp, tau, gt_tau, flag);


        if(gt_tau >= tau || flag == 1 /*|| flag == 3 || flag == 4 */|| pre_face == cur_face )
        {
            en[0] = globalp[0];
            en[1] = globalp[1];
            en[2] = globalp[2];
            en_tri = cur_face;    // 02/11/2010 modified by Guoning
            //en_tri = pre_face; /* the tracing should be stopped at previous triangle */
            delete temp;
            return;
        }
    }

    en[0] = globalp[0];
    en[1] = globalp[1];
    en[2] = globalp[2];
    en_tri = cur_face;

    delete temp;

}



void
MorseDecomp::trace_Ver_f(int tri,
                         double st[3],
                         double en[3],
                         int &en_tri,
                         double tau)
{
    int i;
    int flag = 0;
    double globalp[3];
    int pre_face, cur_face;

    int steps = (int)(tau*20);

    pre_face = cur_face = tri;


    ////if the beginning point is quite close to a vertex of the triangle, move it a little bit away from the vertex

    globalp[0] = st[0];
    globalp[1] = st[1];
    globalp[2] = st[2];

    gt_tau = 0;
    Trajectory *temp = new Trajectory(-2, 1);

    for(i = 0; i < NUMTRACINGTRIANGLE /*&& gt_tau < tau*/; i++)
    {

        if(cur_face < 0)
        {
            en[0] = globalp[0];
            en[1] = globalp[1];
            en[2] = globalp[2];
            en_tri = cur_face;

            delete temp;
            return;
        }

        //if (object->tlist.tris[cur_face]->has_zero_vec)
        //{
        //	en[0] = globalp[0];
        //	en[1] = globalp[1];
        //	en[2] = globalp[2];
        //	en_tri = cur_face;

        //	delete temp;
        //	return;
        //}

        pre_face = cur_face;
        //cur_face = temp->trace_in_triangle_tau(cur_face, globalp, backward, tau, gt_tau, flag);

        cur_face = temp->trace_in_triangle_tau_f(cur_face, globalp, tau, gt_tau, flag);


        if(gt_tau >= tau || flag == 1 /*|| flag == 3 || flag == 4 || pre_face == cur_face*/ )
        {
            en[0] = globalp[0];
            en[1] = globalp[1];
            en[2] = globalp[2];
            en_tri = cur_face;    // 02/11/2010 modified by Guoning
            //en_tri = pre_face; /* the tracing should be stopped at previous triangle */
            delete temp;
            return;
        }
    }

    en[0] = globalp[0];
    en[1] = globalp[1];
    en[2] = globalp[2];
    en_tri = cur_face;

    delete temp;
}


void
MorseDecomp::trace_Ver_b(int tri,
                         double st[3],
                         double en[3],
                         int &en_tri,
                         double tau)
{
    int i;
    int flag = 0;
    double globalp[3];
    int pre_face, cur_face;

    int steps = (int)(tau*20);

    pre_face = cur_face = tri;


    ////if the beginning point is quite close to a vertex of the triangle, move it a little bit away from the vertex

    globalp[0] = st[0];
    globalp[1] = st[1];
    globalp[2] = st[2];

    gt_tau = 0;
    Trajectory *temp = new Trajectory(-2, 1);

    for(i = 0; i < NUMTRACINGTRIANGLE /*&& gt_tau < tau*/; i++)
    {

        if(cur_face < 0)
        {
            en[0] = globalp[0];
            en[1] = globalp[1];
            en[2] = globalp[2];
            en_tri = cur_face;

            delete temp;
            return;
        }

        //if (object->tlist.tris[cur_face]->has_zero_vec)
        //{
        //	en[0] = globalp[0];
        //	en[1] = globalp[1];
        //	en[2] = globalp[2];
        //	en_tri = cur_face;

        //	delete temp;
        //	return;
        //}

        pre_face = cur_face;
        //cur_face = temp->trace_in_triangle_tau(cur_face, globalp, backward, tau, gt_tau, flag);

        cur_face = temp->trace_in_triangle_tau_b(cur_face, globalp, tau, gt_tau, flag);


        if(gt_tau >= tau || flag == 1 /*|| flag == 3 || flag == 4 || pre_face == cur_face*/ )
        {
            en[0] = globalp[0];
            en[1] = globalp[1];
            en[2] = globalp[2];
            en_tri = cur_face;    // 02/11/2010 modified by Guoning
            //en_tri = pre_face; /* the tracing should be stopped at previous triangle */
            delete temp;
            return;
        }
    }

    en[0] = globalp[0];
    en[1] = globalp[1];
    en[2] = globalp[2];
    en_tri = cur_face;

    delete temp;

}


void
MorseDecomp::trace_Ver_local(int tri, double st[3], double en[3], int &en_tri, double tau, int backward, int scc_id)
{
    int i;
    int flag = 0;
    double globalp[3];
    int pre_face, cur_face;

    pre_face = cur_face = tri;


    ////if the beginning point is quite close to a vertex of the triangle, move it a little bit away from the vertex

    globalp[0] = st[0];
    globalp[1] = st[1];
    globalp[2] = st[2];

    gt_tau = 0;
    Trajectory *temp = new Trajectory(-2, 1);

    for(i = 0; i < NUMTRACINGTRIANGLE /*&& gt_tau < tau*/; i++)
    {

        if(cur_face < 0)
        {
            en[0] = globalp[0];
            en[1] = globalp[1];
            en[2] = globalp[2];
            en_tri = cur_face;

            delete temp;
            return;
        }

        //if (object->tlist.tris[cur_face]->has_zero_vec)
        //{
        //	en[0] = globalp[0];
        //	en[1] = globalp[1];
        //	en[2] = globalp[2];
        //	en_tri = cur_face;

        //	delete temp;
        //	return;
        //}

        if(dg->nlist->dirnodes[cur_face]->sscomp_index != scc_id)
        {
            en[0] = globalp[0];
            en[1] = globalp[1];
            en[2] = globalp[2];
            en_tri = cur_face;
            delete temp;
            return;
        }

        pre_face = cur_face;
        //cur_face = temp->trace_in_triangle_tau(cur_face, globalp, backward, tau, gt_tau, flag);
        if (backward == 0)
            cur_face = temp->trace_in_triangle_tau_f(cur_face, globalp, tau, gt_tau, flag);
        else
            cur_face = temp->trace_in_triangle_tau_b(cur_face, globalp, tau, gt_tau, flag);

        if(gt_tau >= tau || flag == 1 || flag == 3 || flag == 4 || pre_face == cur_face )
        {
            en[0] = globalp[0];
            en[1] = globalp[1];
            en[2] = globalp[2];
            en_tri = cur_face;    // 02/11/2010 modified by Guoning
            //en_tri = pre_face; /* the tracing should be stopped at previous triangle */
            delete temp;
            return;
        }
    }

    en[0] = globalp[0];
    en[1] = globalp[1];
    en[2] = globalp[2];
    en_tri = cur_face;

    delete temp;

}

void MorseDecomp::build_edges_all_Vers(int backward)
{
    int i;

    for(i = 0; i < object->vlist.nverts; i++)
    {
        build_edges_Ver(i, backward);
    }
}


void MorseDecomp::build_edges_Ver(int vertid, int backward)
{
    Vertex *v = object->vlist.verts[vertid];

    int i;

    for(i = 0; i < v->ncorners; i++)
    {
        if(v->corners[i]->t < 0 || v->imgtri < 0)
            continue;

        if(backward == 0) /*forward tracing*/
        {
            //if(dg->is_repeated_edge(v->corners[i]->t, v->imgtri))
            if(dg->is_repeated_edge_2(v->corners[i]->t, v->imgtri))
                continue;

            build_one_edge(v->corners[i]->t, v->imgtri);
        }

        else /*backward tracing*/
        {
            //if(dg->is_repeated_edge(v->imgtri, v->corners[i]->t))
            if(dg->is_repeated_edge_2(v->imgtri, v->corners[i]->t))
                continue;

            build_one_edge(v->imgtri, v->corners[i]->t);
        }
    }
}

void MorseDecomp::build_one_edge(int node1, int node2)
{
    dg->add_to_edgelist(node1, node2, dg->elist->nedges);
    ////add the edge to the nodes
    dg->add_edge_to_node(node1, dg->elist->nedges-1);
    dg->add_edge_to_node(node2, dg->elist->nedges-1);
}


void MorseDecomp::trace_center_tri_build_edge(int tri, double tau, int backward)
{
    double center[3] = {0., 0., 0.};

    Triangle *face = object->tlist.tris[tri];

    int i, endtri;

    for(i = 0; i < 3; i++)
    {
        center[0] += face->verts[i]->x;
        center[1] += face->verts[i]->y;
        center[2] += face->verts[i]->z;
    }

    center[0] /= 3.;
    center[1] /= 3.;
    center[2] /= 3.;

    /*start tracing here*/
    //trace_Ver(tri, center, center, endtri, tau, backward);

    if (backward == 0)
        trace_Ver_f(tri, center, center, endtri, tau);
    else
        trace_Ver_b(tri, center, center, endtri, tau);



    if(endtri < 0 || endtri >= object->tlist.ntris)
        return;

    if(backward == 0)  /*forward tracing*/
    {
        //if(!dg->is_repeated_edge(tri, endtri))
        if(!dg->is_repeated_edge_2(tri, endtri))
        {
            build_one_edge(tri, endtri);
        }
    }

    else /*backward tracing*/
    {
        //if(!dg->is_repeated_edge(endtri, tri))
        if(!dg->is_repeated_edge_2(endtri, tri))
        {
            build_one_edge(endtri, tri);
        }
    }
}

void MorseDecomp::trace_all_centers_tris_build_edges(double tau, int backward)
{
    int i;

    for(i = 0; i < object->tlist.ntris; i++)
    {
        /*  Modified by Guoning to avoid the region of zero vectors 03/09/2010  */
        if (object->tlist.tris[i]->has_zero_vec)
            continue;

        trace_center_tri_build_edge(i, tau, backward);
    }
}


//void MorseDecomp::trace_all_edges_build_edges(double tau, int backward)
//{
//}

int current_sample_edge = 0;  // 02/09/2010
bool save_edge_samples = false;
int  investigate_tri = -1;

void MorseDecomp::trace_all_edges_build_di_edges_adp(double tau, int backward)
{
    int i, j;
    Triangle *face;
    Edge *e;
    Vertex *v1, *v2;
    double st1[3], st2[3];
    int endtri;

    init_all_edges();

    for(i = 0; i < object->tlist.ntris; i++)
    {
        face = object->tlist.tris[i];

        //if (i == investigate_tri)
        //	save_edge_samples = true;
        //else
        save_edge_samples = false;

        //save_edge_samples = true;

        /*  Modified by Guoning to avoid the region of zero vectors 03/09/2010  */

        if (face->has_zero_vec)
            continue;

        for(j = 0; j < 3; j++)
        {
            e = face->edges[j];
            if(e->visited )/*we process each edge only once*/
                continue;

            e->visited = true;

            current_sample_edge = e->index;  // 02/09/2010

            v1 = e->verts[0];
            v2 = e->verts[1];

            /*currently, we don't deal with the boundary edges!!!!*/
            if(v1->imgtri < 0 || v2->imgtri < 0)
                continue;

            /*if the images of the two vertices are already continuous, need not
            consider this edge any more*/
            if(v1->imgtri == v2->imgtri
                || are_close_neighbors(v1->imgtri, v2->imgtri))
                continue;

            st1[0] = v1->x;
            st1[1] = v1->y;
            st1[2] = v1->z;

            st2[0] = v2->x;
            st2[1] = v2->y;
            st2[2] = v2->z;

            int neighbor_tri; /*= e->tris[0]->index;
            if(neighbor_tri == i)
                neighbor_tri = e->tris[1]->index*/;

            if(e->tris[0]!=nullptr)
            {
                neighbor_tri = e->tris[0]->index;
                if(neighbor_tri == i && e->tris[1] == nullptr)
                    continue;
                else if(neighbor_tri == i && e->tris[1] !=nullptr)
                    neighbor_tri = e->tris[1]->index;
            }
            else if(e->tris[1] !=nullptr && e->tris[1]->index != i)
                neighbor_tri=e->tris[1]->index;
            else
                continue;


            if (enFastRefinement)
                trace_an_edge_build_di_edges_adp(st1, st2, v1->imgtri, v2->imgtri, i, neighbor_tri,
                                                 tau, backward);

            else
                adp_edge_sampling(st1, st2, v1->imgtri, v2->imgtri, i, neighbor_tri,
                                  tau, backward);

        }
    }
}


void MorseDecomp::init_all_edges()
{
    int i, j;
    Triangle *face;
    Edge *cur_edge;
    for(i = 0; i < object->tlist.ntris; i++)
    {
        face = object->tlist.tris[i];

        for(int j = 0; j < 3; j++)
        {
            cur_edge = face->edges[j];
            cur_edge->OnBoundary = false;
            //cur_edge->BoundaryVisited = false;
            cur_edge->visited = false;
        }
    }
}

bool MorseDecomp::are_close_neighbors(int t1, int t2)
{

    if (t1 < 0 || t2 < 0) return false;

    Triangle *f1 = object->tlist.tris[t1];
    Triangle *f2 = object->tlist.tris[t2];

    int i, j;

    /*if the two triangles share an edge*/
    for(i = 0; i < 3; i++)
    {
        for(j = 0; j < 3; j++)
        {
            if(f1->edges[i] == f2->edges[j])
                return true;
        }
    }

    return false;
}

void MorseDecomp::trace_an_edge_build_di_edges_adp(double st1[3], double st2[3],
                                                   int t1, int t2, int tri,
                                                   int neighbor_tri, double tau, int backward)
{
    /*do one more recursive here*/
    int level = 1;
    double stack_st[3] = {st2[0], st2[1], st2[2]}; /*save vertex v2*/

    st2[0] = (st1[0]+st2[0])/2;
    st2[1] = (st1[1]+st2[1])/2;
    st2[2] = (st1[2]+st2[2])/2;

    double middle_p[3] = {st2[0], st2[1], st2[2]};/*save the middle point*/

    level = 1;

    /*call the recursive adaptive edge sampling for the first half of the edge*/
    trace_recursive_an_edge(st1, st2, t1, tri, neighbor_tri, tau, backward, level);

    /*call the recursive adaptive edge sampling for the second half of the edge*/
    trace_recursive_an_edge(stack_st, middle_p, t2, tri, neighbor_tri, tau, backward, level);
}

/****************************************************************/
/*  02/10/2010
*/
void
MorseDecomp::adp_edge_sampling(double st1[3],
                               double st2[3],
                               int t1, int t2,
                               int tri,
                               int neighbor_tri,
                               double tau,
                               int backward)
{
    if (tri < 0 || neighbor_tri < 0) return;

    /*do one more recursive here*/
    point3 p2;
    p2.p[0] = st2[0]; p2.p[1] = st2[1]; p2.p[2] = st2[2]; /*save vertex v2*/

    //stack
    //stack <point3> stack_pts;
    stack_pts.push(p2);


    // compute the middle point of this edge and push to the stack
    st2[0] = (st1[0]+st2[0])/2;
    st2[1] = (st1[1]+st2[1])/2;
    st2[2] = (st1[2]+st2[2])/2;

    point3 mid_p;
    mid_p.p[0] = st2[0]; mid_p.p[1] = st2[1]; mid_p.p[2] = st2[2]; /*save middle point*/
    stack_pts.push(mid_p);

    // trace from the first point
    double v1_end[3]={0.};
    //trace_Ver(tri, st1, v1_end, t1, tau, backward);

    if(backward == 0)
        trace_Ver_f(tri, st1, v1_end, t1, tau);
    else
        trace_Ver_b(tri, st1, v1_end, t1, tau);

    while(!stack_pts.empty())
    {
        point3 top_p = stack_pts.top();
        stack_pts.pop();

        st2[0] = top_p.p[0];
        st2[1] = top_p.p[1];
        st2[2] = top_p.p[2];

        /*   conduct a tracing from this point  */
        double v2_end[3];
        //trace_Ver(tri, st2, v2_end, t2, tau, backward);

        if (backward == 0)
            trace_Ver_f(tri, st2, v2_end, t2, tau);
        else
            trace_Ver_b(tri, st2, v2_end, t2, tau);

        /************************************************************/
        /*
          Add a sample point here  02/09/2010
        */
        //if (save_edge_samples)
        //{
        //	Edge *cur_e = object->elist.edges[current_sample_edge];
        //	icVector3 len_dir(st2[0]-cur_e->verts[0]->x,
        //		st2[1]-cur_e->verts[0]->y, st2[2]-cur_e->verts[0]->z);
        //	double alpha = length(len_dir)/cur_e->length;
        //	EdgeSamplePt *oneSample = new EdgeSamplePt(current_sample_edge, alpha);
        //	oneSample->end_tri = t2;   // modified at 02/10/2010
        ////	oneSample->backward=backward;
        //	edge_samples->append(oneSample);
        //}

        icVector3 dis;
        dis.entry[0] = st1[0]-st2[0];
        dis.entry[1] = st1[1]-st2[1];
        dis.entry[2] = st1[2]-st2[2];

        if(length(dis) < edge_sample_error)  // try 1.e-15 for fast computation
        //if(length(dis) < 1e-15)  // try 1.e-15 for fast computation
        {
            st1[0] = st2[0];
            st1[1] = st2[1];
            st1[2] = st2[2];
            t1 = t2;
            /*need to build edges here*/
            if (t1 >= 0)  // we consider the planar case here
            {
                if(backward == 0) /*forward tracing*/
                {
                    //if(!dg->is_repeated_edge(tri, t1))
                    if(!dg->is_repeated_edge_2(tri, t1))
                        build_one_edge(tri, t1);
                    //if(neighbor_tri>=0&&!dg->is_repeated_edge(neighbor_tri, t1))
                    if(neighbor_tri>=0&&!dg->is_repeated_edge_2(neighbor_tri, t1))
                        build_one_edge(neighbor_tri, t1);
                }
                else /*backward tracing*/
                {
                    //if(!dg->is_repeated_edge(t1, tri))
                    if(!dg->is_repeated_edge_2(t1, tri))
                        build_one_edge(t1, tri);
                    //if(neighbor_tri>=0&&!dg->is_repeated_edge(t1, neighbor_tri))
                    if(neighbor_tri>=0&&!dg->is_repeated_edge_2(t1, neighbor_tri))
                        build_one_edge(t1, neighbor_tri);
                }
            }
        }

        /*   check connectivity with the image of the previous point  */

        /*   if they intersect, */
        if (t1 == t2 || are_close_neighbors(t1, t2))
        {
            /*need to build edges here*/
            if (t1 >= 0)  // we consider the planar case here
            {
                if(backward == 0) /*forward tracing*/
                {
                    if(!dg->is_repeated_edge_2(tri, t1))
                        build_one_edge(tri, t1);
                    if(neighbor_tri>=0&&!dg->is_repeated_edge_2(neighbor_tri, t1))
                        build_one_edge(neighbor_tri, t1);
                }
                else /*backward tracing*/
                {
                    if(!dg->is_repeated_edge_2(t1, tri))
                        build_one_edge(t1, tri);
                    if(neighbor_tri>=0&&!dg->is_repeated_edge_2(t1, neighbor_tri))
                        build_one_edge(t1, neighbor_tri);
                }
            }

            t1 = t2;
            st1[0] = st2[0];
            st1[1] = st2[1];
            st1[2] = st2[2];
            v1_end[0] = v2_end[0];
            v1_end[1] = v2_end[1];
            v1_end[2] = v2_end[2];
        }

        /*   if they do not intersect, compute the middle point of these two points  */
        else
        {
            mid_p.p[0] = (st1[0]+st2[0])/2;
            mid_p.p[1] = (st1[1]+st2[1])/2;
            mid_p.p[2] = (st1[2]+st2[2])/2;

            stack_pts.push(top_p);
            stack_pts.push(mid_p);
        }
    }

}



void MorseDecomp::trace_recursive_an_edge(double v1[3],
                                          double v2[3],
                                          int &t1,
                                          int tri,
                                          int neighbor_tri,
                                          double tau,
                                          int backward,
                                          int &level)
{
    double stack_v[3];  /*this could only have one element in one level*/
    int top = 0;
    int t2;
    double v2_end[3];

    level++;

    int share_ver;

    /*boundary!!*/
    if(t1 < 0)
        return;

    //trace_Ver(tri, v2, v2_end, t2, tau, backward);
    if (backward == 0)
        trace_Ver_f(tri, v2, v2_end, t2, tau);
    else
        trace_Ver_b(tri, v2, v2_end, t2, tau);

    /************************************************************/
    /*
      Add a sample point here  02/09/2010
    */
    //Edge *cur_e = object->elist.edges[current_sample_edge];
    //icVector3 len_dir(v2[0]-cur_e->verts[0]->x,
    //	v2[1]-cur_e->verts[0]->y, v2[2]-cur_e->verts[0]->z);
    //double alpha = length(len_dir)/cur_e->length;
    //EdgeSamplePt *oneSample = new EdgeSamplePt(current_sample_edge, alpha);
    //oneSample->end_tri = t2;   // modified at 02/10/2010
    //edge_samples->append(oneSample);

    ///////////////////////////////////////////////////////////////

    /*boundary!!*/
    if(t2 < 0)
        return;

    icVector3 dis;
    dis.entry[0] = v1[0]-v2[0];
    dis.entry[1] = v1[1]-v2[1];
    dis.entry[2] = v1[2]-v2[2];

    if(length(dis) < edge_sample_error/*1e-15*/)
    {
        v1[0] = v2[0];
        v1[1] = v2[1];
        v1[2] = v2[2];
        t1 = t2;
        /*need to build edges here*/
        if(backward == 0) /*forward tracing*/
        {
            //if(!dg->is_repeated_edge(tri, t1))
            if(!dg->is_repeated_edge_2(tri, t1))
                build_one_edge(tri, t1);
            //if(neighbor_tri>=0&&!dg->is_repeated_edge(neighbor_tri, t1))
            if(neighbor_tri>=0&&!dg->is_repeated_edge_2(neighbor_tri, t1))
                build_one_edge(neighbor_tri, t1);
        }
        else /*backward tracing*/
        {
            //if(!dg->is_repeated_edge(t1, tri))
            if(!dg->is_repeated_edge_2(t1, tri))
                build_one_edge(t1, tri);
            //if(neighbor_tri>=0&&!dg->is_repeated_edge(t1, neighbor_tri))
            if(neighbor_tri>=0&&!dg->is_repeated_edge_2(t1, neighbor_tri))
                build_one_edge(t1, neighbor_tri);
        }
        return;
    }


    /*if the two points are too close, stop!*/  /*  DO NOT set this limit now 2/9/2010  */
    if(level > 10)
    {
        /*trace v2, get t2 */
        v1[0] = v2[0];
        v1[1] = v2[1];
        v1[2] = v2[2];
        t1 = t2;
        /*need to build edges here*/
        if(backward == 0) /*forward tracing*/
        {
            //if(!dg->is_repeated_edge(tri, t1))
            if(!dg->is_repeated_edge_2(tri, t1))
                build_one_edge(tri, t1);
            //if(neighbor_tri>=0&&!dg->is_repeated_edge(neighbor_tri, t1))
            if(neighbor_tri>=0&&!dg->is_repeated_edge_2(neighbor_tri, t1))
                build_one_edge(neighbor_tri, t1);
        }
        else /*backward tracing*/
        {
            //if(!dg->is_repeated_edge(t1, tri))
            if(!dg->is_repeated_edge_2(t1, tri))
                build_one_edge(t1, tri);
            //if(neighbor_tri>=0&&!dg->is_repeated_edge(t1, neighbor_tri))
            if(neighbor_tri>=0&&!dg->is_repeated_edge_2(t1, neighbor_tri))
                build_one_edge(t1, neighbor_tri);
        }
        return;
    }

    if(t1 == t2 || are_close_neighbors(t1, t2))
    {
        v1[0] = v2[0];
        v1[1] = v2[1];
        v1[2] = v2[2];
        t1 = t2;

        /*need to build edges here*/
        if(backward == 0) /*forward tracing*/
        {
            //if(!dg->is_repeated_edge(tri, t1))
            if(!dg->is_repeated_edge_2(tri, t1))
                build_one_edge(tri, t1);
            //if(neighbor_tri>=0&&!dg->is_repeated_edge(neighbor_tri, t1))
            if(neighbor_tri>=0&&!dg->is_repeated_edge_2(neighbor_tri, t1))
                build_one_edge(neighbor_tri, t1);
        }
        else /*backward tracing*/
        {
            //if(!dg->is_repeated_edge(t1, tri))
            if(!dg->is_repeated_edge_2(t1, tri))
                build_one_edge(t1, tri);
            //if(neighbor_tri>=0&&!dg->is_repeated_edge(t1, neighbor_tri))
            if(neighbor_tri>=0&&!dg->is_repeated_edge_2(t1, neighbor_tri))
                build_one_edge(t1, neighbor_tri);
        }
        return;
    }

    /*consider further approximation of the closure of the image 07/24/07*/
    //else if(are_pseudo_close_neighbors(t1, t2, share_ver))
    //{
    //	v1[0] = v2[0];
    //	v1[1] = v2[1];
    //	v1[2] = v2[2];
    //	t1 = t2;

    //	/*add the one-ring neighborhood of the share_ver to the edge list*/
    //	Vertex* v = Object.vlist[share_ver];
    //	Corner *c;
    //	for(int i=0; i<v->Num_corners; i++)
    //	{
    //		c = Object.clist[v->Corners[i]];
    //		/*need to build edges here*/
    //		if(backward == 0) /*forward tracing*/
    //		{
    //			if(!has_Edge_From_To(tri, c->t))
    //				build_one_edge(tri, c->t);
    //			if(neighbor_tri>=0&&!has_Edge_From_To(neighbor_tri, c->t))
    //				build_one_edge(neighbor_tri, c->t);
    //		}
    //		else /*backward tracing*/
    //		{
    //			if(!has_Edge_From_To(c->t, tri))
    //				build_one_edge(c->t, tri);
    //			if(neighbor_tri>=0&&!has_Edge_From_To(c->t, neighbor_tri))
    //				build_one_edge(c->t, neighbor_tri);
    //		}
    //	}
    //	return;
    //}

    else
    {
        /*push v2 into the stack and indicate that the stack is not empty any more*/
        stack_v[0] = v2[0];
        stack_v[1] = v2[1];
        stack_v[2] = v2[2];
        top++;

        /*use the middle point of this line segment as current v2*/
        v2[0] = (v1[0]+v2[0])/2.;
        v2[1] = (v1[1]+v2[1])/2.;
        v2[2] = (v1[2]+v2[2])/2.;

        /*recursively call the routine itself*/
        trace_recursive_an_edge(v1, v2, t1, tri, neighbor_tri, tau, backward, level);
        level --;
    }

    /*if the stack is not empty*/
    if(top > 0)
    {
        /*pop up the point*/
        v2[0] = stack_v[0];
        v2[1] = stack_v[1];
        v2[2] = stack_v[2];

        top--;
        /*recursively call the routine itself*/
        trace_recursive_an_edge(v1, v2, t1, tri, neighbor_tri, tau, backward, level);
        level--;
    }
}

void
MorseDecomp::show_tri_mapping (int tri)
{
    if (tri<0 || tri>=object->tlist.ntris) return;

    if (dg == nullptr) return;
    Triangle *t = object->tlist.tris[tri];

    DirGraph_Node *cur_n = dg->nlist->dirnodes[tri];

    int i;

    for (i=0; i<cur_n->nedges; i++)
    {
        Graph_Edge *cur_e = dg->elist->edges[cur_n->edges[i]];

        int adj_n = cur_e->node_index1;
        glColor3f (0, 1, 0.2);

        if (adj_n == tri)
        {
            adj_n = cur_e->node_index2;
            glColor3f (1, 0, 0.3);
        }

        /*
            Display this triangle
        */
        Triangle *adj_t = object->tlist.tris[adj_n];
        glBegin(GL_TRIANGLES);
        glVertex3f(adj_t->verts[0]->x, adj_t->verts[0]->y, adj_t->verts[0]->z);
        glVertex3f(adj_t->verts[1]->x, adj_t->verts[1]->y, adj_t->verts[1]->z);
        glVertex3f(adj_t->verts[2]->x, adj_t->verts[2]->y, adj_t->verts[2]->z);
        glEnd();
    }
}

bool
is_connected_reg(int *tris, int ntris)
{
    int i, j;

    for (i=0; i<ntris; i++)
        object->tlist.tris[tris[i]]->visited = false;

    for (i=0; i<ntris; i++)
        object->tlist.tris[tris[i]]->in_img = true;


    DynList_Int *c = new DynList_Int(ntris);
    c->add_New (tris[0]);
    object->tlist.tris[tris[0]]->visited = true;

    for (j=0; j<c->nelems; j++)
    {
        Triangle *t = object->tlist.tris[c->elems[j]];

        for (int k=0; k<3; k++)
        {
            Vertex *v = t->verts[k];

            for (int l=0; l<v->ncorners; l++)
            {
                if (v->tris[l]->in_img && !v->tris[l]->visited)
                {
                    c->add_New (v->tris[l]->index);
                    v->tris[l]->visited = true;
                }
            }
        }
    }

    int nelems_c = c->nelems;

    delete c;

    for (i=0; i<ntris; i++)
    {
        object->tlist.tris[tris[i]]->in_img = false;
        object->tlist.tris[tris[i]]->visited= false;
    }

    if (nelems_c == ntris) return true;
    return false;
}

bool
MorseDecomp::obtain_connected_imgs(int tri)     // judge whether the image and preimage of a given triangle is connected
{
    if (tri<0 || tri>=object->tlist.ntris) return false;

    if (dg == nullptr) return false;
    Triangle *t = object->tlist.tris[tri];

    DirGraph_Node *cur_n = dg->nlist->dirnodes[tri];

    DynList_Int *pre_img = new DynList_Int(cur_n->nedges);
    DynList_Int *img = new DynList_Int(cur_n->nedges);

    int i;

    for (i=0; i<cur_n->nedges; i++)
    {
        Graph_Edge *cur_e = dg->elist->edges[cur_n->edges[i]];

        int adj_n = cur_e->node_index1;

        if (adj_n == tri)
        {
            adj_n = cur_e->node_index2;
            img->add_New(adj_n);         //image
        }
        else
        {
            pre_img->add_New(adj_n);     //pre-image
        }
    }

    /*
       determine the connectivity of these two images
    */
    bool con_img = is_connected_reg (img->elems, img->nelems);
    bool con_pre_img = is_connected_reg(pre_img->elems, pre_img->nelems);

    delete pre_img;
    delete img;

    if (con_img && con_pre_img)
        return true;
    return false;
}

/*
This function will check whether the local updated graph will provide connected image for each triangle inside the refine Morse set or not
   IN: reg_tris -- the list of the indices of the triangles inside the refined Morse set
       ntris    -- the number of triangles in the region
*/
bool
MorseDecomp::obtain_connected_imgs_reg(int *reg_tris, int ntris)
{
    int i;

    for (i=0; i<object->tlist.ntris; i++)
    {
        object->tlist.tris[i]->in_img = false;
        object->tlist.tris[i]->visited = false;
    }

    for (i=0; i<ntris; i++)
    {
        if (!obtain_connected_imgs(reg_tris[i]))
            return false;
    }

    return true;
}


/********************************************************
02/28/2010
 Save the computation results of Morse decompositions for better visualization
*/

void
MorseDecomp::save_dirGraph_FG(char *filename)
{
    if (dg == nullptr) return;

    FILE *fp = fopen(filename, "w");  // write in the ASCII format

    int i;

    fprintf(fp, "#nodes %d\n", dg->nlist->ndirnodes);
    fprintf(fp, "#edges %d\n", dg->elist->nedges);

    /*
       We only need to store the edge list
    */
    for (i=0; i<dg->elist->nedges; i++)
    {
        fprintf(fp, "%d %d\n", dg->elist->edges[i]->node_index1, dg->elist->edges[i]->node_index2);
    }

    fclose(fp);
}


void
MorseDecomp::load_dirGraph_FG(char *filename)
{
    FILE *fp = fopen(filename, "r");

    if (fp == nullptr) return;

    int nodes = 0;

    fscanf(fp, "#nodes %d\n", &nodes);

    if (nodes != object->tlist.ntris) return;  // not match

    int nedges = 0;
    fscanf(fp, "#edges %d\n", &nedges);

    /*
       allocate memory for the directed graph
    */

    if (dg != nullptr) delete dg;

    dg = new DirGraph(nodes, nedges);

    int i;

    for (i=0; i<nodes; i++)
    {
        dg->nlist->dirnodes[i] = new DirGraph_Node();
        dg->nlist->dirnodes[i]->edges = nullptr;
        dg->nlist->dirnodes[i]->nedges = 0;
        dg->nlist->dirnodes[i]->node_index = i;
    }

    dg->nlist->ndirnodes = nodes;

    for(i = 0; i < nedges; i++)
    {
        dg->elist->edges[i] = new Graph_Edge();
    }

    for (i=0; i<nedges; i++)
    {
        int node1, node2;

        fscanf(fp, "%d %d\n", &node1, &node2);

        //add a new edge

        build_one_edge(node1, node2);
    }

    fclose(fp);
}

/*************************************/
void
MorseDecomp::out_sub_FG(char *filename)
{
    if (dg == nullptr) return;

    FILE *fp = fopen(filename, "w");  // write in the ASCII format

    int i;

    //fprintf(fp, "#nodes %d\n", dg->nlist->ndirnodes);
    //fprintf(fp, "#edges %d\n", dg->elist->nedges);

    int num_nodes = 0;
    int num_edges = 0;
    for (i=0; i<object->tlist.ntris; i++)
    {
        if (!object->tlist.tris[i]->exclude)
            continue;

        num_nodes ++;

        DirGraph_Node *node = dg->nlist->dirnodes[i];

        for (int j=0; j<node->nedges; j++)
        {
            Graph_Edge *edge=dg->elist->edges[node->edges[j]];

            int adj_n = edge->node_index1;

            if (adj_n==i)
                adj_n = edge->node_index2;

            if (object->tlist.tris[adj_n]->exclude
                && object->tlist.tris[i]->local_index == object->tlist.tris[adj_n]->local_index)
                num_edges ++;
        }
    }

    fprintf(fp, "#nodes %d\n", num_nodes);
    fprintf(fp, "#edges %d\n", num_edges);

    /*
       We need to store the both node and edge lists
    */
    for (i=0; i<object->tlist.ntris; i++)
    {
        if (!object->tlist.tris[i]->exclude)
            continue;
        fprintf(fp, "%d\n", i);
    }

    for (i=0; i<object->tlist.ntris; i++)
    {
        if (!object->tlist.tris[i]->exclude)
            continue;

        DirGraph_Node *node = dg->nlist->dirnodes[i];

        for (int j=0; j<node->nedges; j++)
        {
            Graph_Edge *edge=dg->elist->edges[node->edges[j]];

            int adj_n = edge->node_index1;

            if (adj_n==i)
                adj_n = edge->node_index2;

            if (object->tlist.tris[adj_n]->exclude
                && object->tlist.tris[i]->local_index == object->tlist.tris[adj_n]->local_index)
            {
                fprintf(fp, "%d %d\n", edge->node_index1, edge->node_index2);
            }
        }
    }
}


void
MorseDecomp::load_sub_FG(char *filename)
{
    FILE *fp = fopen(filename, "r");

    if (fp == nullptr) return;

    int nodes = 0;

    fscanf(fp, "#nodes %d\n", &nodes);

    int nedges = 0;
    fscanf(fp, "#edges %d\n", &nedges);

    /*
       allocate memory for the directed graph
    */

    if (dg != nullptr) delete dg;

    dg = new DirGraph(nodes, nedges);

    int i;

    for (i=0; i<nodes; i++)
    {
        int node_id=-1;
        fscanf(fp, "%d\n", &node_id);
        dg->nlist->dirnodes[i] = new DirGraph_Node();
        dg->nlist->dirnodes[i]->edges = nullptr;
        dg->nlist->dirnodes[i]->nedges = 0;
        dg->nlist->dirnodes[i]->node_index = node_id;
    }

    dg->nlist->ndirnodes = nodes;

    for(i = 0; i < nedges; i++)
    {
        dg->elist->edges[i] = new Graph_Edge();
    }

    for (i=0; i<nedges; i++)
    {
        int node1, node2;

        fscanf(fp, "%d %d\n", &node1, &node2);

        //add a new edge

        // search the corresponding node index in the node list

        int node1_id, node2_id;
        for (int j=0; j<nodes; j++)
        {
            if (dg->nlist->dirnodes[j]->node_index == node1)
                node1_id = j;

            if (dg->nlist->dirnodes[j]->node_index == node2)
                node2_id = j;
        }

        build_one_edge(node1_id, node2_id);
    }

    fclose(fp);
}

//////////////////////////////////////////////////////////////////////


void
MorseDecomp::get_t_max_min()
{
    int i;
    tau_max = -1;
    tau_min = 100;

    for (i=0; i<object->tlist.ntris; i++)
    {
        if (object->tlist.tris[i]->used_tau > tau_max)
            tau_max = object->tlist.tris[i]->used_tau;
        if (object->tlist.tris[i]->used_tau < tau_min)
            tau_min = object->tlist.tris[i]->used_tau;
    }
}


void
MorseDecomp::compare_with_ECG ()
{
    /*  We first mark all the triangles that are included in either Morse sets or connection regions  */
    int i, j;

    /*
       First, initialize all the triangle flags
    */
    //for (i=0; i<object->tlist.ntris; i++)
    //{
    //	object->tlist.tris[i]->is_part_of_MCG = false;
    //	object->tlist.tris[i]->diff_ECG_MCG = false;
    //}

    //for (i=0; i<mcg->nlist->nmnodes; i++)
    //{
    //	int scc_id = mcg->nlist->mnodes[i]->scc_index;
    //	SCComponent *scc_comp = scclist->scccomponents[scc_id];
    //	for (j=0; j<scc_comp->nnodes; j++)
    //	{
    //		object->tlist.tris[scc_comp->nodes[j]]->is_part_of_MCG = true;
    //	}
    //}

    //for (i=0; i<mcg->elist->nedges; i++)
    //{
    //	Graph_Edge *cur_e = mcg->elist->edges[i];

    //	for (j=0; j<cur_e->ntris; j++)
    //	{
    //		object->tlist.tris[cur_e->triangles[j]]->is_part_of_MCG = true;
    //	}
    //}

    //for (i=0; i<object->tlist.ntris; i++)
    //	object->tlist.tris[i]->counter ++;

    mark_cur_MorseSets();
    for (i=0; i<mcg->elist->nedges; i++)
    {
        Graph_Edge *cur_e = mcg->elist->edges[i];

        for (j=0; j<cur_e->ntris; j++)
        {
            object->tlist.tris[cur_e->triangles[j]]->is_part_of_MCG = true;
        }
    }

    /*
       Do the ECG computation and mark those triangles that are not covered by the ECG
    */
    //object->capture_Singularities();
    //object->detect_PeriodicOrbit();
    //object->cal_separatrices();

    /*   Morse sets   */
    for (i=0; i<object->slist.nsingularities; i++)
    {
        if (!object->tlist.tris[object->slist.slist[i]->TriangleID]->is_part_of_MCG)
            object->tlist.tris[object->slist.slist[i]->TriangleID]->diff_ECG_MCG = true;
    }

    for (i=0; i<ndisplay_POs; i++)
    {
        for (j=0; j<periodic_orbits->polist[i]->traj->nlinesegs; j++)
        {
            LineSeg *cur_l = &periodic_orbits->polist[i]->traj->linesegs[j];

            if (!object->tlist.tris[cur_l->Triangle_ID]->is_part_of_MCG)
                object->tlist.tris[cur_l->Triangle_ID]->diff_ECG_MCG = true;
        }
    }

    /*  Connection regions  */

    //for (i=0; i
}




void
MorseDecomp::mark_cur_MorseSets()
{
    int i, j;

    /*
       First, initialize all the triangle flags
    */
    for (i=0; i<object->tlist.ntris; i++)
    {
        object->tlist.tris[i]->is_part_of_MCG = false;
        object->tlist.tris[i]->diff_ECG_MCG = false;
    }

    for (i=0; i<mcg->nlist->nmnodes; i++)
    {
        int scc_id = mcg->nlist->mnodes[i]->scc_index;
        SCComponent *scc_comp = scclist->scccomponents[scc_id];
        for (j=0; j<scc_comp->nnodes; j++)
        {
            object->tlist.tris[scc_comp->nodes[j]]->is_part_of_MCG = true;
        }
    }


    for (i=0; i<object->tlist.ntris; i++)
    {
        if (object->tlist.tris[i]->is_part_of_MCG)
            object->tlist.tris[i]->counter ++;
    }
}









/**********************************************************************************/

/////////////////////////////////////////////////////////////////////////////////////
/*   Implementation of the functions in class Trajectory for Morse decomposition   */
/////////////////////////////////////////////////////////////////////////////////////

int Trajectory::trace_in_triangle_tau(int &face_id, double globalp[3],
                                      int type, double tau, double &gt_tau, int &flag)
{
    int i;
    double alpha[3];
    double cur_point[2], pre_point[2];
    double vert0[3];
    icVector3 VP, globalv;

    Triangle *face = object->tlist.tris[face_id];

    /************************************************************/
    /*
         Trick: Do NOT trace if the triangle contains a fixed points (02/11/2010)
    */
    //if (face->singularityID>=0)
    //{
    //	flag = 1;
    //	return face_id;
    //}

    Triangle *pre_f = face;

    icVector2 dis;

    ////Temporary curve point array

    ////initialize
    VP.entry[0] = globalp[0] - face->verts[0]->x;
    VP.entry[1] = globalp[1] - face->verts[0]->y;
    VP.entry[2] = globalp[2] - face->verts[0]->z;

    pre_point[0] = cur_point[0] = dot(VP, face->LX);
    pre_point[1] = cur_point[1] = dot(VP, face->LY);

    vert0[0] = face->verts[0]->x;   ////for update the global point
    vert0[1] = face->verts[0]->y;
    vert0[2] = face->verts[0]->z;

    //globalface = face_id;

    ////////////////////////////////////////////////////
    for(i = 0; i < 100; i++)
    {
        ////1. calculate the barycentric coordinates for current point

        object->get_2D_Barycentric_Facters(face_id, cur_point[0], cur_point[1], alpha);

        ////2. if current point is inside current triangle
        if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1
            && (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1
            && (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
        {
            ////store the point into the temp curve points list

            pre_point[0] = cur_point[0];
            pre_point[1] = cur_point[1];

            //if(cal_next_point_euler1_tau(pre_point, cur_point, face_id, alpha, type))
            //if(cal_next_point_euler1_tau_nonNorm(pre_point, cur_point, face_id, alpha, type))
            //if(cal_next_point_2ndEuler_tau_nonNorm(pre_point, cur_point, face_id, alpha, type))
            //if(cal_next_point_2ndEuler_tau_Norm(pre_point, cur_point, face_id, alpha, type))
            //if(cal_next_pt_RK4_tau(pre_point, cur_point, face_id, alpha, type))
            //if (get_next_pt_tau(pre_point, cur_point, face_id, alpha, type, unsigned char(Integrator_opt)))

            bool trace_result=false;

            if (type == 0)
                trace_result = get_next_pt_tau_f(pre_point, cur_point, face_id, alpha, (unsigned char)(Integrator_opt));
            else
                trace_result = get_next_pt_tau_b(pre_point, cur_point, face_id, alpha, unsigned char(Integrator_opt));

            if (trace_result)
            {
                ////update the global point

                globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

                //Get global coordinates of the point
                globalp[0] = vert0[0] + globalv.entry[0];
                globalp[1] = vert0[1] + globalv.entry[1];
                globalp[2] = vert0[2] + globalv.entry[2];

                /*This is more accurate method to accumlate time*/
                dis.entry[0] = cur_point[0] - pre_point[0];
                dis.entry[1] = cur_point[1] - pre_point[1];

                /*temporal tau*/
                gt_tau+=length(dis)/cur_vec_mag;
                //////////////g_dt += len*len/fabs(dot(dis, glob_loc_v));

                /////////////*spatial tau*/
                //////////////double len = length(dis);
                //////////////g_dt += len;

                /*   now the tau increment should be constant! 2/8/2010 */
                //gt_tau += 0.01;

                if(gt_tau >= tau)
                {
                    flag = 1;
                    return face_id;
                }
            }

            else{  ////the curve reach a singularity
                flag = 1;

                ////Store the record into global line segment array

                return face_id;
            }
        }

        ////3. if the point is out of current triangle
        else{
            double t[2] = {0.};

            int PassVertornot = 0;

            get_next_triangle(face_id, pre_point, cur_point, t, type, PassVertornot, alpha);

            ////update the globalpoint here (Modified on 01/30/07)
            if(PassVertornot > 0)
            {
                if (face_id < 0)
                    return -1;

                //we first need to know which vertex it is in the new triangle 01/30/07
                Vertex* vertid = pre_f->verts[PassVertornot-1];
                Triangle *cur_f = object->tlist.tris[face_id];
                int vert_new = 0;
                for(int k = 0; k < 3; k++)
                {
                    if(cur_f->verts[k] == vertid)
                    {
                        vert_new = k;
                        break;
                    }
                }

                alpha[vert_new]=1-0.0001;
                alpha[(vert_new+1)%3]=0.00005;
                alpha[(vert_new+2)%3]=0.00005;


                /* Get the new cur_point */
                cur_point[0] = alpha[1]*cur_f->x1+alpha[2]*cur_f->x2;
                cur_point[1] = alpha[2]*cur_f->y2;

                globalv = cur_point[0] * cur_f->LX + cur_point[1] * cur_f->LY;

                globalp[0] = cur_f->verts[0]->x + globalv.entry[0];
                globalp[1] = cur_f->verts[0]->y + globalv.entry[1];
                globalp[2] = cur_f->verts[0]->z + globalv.entry[2];
                face=cur_f;
            }

            else{
                //// transfer it to the global coordinates
                globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;

                globalp[0] = vert0[0] + globalv.entry[0];
                globalp[1] = vert0[1] + globalv.entry[1];
                globalp[2] = vert0[2] + globalv.entry[2];
            }

            dis.entry[0] = cur_point[0] - pre_point[0];
            dis.entry[1] = cur_point[1] - pre_point[1];

            /*temporal tau*/
            double temp_ratio=length(dis)/cur_vec_mag;
            gt_tau+=length(dis)/cur_vec_mag;  // This computation is problematic  02/08/2010
            ////////g_dt += len*len/fabs(dot(dis, glob_loc_v));

            ///////*spatial tau*/
            ////////double len = length(dis);
            ////////g_dt += len;

            //gt_tau += 0.01;

            if (gt_tau > tau)                 // Added by Guoning at 02/08/2010
            {
                flag = 1;
                //globalv = pre_point[0] * pre_f->LX + pre_point[1] * pre_f->LY;
                //globalp[0] = pre_f->verts[0]->x + globalv.entry[0];
                //globalp[1] = pre_f->verts[0]->y + globalv.entry[1];
                //globalp[2] = pre_f->verts[0]->z + globalv.entry[2];

                //gt_tau-=temp_ratio;
                return pre_f->index;
            }

            return face_id;
        }

    }

    return face_id;
}



int Trajectory::trace_in_triangle_tau_f(int &face_id, double globalp[3],
                                        double tau, double &gt_tau, int &flag)
{
    int i;
    double alpha[3];
    double cur_point[2], pre_point[2];
    double vert0[3];
    icVector3 VP, globalv;

    Triangle *face = object->tlist.tris[face_id];

    //eulerstep_scalar = face->euler_step;

    /************************************************************/
    /*
         Trick: Do NOT trace if the triangle contains a fixed points (02/11/2010)
    */
    //if (face->singularityID>=0)
    //{
    //	flag = 1;
    //	return face_id;
    //}

    Triangle *pre_f = face;

    icVector2 dis;

    ////Temporary curve point array

    ////initialize
    VP.entry[0] = globalp[0] - face->verts[0]->x;
    VP.entry[1] = globalp[1] - face->verts[0]->y;
    VP.entry[2] = globalp[2] - face->verts[0]->z;

    pre_point[0] = cur_point[0] = dot(VP, face->LX);
    pre_point[1] = cur_point[1] = dot(VP, face->LY);

    vert0[0] = face->verts[0]->x;   ////for update the global point
    vert0[1] = face->verts[0]->y;
    vert0[2] = face->verts[0]->z;

    //globalface = face_id;

    ////////////////////////////////////////////////////
    for(i = 0; i < 100; i++)
    {
        ////1. calculate the barycentric coordinates for current point

        object->get_2D_Barycentric_Facters(face_id, cur_point[0], cur_point[1], alpha);

        ////2. if current point is inside current triangle
        if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1
            && (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1
            && (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
        {
            ////store the point into the temp curve points list

            pre_point[0] = cur_point[0];
            pre_point[1] = cur_point[1];

            //if(cal_next_point_euler1_tau(pre_point, cur_point, face_id, alpha, type))
            //if(cal_next_point_euler1_tau_nonNorm(pre_point, cur_point, face_id, alpha, type))
            //if(cal_next_point_2ndEuler_tau_nonNorm(pre_point, cur_point, face_id, alpha, type))
            //if(cal_next_point_2ndEuler_tau_Norm(pre_point, cur_point, face_id, alpha, type))
            //if(cal_next_pt_RK4_tau(pre_point, cur_point, face_id, alpha, type))
            //if (get_next_pt_tau(pre_point, cur_point, face_id, alpha, type, unsigned char(Integrator_opt)))


            if (get_next_pt_tau_f(pre_point, cur_point, face_id, alpha, unsigned char(Integrator_opt)))
            {
                ////update the global point

                globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

                //Get global coordinates of the point
                globalp[0] = vert0[0] + globalv.entry[0];
                globalp[1] = vert0[1] + globalv.entry[1];
                globalp[2] = vert0[2] + globalv.entry[2];

                /*This is more accurate method to accumlate time*/
                //dis.entry[0] = cur_point[0] - pre_point[0];
                //dis.entry[1] = cur_point[1] - pre_point[1];

                /*temporal tau*/
                gt_tau+=/*length(dis)*/move_dist/cur_vec_mag;
                //////////////g_dt += len*len/fabs(dot(dis, glob_loc_v));

                /////////////*spatial tau*/
                //////////////double len = length(dis);
                //////////////g_dt += len;

                /*   now the tau increment should be constant! 2/8/2010 */
                //gt_tau += 0.01;

                if(gt_tau >= tau)
                {
                    flag = 1;
                    return face_id;
                }
            }

            else{  ////the curve reach a singularity
                flag = 1;

                ////Store the record into global line segment array

                return face_id;
            }
        }

        ////3. if the point is out of current triangle
        else{
            double t[2] = {0.};

            int PassVertornot = 0;

            get_next_triangle(face_id, pre_point, cur_point, t, 0, PassVertornot, alpha);

            ////update the globalpoint here (Modified on 01/30/07)
            if(PassVertornot > 0)
            {

                if (face_id < 0)
                    return -1;

                //we first need to know which vertex it is in the new triangle 01/30/07
                Vertex* vertid = pre_f->verts[PassVertornot-1];
                Triangle *cur_f = object->tlist.tris[face_id];
                int vert_new = 0;
                for(int k = 0; k < 3; k++)
                {
                    if(cur_f->verts[k] == vertid)
                    {
                        vert_new = k;
                        break;
                    }
                }

                alpha[vert_new]=1-0.0001;
                alpha[(vert_new+1)%3]=0.00005;
                alpha[(vert_new+2)%3]=0.00005;


                /* Get the new cur_point */
                cur_point[0] = alpha[1]*cur_f->x1+alpha[2]*cur_f->x2;
                cur_point[1] = alpha[2]*cur_f->y2;

                globalv = cur_point[0] * cur_f->LX + cur_point[1] * cur_f->LY;

                globalp[0] = cur_f->verts[0]->x + globalv.entry[0];
                globalp[1] = cur_f->verts[0]->y + globalv.entry[1];
                globalp[2] = cur_f->verts[0]->z + globalv.entry[2];
                face=cur_f;
            }

            else{
                //// transfer it to the global coordinates
                globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;

                globalp[0] = vert0[0] + globalv.entry[0];
                globalp[1] = vert0[1] + globalv.entry[1];
                globalp[2] = vert0[2] + globalv.entry[2];
            }

            dis.entry[0] = cur_point[0] - pre_point[0];
            dis.entry[1] = cur_point[1] - pre_point[1];

            /*temporal tau*/
            //double temp_ratio=length(dis)/cur_vec_mag;
            gt_tau+=length(dis)/cur_vec_mag;  // This computation is problematic  02/08/2010
            ////////g_dt += len*len/fabs(dot(dis, glob_loc_v));

            ///////*spatial tau*/
            ////////double len = length(dis);
            ////////g_dt += len;

            //gt_tau += 0.01;

            if (gt_tau > tau)                 // Added by Guoning at 02/08/2010
            {
                flag = 1;
                //globalv = pre_point[0] * pre_f->LX + pre_point[1] * pre_f->LY;
                //globalp[0] = pre_f->verts[0]->x + globalv.entry[0];
                //globalp[1] = pre_f->verts[0]->y + globalv.entry[1];
                //globalp[2] = pre_f->verts[0]->z + globalv.entry[2];

                //gt_tau-=temp_ratio;
                return pre_f->index;
            }

            return face_id;
        }

    }

    return face_id;
}


int Trajectory::trace_in_triangle_tau_b(int &face_id, double globalp[3],
                                        double tau, double &gt_tau, int &flag)
{
    int i;
    double alpha[3];
    double cur_point[2], pre_point[2];
    double vert0[3];
    icVector3 VP, globalv;

    Triangle *face = object->tlist.tris[face_id];

    //eulerstep_scalar = face->euler_step;

    /************************************************************/
    /*
         Trick: Do NOT trace if the triangle contains a fixed points (02/11/2010)
    */
    //if (face->singularityID>=0)
    //{
    //	flag = 1;
    //	return face_id;
    //}

    Triangle *pre_f = face;

    icVector2 dis;

    ////Temporary curve point array

    ////initialize
    VP.entry[0] = globalp[0] - face->verts[0]->x;
    VP.entry[1] = globalp[1] - face->verts[0]->y;
    VP.entry[2] = globalp[2] - face->verts[0]->z;

    pre_point[0] = cur_point[0] = dot(VP, face->LX);
    pre_point[1] = cur_point[1] = dot(VP, face->LY);

    vert0[0] = face->verts[0]->x;   ////for update the global point
    vert0[1] = face->verts[0]->y;
    vert0[2] = face->verts[0]->z;

    //globalface = face_id;

    ////////////////////////////////////////////////////
    for(i = 0; i < 100; i++)
    {
        ////1. calculate the barycentric coordinates for current point

        object->get_2D_Barycentric_Facters(face_id, cur_point[0], cur_point[1], alpha);

        ////2. if current point is inside current triangle
        if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1
            && (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1
            && (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
        {
            ////store the point into the temp curve points list

            pre_point[0] = cur_point[0];
            pre_point[1] = cur_point[1];

            //if(cal_next_point_euler1_tau(pre_point, cur_point, face_id, alpha, type))
            //if(cal_next_point_euler1_tau_nonNorm(pre_point, cur_point, face_id, alpha, type))
            //if(cal_next_point_2ndEuler_tau_nonNorm(pre_point, cur_point, face_id, alpha, type))
            //if(cal_next_point_2ndEuler_tau_Norm(pre_point, cur_point, face_id, alpha, type))
            //if(cal_next_pt_RK4_tau(pre_point, cur_point, face_id, alpha, type))
            //if (get_next_pt_tau(pre_point, cur_point, face_id, alpha, type, unsigned char(Integrator_opt)))

            if (get_next_pt_tau_b(pre_point, cur_point, face_id, alpha, unsigned char(Integrator_opt)))
            {
                ////update the global point

                globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

                //Get global coordinates of the point
                globalp[0] = vert0[0] + globalv.entry[0];
                globalp[1] = vert0[1] + globalv.entry[1];
                globalp[2] = vert0[2] + globalv.entry[2];

                /*This is more accurate method to accumlate time*/
                //dis.entry[0] = cur_point[0] - pre_point[0];
                //dis.entry[1] = cur_point[1] - pre_point[1];

                /*temporal tau*/
                gt_tau+=/*length(dis)*/move_dist/cur_vec_mag;
                //////////////g_dt += len*len/fabs(dot(dis, glob_loc_v));

                /////////////*spatial tau*/
                //////////////double len = length(dis);
                //////////////g_dt += len;

                /*   now the tau increment should be constant! 2/8/2010 */
                //gt_tau += 0.01;

                if(gt_tau >= tau)
                {
                    flag = 1;
                    return face_id;
                }
            }

            else{  ////the curve reach a singularity
                flag = 1;

                ////Store the record into global line segment array

                return face_id;
            }
        }

        ////3. if the point is out of current triangle
        else{
            double t[2] = {0.};

            int PassVertornot = 0;

            get_next_triangle(face_id, pre_point, cur_point, t, 1, PassVertornot, alpha);

            ////update the globalpoint here (Modified on 01/30/07)
            if(PassVertornot > 0)
            {
                if (face_id < 0)
                    return -1;

                //we first need to know which vertex it is in the new triangle 01/30/07
                Vertex* vertid = pre_f->verts[PassVertornot-1];
                Triangle *cur_f = object->tlist.tris[face_id];
                int vert_new = 0;
                for(int k = 0; k < 3; k++)
                {
                    if(cur_f->verts[k] == vertid)
                    {
                        vert_new = k;
                        break;
                    }
                }

                alpha[vert_new]=1-0.0001;
                alpha[(vert_new+1)%3]=0.00005;
                alpha[(vert_new+2)%3]=0.00005;


                /* Get the new cur_point */
                cur_point[0] = alpha[1]*cur_f->x1+alpha[2]*cur_f->x2;
                cur_point[1] = alpha[2]*cur_f->y2;

                globalv = cur_point[0] * cur_f->LX + cur_point[1] * cur_f->LY;

                globalp[0] = cur_f->verts[0]->x + globalv.entry[0];
                globalp[1] = cur_f->verts[0]->y + globalv.entry[1];
                globalp[2] = cur_f->verts[0]->z + globalv.entry[2];
                face=cur_f;
            }

            else{
                //// transfer it to the global coordinates
                globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;

                globalp[0] = vert0[0] + globalv.entry[0];
                globalp[1] = vert0[1] + globalv.entry[1];
                globalp[2] = vert0[2] + globalv.entry[2];
            }

            dis.entry[0] = cur_point[0] - pre_point[0];
            dis.entry[1] = cur_point[1] - pre_point[1];

            /*temporal tau*/
            //double temp_ratio=length(dis)/cur_vec_mag;
            gt_tau+=length(dis)/cur_vec_mag;  // This computation is problematic  02/08/2010
            ////////g_dt += len*len/fabs(dot(dis, glob_loc_v));

            ///////*spatial tau*/
            ////////double len = length(dis);
            ////////g_dt += len;

            //gt_tau += 0.01;

            if (gt_tau > tau)                 // Added by Guoning at 02/08/2010
            {
                flag = 1;
                //globalv = pre_point[0] * pre_f->LX + pre_point[1] * pre_f->LY;
                //globalp[0] = pre_f->verts[0]->x + globalv.entry[0];
                //globalp[1] = pre_f->verts[0]->y + globalv.entry[1];
                //globalp[2] = pre_f->verts[0]->z + globalv.entry[2];

                //gt_tau-=temp_ratio;
                return pre_f->index;
            }

            return face_id;
        }

    }

    return face_id;
}


//-----------------------------------------------------------------------------------

bool Trajectory::cal_next_point_euler1_tau(double first[2], double second[2], int &face_id,
                                           double alpha[3], int type)
{
    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    icVector2 VecAtPoint = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if(length(VecAtPoint) < 1e-15 ) return false;

    //scale the vector??? 01/23/07 !!!!!!!!!!!!!Important modification
    normalize(VecAtPoint); // removed by Guoning on 07/13/2010

    Triangle *t = object->tlist.tris[face_id];
    VecAtPoint = 0.4*length(t->verts[0]->t_vec)*VecAtPoint;

    ////Using first order Euler method to test 12/11/05
    if(type == 0)
    {
        second[0] = first[0] + 2*eulerstep_scalar * VecAtPoint.entry[0];
        second[1] = first[1] + 2*eulerstep_scalar * VecAtPoint.entry[1];
    }
    else
    {
        second[0] = first[0] - 2*eulerstep_scalar * VecAtPoint.entry[0];
        second[1] = first[1] - 2*eulerstep_scalar * VecAtPoint.entry[1];
    }

    /*ready for the temporary tau  08/30/07*/
    {
        icVector2 line_v;
        line_v.entry[0] = second[0]-first[0];
        line_v.entry[1] = second[1]-first[1];
        normalize(line_v);
        double proj_len = dot(line_v, VecAtPoint); /*here we use the dot product in the local frame*/
        cur_vec_mag = fabs(proj_len);
    }

    return true;
}



bool
Trajectory::cal_next_point_euler1_tau_f(double first[2], double second[2], int &face_id,
                                        double alpha[3])
{
    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    icVector2 VecAtPoint = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if(length(VecAtPoint) < 1e-10 ) return false;

    //scale the vector??? 01/23/07 !!!!!!!!!!!!!Important modification
    normalize(VecAtPoint); // removed by Guoning on 07/13/2010

    Triangle *t = object->tlist.tris[face_id];
    //VecAtPoint = 0.4/**length(t->verts[0]->t_vec)*/*VecAtPoint;
    //eulerstep_scalar = sqrt(t->area)/20.;  // this is applied whether the vector is normalized

    //double eulerstep_scalar = t->euler_step;

    ////Using first order Euler method to test 12/11/05
    second[0] = first[0] + /*2**/eulerstep_scalar * VecAtPoint.entry[0];
    second[1] = first[1] + /*2**/eulerstep_scalar * VecAtPoint.entry[1];

    /*ready for the temporary tau  08/30/07*/
    {
        icVector2 line_v;
        line_v.entry[0] = second[0]-first[0];
        line_v.entry[1] = second[1]-first[1];
        move_dist = length(line_v);
        //normalize(line_v);

        if (move_dist>0)
            line_v = 1./move_dist * line_v;
        double proj_len = dot(line_v, VecAtPoint); /*here we use the dot product in the local frame*/
        cur_vec_mag = fabs(proj_len);
    }

    return true;
}


bool
Trajectory::cal_next_point_euler1_tau_b(double first[2], double second[2], int &face_id,
                                        double alpha[3])
{
    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    icVector2 VecAtPoint = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if(length(VecAtPoint) < 1e-10 ) return false;

    //scale the vector??? 01/23/07 !!!!!!!!!!!!!Important modification
    normalize(VecAtPoint); // removed by Guoning on 07/13/2010

    Triangle *t = object->tlist.tris[face_id];
    //VecAtPoint = 0.4*length(t->verts[0]->t_vec)*VecAtPoint;
    //VecAtPoint = 0.4*VecAtPoint;
    //eulerstep_scalar = sqrt(t->area)/20.;  // this is applied whether the vector is normalized
    //double eulerstep_scalar = t->euler_step;

    ////Using first order Euler method to test 12/11/05
    second[0] = first[0] - /*2**/eulerstep_scalar * VecAtPoint.entry[0];
    second[1] = first[1] - /*2**/eulerstep_scalar * VecAtPoint.entry[1];

    /*ready for the temporary tau  08/30/07*/
    {
        icVector2 line_v;
        line_v.entry[0] = second[0]-first[0];
        line_v.entry[1] = second[1]-first[1];
        move_dist = length(line_v);
        //normalize(line_v);

        if (move_dist>0)
            line_v = 1./move_dist * line_v;
        double proj_len = dot(line_v, VecAtPoint); /*here we use the dot product in the local frame*/
        cur_vec_mag = fabs(proj_len);
    }

    return true;
}


/*****************************************************************
    The non normalized tracing
*/
bool Trajectory::cal_next_point_euler1_tau_nonNorm(double first[2], double second[2], int &face_id,
                                                   double alpha[3], int type)
{
    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    icVector2 VecAtPoint = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if (object->tlist.tris[face_id]->singularityID >= 0)
    {
        double scale = 1./(1+length(VecAtPoint));
        VecAtPoint = scale *VecAtPoint;
    }


    if(length(VecAtPoint) < 1.e-10 ) return false;  // modified by Guoning 02/18/2010


    //scale the vector??? 01/23/07 !!!!!!!!!!!!!Important modification
    //normalize(VecAtPoint);

    //Triangle *t = object->tlist.tris[face_id];
    //VecAtPoint = 0.4*length(t->verts[0]->t_vec)*VecAtPoint;

    ////Using first order Euler method to test 12/11/05
    if(type == 0)
    {
        second[0] = first[0] + 2*eulerstep_scalar * VecAtPoint.entry[0];
        second[1] = first[1] + 2*eulerstep_scalar * VecAtPoint.entry[1];
    }
    else
    {
        second[0] = first[0] - 2*eulerstep_scalar * VecAtPoint.entry[0];
        second[1] = first[1] - 2*eulerstep_scalar * VecAtPoint.entry[1];
    }

    /*ready for the temporary tau  08/30/07*/
    {
        icVector2 line_v;
        line_v.entry[0] = second[0]-first[0];
        line_v.entry[1] = second[1]-first[1];
        normalize(line_v);
        double proj_len = dot(line_v, VecAtPoint); /*here we use the dot product in the local frame*/
        cur_vec_mag = fabs(proj_len);
    }

    return true;
}


bool Trajectory::cal_next_point_euler1_tau_nonNorm_f(double first[2], double second[2], int &face_id,
                                                     double alpha[3])
{
    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    icVector2 VecAtPoint = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if (object->tlist.tris[face_id]->singularityID >= 0)
    {
        double scale = 1./(1+length(VecAtPoint));
        VecAtPoint = scale *VecAtPoint;
    }


    if(length(VecAtPoint) < 1.e-10 ) return false;  // modified by Guoning 02/18/2010


    //scale the vector??? 01/23/07 !!!!!!!!!!!!!Important modification
    //normalize(VecAtPoint);

    //Triangle *t = object->tlist.tris[face_id];
    //VecAtPoint = 0.4*length(t->verts[0]->t_vec)*VecAtPoint;

    ////Using first order Euler method to test 12/11/05
    second[0] = first[0] + 2*eulerstep_scalar * VecAtPoint.entry[0];
    second[1] = first[1] + 2*eulerstep_scalar * VecAtPoint.entry[1];

    /*ready for the temporary tau  08/30/07*/
    {
        icVector2 line_v;
        line_v.entry[0] = second[0]-first[0];
        line_v.entry[1] = second[1]-first[1];
        normalize(line_v);
        double proj_len = dot(line_v, VecAtPoint); /*here we use the dot product in the local frame*/
        cur_vec_mag = fabs(proj_len);
    }

    return true;
}

bool Trajectory::cal_next_point_euler1_tau_nonNorm_b(double first[2], double second[2], int &face_id,
                                                     double alpha[3])
{
    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    icVector2 VecAtPoint = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if (object->tlist.tris[face_id]->singularityID >= 0)
    {
        double scale = 1./(1+length(VecAtPoint));
        VecAtPoint = scale *VecAtPoint;
    }


    if(length(VecAtPoint) < 1.e-10 ) return false;  // modified by Guoning 02/18/2010


    //scale the vector??? 01/23/07 !!!!!!!!!!!!!Important modification
    //normalize(VecAtPoint);

    //Triangle *t = object->tlist.tris[face_id];
    //VecAtPoint = 0.4*length(t->verts[0]->t_vec)*VecAtPoint;

    ////Using first order Euler method to test 12/11/05

    second[0] = first[0] - 2*eulerstep_scalar * VecAtPoint.entry[0];
    second[1] = first[1] - 2*eulerstep_scalar * VecAtPoint.entry[1];


    /*ready for the temporary tau  08/30/07*/
    {
        icVector2 line_v;
        line_v.entry[0] = second[0]-first[0];
        line_v.entry[1] = second[1]-first[1];
        normalize(line_v);
        double proj_len = dot(line_v, VecAtPoint); /*here we use the dot product in the local frame*/
        cur_vec_mag = fabs(proj_len);
    }

    return true;
}

/************************************************************************/
/*
*/
bool Trajectory::cal_next_point_2ndEuler_tau_nonNorm(double first[2],
                                                     double second[2],
                                                     int &face_id,
                                                     double alpha[3],
                                                     int type)
{
    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    icVector2 VecAtPoint = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if(length(VecAtPoint) < 1e-10 ) return false;

    //scale the vector??? 01/23/07 !!!!!!!!!!!!!Important modification
    //normalize(VecAtPoint);

    //Triangle *t = object->tlist.tris[face_id];
    //VecAtPoint = 0.4*length(t->verts[0]->t_vec)*VecAtPoint;

    ////Using first order Euler method to test 12/11/05

    double euler2nd_step = eulerstep_scalar/*/2*/;
    //double euler2nd_step = 2*eulerstep_scalar/*/2*/;

    if(type == 0)
    {
        second[0] = first[0] + euler2nd_step/2. * VecAtPoint.entry[0];
        second[1] = first[1] + euler2nd_step/2. * VecAtPoint.entry[1];
    }
    else
    {
        second[0] = first[0] - euler2nd_step/2. * VecAtPoint.entry[0];
        second[1] = first[1] - euler2nd_step/2. * VecAtPoint.entry[1];
    }

    /*compute K2*/
    ////get the vector at next point
    double alpha1[3] = {0.};
    local_To_global(face_id, second, t_gp);
    object->get_2D_Barycentric_Facters(face_id, second[0], second[1], alpha1);

    if ( (alpha1[0] >= 0 ||fabs(alpha1[0]) <= 1e-8) && alpha1[0] <= 1
        && (alpha1[1] >= 0 || fabs(alpha1[1]) <= 1e-8) && alpha1[1] <= 1
        && (alpha1[2] >= 0 || fabs(alpha1[2]) <= 1e-8) && alpha1[2] <= 1)
    {
        icVector2 VecAtPoint2 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha1, second[0], second[1]);

        //Triangle *face = object->tlist.tris[face_id];
        //icVector3 gp = second[0]*face->LX + second[1]*face->LY;
        //icVector2 VecAtPoint2 = object->get_Vector_At_Point(face_id, gp.entry, alpha, second[0], second[1]);

        VecAtPoint = 0.5*(VecAtPoint + VecAtPoint2);
        if(type == 0)
        {
            second[0] = first[0] + euler2nd_step * VecAtPoint.entry[0];
            second[1] = first[1] + euler2nd_step * VecAtPoint.entry[1];
        }
        else
        {
            second[0] = first[0] - euler2nd_step * VecAtPoint.entry[0];
            second[1] = first[1] - euler2nd_step * VecAtPoint.entry[1];
        }
    }

    /*ready for the temporary tau  08/30/07*/
    icVector2 line_v;
    line_v.entry[0] = second[0]-first[0];
    line_v.entry[1] = second[1]-first[1];
    normalize(line_v);
    double proj_len = dot(line_v, VecAtPoint); /*here we use the dot product in the local frame*/
    cur_vec_mag = fabs(proj_len);

    return true;
}


bool
Trajectory::cal_next_point_2ndEuler_tau_Norm(double first[2],
                                             double second[2],
                                             int &face_id,
                                             double alpha[3],
                                             int type)
{
    ////Using first order Euler method to get the next point
    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    icVector2 VecAtPoint = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if(length(VecAtPoint) < 1e-10) return false;

    //scale the vector??? 01/23/07 !!!!!!!!!!!!!Important modification
    normalize(VecAtPoint);

    Triangle *t = object->tlist.tris[face_id];
    VecAtPoint = 0.08*length(t->verts[0]->t_vec)*VecAtPoint; /*evaluate in very small step now*/

    double temp[2] = {0.};
    double euler2nd_step = eulerstep_scalar;

    if(type == 0)
    {
        temp[0] = first[0] + /*2**/euler2nd_step/2.*VecAtPoint.entry[0];
        temp[1] = first[1] + /*2**/euler2nd_step/2*VecAtPoint.entry[1];
    }
    else
    {
        temp[0] = first[0] - /*2**/euler2nd_step/2*VecAtPoint.entry[0];
        temp[1] = first[1] - /*2**/euler2nd_step/2*VecAtPoint.entry[1];
    }

    ////get the vector at next point
    double alpha1[3] = {0.};
    local_To_global(face_id, temp, t_gp);
    object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha1);

    if ( (alpha1[0] >= 0 ||fabs(alpha1[0]) <= 1e-8) && alpha1[0] <= 1
        && (alpha1[1] >= 0 || fabs(alpha1[1]) <= 1e-8) && alpha1[1] <= 1
        && (alpha1[2] >= 0 || fabs(alpha1[2]) <= 1e-8) && alpha1[2] <= 1)

    {
        icVector2 VecAtPoint2 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha1, temp[0], temp[1]);
        normalize(VecAtPoint);
        normalize(VecAtPoint2);

        icVector2 total_v;
        total_v.entry[0] = (VecAtPoint.entry[0] +  VecAtPoint2.entry[0])/2;
        total_v.entry[1] = (VecAtPoint.entry[1] +  VecAtPoint2.entry[1])/2;
        total_v = length(t->verts[0]->t_vec)*VecAtPoint;

        //if(type == 0)
        //{
        //	second[0] = first[0] + 2*eulerstep_scalar*total_v.entry[0];
        //	second[1] = first[1] + 2*eulerstep_scalar*total_v.entry[1];
        //}
        //else
        //{
        //	second[0] = first[0] - 2*eulerstep_scalar*total_v.entry[0];
        //	second[1] = first[1] - 2*eulerstep_scalar*total_v.entry[1];
        //}
        if(type == 0)
        {
            second[0] = first[0] + euler2nd_step/*2*eulerstep_scalar*/*total_v.entry[0];
            second[1] = first[1] + euler2nd_step/*2*eulerstep_scalar*/*total_v.entry[1];
        }
        else
        {
            second[0] = first[0] - euler2nd_step/*2*eulerstep_scalar*/*total_v.entry[0];
            second[1] = first[1] - euler2nd_step/*2*eulerstep_scalar*/*total_v.entry[1];
        }

        icVector2 line_v;
        line_v.entry[0] = second[0]-first[0];
        line_v.entry[1] = second[1]-first[1];
        normalize(line_v);
        double proj_len = dot(line_v, total_v); /*here we use the dot product in the local frame*/
        cur_vec_mag = fabs(proj_len);
    }

    else
    {
        if(type == 0)
        {
            second[0] = first[0] + euler2nd_step/*2*eulerstep_scalar*/*VecAtPoint.entry[0];
            second[1] = first[1] + euler2nd_step/*2*eulerstep_scalar*/*VecAtPoint.entry[1];
        }
        else
        {
            second[0] = first[0] - euler2nd_step/*2*eulerstep_scalar*/*VecAtPoint.entry[0];
            second[1] = first[1] - euler2nd_step/*2*eulerstep_scalar*/*VecAtPoint.entry[1];
        }
        icVector2 line_v;
        line_v.entry[0] = second[0]-first[0];
        line_v.entry[1] = second[1]-first[1];
        normalize(line_v);
        double proj_len = dot(line_v, VecAtPoint); /*here we use the dot product in the local frame*/
        cur_vec_mag = fabs(proj_len);
    }

    return true;

}


/////////////////////////////////////////////////////////////////////////////////////

bool
Trajectory::cal_next_point_2ndEuler_Norm_b(  double first[2],
                                           double second[2],
                                           int &face_id,
                                           double alpha[3])
{
    ////Using first order Euler method to get the next point
    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    icVector2 VecAtPoint = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if(length(VecAtPoint) < 1e-10) return false;

    //scale the vector??? 01/23/07 !!!!!!!!!!!!!Important modification
    normalize(VecAtPoint);

    Triangle *t = object->tlist.tris[face_id];
    //VecAtPoint = 0.08*/*length(t->verts[0]->t_vec)**/VecAtPoint; /*evaluate in very small step now*/
    //eulerstep_scalar = sqrt(t->area)/30.;  // this is applied whether the vector is normalized

    double temp[2] = {0.};
    double euler2nd_step = t->euler_step/*eulerstep_scalar*/;

    temp[0] = first[0] - euler2nd_step/2*VecAtPoint.entry[0];
    temp[1] = first[1] - euler2nd_step/2*VecAtPoint.entry[1];

    ////get the vector at next point
    double alpha1[3] = {0.};
    local_To_global(face_id, temp, t_gp);
    object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha1);

    if ( (alpha1[0] >= 0 ||fabs(alpha1[0]) <= 1e-8) && alpha1[0] <= 1
        && (alpha1[1] >= 0 || fabs(alpha1[1]) <= 1e-8) && alpha1[1] <= 1
        && (alpha1[2] >= 0 || fabs(alpha1[2]) <= 1e-8) && alpha1[2] <= 1)

    {
        icVector2 VecAtPoint2 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha1, temp[0], temp[1]);
        normalize(VecAtPoint2);

        icVector2 total_v = 0.5*(VecAtPoint + VecAtPoint2);

        normalize(total_v);

        second[0] = first[0] - euler2nd_step*total_v.entry[0];
        second[1] = first[1] - euler2nd_step*total_v.entry[1];
    }

    else
    {
        second[0] = first[0] - euler2nd_step*VecAtPoint.entry[0];
        second[1] = first[1] - euler2nd_step*VecAtPoint.entry[1];

    }

    return true;
}

bool
Trajectory::cal_next_point_2ndEuler_Norm_f(  double first[2],
                                           double second[2],
                                           int &face_id,
                                           double alpha[3])
{
    ////Using first order Euler method to get the next point
    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    icVector2 VecAtPoint = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if(length(VecAtPoint) < 1e-10) return false;

    //scale the vector??? 01/23/07 !!!!!!!!!!!!!Important modification
    normalize(VecAtPoint);

    Triangle *t = object->tlist.tris[face_id];
    //VecAtPoint = 0.08*/*length(t->verts[0]->t_vec)**/VecAtPoint; /*evaluate in very small step now*/
    //eulerstep_scalar = sqrt(t->area)/30.;  // this is applied whether the vector is normalized

    double temp[2] = {0.};
    double euler2nd_step = t->euler_step/*eulerstep_scalar*/;

    temp[0] = first[0] - euler2nd_step/2*VecAtPoint.entry[0];
    temp[1] = first[1] - euler2nd_step/2*VecAtPoint.entry[1];

    ////get the vector at next point
    double alpha1[3] = {0.};
    local_To_global(face_id, temp, t_gp);
    object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha1);

    if ( (alpha1[0] >= 0 ||fabs(alpha1[0]) <= 1e-8) && alpha1[0] <= 1
        && (alpha1[1] >= 0 || fabs(alpha1[1]) <= 1e-8) && alpha1[1] <= 1
        && (alpha1[2] >= 0 || fabs(alpha1[2]) <= 1e-8) && alpha1[2] <= 1)

    {
        icVector2 VecAtPoint2 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha1, temp[0], temp[1]);
        normalize(VecAtPoint2);

        icVector2 total_v = 0.5*(VecAtPoint + VecAtPoint2);

        normalize(total_v);

        second[0] = first[0] - euler2nd_step*total_v.entry[0];
        second[1] = first[1] - euler2nd_step*total_v.entry[1];
    }

    else
    {
        second[0] = first[0] - euler2nd_step*VecAtPoint.entry[0];
        second[1] = first[1] - euler2nd_step*VecAtPoint.entry[1];

    }

    return true;
}

bool
Trajectory::cal_next_point_euler1_f(double first[2], double second[2], int &face_id,
                                    double alpha[3])
{
    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    icVector2 VecAtPoint = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if(length(VecAtPoint) < 1e-10 ) return false;

    //scale the vector??? 01/23/07 !!!!!!!!!!!!!Important modification
    normalize(VecAtPoint); // removed by Guoning on 07/13/2010

    Triangle *t = object->tlist.tris[face_id];
    //VecAtPoint = 0.4/**length(t->verts[0]->t_vec)*/*VecAtPoint;
    //eulerstep_scalar = sqrt(t->area)/20.;  // this is applied whether the vector is normalized

    double eulerstep_scalar = t->euler_step;

    ////Using first order Euler method to test 12/11/05
    second[0] = first[0] + /*2**/eulerstep_scalar * VecAtPoint.entry[0];
    second[1] = first[1] + /*2**/eulerstep_scalar * VecAtPoint.entry[1];

    return true;
}


bool
Trajectory::cal_next_point_euler1_b(double first[2], double second[2], int &face_id,
                                    double alpha[3])
{
    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    icVector2 VecAtPoint = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if(length(VecAtPoint) < 1e-10 ) return false;

    //scale the vector??? 01/23/07 !!!!!!!!!!!!!Important modification
    normalize(VecAtPoint); // removed by Guoning on 07/13/2010

    Triangle *t = object->tlist.tris[face_id];
    //VecAtPoint = 0.4*length(t->verts[0]->t_vec)*VecAtPoint;
    //VecAtPoint = 0.4*VecAtPoint;
    //eulerstep_scalar = sqrt(t->area)/20.;  // this is applied whether the vector is normalized
    double eulerstep_scalar = t->euler_step;

    ////Using first order Euler method to test 12/11/05
    second[0] = first[0] - /*2**/eulerstep_scalar * VecAtPoint.entry[0];
    second[1] = first[1] - /*2**/eulerstep_scalar * VecAtPoint.entry[1];

    return true;
}



bool
Trajectory::get_next_pt_rotsum_f(double first[2],
                                 double second[2],
                                 int &face_id,
                                 double alpha[3],
                                 unsigned char opt)
{
    switch (opt)
    {
    case 0:
        //return cal_next_point_euler1_tau_nonNorm_f(first, second, face_id, alpha);
        //return cal_next_point_euler1_tau_f(first, second, face_id, alpha);
        return cal_next_point_euler1_f(first, second, face_id, alpha);

    case 1:
        //return cal_next_point_2ndEuler_tau_nonNorm(first, second, face_id, alpha, type);
        return cal_next_point_2ndEuler_Norm_f(first, second, face_id, alpha);
    case 2:
        return cal_next_pt_RK4_tau_f(first, second, face_id, alpha);

    default:
        return false;
    }
}


bool
Trajectory::get_next_pt_rotsum_b(double first[2],
                                 double second[2],
                                 int &face_id,
                                 double alpha[3],
                                 unsigned char opt)
{
    switch (opt)
    {
    case 0:
        //return cal_next_point_euler1_tau_nonNorm_b(first, second, face_id, alpha);
        //return cal_next_point_euler1_tau_b(first, second, face_id, alpha);
        return cal_next_point_euler1_b(first, second, face_id, alpha);

    case 1:
        return cal_next_point_2ndEuler_Norm_b(first, second, face_id, alpha);
    case 2:
        return cal_next_pt_RK4_tau_b(first, second, face_id, alpha);

    default:
        return false;
    }
}

///////////////////////////////////////////////////////////////////////////////

bool
Trajectory::cal_next_point_2ndEuler_tau_Norm_f(double first[2],
                                               double second[2],
                                               int &face_id,
                                               double alpha[3])
{
    ////Using first order Euler method to get the next point
    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    icVector2 VecAtPoint = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if(length(VecAtPoint) < 1e-10) return false;

    //scale the vector??? 01/23/07 !!!!!!!!!!!!!Important modification
    normalize(VecAtPoint);

    Triangle *t = object->tlist.tris[face_id];
    //VecAtPoint = 0.08*/*length(t->verts[0]->t_vec)*/VecAtPoint; /*evaluate in very small step now*/
    //eulerstep_scalar = sqrt(t->area)/30.;  // this is applied whether the vector is normalized
    //double eulerstep_scalar = t->euler_step;

    double temp[2] = {0.};
    double euler2nd_step = t->euler_step;

    temp[0] = first[0] + /*2**/euler2nd_step/2.*VecAtPoint.entry[0];
    temp[1] = first[1] + /*2**/euler2nd_step/2*VecAtPoint.entry[1];

    ////get the vector at next point
    double alpha1[3] = {0.};
    local_To_global(face_id, temp, t_gp);
    object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha1);

    if ( (alpha1[0] >= 0 ||fabs(alpha1[0]) <= 1e-8) && alpha1[0] <= 1
        && (alpha1[1] >= 0 || fabs(alpha1[1]) <= 1e-8) && alpha1[1] <= 1
        && (alpha1[2] >= 0 || fabs(alpha1[2]) <= 1e-8) && alpha1[2] <= 1)

    {
        icVector2 VecAtPoint2 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha1, temp[0], temp[1]);
        normalize(VecAtPoint);
        normalize(VecAtPoint2);

        icVector2 total_v;
        total_v.entry[0] = (VecAtPoint.entry[0] +  VecAtPoint2.entry[0])/2;
        total_v.entry[1] = (VecAtPoint.entry[1] +  VecAtPoint2.entry[1])/2;
        //total_v = /*length(t->verts[0]->t_vec)**/VecAtPoint;

        second[0] = first[0] + euler2nd_step/*2*eulerstep_scalar*/*total_v.entry[0];
        second[1] = first[1] + euler2nd_step/*2*eulerstep_scalar*/*total_v.entry[1];

        icVector2 line_v;
        line_v.entry[0] = second[0]-first[0];
        line_v.entry[1] = second[1]-first[1];
        move_dist = length(line_v);
        //normalize(line_v);

        if (move_dist>0)
            line_v = 1./move_dist * line_v;
        double proj_len = dot(line_v, total_v); /*here we use the dot product in the local frame*/
        cur_vec_mag = fabs(proj_len);
    }

    else
    {
        second[0] = first[0] + euler2nd_step/*2*eulerstep_scalar*/*VecAtPoint.entry[0];
        second[1] = first[1] + euler2nd_step/*2*eulerstep_scalar*/*VecAtPoint.entry[1];

        icVector2 line_v;
        line_v.entry[0] = second[0]-first[0];
        line_v.entry[1] = second[1]-first[1];
        move_dist = length(line_v);
        //normalize(line_v);

        if (move_dist>0)
            line_v = 1./move_dist * line_v;
        double proj_len = dot(line_v, VecAtPoint); /*here we use the dot product in the local frame*/
        cur_vec_mag = fabs(proj_len);
    }

    return true;

}


bool
Trajectory::cal_next_point_2ndEuler_tau_Norm_b(double first[2],
                                               double second[2],
                                               int &face_id,
                                               double alpha[3])
{
    ////Using first order Euler method to get the next point
    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    icVector2 VecAtPoint = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if(length(VecAtPoint) < 1e-10) return false;

    //scale the vector??? 01/23/07 !!!!!!!!!!!!!Important modification
    normalize(VecAtPoint);

    Triangle *t = object->tlist.tris[face_id];
    //VecAtPoint = 0.08*/*length(t->verts[0]->t_vec)**/VecAtPoint; /*evaluate in very small step now*/
    //eulerstep_scalar = sqrt(t->area)/30.;  // this is applied whether the vector is normalized

    //double eulerstep_scalar = t->euler_step;

    double temp[2] = {0.};
    double euler2nd_step = t->euler_step/*eulerstep_scalar*/;

    temp[0] = first[0] - /*2**/euler2nd_step/2*VecAtPoint.entry[0];
    temp[1] = first[1] - /*2**/euler2nd_step/2*VecAtPoint.entry[1];

    ////get the vector at next point
    double alpha1[3] = {0.};
    local_To_global(face_id, temp, t_gp);
    object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha1);

    if ( (alpha1[0] >= 0 ||fabs(alpha1[0]) <= 1e-8) && alpha1[0] <= 1
        && (alpha1[1] >= 0 || fabs(alpha1[1]) <= 1e-8) && alpha1[1] <= 1
        && (alpha1[2] >= 0 || fabs(alpha1[2]) <= 1e-8) && alpha1[2] <= 1)

    {
        icVector2 VecAtPoint2 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha1, temp[0], temp[1]);
        normalize(VecAtPoint);
        normalize(VecAtPoint2);

        icVector2 total_v;
        total_v.entry[0] = (VecAtPoint.entry[0] +  VecAtPoint2.entry[0])/2;
        total_v.entry[1] = (VecAtPoint.entry[1] +  VecAtPoint2.entry[1])/2;
        //total_v = length(t->verts[0]->t_vec)*VecAtPoint;

        second[0] = first[0] - euler2nd_step/*2*eulerstep_scalar*/*total_v.entry[0];
        second[1] = first[1] - euler2nd_step/*2*eulerstep_scalar*/*total_v.entry[1];

        icVector2 line_v;
        line_v.entry[0] = second[0]-first[0];
        line_v.entry[1] = second[1]-first[1];
        move_dist = length(line_v);
        //normalize(line_v);

        if (move_dist>0)
            line_v = 1./move_dist * line_v;
        double proj_len = dot(line_v, total_v); /*here we use the dot product in the local frame*/
        cur_vec_mag = fabs(proj_len);
    }

    else
    {
        second[0] = first[0] - euler2nd_step/*2*eulerstep_scalar*/*VecAtPoint.entry[0];
        second[1] = first[1] - euler2nd_step/*2*eulerstep_scalar*/*VecAtPoint.entry[1];
        icVector2 line_v;
        line_v.entry[0] = second[0]-first[0];
        line_v.entry[1] = second[1]-first[1];
        move_dist = length(line_v);
        //normalize(line_v);

        if (move_dist>0)
            line_v = 1./move_dist * line_v;
        double proj_len = dot(line_v, VecAtPoint); /*here we use the dot product in the local frame*/
        cur_vec_mag = fabs(proj_len);
    }

    return true;

}

//-------------------------------------------------------

bool
Trajectory::get_next_pt_tau(double first[2],
                            double second[2],
                            int &face_id,
                            double alpha[3],
                            int type,
                            unsigned char opt)
{
    switch (opt)
    {
    case 0:
        return cal_next_point_euler1_tau_nonNorm(first, second, face_id, alpha, type);
    case 1:
        //return cal_next_point_2ndEuler_tau_nonNorm(first, second, face_id, alpha, type);
        return cal_next_point_2ndEuler_tau_Norm(first, second, face_id, alpha, type);
    case 2:
        return cal_next_pt_RK4_tau(first, second, face_id, alpha, type);

    default:
        return false;
    }
}


bool
Trajectory::get_next_pt_tau_f(double first[2],
                              double second[2],
                              int &face_id,
                              double alpha[3],
                              unsigned char opt)
{
    switch (opt)
    {
    case 0:
        //return cal_next_point_euler1_tau_nonNorm_f(first, second, face_id, alpha);
        return cal_next_point_euler1_tau_f(first, second, face_id, alpha);
    case 1:
        //return cal_next_point_2ndEuler_tau_nonNorm(first, second, face_id, alpha, type);
        return cal_next_point_2ndEuler_tau_Norm_f(first, second, face_id, alpha);
        //return cal_next_point_2ndEuler_Norm_f(first, second, face_id, alpha);
    case 2:
        return cal_next_pt_RK4_tau_f(first, second, face_id, alpha);

    default:
        return false;
    }
}


bool
Trajectory::get_next_pt_tau_b(double first[2],
                              double second[2],
                              int &face_id,
                              double alpha[3],
                              unsigned char opt)
{
    switch (opt)
    {
    case 0:
        //return cal_next_point_euler1_tau_nonNorm_b(first, second, face_id, alpha);
        return cal_next_point_euler1_tau_b(first, second, face_id, alpha);
    case 1:
        //return cal_next_point_2ndEuler_tau_nonNorm(first, second, face_id, alpha, type);
        return cal_next_point_2ndEuler_tau_Norm_b(first, second, face_id, alpha);
        //return cal_next_point_2ndEuler_Norm_b(first, second, face_id, alpha);
    case 2:
        return cal_next_pt_RK4_tau_b(first, second, face_id, alpha);

    default:
        return false;
    }
}


/*Find the smallest edge distance in the triangle*/
double get_shortestedge_tri(int tri)
{
    int i;
    Triangle *face = object->tlist.tris[tri];
    double shortest = face->edges[0]->length;

    for(i=1; i<3; i++)
    {
        if(face->edges[i]->length < shortest)
            shortest = face->edges[i]->length;
    }

    return shortest;
}


//bool
//Trajectory::cal_next_pt_RK4(double first[2], double second[2], int &face_id, double alpha[3], int type)
//{
//	////Using 4th Runge-Kutta method to get the next point
//	icVector2 t_v;
//	double temp[2] = {0.};
//
//	double RK4_step;
//
//	/*compute K1*/
//	icVector3 t_gp;
//	local_To_global(face_id, first, t_gp);
//	icVector2 V1 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);
//
//	if(length(V1) < 1e-14) return false;
//
//	double shortest = get_shortestedge_tri(face_id);
//	double scalar = 0.06;
//
//	//if(Mag_Scheme == 0)
//	//{
//	//	RK4_step = scalar*shortest*0.5;
//	//}
//	//else if(Mag_Scheme == 1)
//	//{
//		RK4_step = scalar*shortest;
//	//}
//	//else
//	//{
//	//	RK4_step = scalar*shortest*2;
//	//}
//
//	t_v = V1;
//	normalize(t_v);
//
//	if(type == 0)
//	{
//		temp[0] = first[0] + RK4_step/2*t_v.entry[0];
//		temp[1] = first[1] + RK4_step/2*t_v.entry[1];
//	}
//	else
//	{
//		temp[0] = first[0] - RK4_step/2*t_v.entry[0];
//		temp[1] = first[1] - RK4_step/2*t_v.entry[1];
//	}
//
//	/*compute K2*/
//	local_To_global(face_id, temp, t_gp);
//	object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha);
//	icVector2 V2 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, temp[0], temp[1]);
//	t_v = V2;
//	normalize(t_v);
//
//	if(type == 0)
//	{
//		temp[0] = first[0] + RK4_step/2*t_v.entry[0];
//		temp[1] = first[1] + RK4_step/2*t_v.entry[1];
//	}
//	else
//	{
//		temp[0] = first[0] - RK4_step/2*t_v.entry[0];
//		temp[1] = first[1] - RK4_step/2*t_v.entry[1];
//	}
//
//	/*compute K3*/
//	local_To_global(face_id, temp, t_gp);
//	object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha);
//	icVector2 V3 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, temp[0], temp[1]);
//	t_v = V3;
//	normalize(t_v);
//
//	if(type == 0)
//	{
//		temp[0] = first[0] + RK4_step*t_v.entry[0];
//		temp[1] = first[1] + RK4_step*t_v.entry[1];
//	}
//	else
//	{
//		temp[0] = first[0] - RK4_step*t_v.entry[0];
//		temp[1] = first[1] - RK4_step*t_v.entry[1];
//	}
//
//	/*compute K4*/
//	local_To_global(face_id, temp, t_gp);
//	object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha);
//	icVector2 V4 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, temp[0], temp[1]);
//
//	icVector2 total_v = 1./6.*(V1+2*V2+2*V3+V4);
//	t_v = total_v;
//	normalize(t_v);
//
//	if(type == 0)
//	{
//		second[0] = first[0] + RK4_step*t_v.entry[0];
//		second[1] = first[1] + RK4_step*t_v.entry[1];
//	}
//
//	else
//	{
//		second[0] = first[0] - RK4_step*t_v.entry[0];
//		second[1] = first[1] - RK4_step*t_v.entry[1];
//	}
//
//	//hstep = euler_stepsize*length(total_v);
//
//	/*********07/05/2007**********/
//	//if(time_or_patial == 0)  /*use spatial tau*/
//		//g_dt += hstep;
//
//
//	/*we need to accumulate the time rather than the distance */
//	/*calculate the t according to the step and the magnitude of the vector,
//	but we need to first project the vector onto the line segment to get the real length*/
//	icVector2 line_v;
//	line_v.entry[0] = second[0]-first[0];
//	line_v.entry[1] = second[1]-first[1];
//	normalize(line_v);
//	double proj_len = dot(line_v, total_v); /*here we use the dot product in the local frame*/
//
//	/*for temporal tau*/
//	//g_cur_vec = total_v;
//	//g_vec_mag = fabs(proj_len);
//}


/////////////////////////////////////////////////////////////////////////
/***********************************************************************/
/* we implement part of the routines for the Region_Tri_Index class here*/

/*calculate the Euler characteristics value*/
int Region_Tri_Index::cal_euler_value()
{
    int i, j;
    Triangle *face;
    Edge *cur_edge, **temp_edges;
    //Vertex *vert;
    int *temp_verts;
    int num_edges, num_verts;
    int MaxEdges = 3.5*ntris;
    int MaxVerts = 2*ntris;

    ////Initialize part
    num_edges = num_verts = 0;

    ////allocate memory for the edges' list and vertices' list
    temp_edges = (Edge **) malloc(sizeof(Edge *) * MaxEdges);
    temp_verts = (int *)malloc(sizeof(int) * MaxVerts);

    for(i = 0; i < ntris; i++)
    {
        if(tris[i]<0)
            continue;

        face = object->tlist.tris[tris[i]];
        ////count the number of the edges and vertices
        for(j = 0; j < 3; j++)
        {
            cur_edge = face->edges[j];
            cur_edge->visited = false;

            face->verts[j]->visited = false;

        }
    }

    int num_tris = ntris;


    for(i = 0; i < ntris; i++)
    {
        //if(tris[i]<0)  // removed by Guoning on 07/07/2010 if everything is correct, we do not need this.
        //{
        //	num_tris --;
        //	continue;
        //}

        face = object->tlist.tris[tris[i]];
        ////count the number of the edges and vertices
        for(j = 0; j < 3; j++)
        {
            cur_edge = face->edges[j];

            if(!cur_edge->visited)
            {
                ////judge and count the edge
                if(num_edges >= MaxEdges)
                {
                    MaxEdges += 100;
                    temp_edges = (Edge **)realloc(temp_edges, sizeof(Edge*) * MaxEdges);

                    if(temp_edges == nullptr)
                    {
                        //MessageBox(nullptr, "failed to reallocate memory for temp_edges!", "Error", MB_OK);
                        exit(-1);
                    }
                }
                //if(!in_edge_list(temp_edges, num_edges, cur_edge))
                //{
                ////add to the list and count the number of edges
                temp_edges[num_edges] = cur_edge;
                num_edges ++;
                cur_edge->visited=true;
                //}
            }

            ////judge and count the vertices;
            if(!face->verts[j]->visited)
            {
                if(num_verts >= MaxVerts)
                {
                    MaxVerts += 200;
                    temp_verts = (int *)realloc(temp_verts, sizeof(int ) *MaxVerts);

                    if(temp_verts == nullptr)
                    {
                        //MessageBox(nullptr, "failed to reallocate memory for temp_verts!", "Error", MB_OK);
                        exit(-1);
                    }
                }
                //if(!in_vert_list(temp_verts, num_verts, face->verts[j]->index))
                //{
                temp_verts[num_verts] = face->verts[j]->index;
                num_verts ++;
                face->verts[j]->visited=true;
                //}
            }
        }
    }

    free(temp_edges);
    free(temp_verts);

    return (num_verts + num_tris - num_edges);
}

/*judge whether this region contains fixed points or not*/
bool Region_Tri_Index::contain_fixedpts()
{
    int i;
    for(i = 0; i < ntris; i++)
    {
        if(object->tlist.tris[tris[i]]->contain_fixedpt())
            return true;
    }

    return false;
}

/*count the number of fixed points in this region*/
int Region_Tri_Index::count_num_fixedpts()
{
    int count = 0;
    int i;
    for(i = 0; i < ntris; i++)
    {
        if(object->tlist.tris[tris[i]]->contain_fixedpt())
            count++;
    }
    return count;
}


bool Region_Tri_Index::in_edge_list(Edge **edges, int num, Edge *anedge)
{
    int i;
    for(i = 0; i < num; i++)
    {
        if(edges[i] == anedge)
            return true;
    }
    return false;
}


bool Region_Tri_Index::in_vert_list(int *verts, int num, int vert_index)
{
    int i;
    for(i = 0; i < num; i++)
    {
        if(verts[i] == vert_index)
            return true;
    }
    return false;
}



void MorseDecomp::init_rot_sum()
{
    int i;

    for (i=0; i<object->vlist.nverts; i++)
    {
        object->vlist.verts[i]->total_rot_sum = 0;
    }
}


void MorseDecomp::trace_all_Verts_rot_sum(double tau, int backward)
{
    int i, k;
    int tri_id, end_tris;
    Vertex *v;
    icVector3 stP, newP;

    double lp[2];
    double alpha[3];
    Triangle *face;

    Trajectory *temp= new Trajectory(-1, 1);

    for(i = 0; i < object->vlist.nverts; i++)
    {
        v = object->vlist.verts[i];
        //v->total_rot_sum = 0;

        temp->pass_vertex(v->index, tri_id, backward);


        if(tri_id < 0)
        {
            tri_id = v->corners[0]->t;

            // Here, we are going to develop a new way to determine the starting triangle!
            temp->pass_vertex_2 (v->index, tri_id, backward);

            if (tri_id < 0)
                continue;
        }

        stP.entry[0] = v->x;
        stP.entry[1] = v->y;
        stP.entry[2] = v->z;

        if (backward == 0)
            trace_Ver_f_rot_sum(tri_id, stP.entry, newP.entry, end_tris, tau, v->total_rot_sum);
        else
            trace_Ver_b_rot_sum(tri_id, stP.entry, newP.entry, end_tris, tau, v->total_rot_sum);

        if (end_tris > -1)
        {
            int stp = 0;
        }

        v->img_tau[0] = newP.entry[0];
        v->img_tau[1] = newP.entry[1];
        v->img_tau[2] = newP.entry[2];

        v->imgtri = end_tris;

    }

    delete temp;
}


void
MorseDecomp::trace_Ver_f_rot_sum(int tri,
                                 double st[3],
                                 double en[3],
                                 int &en_tri,
                                 double tau,
                                 float &rot_sum)
{
    int i;
    int flag = 0;
    double globalp[3];
    int pre_face, cur_face;

    int steps = (int)(tau*20);

    pre_face = cur_face = tri;


    ////if the beginning point is quite close to a vertex of the triangle, move it a little bit away from the vertex

    globalp[0] = st[0];
    globalp[1] = st[1];
    globalp[2] = st[2];

    gt_tau = 0;
    Trajectory *temp = new Trajectory(-2, 1);

    for(i = 0; i < NUMTRACINGTRIANGLE*tau/*&& gt_tau < tau*/; i++)
    {

        if(cur_face < 0)
        {
            en[0] = globalp[0];
            en[1] = globalp[1];
            en[2] = globalp[2];
            en_tri = cur_face;

            delete temp;
            return;
        }

        pre_face = cur_face;

        cur_face = temp->trace_in_triangle_tau_f_rot_sum(cur_face, globalp, tau, gt_tau, rot_sum, flag);


        if(/*gt_tau >= tau ||*/ flag == 1 )
        {
            en[0] = globalp[0];
            en[1] = globalp[1];
            en[2] = globalp[2];
            en_tri = cur_face;    // 02/11/2010 modified by Guoning
            //en_tri = pre_face; /* the tracing should be stopped at previous triangle */
            delete temp;
            return;
        }
    }

    en[0] = globalp[0];
    en[1] = globalp[1];
    en[2] = globalp[2];
    en_tri = cur_face;

    delete temp;
}


void
MorseDecomp::trace_Ver_b_rot_sum(int tri,
                                 double st[3],
                                 double en[3],
                                 int &en_tri,
                                 double tau,
                                 float &rot_sum)
{
    int i;
    int flag = 0;
    double globalp[3];
    int pre_face, cur_face;

    int steps = (int)(tau*20);

    pre_face = cur_face = tri;


    ////if the beginning point is quite close to a vertex of the triangle, move it a little bit away from the vertex

    globalp[0] = st[0];
    globalp[1] = st[1];
    globalp[2] = st[2];

    gt_tau = 0;
    Trajectory *temp = new Trajectory(-2, 1);

    for(i = 0; i < NUMTRACINGTRIANGLE*tau; i++)
    {

        if(cur_face < 0)
        {
            en[0] = globalp[0];
            en[1] = globalp[1];
            en[2] = globalp[2];
            en_tri = cur_face;

            delete temp;
            return;
        }

        pre_face = cur_face;
        //cur_face = temp->trace_in_triangle_tau(cur_face, globalp, backward, tau, gt_tau, flag);

        cur_face = temp->trace_in_triangle_tau_b_rot_sum(cur_face, globalp, tau, gt_tau, rot_sum, flag);


        if(/*gt_tau >= tau ||*/ flag == 1 )
        {
            en[0] = globalp[0];
            en[1] = globalp[1];
            en[2] = globalp[2];
            en_tri = cur_face;    // 02/11/2010 modified by Guoning
            //en_tri = pre_face; /* the tracing should be stopped at previous triangle */
            delete temp;
            return;
        }
    }

    en[0] = globalp[0];
    en[1] = globalp[1];
    en[2] = globalp[2];
    en_tri = cur_face;

    delete temp;

}

//////

void MorseDecomp::trace_all_Verts_smooth(double tau, int backward, int step)
{
    int i, k;
    int tri_id, end_tris;
    Vertex *v;
    icVector3 stP, newP;

    double lp[2];
    double alpha[3];
    Triangle *face;

    Trajectory *temp= new Trajectory(-1, 1);

    for(i = 0; i < object->vlist.nverts; i++)
    {
        v = object->vlist.verts[i];
        //v->total_rot_sum = 0;

        temp->pass_vertex(v->index, tri_id, backward);


        if(tri_id < 0)
        {
            tri_id = v->corners[0]->t;

            // Here, we are going to develop a new way to determine the starting triangle!
            temp->pass_vertex_2 (v->index, tri_id, backward);

            if (tri_id < 0)
                continue;
        }

        stP.entry[0] = v->x;
        stP.entry[1] = v->y;
        stP.entry[2] = v->z;

        if (backward == 0)
            trace_Ver_f_smooth(tri_id, stP.entry, newP.entry, end_tris, tau, v->ave_sum, step);
        else
            trace_Ver_b_smooth(tri_id, stP.entry, newP.entry, end_tris, tau, v->ave_sum, step);

    }

    delete temp;
}

void
MorseDecomp::trace_Ver_f_smooth(int tri,
                                double st[3],
                                double en[3],
                                int &en_tri,
                                double tau,
                                float &rot_sum,
                                int step)
{
    int i;
    int flag = 0;
    double globalp[3];
    int pre_face, cur_face;

    pre_face = cur_face = tri;


    ////if the beginning point is quite close to a vertex of the triangle, move it a little bit away from the vertex

    globalp[0] = st[0];
    globalp[1] = st[1];
    globalp[2] = st[2];

    gt_tau = 0;
    Trajectory *temp = new Trajectory(-2, 1);

    for(i = 0; i < step; i++)
    {

        if(cur_face < 0)
        {
            en[0] = globalp[0];
            en[1] = globalp[1];
            en[2] = globalp[2];
            en_tri = cur_face;

            if (i>0)
                rot_sum /= i;
            delete temp;
            return;
        }

        pre_face = cur_face;

        cur_face = temp->trace_in_triangle_tau_f_smooth(cur_face, globalp, tau, gt_tau, rot_sum, flag);


        if(flag == 1 )
        {
            en[0] = globalp[0];
            en[1] = globalp[1];
            en[2] = globalp[2];
            en_tri = cur_face;    // 02/11/2010 modified by Guoning
            //en_tri = pre_face; /* the tracing should be stopped at previous triangle */
            delete temp;

            i++;
            rot_sum /= i;
            return;
        }
    }

    en[0] = globalp[0];
    en[1] = globalp[1];
    en[2] = globalp[2];
    en_tri = cur_face;


    rot_sum /= i;

    delete temp;
}

void
MorseDecomp::trace_Ver_b_smooth(int tri,
                                double st[3],
                                double en[3],
                                int &en_tri,
                                double tau,
                                float &rot_sum,
                                int step)
{
    int i;
    int flag = 0;
    double globalp[3];
    int pre_face, cur_face;

    int steps = (int)(tau*20);

    pre_face = cur_face = tri;


    ////if the beginning point is quite close to a vertex of the triangle, move it a little bit away from the vertex

    globalp[0] = st[0];
    globalp[1] = st[1];
    globalp[2] = st[2];

    gt_tau = 0;
    Trajectory *temp = new Trajectory(-2, 1);

    for(i = 0; i < step; i++)
    {

        if(cur_face < 0)
        {
            en[0] = globalp[0];
            en[1] = globalp[1];
            en[2] = globalp[2];
            en_tri = cur_face;

            if (i>0)
                rot_sum /= i;

            delete temp;
            return;
        }

        pre_face = cur_face;
        //cur_face = temp->trace_in_triangle_tau(cur_face, globalp, backward, tau, gt_tau, flag);

        cur_face = temp->trace_in_triangle_tau_b_smooth(cur_face, globalp, tau, gt_tau, rot_sum, flag);


        if(flag == 1 )
        {
            en[0] = globalp[0];
            en[1] = globalp[1];
            en[2] = globalp[2];
            en_tri = cur_face;    // 02/11/2010 modified by Guoning
            //en_tri = pre_face; /* the tracing should be stopped at previous triangle */
            delete temp;
            i++;
            rot_sum /= i;
            return;
        }
    }

    en[0] = globalp[0];
    en[1] = globalp[1];
    en[2] = globalp[2];
    en_tri = cur_face;

    rot_sum /= i;

    delete temp;

}

///////////////

int Trajectory::trace_in_triangle_tau_f_rot_sum(int &face_id, double globalp[3],
                                                double tau, double &gt_tau, float &rot_sum, int &flag)
{
    int i;
    double alpha[3];
    double cur_point[2], pre_point[2];
    double vert0[3];
    icVector3 VP, globalv;

    Triangle *face = object->tlist.tris[face_id];

    //eulerstep_scalar = face->euler_step;

    /************************************************************/
    /*
         Trick: Do NOT trace if the triangle contains a fixed points (02/11/2010)
    */
    //if (face->singularityID>=0)
    //{
    //	flag = 1;
    //	return face_id;
    //}

    Triangle *pre_f = face;

    icVector2 dis;

    ////Temporary curve point array

    ////initialize
    VP.entry[0] = globalp[0] - face->verts[0]->x;
    VP.entry[1] = globalp[1] - face->verts[0]->y;
    VP.entry[2] = globalp[2] - face->verts[0]->z;

    pre_point[0] = cur_point[0] = dot(VP, face->LX);
    pre_point[1] = cur_point[1] = dot(VP, face->LY);

    icVector2 pre_vec, cur_vec;

    float pre_ang = 0;
    float cur_ang = 0;

    std::vector<float> ang_list;

    vert0[0] = face->verts[0]->x;   ////for update the global point
    vert0[1] = face->verts[0]->y;
    vert0[2] = face->verts[0]->z;

    //globalface = face_id;

    ////////////////////////////////////////////////////
    for(i = 0; i < 100; i++)
    {
        ////1. calculate the barycentric coordinates for current point

        object->get_2D_Barycentric_Facters(face_id, cur_point[0], cur_point[1], alpha);

        ////2. if current point is inside current triangle
        if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1
            && (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1
            && (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
        {
            ////store the point into the temp curve points list

            ////////////////////////////////////////////////////////////////////
            // compute the total rotation sum
            // we can do this locally is because we assume the mesh has been properly oriented
            //if (i>0)
            //{
            cur_vec.set (cur_point[0]-pre_point[0], cur_point[1]-pre_point[1]);
            cur_ang = atan2 (cur_vec.entry[1], cur_vec.entry[0]);

            ang_list.push_back(cur_ang);

            //	if (i>1)
            //	{
            //		float ang_diff = cur_ang - pre_ang;
            //		if (ang_diff > M_PI)
            //			ang_diff -= (2*M_PI);
            //		else if (ang_diff < -M_PI)
            //			ang_diff += (2*M_PI);
            //		rot_sum += ang_diff;
            //	}
            //	pre_ang = cur_ang;
            //	pre_vec = cur_vec;
            //}
            ////////////////////////////////////////////////////////////////////

            pre_point[0] = cur_point[0];
            pre_point[1] = cur_point[1];

            //if(cal_next_point_euler1_tau(pre_point, cur_point, face_id, alpha, type))
            //if(cal_next_point_euler1_tau_nonNorm(pre_point, cur_point, face_id, alpha, type))
            //if(cal_next_point_2ndEuler_tau_nonNorm(pre_point, cur_point, face_id, alpha, type))
            //if(cal_next_point_2ndEuler_tau_Norm(pre_point, cur_point, face_id, alpha, type))
            //if(cal_next_pt_RK4_tau(pre_point, cur_point, face_id, alpha, type))
            //if (get_next_pt_tau(pre_point, cur_point, face_id, alpha, type, unsigned char(Integrator_opt)))


            //if (get_next_pt_tau_f(pre_point, cur_point, face_id, alpha, unsigned char(Integrator_opt)))
            if (get_next_pt_rotsum_f(pre_point, cur_point, face_id, alpha, unsigned char(Integrator_opt)))
            {
                ////update the global point

                globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

                //Get global coordinates of the point
                globalp[0] = vert0[0] + globalv.entry[0];
                globalp[1] = vert0[1] + globalv.entry[1];
                globalp[2] = vert0[2] + globalv.entry[2];

                /*This is more accurate method to accumlate time*/
                dis.entry[0] = cur_point[0] - pre_point[0];
                dis.entry[1] = cur_point[1] - pre_point[1];

                /*temporal tau*/
                //gt_tau+=/*length(dis)*/move_dist/cur_vec_mag;
                //gt_tau+=length(dis);//move_dist;
                //////////////g_dt += len*len/fabs(dot(dis, glob_loc_v));

                /////////////*spatial tau*/
                //////////////double len = length(dis);
                //////////////g_dt += len;

                /*   now the tau increment should be constant! 2/8/2010 */
                //gt_tau += 0.01;

                //if(gt_tau >= tau)
                //{
                //	flag = 1;
                //	//return face_id;
                //	break;
                //}

            }

            else{  ////the curve reach a singularity
                flag = 1;

                ////Store the record into global line segment array

                //return face_id;
                break;
            }
        }

        ////3. if the point is out of current triangle
        else{
            double t[2] = {0.};

            int PassVertornot = 0;

            get_next_triangle(face_id, pre_point, cur_point, t, 0, PassVertornot, alpha);

            ////update the globalpoint here (Modified on 01/30/07)
            if(PassVertornot > 0)
            {

                if (face_id < 0)
                    //return -1;
                    break;

                //we first need to know which vertex it is in the new triangle 01/30/07
                Vertex* vertid = pre_f->verts[PassVertornot-1];
                Triangle *cur_f = object->tlist.tris[face_id];
                int vert_new = 0;
                for(int k = 0; k < 3; k++)
                {
                    if(cur_f->verts[k] == vertid)
                    {
                        vert_new = k;
                        break;
                    }
                }

                alpha[vert_new]=1-0.0001;
                alpha[(vert_new+1)%3]=0.00005;
                alpha[(vert_new+2)%3]=0.00005;


                /* Get the new cur_point */
                cur_point[0] = alpha[1]*cur_f->x1+alpha[2]*cur_f->x2;
                cur_point[1] = alpha[2]*cur_f->y2;

                globalv = cur_point[0] * cur_f->LX + cur_point[1] * cur_f->LY;

                globalp[0] = cur_f->verts[0]->x + globalv.entry[0];
                globalp[1] = cur_f->verts[0]->y + globalv.entry[1];
                globalp[2] = cur_f->verts[0]->z + globalv.entry[2];
                face=cur_f;
            }

            else{
                //// transfer it to the global coordinates
                globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;

                globalp[0] = vert0[0] + globalv.entry[0];
                globalp[1] = vert0[1] + globalv.entry[1];
                globalp[2] = vert0[2] + globalv.entry[2];
            }

            ////////////////////////////////////////////////////////////////////
            // compute the total rotation sum
            // we can do this locally is because we assume the mesh has been properly oriented
            //cur_ang = atan2(cur_point[1], cur_point[0]);

            //float ang_diff = cur_ang - pre_ang;
            //if (ang_diff > M_PI)
            //	ang_diff -= (2*M_PI);
            //else if (ang_diff < -M_PI)
            //	ang_diff += (2*M_PI);
            //rot_sum += ang_diff;
            //if (i>0)
            //{
            //	cur_vec.set (cur_point[0]-pre_point[0], cur_point[1]-pre_point[1]);
            //	cur_ang = atan2 (cur_vec.entry[1], cur_vec.entry[0]);

            //	if (i>1)
            //	{
            //		float ang_diff = cur_ang - pre_ang;
            //		if (ang_diff > M_PI)
            //			ang_diff -= (2*M_PI);
            //		else if (ang_diff < -M_PI)
            //			ang_diff += (2*M_PI);
            //		rot_sum += ang_diff;
            //	}
            //	pre_ang = cur_ang;
            //	pre_vec = cur_vec;
            //}

            cur_vec.set (cur_point[0]-pre_point[0], cur_point[1]-pre_point[1]);
            cur_ang = atan2 (cur_vec.entry[1], cur_vec.entry[0]);

            ang_list.push_back(cur_ang);

            ////////////////////////////////////////////////////////////////////

            //dis.entry[0] = cur_point[0] - pre_point[0];
            //dis.entry[1] = cur_point[1] - pre_point[1];

            /*temporal tau*/
            //double temp_ratio=length(dis)/cur_vec_mag;
            //gt_tau+=length(dis)/cur_vec_mag;  // This computation is problematic  02/08/2010
            //gt_tau=length(dis);  // This computation is problematic  02/08/2010
            ////////g_dt += len*len/fabs(dot(dis, glob_loc_v));

            ///////*spatial tau*/
            ////////double len = length(dis);
            ////////g_dt += len;

            //gt_tau += 0.01;

            //if (gt_tau > tau)                 // Added by Guoning at 02/08/2010
            //{
            //	flag = 1;
            //	//globalv = pre_point[0] * pre_f->LX + pre_point[1] * pre_f->LY;
            //	//globalp[0] = pre_f->verts[0]->x + globalv.entry[0];
            //	//globalp[1] = pre_f->verts[0]->y + globalv.entry[1];
            //	//globalp[2] = pre_f->verts[0]->z + globalv.entry[2];

            //	//gt_tau-=temp_ratio;
            //	//return pre_f->index;
            //	break;
            //}

            //return face_id;
            break;
        }

    }

    // sum up the angles
    for (i=1; i<ang_list.size()-1; i++)
    {
        float ang_diff = ang_list[i+1] - ang_list[i];
        if (ang_diff > M_PI)
            ang_diff -= (2*M_PI);
        else if (ang_diff < -M_PI)
            ang_diff += (2*M_PI);
        rot_sum += ang_diff;
    }

    return face_id;
}














int Trajectory::trace_in_triangle_tau_f_smooth(int &face_id, double globalp[3],
                                               double tau, double &gt_tau, float &rot_sum, int &flag)
{
    int i;
    double alpha[3];
    double cur_point[2], pre_point[2];
    double vert0[3];
    icVector3 VP, globalv;

    Triangle *face = object->tlist.tris[face_id];


    Triangle *pre_f = face;

    icVector2 dis;

    ////Temporary curve point array

    ////initialize
    VP.entry[0] = globalp[0] - face->verts[0]->x;
    VP.entry[1] = globalp[1] - face->verts[0]->y;
    VP.entry[2] = globalp[2] - face->verts[0]->z;

    pre_point[0] = cur_point[0] = dot(VP, face->LX);
    pre_point[1] = cur_point[1] = dot(VP, face->LY);

    icVector2 pre_vec, cur_vec;

    float pre_ang = 0;
    float cur_ang = 0;

    std::vector<float> ang_list;

    vert0[0] = face->verts[0]->x;   ////for update the global point
    vert0[1] = face->verts[0]->y;
    vert0[2] = face->verts[0]->z;

    float sum_rot = 0;

    ////////////////////////////////////////////////////
    for(i = 0; i < 100; i++)
    {
        ////1. calculate the barycentric coordinates for current point

        object->get_2D_Barycentric_Facters(face_id, cur_point[0], cur_point[1], alpha);

        ////2. if current point is inside current triangle
        if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1
            && (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1
            && (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
        {

            pre_point[0] = cur_point[0];
            pre_point[1] = cur_point[1];

            // collect the rotation value value
            sum_rot += (alpha[0]*face->verts[0]->total_rot_sum
                        +alpha[1]*face->verts[1]->total_rot_sum
                        +alpha[2]*face->verts[2]->total_rot_sum);

            //if (get_next_pt_tau_f(pre_point, cur_point, face_id, alpha, unsigned char(Integrator_opt)))
            if (get_next_pt_rotsum_f(pre_point, cur_point, face_id, alpha, unsigned char(Integrator_opt)))
            {
                ////update the global point

                globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

                //Get global coordinates of the point
                globalp[0] = vert0[0] + globalv.entry[0];
                globalp[1] = vert0[1] + globalv.entry[1];
                globalp[2] = vert0[2] + globalv.entry[2];

                /*This is more accurate method to accumlate time*/
            }

            else{  ////the curve reach a singularity
                flag = 1;

                ////Store the record into global line segment array

                //return face_id;
                break;
            }
        }

        ////3. if the point is out of current triangle
        else{
            double t[2] = {0.};

            int PassVertornot = 0;

            get_next_triangle(face_id, pre_point, cur_point, t, 0, PassVertornot, alpha);

            ////update the globalpoint here (Modified on 01/30/07)
            if(PassVertornot > 0)
            {

                if (face_id < 0)
                    //return -1;
                    break;

                //we first need to know which vertex it is in the new triangle 01/30/07
                Vertex* vertid = pre_f->verts[PassVertornot-1];
                Triangle *cur_f = object->tlist.tris[face_id];
                int vert_new = 0;
                for(int k = 0; k < 3; k++)
                {
                    if(cur_f->verts[k] == vertid)
                    {
                        vert_new = k;
                        break;
                    }
                }

                alpha[vert_new]=1-0.0001;
                alpha[(vert_new+1)%3]=0.00005;
                alpha[(vert_new+2)%3]=0.00005;


                /* Get the new cur_point */
                cur_point[0] = alpha[1]*cur_f->x1+alpha[2]*cur_f->x2;
                cur_point[1] = alpha[2]*cur_f->y2;

                globalv = cur_point[0] * cur_f->LX + cur_point[1] * cur_f->LY;

                globalp[0] = cur_f->verts[0]->x + globalv.entry[0];
                globalp[1] = cur_f->verts[0]->y + globalv.entry[1];
                globalp[2] = cur_f->verts[0]->z + globalv.entry[2];
                face=cur_f;
            }

            else{
                //// transfer it to the global coordinates
                globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;

                globalp[0] = vert0[0] + globalv.entry[0];
                globalp[1] = vert0[1] + globalv.entry[1];
                globalp[2] = vert0[2] + globalv.entry[2];
            }

            break;
        }

    }

    if (i>0)
        sum_rot /= i;

    rot_sum += sum_rot;

    return face_id;
}



int Trajectory::trace_in_triangle_tau_b_smooth(int &face_id, double globalp[3],
                                               double tau, double &gt_tau, float &rot_sum, int &flag)
{
    int i;
    double alpha[3];
    double cur_point[2], pre_point[2];
    double vert0[3];
    icVector3 VP, globalv;

    Triangle *face = object->tlist.tris[face_id];


    Triangle *pre_f = face;

    icVector2 dis;

    ////Temporary curve point array

    ////initialize
    VP.entry[0] = globalp[0] - face->verts[0]->x;
    VP.entry[1] = globalp[1] - face->verts[0]->y;
    VP.entry[2] = globalp[2] - face->verts[0]->z;

    pre_point[0] = cur_point[0] = dot(VP, face->LX);
    pre_point[1] = cur_point[1] = dot(VP, face->LY);

    icVector2 pre_vec, cur_vec;

    float pre_ang = 0;
    float cur_ang = 0;

    std::vector<float> ang_list;

    vert0[0] = face->verts[0]->x;   ////for update the global point
    vert0[1] = face->verts[0]->y;
    vert0[2] = face->verts[0]->z;

    float sum_rot = 0;

    ////////////////////////////////////////////////////
    for(i = 0; i < 100; i++)
    {
        ////1. calculate the barycentric coordinates for current point

        object->get_2D_Barycentric_Facters(face_id, cur_point[0], cur_point[1], alpha);

        ////2. if current point is inside current triangle
        if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1
            && (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1
            && (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
        {

            pre_point[0] = cur_point[0];
            pre_point[1] = cur_point[1];

            // collect the rotation value value
            sum_rot += (alpha[0]*face->verts[0]->total_rot_sum
                        +alpha[1]*face->verts[1]->total_rot_sum
                        +alpha[2]*face->verts[2]->total_rot_sum);

            //if (get_next_pt_tau_b(pre_point, cur_point, face_id, alpha, unsigned char(Integrator_opt)))
            if (get_next_pt_rotsum_b(pre_point, cur_point, face_id, alpha, unsigned char(Integrator_opt)))
            {
                ////update the global point

                globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

                //Get global coordinates of the point
                globalp[0] = vert0[0] + globalv.entry[0];
                globalp[1] = vert0[1] + globalv.entry[1];
                globalp[2] = vert0[2] + globalv.entry[2];

                /*This is more accurate method to accumlate time*/
            }

            else{  ////the curve reach a singularity
                flag = 1;

                ////Store the record into global line segment array

                //return face_id;
                break;
            }
        }

        ////3. if the point is out of current triangle
        else{
            double t[2] = {0.};

            int PassVertornot = 0;

            get_next_triangle(face_id, pre_point, cur_point, t, 1, PassVertornot, alpha);

            ////update the globalpoint here (Modified on 01/30/07)
            if(PassVertornot > 0)
            {

                if (face_id < 0)
                    //return -1;
                    break;

                //we first need to know which vertex it is in the new triangle 01/30/07
                Vertex* vertid = pre_f->verts[PassVertornot-1];
                Triangle *cur_f = object->tlist.tris[face_id];
                int vert_new = 0;
                for(int k = 0; k < 3; k++)
                {
                    if(cur_f->verts[k] == vertid)
                    {
                        vert_new = k;
                        break;
                    }
                }

                alpha[vert_new]=1-0.0001;
                alpha[(vert_new+1)%3]=0.00005;
                alpha[(vert_new+2)%3]=0.00005;


                /* Get the new cur_point */
                cur_point[0] = alpha[1]*cur_f->x1+alpha[2]*cur_f->x2;
                cur_point[1] = alpha[2]*cur_f->y2;

                globalv = cur_point[0] * cur_f->LX + cur_point[1] * cur_f->LY;

                globalp[0] = cur_f->verts[0]->x + globalv.entry[0];
                globalp[1] = cur_f->verts[0]->y + globalv.entry[1];
                globalp[2] = cur_f->verts[0]->z + globalv.entry[2];
                face=cur_f;
            }

            else{
                //// transfer it to the global coordinates
                globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;

                globalp[0] = vert0[0] + globalv.entry[0];
                globalp[1] = vert0[1] + globalv.entry[1];
                globalp[2] = vert0[2] + globalv.entry[2];
            }

            break;
        }

    }

    if (i>0)
        sum_rot /= i;

    rot_sum += sum_rot;

    return face_id;
}


int Trajectory::trace_in_triangle_tau_b_rot_sum(int &face_id, double globalp[3],
                                                double tau, double &gt_tau, float &rot_sum, int &flag)
{
    int i;
    double alpha[3];
    double cur_point[2], pre_point[2];
    double vert0[3];
    icVector3 VP, globalv;

    Triangle *face = object->tlist.tris[face_id];

    //eulerstep_scalar = face->euler_step;

    /************************************************************/
    /*
         Trick: Do NOT trace if the triangle contains a fixed points (02/11/2010)
    */
    //if (face->singularityID>=0)
    //{
    //	flag = 1;
    //	return face_id;
    //}

    Triangle *pre_f = face;

    icVector2 dis;

    ////Temporary curve point array

    ////initialize
    VP.entry[0] = globalp[0] - face->verts[0]->x;
    VP.entry[1] = globalp[1] - face->verts[0]->y;
    VP.entry[2] = globalp[2] - face->verts[0]->z;

    pre_point[0] = cur_point[0] = dot(VP, face->LX);
    pre_point[1] = cur_point[1] = dot(VP, face->LY);

    icVector2 pre_vec, cur_vec;
    float pre_ang =0;
    float cur_ang =0;
    std::vector<float> ang_list;

    vert0[0] = face->verts[0]->x;   ////for update the global point
    vert0[1] = face->verts[0]->y;
    vert0[2] = face->verts[0]->z;

    //globalface = face_id;

    ////////////////////////////////////////////////////
    for(i = 0; i < 100; i++)
    {
        ////1. calculate the barycentric coordinates for current point

        object->get_2D_Barycentric_Facters(face_id, cur_point[0], cur_point[1], alpha);

        ////2. if current point is inside current triangle
        if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1
            && (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1
            && (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
        {
            ////store the point into the temp curve points list

            ////////////////////////////////////////////////////////////////////
            // compute the total rotation sum
            // we can do this locally is because we assume the mesh has been properly oriented
            //cur_ang = atan2(cur_point[1], cur_point[0]);

            //float ang_diff = cur_ang - pre_ang;
            //if (ang_diff > M_PI)
            //	ang_diff -= (2*M_PI);
            //else if (ang_diff < -M_PI)
            //	ang_diff += (2*M_PI);
            //rot_sum -= ang_diff;
            //if (i>0)
            //{
            cur_vec.set (cur_point[0]-pre_point[0], cur_point[1]-pre_point[1]);
            cur_ang = atan2 (cur_vec.entry[1], cur_vec.entry[0]);

            ang_list.push_back(cur_ang);

            //	if (i>1)
            //	{
            //		float ang_diff = cur_ang - pre_ang;
            //		if (ang_diff > M_PI)
            //			ang_diff -= (2*M_PI);
            //		else if (ang_diff < -M_PI)
            //			ang_diff += (2*M_PI);
            //		rot_sum += ang_diff;
            //	}
            //	pre_ang = cur_ang;
            //	pre_vec = cur_vec;
            //}
            ////////////////////////////////////////////////////////////////////

            pre_point[0] = cur_point[0];
            pre_point[1] = cur_point[1];

            //if(cal_next_point_euler1_tau(pre_point, cur_point, face_id, alpha, type))
            //if(cal_next_point_euler1_tau_nonNorm(pre_point, cur_point, face_id, alpha, type))
            //if(cal_next_point_2ndEuler_tau_nonNorm(pre_point, cur_point, face_id, alpha, type))
            //if(cal_next_point_2ndEuler_tau_Norm(pre_point, cur_point, face_id, alpha, type))
            //if(cal_next_pt_RK4_tau(pre_point, cur_point, face_id, alpha, type))
            //if (get_next_pt_tau(pre_point, cur_point, face_id, alpha, type, unsigned char(Integrator_opt)))

            //if (get_next_pt_tau_b(pre_point, cur_point, face_id, alpha, unsigned char(Integrator_opt)))
            if (get_next_pt_rotsum_b(pre_point, cur_point, face_id, alpha, unsigned char(Integrator_opt)))
            {
                ////update the global point

                globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

                //Get global coordinates of the point
                globalp[0] = vert0[0] + globalv.entry[0];
                globalp[1] = vert0[1] + globalv.entry[1];
                globalp[2] = vert0[2] + globalv.entry[2];

                /*This is more accurate method to accumlate time*/
                dis.entry[0] = cur_point[0] - pre_point[0];
                dis.entry[1] = cur_point[1] - pre_point[1];

                /*temporal tau*/
                //gt_tau+=/*length(dis)*/move_dist/cur_vec_mag;
                //gt_tau+=length(dis)/*move_dist*/;
                //////////////g_dt += len*len/fabs(dot(dis, glob_loc_v));

                /////////////*spatial tau*/
                //////////////double len = length(dis);
                //////////////g_dt += len;

                /*   now the tau increment should be constant! 2/8/2010 */
                //gt_tau += 0.01;

                //if(gt_tau >= tau)
                //{
                //	flag = 1;
                //	//return face_id;
                //	break;
                //}
            }

            else{  ////the curve reach a singularity
                flag = 1;

                ////Store the record into global line segment array

                //return face_id;
                break;
            }
        }

        ////3. if the point is out of current triangle
        else{
            double t[2] = {0.};

            int PassVertornot = 0;

            get_next_triangle(face_id, pre_point, cur_point, t, 1, PassVertornot, alpha);

            ////update the globalpoint here (Modified on 01/30/07)
            if(PassVertornot > 0)
            {
                if (face_id < 0)
                    //return -1;
                    break;

                //we first need to know which vertex it is in the new triangle 01/30/07
                Vertex* vertid = pre_f->verts[PassVertornot-1];
                Triangle *cur_f = object->tlist.tris[face_id];
                int vert_new = 0;
                for(int k = 0; k < 3; k++)
                {
                    if(cur_f->verts[k] == vertid)
                    {
                        vert_new = k;
                        break;
                    }
                }

                alpha[vert_new]=1-0.0001;
                alpha[(vert_new+1)%3]=0.00005;
                alpha[(vert_new+2)%3]=0.00005;


                /* Get the new cur_point */
                cur_point[0] = alpha[1]*cur_f->x1+alpha[2]*cur_f->x2;
                cur_point[1] = alpha[2]*cur_f->y2;

                globalv = cur_point[0] * cur_f->LX + cur_point[1] * cur_f->LY;

                globalp[0] = cur_f->verts[0]->x + globalv.entry[0];
                globalp[1] = cur_f->verts[0]->y + globalv.entry[1];
                globalp[2] = cur_f->verts[0]->z + globalv.entry[2];
                face=cur_f;
            }

            else{
                //// transfer it to the global coordinates
                globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;

                globalp[0] = vert0[0] + globalv.entry[0];
                globalp[1] = vert0[1] + globalv.entry[1];
                globalp[2] = vert0[2] + globalv.entry[2];
            }

            dis.entry[0] = cur_point[0] - pre_point[0];
            dis.entry[1] = cur_point[1] - pre_point[1];

            ////////////////////////////////////////////////////////////////////
            // compute the total rotation sum
            // we can do this locally is because we assume the mesh has been properly oriented
            //cur_ang = atan2(cur_point[1], cur_point[0]);

            //float ang_diff = cur_ang - pre_ang;
            //if (ang_diff > M_PI)
            //	ang_diff -= (2*M_PI);
            //else if (ang_diff < -M_PI)
            //	ang_diff += (2*M_PI);
            //rot_sum -= ang_diff;
            //if (i>0)
            //{
            //cur_vec.set (cur_point[0]-pre_point[0], cur_point[1]-pre_point[1]);
            cur_ang = atan2 (dis.entry[1], dis.entry[0]);

            ang_list.push_back(cur_ang);

            //	if (i>1)
            //	{
            //		float ang_diff = cur_ang - pre_ang;
            //		if (ang_diff > M_PI)
            //			ang_diff -= (2*M_PI);
            //		else if (ang_diff < -M_PI)
            //			ang_diff += (2*M_PI);
            //		rot_sum += ang_diff;
            //	}
            //	pre_ang = cur_ang;
            //	pre_vec = cur_vec;
            //}
            ////////////////////////////////////////////////////////////////////

            /*temporal tau*/
            //double temp_ratio=length(dis)/cur_vec_mag;
            //gt_tau+=length(dis)/cur_vec_mag;  // This computation is problematic  02/08/2010
            //gt_tau+=length(dis);  // This computation is problematic  02/08/2010
            ////////g_dt += len*len/fabs(dot(dis, glob_loc_v));

            ///////*spatial tau*/
            ////////double len = length(dis);
            ////////g_dt += len;

            //gt_tau += 0.01;

            //if (gt_tau > tau)                 // Added by Guoning at 02/08/2010
            //{
            //	flag = 1;
            //	//globalv = pre_point[0] * pre_f->LX + pre_point[1] * pre_f->LY;
            //	//globalp[0] = pre_f->verts[0]->x + globalv.entry[0];
            //	//globalp[1] = pre_f->verts[0]->y + globalv.entry[1];
            //	//globalp[2] = pre_f->verts[0]->z + globalv.entry[2];

            //	//gt_tau-=temp_ratio;
            //	//return pre_f->index;
            //	break;
            //}

            //return face_id;
            break;
        }

    }

    // sum up the angles
    for (i=1; i<ang_list.size()-1; i++)
    {
        float ang_diff = ang_list[i+1] - ang_list[i];
        if (ang_diff > M_PI)
            ang_diff -= (2*M_PI);
        else if (ang_diff < -M_PI)
            ang_diff += (2*M_PI);
        rot_sum -= ang_diff;
    }

    return face_id;
}


void
MorseDecomp::find_out_min_max_rot_sum()
{
    int i;
    object->min_rot_sum = 1.e8;
    object->max_rot_sum = -1.e8;

    for (i=0; i<object->vlist.nverts; i++)
    {
        if (object->vlist.verts[i]->total_rot_sum < object->min_rot_sum)
            object->min_rot_sum = object->vlist.verts[i]->total_rot_sum;
        if (object->vlist.verts[i]->total_rot_sum > object->max_rot_sum)
            object->max_rot_sum = object->vlist.verts[i]->total_rot_sum;
    }
}


void
MorseDecomp::smooth_A_field()
{
    // we try to remove the bubble in the A field

    int i, j;
    for (i=0; i<object->vlist.nverts; i++)
    {
        Vertex *cv = object->vlist.verts[i];

        // see whether this vertex is a local minimum or maximum
        int min_counter =0, max_counter=0;
        float ave_sum = 0;
        for (j=0; j<cv->nedges; j++)
        {
            Edge *e = cv->edges[j];

            Vertex *nv = e->verts[0];

            if (nv == cv)
                nv = e->verts[1];

            if (nv->total_rot_sum > cv->total_rot_sum) min_counter ++;

            if (nv->total_rot_sum < cv->total_rot_sum) max_counter ++;

            ave_sum += nv->total_rot_sum;
        }

        if (min_counter == cv->nedges || max_counter == cv->nedges)
        {
            cv->total_rot_sum = ave_sum/cv->nedges;
        }
    }
}

void
MorseDecomp::LIC_like_smoothing()
{
    int L=2; // we perform a forward and backward integration and accumulate the rotation value via interpolation...
    int i, j;
    for (i=0; i<object->vlist.nverts; i++)
    {
        object->vlist.verts[i]->ave_sum = 0;
    }

    //forward tracing
    trace_all_Verts_smooth(0.1, 0, L);

    //backward tracing
    trace_all_Verts_smooth(0.1, 1, L);

    for (i=0; i<object->vlist.nverts; i++)
    {
        object->vlist.verts[i]->total_rot_sum = object->vlist.verts[i]->ave_sum;
    }
}
