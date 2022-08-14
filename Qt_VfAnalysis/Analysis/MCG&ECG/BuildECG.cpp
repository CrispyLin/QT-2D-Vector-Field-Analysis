/*
This file contains routines of class MCG_Graph that we use to build the ECG graph based on
the obtain Morse decomposition

Created and Modified by Guoning Chen
        copyright @2007
    */


#include "Analysis/MorseDecomp.h"
#include "Others/common_routines.h"
#include "Analysis/PeriodicOrbitExtraction.h"


extern Polyhedron *object;
extern PeriodicOrbitList *periodic_orbits;
extern MorseDecomp *morse_decomp;
extern TrajectoryList *separatrices;
extern PeriodicOrbit_Detector *periodicorbit_detector ;

/*the data structure for ECG*/
ECG_Graph *ecg = NULL;

extern double ConleyGraphWin_x ;
extern double ConleyGraphWin_y ;
extern double types_interval_y ;
extern const int DebugOn;

/*
*/
int Trajectory::trace_in_triangle_connection(double g[3], int &face_id, int type, int &flag)
{
    //similar to the regular tracing, except that you need to mark those special points too close to the cycle
    int i;
    double alpha[3];
    double cur_point[2], pre_point[2];
    double vert0[3];
    icVector3 VP, globalv;
    icVector2 dis;  //to calculate the distance for each line segment  08/10/06

    if(face_id < 0)
        return -1;

    Triangle *face = object->tlist.tris[face_id];

    Triangle *pre_f = face;

    ////initialize
    VP.entry[0] = g[0] - face->verts[0]->x;
    VP.entry[1] = g[1] - face->verts[0]->y;
    VP.entry[2] = g[2] - face->verts[0]->z;

    pre_point[0] = cur_point[0] = dot(VP, face->LX);
    pre_point[1] = cur_point[1] = dot(VP, face->LY);

    vert0[0] = face->verts[0]->x;   ////for update the global point
    vert0[1] = face->verts[0]->y;
    vert0[2] = face->verts[0]->z;

    ////////////////////////////////////////////////////
    for(i = 0; i < 300; i++)
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

            if(cal_next_point_euler1(pre_point, cur_point, face_id, alpha, type))
            //if(cal_nextpt_2ndeuler(pre_point, cur_point, face_id, alpha, type))
            //if(cal_nextpt_RK4(pre_point, cur_point, face_id, alpha, type))
            {
                ////update the global point
                //face = Object.flist[face_id];

                globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

                g[0] = vert0[0] + globalv.entry[0];
                g[1] = vert0[1] + globalv.entry[1];
                g[2] = vert0[2] + globalv.entry[2];

                ////08/10/06
                dis.entry[0] = cur_point[0] - pre_point[0];
                dis.entry[1] = cur_point[1] - pre_point[1];

                //sum_flow_length += length(dis);
            }

            else{  ////the curve reaches a singularity
                flag = 1;
                return face_id;
            }

        }

        ////3. if the point is out of current triangle
        else{
            double t[2] = {0.};

            int PassVertornot = 0;

            int presave_face = face_id;

            //cross_a_vertex(face_id, cur_point, pre_point, type, PassVertornot);
            get_next_triangle(face_id, pre_point, cur_point, t, type, PassVertornot, alpha);

            ////update the global point here (Big modification on 01/30/07)
            if(PassVertornot > 0)
            {
                ////we should not directly use the vertex as next point!!
                if (face_id < 0)
                {
                    return -1;
                }
                ////we may move a little bit along the VF direction, but make sure it is still inside
                ////the new triangle
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

                g[0] = cur_f->verts[0]->x + globalv.entry[0];
                g[1] = cur_f->verts[0]->y + globalv.entry[1];
                g[2] = cur_f->verts[0]->z + globalv.entry[2];
                face=cur_f;

                ////08/10/06 get the length of the last line segment
                dis.entry[0] = cur_point[0] - pre_point[0];
                dis.entry[1] = cur_point[1] - pre_point[1];
                //sum_flow_length += length(dis);


            }


            else{

                //// transfer it to the global coordinates
                globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;

                g[0] = vert0[0] + globalv.entry[0];
                g[1] = vert0[1] + globalv.entry[1];
                g[2] = vert0[2] + globalv.entry[2];

                ////08/10/06 get the length of the last line segment
                dis.entry[0] = cur_point[0] - pre_point[0];
                dis.entry[1] = cur_point[1] - pre_point[1];
                //sum_flow_length += length(dis);
            }

            return face_id;
        }
    }

    return face_id;
}



void Singularity::update_list_to_PO(int limitcycle)
{
    ////1.Extend current list
    ////2.Add to the end of current list
    ////3.Update the counter

    int i;
    int *temp_list = connect_cycles;

    if(num_connect_cycles == 0 || temp_list == NULL)
    {
        connect_cycles = (int*)malloc(sizeof(int));
        connect_cycles[0] = limitcycle;
        num_connect_cycles = 1;
    }

    else{
        connect_cycles = (int*)malloc(sizeof(int)*(num_connect_cycles+1));

        for(i = 0; i < num_connect_cycles; i++)
        {
            connect_cycles[i] = temp_list[i];
        }

        connect_cycles[i] = limitcycle;

        num_connect_cycles += 1;

        free(temp_list);
    }
}


void PeriodicOrbit::update_list_to_fixedpt(int singID)
{
    ////Extend current list	and add to the end of current list

    int i;
    int *temp_list = connected_fixedpts;
    //double *temp_length = flow_length_connectedsing;

    if(temp_list == NULL)
    {
        connected_fixedpts = (int*)malloc(sizeof(int));
        connected_fixedpts[0] = singID;
        num_connect_fixedpts = 1;

        ////Allocate space to save the flow length of the connection
        //flow_length_connectedsing =
        //	(double*)malloc(sizeof(double)*	num_connect_fixedpts); //08/10/06
        //flow_length_connectedsing[0] = sum_flow_length;  //08/10/06
    }

    else{
        connected_fixedpts =
            (int*)malloc(sizeof(int)*(num_connect_fixedpts+1));

        //flow_length_connectedsing =
        //	(double*)malloc(sizeof(double)*(num_connect_fixedpts+1)); //08/10/06

        for(i = 0; i < num_connect_fixedpts; i++)
        {
            connected_fixedpts[i] = temp_list[i];
            //flow_length_connectedsing[i] = temp_length[i]; //08/10/06
        }

        connected_fixedpts[i] = singID;

        //flow_length_connectedsing[i] = sum_flow_length;  //08/10/06

        num_connect_fixedpts += 1;
        free(temp_list);
        //free(temp_length);
    }
}


void PeriodicOrbit::update_list_to_cycle(int cycle)
{
    ////1. Extend current list
    ////2. Add to the end of current list
    ////3. Update the counter

    int i;
    int *temp_list = connected_POs;
    //double *temp_length = flow_length_connectedcycle;

    if(temp_list == NULL)
    {
        connected_POs = (int*)malloc(sizeof(int));
        connected_POs[0] = cycle;
        num_connect_POs = 1;

        ////08/14/06 store the flow length of this connection
        //flow_length_connectedcycle = (double*)malloc(sizeof(double));
        //flow_length_connectedcycle[0] = sum_flow_length;
    }

    else{
        connected_POs =
            (int*)malloc(sizeof(int)*(num_connect_POs+1));

        ////08/10/06 store the flow length of this connection
        //flow_length_connectedcycle =
        //	(double*)malloc(sizeof(double)*(num_connect_POs+1));

        for(i = 0; i < num_connect_POs; i++)
        {
            connected_POs[i] = temp_list[i];
            //flow_length_connectedcycle[i] = temp_length[i]; //08/10/06
        }

        connected_POs[i] = cycle;

        //flow_length_connectedcycle[i] = sum_flow_length;  ////08/10/06
        num_connect_POs += 1;

        free(temp_list);
        //free(temp_length);
    }
}



/*member functions for MCG construction and modification*/
/*build the mcg according to the obtained Morse sets*/
void ECG_Graph::build_ecg()
{
    assign_ecgnodes();
    layout_ecg();

    if(periodic_orbits != NULL)
        build_connections();

    build_ecg_edges();
}


/*
Assign the ECG nodes according to the obtained fixed points and periodic orbits
*/
void ECG_Graph::assign_ecgnodes()
{
    int i;
    ////Assign one node for each singularity and limit cycle
    r_counter = a_counter = s_counter = 0;
    cur_ecgnode_index = 0;

    for(i = 0; i < object->slist.nsingularities; i++)
    {
        if(object->slist.slist[i]->type != CWCENTER && object->slist.slist[i]->type != CCWCENTER)
        {
            nlist->enodes[cur_ecgnode_index]->node_index = cur_ecgnode_index;
            nlist->enodes[cur_ecgnode_index]->fixedPtID = i;
            nlist->enodes[cur_ecgnode_index]->periodicorbitID = -1;
            nlist->enodes[cur_ecgnode_index]->cancelled = false;

            object->slist.slist[i]->node_index = cur_ecgnode_index;

            if(object->slist.slist[i]->type == SOURCE || object->slist.slist[i]->type == RFOCUS)
            {
                nlist->enodes[cur_ecgnode_index]->type = 0;   ////it is a repeller
                nlist->enodes[cur_ecgnode_index]->labelindex = r_counter;
                r_counter++;
            }
            else if(object->slist.slist[i]->type == SINK || object->slist.slist[i]->type == AFOCUS)
            {
                nlist->enodes[cur_ecgnode_index]->type = 1;   ////it is a attractor
                nlist->enodes[cur_ecgnode_index]->labelindex = a_counter;
                a_counter++;
            }
            else
            {
                nlist->enodes[cur_ecgnode_index]->type = 2;   ////it is a saddle
                nlist->enodes[cur_ecgnode_index]->labelindex = s_counter;
                s_counter++;
            }

            nlist->enodes[cur_ecgnode_index]->nedges = 0;
            nlist->enodes[cur_ecgnode_index]->graph_edges = NULL;

            cur_ecgnode_index ++;
            nlist->nenodes ++;
        }
    }

    //FILE *fp;
    //	if(DebugOn == 1)
    //	{
    //		fp = fopen("detect_process_log.txt", "a");
    //		fprintf(fp, "finish assigning fixed points for ECG.\n");
    //		fclose(fp);
    //	}

    if(periodic_orbits == NULL)
        return;

    ////periodic orbits
    for(i = 0; i < periodic_orbits->nporbits; i++)
    {
        nlist->enodes[cur_ecgnode_index]->node_index = cur_ecgnode_index;
        nlist->enodes[cur_ecgnode_index]->periodicorbitID = i;
        nlist->enodes[cur_ecgnode_index]->fixedPtID = -1;
        nlist->enodes[cur_ecgnode_index]->cancelled = false;
        periodic_orbits->polist[i]->node_index = cur_ecgnode_index;

        if(periodic_orbits->polist[i]->type == 0)  ////it is a repeller
        {
            nlist->enodes[cur_ecgnode_index]->type = 0;
            nlist->enodes[cur_ecgnode_index]->labelindex = r_counter;
            r_counter++;
        }
        else                          ////it is an attractor
        {
            nlist->enodes[cur_ecgnode_index]->type = 1;
            nlist->enodes[cur_ecgnode_index]->labelindex = a_counter;
            a_counter++;
        }

        nlist->enodes[cur_ecgnode_index]->nedges = 0;
        nlist->enodes[cur_ecgnode_index]->graph_edges = NULL;

        nlist->nenodes ++;
        cur_ecgnode_index ++;
    }
}


void ECG_Graph::layout_ecg()
{
    int i;
    int num_repellers, num_attractors, num_saddles;
    int cur_repeller_index, cur_attractor_index, cur_saddle_index;
    double repeller_interval, attractor_interval, saddle_interval;
    double repeller_y, attractor_y, saddle_y;


    num_repellers = num_attractors = num_saddles = 0;
    cur_repeller_index = cur_attractor_index = cur_saddle_index = 0;

    ////Actually, they may be counted during the building of the graph
    ////1. First, count the number of the nodes of 3 types respectively
    for(i = 0; i < cur_ecgnode_index; i++)
    {
        if(nlist->enodes[i]->type == 0)
            num_repellers ++;
        else if(nlist->enodes[i]->type == 1)
            num_attractors ++;
        else
            num_saddles ++;
    }

    repeller_interval = ConleyGraphWin_x / (num_repellers+1);
    attractor_interval = ConleyGraphWin_x / (num_attractors+1);
    saddle_interval = ConleyGraphWin_x / (num_saddles+1);

    attractor_y = types_interval_y;
    saddle_y = 2 * types_interval_y;
    repeller_y = 3 * types_interval_y;


    ////2. Then, evenly space the positions of these nodes
    for(i = 0; i < cur_ecgnode_index; i++)
    {
        if(nlist->enodes[i]->type == 0)
        {
            nlist->enodes[i]->pos.entry[0] = (cur_repeller_index+1)*repeller_interval;
            nlist->enodes[i]->pos.entry[1] = repeller_y;
            cur_repeller_index ++;
        }
        else if(nlist->enodes[i]->type == 1)
        {
            nlist->enodes[i]->pos.entry[0] = (cur_attractor_index+1)*attractor_interval;
            nlist->enodes[i]->pos.entry[1] = attractor_y;
            cur_attractor_index ++;
        }
        else{
            nlist->enodes[i]->pos.entry[0] = (cur_saddle_index+1)*saddle_interval;
            nlist->enodes[i]->pos.entry[1] = saddle_y;
            cur_saddle_index ++;
        }
    }
}


void ECG_Graph::add_to_ecg_nodelist(int type, int &index, int fixedorcycle)
{
    if(index >= nlist->curMaxNumENodes)
    {
        int oldnum = nlist->curMaxNumENodes;

        /*extend the memory*/
        if(!nlist->extend(50))
        {
            /*report errors*/
            return;
        }

        for(int i = oldnum; i < nlist->curMaxNumENodes; i++)
            nlist->enodes[i] = new ECG_Node();
    }

    nlist->enodes[index]->node_index = index;
    //nlist->enodes[index]->scc_index = scc_index;
    nlist->enodes[index]->graph_edges = NULL;
    nlist->enodes[index]->nedges = 0;
    nlist->enodes[index]->type = type;
    nlist->enodes[index]->visited = false;
    nlist->enodes[index]->cancelled = false;


    /*------------------------------------------------------------------*/
    /*the inaccurate judgement could happen for larger Morse sets!!!!*/

    /*------------------------------------------------------------------*/

    if( nlist->enodes[index]->type == 0)
    {
        nlist->enodes[index]->labelindex = r_counter;
        r_counter++ ;
    }
    else if( nlist->enodes[index]->type == 1)
    {
        nlist->enodes[index]->labelindex = a_counter;
        a_counter++;
    }
    else
    {
        nlist->enodes[index]->labelindex = s_counter;
        s_counter++;
    }

    nlist->nenodes++;
    index++;
}


void ECG_Graph::add_to_ecg_edgelist(int node1, int node2)
{
    if(elist->isFull())
    {
        int oldnum = elist->curMaxNumGedges;
        if(!elist->extend())
        {
            /*report errors*/
            return;
        }

        /*
        we need to add new space for the new allocate pointers
        */
        for(int i = oldnum; i < elist->curMaxNumGedges; i++)
        {
            elist->edges[i] = new Graph_Edge();
        }
    }

    elist->edges[cur_ecgedge_index]->edge_index = cur_ecgedge_index;
    elist->edges[cur_ecgedge_index]->node_index1 = node1;
    elist->edges[cur_ecgedge_index]->node_index2 = node2;
    elist->edges[cur_ecgedge_index]->cancel = false;
    elist->edges[cur_ecgedge_index]->visited = false;

    elist->nedges++;
    cur_ecgedge_index++;
}


void ECG_Graph::add_edge_to_node(int node, int edgeindex)
{
    if(node < 0)
        return;

    nlist->enodes[node]->graph_edges = extend_link(nlist->enodes[node]->graph_edges,
                                                   nlist->enodes[node]->nedges);

    nlist->enodes[node]->graph_edges[nlist->enodes[node]->nedges] = edgeindex;
    nlist->enodes[node]->nedges++;
}


/*
Here, we assume that we find out all the connections first.
*/
void ECG_Graph::build_ecg_edges()
{
    int i, j;
    Trajectory *cur_traj;
    int lineseg_id, last_triangle;
    int end_sing;
    int node1, node2;
    cur_ecgedge_index = 0;

    //FILE *fp;

    // Build the edges
    ////The graph is a directed graph!!!! 11/19/05, the direction of an edge is node1->node2
    ////i.e. node1 should be a repeller and node2 should be an attractor
    ////Note that, during the building of each edge, we need to check whether it is a previous edge!

    /*(a) Find the edges between saddles and other singularities or limit cycles */
    for(i = 0; i < object->slist.nsingularities; i++)
    {

        ////Begin from each saddle, check its separatrices
        ////if the last triangle of the separatrix contains singularity, build an edge between the singularity and current saddle
        ////Note that we need to allow double connections here 08/14/06
        if(object->slist.slist[i]->type == SADDLE)
        {
            ////if the last triangle of the separatrix does not contain any singularity
            ////it maybe forms a closed streamline (we do not consider the 'center' here now) or reach the boundary
            ////Build an edge between the limit cycle and current saddle

            ////Then, follow the separatrices of the saddle
            for(j = 0; j < 4; j++)
            {
                switch(j)
                {
                case 0:
                    cur_traj = separatrices->trajs[object->slist.slist[i]->separatrices];
                    break;
                case 1:
                    cur_traj = separatrices->trajs[object->slist.slist[i]->separatrices+1];
                    break;
                case 2:
                    cur_traj = separatrices->trajs[object->slist.slist[i]->separatrices+2];
                    break;
                case 3:
                    cur_traj = separatrices->trajs[object->slist.slist[i]->separatrices+3];
                    break;
                }

                lineseg_id = cur_traj->nlinesegs;
                if(lineseg_id <= 0) continue;
                last_triangle = cur_traj->linesegs[lineseg_id-1].Triangle_ID;
                end_sing = object->tlist.tris[last_triangle]->singularityID;
                /*or if its neighboring triangle contains fixed point very close to
                the shared edges of the two triangles*/


                if(last_triangle >= 0 && end_sing >= 0) ////contains singularity
                {
                    ////Get the node index for the other singularity
                    node1 = object->slist.slist[i]->node_index;
                    node2 = object->slist.slist[end_sing]->node_index;

                    /////////////////////////////

                    if(object->slist.slist[end_sing]->type == SOURCE
                        ||object->slist.slist[end_sing]->type == RFOCUS)
                    {
                        ////here, saddle is the attractor
                        node2 = node1;
                        node1 = object->slist.slist[end_sing]->node_index;
                    }

                    ////Add to graph edge array
                    add_to_ecg_edgelist(node1, node2);

                    ////Save the edge length (the length of the separatrix) to the new added edge
                    ////08/10/06, here i happens to be the index of the saddle, j is the index of the sep of the saddle
                    //graphedges[cur_graphedge_index-1].flow_length = GetLength_separatrix(i, j+1);

                    ////Need to add to the edge list of node1 and node2
                    add_edge_to_node(node1, elist->nedges-1);
                    add_edge_to_node(node2, elist->nedges-1);

                    ////Store the separatrix information to the edge  08/16/06
                    //graphedges[cur_graphedge_index-1].trajID = trajectory_id;
                }
            }

            if(periodic_orbits == NULL)
                continue;

            /*(b) build the connections between saddles and limit cycles */
            for(j = 0; j < object->slist.slist[i]->num_connect_cycles; j++)
            {
                node1 = object->slist.slist[i]->node_index;
                node2 = periodic_orbits->polist[object->slist.slist[i]->connect_cycles[j]]->node_index;


                if(periodic_orbits->polist[object->slist.slist[i]->connect_cycles[j]]->type == 0)
                {
                    node2 = node1;
                    node1 = periodic_orbits->polist[object->slist.slist[i]->connect_cycles[j]]->node_index;
                }


                ////If there is no connection between them, connect them
                /* We should allow doubly connected here 02/12/07 */
                //if(!IsConnected(node1, node2))
                //{
                ////Add to graph edge array
                add_to_ecg_edgelist(node1, node2);

                ///*following codes are for flow distance calculation for automatic simplification*/
                //int k;
                //for(k = 0; k < periodic_orbits->polist[object->slist.slist[i]->connect_cycles[j]]->num_connect_fixedpts;
                //	k++)
                //{
                //	if(periodic_orbits->polist[object->slist.slist[i]->connect_cycles[j]]->connected_fixedpts[k]==i)
                //		break;
                //}
                //graphedges[cur_graphedge_index-1].flow_length =
                //	periodic_orbits->polist[object->slist.slist[i]->connect_cycles[j]]->flow_length_connectedsing[k];

                ////Need to add to the edge list of node1 and node2
                add_edge_to_node(node1, cur_ecgedge_index-1);
                add_edge_to_node(node2, cur_ecgedge_index-1);
                //}
            }
        }
    }

    /*(c) Find the edges between limit cycles and other singularities */
    //////Here we build the saddle-limit cycle connections first
    //////which can be fulfilled in previous step 2/15/06

    //////Here we begin the searching from the connected list in limit cycle data structure
    //////Note that we also need to consider the directions of the connections
    if(periodic_orbits == NULL)
        return;

    for(i = 0; i < periodic_orbits->nporbits; i++)
    {
        node1 = periodic_orbits->polist[i]->node_index;

        for(j = 0; j < periodic_orbits->polist[i]->num_connect_fixedpts; j++)
        {
            node2 = object->slist.slist[periodic_orbits->polist[i]->connected_fixedpts[j]]->node_index;

            if(periodic_orbits->polist[i]->type == 1)
            {
                ////this is an attractive limit cycle, hence saddle becomes a repeller here
                node1 = node2;
                node2 = periodic_orbits->polist[i]->node_index;
            }

            if(has_interval(node1, node2))  ////there is saddel connecting them
                continue;
            if(is_connected(node1, node2)) ////avoid repeating  08/10/06
                continue;

            ////Add to graph edge array
            add_to_ecg_edgelist(node1, node2);

            ////Need to add to the edge list of node1 and node2
            add_edge_to_node(node1, cur_ecgedge_index-1);
            add_edge_to_node(node2, cur_ecgedge_index-1);

            ////we need to save the edge length here
            ////We need to first find out the location of the singularity in the list
            //elist->edges[cur_ecgedge_index-1]->flow_length =
            //	periodic_orbits->polist[i]->flow_length_connectedsing[j];
        }
    }

    /*(c) Find the edges between embeded limit cycles */
    for(i = 0; i < periodic_orbits->nporbits; i++)
    {
        node1 = periodic_orbits->polist[i]->node_index;

        for(j = 0; j < periodic_orbits->polist[i]->num_connect_POs; j++)
        {
            ////if there has already one edge between them, continue
            node2 = periodic_orbits->polist[periodic_orbits->polist[i]->connected_POs[j]]->node_index;

            if(periodic_orbits->polist[periodic_orbits->polist[i]->connected_POs[j]]->type ==
                periodic_orbits->polist[i]->type)
                continue;

            if(periodic_orbits->polist[i]->type == 1)
            {
                ////this is an attractive limit cycle, hence the other becomes a repeller here
                node1 = node2;
                node2 = periodic_orbits->polist[i]->node_index;
            }

            if(is_connected(node1, node2)) ////avoid repeating  08/10/06
                continue;

            if(has_interval(node1, node2))  ////there is saddel connecting them
                continue;

            add_to_ecg_edgelist(node1, node2);

            add_edge_to_node(node1, cur_ecgedge_index-1);
            add_edge_to_node(node2, cur_ecgedge_index-1);

            ////store the edge length
            //elist->edges[cur_ecgedge_index-1]->flow_length =
            //	periodic_orbits->polist[i]->flow_length_connectedsing[j];
        }
    }
}



bool ECG_Graph::has_interval(int node1, int node2)
{
    int i, j;
    ////we need to make sure node1 is the repeller
    if(nlist->enodes[node1]->type == 1)
    {
        int temp = node1;
        node1 = node2;
        node2 = temp;
    }


    ////Begin to judge whether there is at least one interval connecting them
    int *temp_connect_saddlelist = new int[nlist->enodes[node1]->nedges];
    int num_connect_saddlelist = 0;
    int cur_edge, n2;

    ////Begin searching from this repeller to search all the connected saddles
    for(j = 0; j < nlist->enodes[node1]->nedges; j++)
    {
        cur_edge = nlist->enodes[node1]->graph_edges[j];
        n2 = elist->edges[cur_edge]->node_index2;

        if(nlist->enodes[n2]->type == 2) //if this is a saddle
        {
            ////we need to perform repeating test if we allow multiple connections between two nodes!!
            temp_connect_saddlelist[num_connect_saddlelist] = n2;
            num_connect_saddlelist++;
        }
    }

    for(j = 0; j < num_connect_saddlelist; j++)
    {
        for(i = 0; i < nlist->enodes[temp_connect_saddlelist[j]]->nedges; i++)
        {
            cur_edge = nlist->enodes[temp_connect_saddlelist[j]]->graph_edges[i];
            n2 = elist->edges[cur_edge]->node_index2;

            if(n2 == node2) //we find an interval
            {
                free(temp_connect_saddlelist);
                return true;
            }
        }
    }

    free(temp_connect_saddlelist);
    return false;
}


bool ECG_Graph::is_connected(int node1, int node2)
{
    int i, cur_edge;

    for(i = 0; i < nlist->enodes[node1]->nedges; i++)
    {
        cur_edge = nlist->enodes[node1]->graph_edges[i];

        if(elist->edges[cur_edge]->node_index1 == node1)
        {
            if(elist->edges[cur_edge]->node_index2 == node2)
                return true;
        }

        else{
            if(elist->edges[cur_edge]->node_index1 == node2)
                return true;
        }
    }

    return false;
}


void ECG_Graph::trace_from_saddle_to_PO(int saddle)
{
    double sing_center[3], newpos[3];
    icVector3 sep_vector;
    int start_triangle;
    Triangle *face = object->tlist.tris[object->slist.slist[saddle]->TriangleID];

    sing_center[0] = object->slist.slist[saddle]->gpos.entry[0];
    sing_center[1] = object->slist.slist[saddle]->gpos.entry[1];
    sing_center[2] = object->slist.slist[saddle]->gpos.entry[2];

    //follow the four separatrices respectively

    //positive outgoing
    sep_vector = object->slist.slist[saddle]->loutgoing.entry[0]*face->LX
                 + object->slist.slist[saddle]->loutgoing.entry[1]*face->LY;

    start_triangle = object->slist.slist[saddle]->TriangleID;
    separatrices->cal_startpt_sep(start_triangle, sep_vector, sing_center, newpos);

    trace_and_build_connection(newpos, start_triangle, saddle, 0);

    //positive incoming
    sep_vector = object->slist.slist[saddle]->lincoming.entry[0]*face->LX
                 + object->slist.slist[saddle]->lincoming.entry[1]*face->LY;

    start_triangle = object->slist.slist[saddle]->TriangleID;
    separatrices->cal_startpt_sep(start_triangle, sep_vector, sing_center, newpos);

    trace_and_build_connection(newpos, start_triangle, saddle, 1);

    //negative outgoing
    sep_vector = - object->slist.slist[saddle]->loutgoing.entry[0]*face->LX
                 - object->slist.slist[saddle]->loutgoing.entry[1]*face->LY;

    start_triangle = object->slist.slist[saddle]->TriangleID;
    separatrices->cal_startpt_sep(start_triangle, sep_vector, sing_center, newpos);

    trace_and_build_connection(newpos, start_triangle, saddle, 0);

    //negative incoming
    sep_vector = - object->slist.slist[saddle]->lincoming.entry[0]*face->LX
                 - object->slist.slist[saddle]->lincoming.entry[1]*face->LY;

    start_triangle = object->slist.slist[saddle]->TriangleID;
    separatrices->cal_startpt_sep(start_triangle, sep_vector, sing_center, newpos);

    trace_and_build_connection(newpos, start_triangle, saddle, 1);
}



void ECG_Graph::trace_from_PO_to_PO(double s[3], int triangle, int cycleID, int type)
{
    int i;
    int flag = 0;
    double globalp[3];
    int pre_face, cur_face;

    pre_face = cur_face = triangle;
    //sum_flow_length = 0.;

    globalp[0] = s[0];   globalp[1] = s[1];  globalp[2] = s[2];
    int NUMTRACETRIS = (int)sqrt((double)object->tlist.ntris);
    Trajectory *temp_traj = new Trajectory(-1);

    for(i = 0; i < 10*NUMTRACETRIS; i++)
    {

        if(cur_face == -1)
        {
            return ;
        }

        pre_face = cur_face;

        ////Here we need a sub routine that can return the intersection and the corresponding edge
        ////when it enters a new triangle

        cur_face = temp_traj->trace_in_triangle_connection(globalp, cur_face, type, flag);

        ////Move the repeated cycle judgement here 07/25/06
        int pre_limit_index = -1;
        if(periodicorbit_detector->is_existing_porbit(cur_face, pre_limit_index, 1-type))
        {
            //build connection here
            if(pre_limit_index < 0)
                continue;

            periodic_orbits->polist[pre_limit_index]->update_list_to_cycle(cycleID);
            periodic_orbits->polist[cycleID]->update_list_to_cycle(pre_limit_index);
            delete temp_traj;
            return;
        }

        if(flag == 1 /*|| pre_face == cur_face*/)
        {
            ////build the connection with the center singularity

            int singID = object->tlist.tris[cur_face]->singularityID;

            /*if singID < 0, we need to find the neighboring fixed point 05/11/07*/

            if(singID >=0 && object->slist.slist[singID]->type!= SADDLE)
            {
                //build the connection between the singularity and the limit cycle
                object->slist.slist[singID]->update_list_to_PO(cycleID);
                periodic_orbits->polist[cycleID]->update_list_to_fixedpt(singID);
                object->slist.slist[singID]->connected = true;

                /*do we need to set the central singularity of a periodic orbit*/
                //limitcycles[cycleID].singularID = singID;  ////if we set the flag correctly, this will be fine!
                flag = 1;      //fail, reach singularity or other error, no cycle has been found on this path
                delete temp_traj;
                return;
            }
        }
    }
}


void ECG_Graph::trace_and_build_connection(double s[3], int triangle, int singID, int type)
{
    int i;
    int flag = 0;
    double globalp[3];
    int pre_face, cur_face;

    pre_face = cur_face = triangle;
    //sum_flow_length = 0.;

    Trajectory *temp_traj = new Trajectory(-1);

    globalp[0] = s[0];   globalp[1] = s[1];  globalp[2] = s[2];
    int NUMTRACETRIS = (int)sqrt((double)object->tlist.ntris);

    for(i = 0; i < 6*NUMTRACETRIS; i++)
    {

        if(cur_face == -1)
        {
            delete temp_traj;
            return ;
        }

        pre_face = cur_face;

        ////Here we need a sub routine that can return the intersection and the corresponding edge
        ////when it enters a new triangle

        cur_face = temp_traj->trace_in_triangle_connection(globalp, cur_face, type, flag);

        ////Move the repeated cycle judgement here 07/25/06
        int pre_limit_index = -1;
        if(periodicorbit_detector->is_existing_porbit(cur_face, pre_limit_index, 1-type))
        {
            //build connection here
            if(pre_limit_index < 0)
                continue;

            object->slist.slist[singID]->update_list_to_PO(pre_limit_index);
            periodic_orbits->polist[pre_limit_index]->update_list_to_fixedpt(singID);

            ////Update the singularity "connected" component
            object->slist.slist[singID]->connected = true;

            delete temp_traj;
            return;
        }

        if(flag == 1 /*|| pre_face == cur_face*/)
        {
            flag = 1;      //fail, reach singularity or other error, no cycle has been found on this path

            delete temp_traj;
            return;
        }
    }

    delete temp_traj;
}



void ECG_Graph::trace_connection_cycle(int cycle)
{
    //find out the edge contain the fixed point
    double fixed_p[2], gfixed_p[3], v1[3], v2[3], out[3], out2[3], alpha[3];
    int triangle = periodic_orbits->polist[cycle]->traj->linesegs[0].Triangle_ID;
    Triangle *face = object->tlist.tris[triangle];

    //get the local point to calculate the barycentric coordinates
    fixed_p[0] = periodic_orbits->polist[cycle]->traj->linesegs[0].start.entry[0];
    fixed_p[1] = periodic_orbits->polist[cycle]->traj->linesegs[0].start.entry[1];

    object->get_2D_Barycentric_Facters(triangle, fixed_p[0], fixed_p[1], alpha);

    //get the global coordinates of the fixed point
    gfixed_p[0] = periodic_orbits->polist[cycle]->traj->linesegs[0].gstart.entry[0];
    gfixed_p[1] = periodic_orbits->polist[cycle]->traj->linesegs[0].gstart.entry[1];
    gfixed_p[2] = periodic_orbits->polist[cycle]->traj->linesegs[0].gstart.entry[2];


    //Do we need to consider the case that fixed point locates at a vertex ?
    if(abs(alpha[0])<1e-8) //fixed point is on the edge v1v2
    {
        v1[0] = face->verts[1]->x;
        v1[1] = face->verts[1]->y;
        v1[2] = face->verts[1]->z;

        v2[0] = face->verts[2]->x;
        v2[1] = face->verts[2]->y;
        v2[2] = face->verts[2]->z;
    }

    else if(abs(alpha[1])<1e-8) //fixed point iss on the edge v2v0
    {
        v1[0] = face->verts[2]->x;
        v1[1] = face->verts[2]->y;
        v1[2] = face->verts[2]->z;

        v2[0] = face->verts[0]->x;
        v2[1] = face->verts[0]->y;
        v2[2] = face->verts[0]->z;
    }

    else{  //fixed point iss on the edge v0v1
        v1[0] = face->verts[0]->x;
        v1[1] = face->verts[0]->y;
        v1[2] = face->verts[0]->z;

        v2[0] = face->verts[1]->x;
        v2[1] = face->verts[1]->y;
        v2[2] = face->verts[1]->z;
    }

    //get the first point
    out[0] = (v1[0]+gfixed_p[0])/2.;
    out[1] = (v1[1]+gfixed_p[1])/2.;
    out[2] = (v1[2]+gfixed_p[2])/2.;

    //get the second point
    out2[0] = (v2[0]+gfixed_p[0])/2.;
    out2[1] = (v2[1]+gfixed_p[1])/2.;
    out2[2] = (v2[2]+gfixed_p[2])/2.;

    ////Adjust the point selection
    icVector3 distvec1, distvec2;
    double dist1, dist2;
    distvec1.entry[0] = v1[0] - gfixed_p[0];
    distvec1.entry[1] = v1[1] - gfixed_p[1];
    distvec1.entry[2] = v1[2] - gfixed_p[2];
    dist1 = length(distvec1);
    distvec2.entry[0] = v2[0] - gfixed_p[0];
    distvec2.entry[1] = v2[1] - gfixed_p[1];
    distvec2.entry[2] = v2[2] - gfixed_p[2];
    dist2 = length(distvec2);

    if(dist1 > dist2)
    {
        out[0] = gfixed_p[0] - distvec2.entry[0]*0.5;
        out[1] = gfixed_p[1] - distvec2.entry[1]*0.5;
        out[2] = gfixed_p[2] - distvec2.entry[2]*0.5;
    }

    else
    {
        out2[0] = gfixed_p[0] - distvec1.entry[0]*0.5;
        out2[1] = gfixed_p[1] - distvec1.entry[1]*0.5;
        out2[2] = gfixed_p[2] - distvec1.entry[2]*0.5;
    }

    trace_from_PO_to_PO(out, triangle, cycle, periodic_orbits->polist[cycle]->type);

    trace_from_PO_to_PO(out2, triangle, cycle, periodic_orbits->polist[cycle]->type);
}


void ECG_Graph::build_connections()
{

    time_t rawtime;
    clock_t start, finish;
    FILE *fp ;

    time ( &rawtime );

    start = clock();

    //fp = fopen("detect_porbit_cooling.txt", "a");
    //fp = fopen("detect_porbit_swirl.txt", "a");
    if(DebugOn == 1){
        fp = fopen("detect_porbit.txt", "a");
        fprintf(fp, "start finding the other connections. \n");
        fprintf(fp, "current date and time are :  %s. \n", ctime (&rawtime) );
        fclose(fp);
    }
    ///////////////////////////////////////////////

    ///////////////////////////////////////

    ///////////////////////////////////////////

    ////////////////////////////////////////////

    int i;

    for(i = 0; i < periodic_orbits->nporbits; i++)
    {
        trace_connection_cycle(i);
    }

    int limitcycle;

    for(i = 0; i < object->slist.nsingularities; i++)
    {
        if(object->slist.slist[i]->type == SOURCE || object->slist.slist[i]->type == RFOCUS
            ||object->slist.slist[i]->type == SINK || object->slist.slist[i]->type == AFOCUS)
        {
            trace_from_nonsaddle_to_PO(i);
        }

        else if(object->slist.slist[i]->type == SADDLE)
        {
            trace_from_saddle_to_PO(i);
        }
    }
    ///////////////////////////////////////////////

    time ( &rawtime );

    finish = clock();

    //fp = fopen("detect_porbit_cooling.txt", "a");
    //fp = fopen("detect_porbit_swirl.txt", "w");
    if(DebugOn == 1){
        fp = fopen("detect_porbit.txt", "a");
        fprintf(fp, "finish building the remaining connections.\n");
        fprintf(fp, "current date and time are :  %s. \n", ctime (&rawtime) );
        fprintf(fp, "time spent on building all connections is %f. \n", (double)(finish - start)/CLOCKS_PER_SEC);
        fclose(fp);
    }
}



void ECG_Graph::trace_from_nonsaddle_to_PO(int singID)
{
    double sing_c[3], out[3];

    if(object->slist.slist[singID]->connected)
        return;

    //first, find a point close to the center of the singularity
    sing_c[0] = object->slist.slist[singID]->gpos.entry[0];
    sing_c[1] = object->slist.slist[singID]->gpos.entry[1];
    sing_c[2] = object->slist.slist[singID]->gpos.entry[2];
    Triangle *t = object->tlist.tris[object->slist.slist[singID]->TriangleID];

    t->cal_a_close_pt(sing_c, out);

    if(object->slist.slist[singID]->type == SOURCE || object->slist.slist[singID]->type == RFOCUS)
    {
        //we need to trace forward
        trace_and_build_connection(out, object->slist.slist[singID]->TriangleID, singID, 0);
    }
    else if(object->slist.slist[singID]->type == SINK || object->slist.slist[singID]->type == AFOCUS)
    {
        //we need to trace backward
        trace_and_build_connection(out, object->slist.slist[singID]->TriangleID, singID, 1);
    }
}
