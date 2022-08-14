/*
This file contains routines for periodic orbit extraction
Classes involved in this file are PeriodicOrbit, PeriodicOrbit_Detector, Polyhedron, Edge,
SCComponent
Created and Modified by Guoning Chen
copyright @2007
*/

#include "Analysis/PeriodicOrbitExtraction.h"


extern MorseDecomp *morse_decomp;
extern Polyhedron *object;
extern const int DebugOn;
extern int ndisplay_POs;


/*We now separate the periodic orbit list from the geometry surface*/
PeriodicOrbitList *periodic_orbits = nullptr;
PeriodicOrbit_Detector *periodicorbit_detector = nullptr;

/**/
extern time_t g_rawtime;
extern clock_t g_start, g_finish;

extern int Integrator_opt;


Edge *g_chosen_edge;  /*for some reason, I have to use global edge to save the selected edge*/

/*
The main routine for periodic orbit extraction
Need to record some information for the performance evaluation
*/
void Polyhedron::detect_PeriodicOrbit()
{
    /*first, we need to perform Morse decomposition*/
    if(morse_decomp == nullptr)
        morse_decomp = new MorseDecomp();

    FILE *fp;
    clock_t l_start, l_finish;

    l_start = clock();
    morse_decomp->morse_decomp();
    l_finish = clock();

    //fp = fopen("detect_porbit_cooling.txt", "a");
    //fp = fopen("detect_porbit_swirl.txt", "a");
    //fp = fopen("detect_porbit_tumble.txt", "a");

    if(DebugOn == 1){
        fp = fopen("detect_porbit.txt", "a");
        fprintf(fp, "time for computing Morse Decomposition is %f...\n",
                (double)(l_finish - l_start)/CLOCKS_PER_SEC);
        fclose(fp);

        _cprintf("time for computing Morse Decomposition is %f...\n",
                 (double)(l_finish - l_start)/CLOCKS_PER_SEC);
    }

    /*second, calculate the separation and attachment points for the valid SCCs*/
    l_start = clock();
    morse_decomp->cal_sepsandatts_valid_SCCs();
    l_finish = clock();

    //fp = fopen("detect_porbit_cooling.txt", "a");
    //fp = fopen("detect_porbit_swirl.txt", "a");
    //fp = fopen("detect_porbit_tumble.txt", "a");
    if(DebugOn == 1){
        fp = fopen("detect_porbit.txt", "a");
        fprintf(fp, "time for computing special points is %f...\n",
                (double)(l_finish - l_start)/CLOCKS_PER_SEC);
        fclose(fp);

        _cprintf("time for computing special points is %f...\n",
                 (double)(l_finish - l_start)/CLOCKS_PER_SEC);
    }

    /*third, start periodic orbit extraction*/
    //object->p
    //if(periodic_orbits != nullptr)
    //	delete periodic_orbits;

    /*initialize the periodic orbit list*/
    if(periodic_orbits == nullptr)
    {
        periodic_orbits = new PeriodicOrbitList(100);

        for(int i = 0; i < periodic_orbits->curMaxNumPOrbits; i++)
            periodic_orbits->polist[i] = new PeriodicOrbit();
    }

    periodic_orbits->reset_all_porbits();

    periodic_orbits->nporbits = 0;

    if(periodicorbit_detector == nullptr)
        periodicorbit_detector = new PeriodicOrbit_Detector();


    //_cprintf("start detecting periodic orbits.\n");

    l_start = clock();
    periodicorbit_detector->detect_periodicorbit_morsedecomp();
    l_finish = clock();

    //fp = fopen("detect_porbit_cooling.txt", "a");
    //fp = fopen("detect_porbit_swirl.txt", "a");
    //fp = fopen("detect_porbit_tumble.txt", "a");

    if(DebugOn == 1){
        fp = fopen("detect_porbit.txt", "a");
        fprintf(fp, "time for extracting all the periodic orbit is %f...\n",
                (double)(l_finish - l_start)/CLOCKS_PER_SEC);
        fclose(fp);
    }
}


/*
Judge whether a separation/attachment point on an edge is a valid point or not. It is valid, If it satisfies:
1) not close to the being visited separation point in the same SCC
2) not close to the previously detected limit cycle
We assume that the "cur_e" contain a separation point
Note: we need to set a threshold here. Currently, the threshold is relative to the average edge length
or the Object.radius, we need to make sure it is less than one triangle width
*/

bool Edge::is_valid_sep(int tri, int scc_index)
{
    // The following setting seems should be put in the tracing code as well

    int *NearbyTriangles = nullptr;
    int num_triangles = 0;
    int i, j;
    Triangle *face;
    Edge *cur_e;
    double t_s[3] = {sep->x, sep->y, sep->z};

    NearbyTriangles = get_disc_coarse(t_s, triangle, NearbyTriangles, num_triangles);

    for(i = 0; i < num_triangles; i++)
    {
        //Do not change the states of those special points not belonging to current SCC
        if(!is_repeated_elem(morse_decomp->scclist->scccomponents[scc_index]->nodes, NearbyTriangles[i],
                              morse_decomp->scclist->scccomponents[scc_index]->nnodes))
            continue;

        face = object->tlist.tris[NearbyTriangles[i]];

        for(j = 0; j < 3; j++)
        {
            cur_e = face->edges[j];

            if(!cur_e->valid || !cur_e->find_sep || cur_e == this)
                continue;

            if(cur_e->sep_visit == 1)  //it is a dead point
            {
                this->sep_visit = 2;
                return false;
            }

        }
    }

    free(NearbyTriangles);
    return true;
}


bool Edge::is_valid_attp(int tri, int scc_index)
{
    int *NearbyTriangles = nullptr;
    int num_triangles = 0;
    int i, j;
    Triangle *face;
    Edge *cur_e;

    //NearbyTriangles = GetDisc(the_e->attp.entry, triangle, ave_length, 1, NearbyTriangles, num_triangles);

    for(i = 0; i < num_triangles; i++)
    {
        //Do not change the states of those special points not belonging to current SCC
        if(!is_repeated_elem(morse_decomp->scclist->scccomponents[scc_index]->nodes, NearbyTriangles[i],
                              morse_decomp->scclist->scccomponents[scc_index]->nnodes))
            continue;

        face = object->tlist.tris[NearbyTriangles[i]];

        for(j = 0; j < 3; j++)
        {
            cur_e = face->edges[j];

            if(!cur_e->valid || !cur_e->find_attp  || cur_e == this)
                continue;

            if(cur_e->att_visit == 1)  //it is a dead point
            {
                this->att_visit = 2;
                return false;
            }

        }
    }

    free(NearbyTriangles);
    return true;
}


void Edge::reset_intersections()
{
    num_intersections = 0;
    intersections[0].set(0, 0, 0);
    intersections[1].set(0, 0, 0);
}

/*
It will be called if the same edge has been intersected consectively, which may corresponding
to a tangent point on the edge.
This routine is used in periodic orbit extraction
*/
bool Edge::approx_tangent(double seg_length, icVector3 new_intersect)
{
    icVector3 approx_dis;
    if(num_intersections > 0)
    {
        approx_dis = new_intersect - intersections[0];

        if(abs(::length(approx_dis)-seg_length) < 1e-4)
            return true;
    }

    return false;
}

/*
Judge whether the two intersections on the same edge means finding a cycle
*/
bool Edge::is_good_intersections(int triangle, int scc_index)
{
    //if the line segments between the two intersections is small

    //if(morse_decomp->scclist->scccomponents[scc_index]->num_boundaries > 1 &&
    //	IsMixedEdge(cur_e, triangle))
    //	return false;               //we probably may miss the tiny cycle in a ring shaped SCC 07/22/06

    ////if the vector values on the two intersections are not almost parallel to each other

    ////first, we need to get the interpolated vector values on the two intersections of the edge
    //icVector2 vec1, vec2;

    //////Get the local flow vector under the local frame of "triangle"
    //vec1 = GetVecOnIntersection_local(cur_e->intersections[0].entry, triangle);
    //vec2 = GetVecOnIntersection_local(cur_e->intersections[1].entry, triangle);
    //

    //double cross_r = vec1.entry[0]*vec2.entry[1] - vec1.entry[1]*vec2.entry[0];

    //if(abs(cross_r) > 1e-8)   //modified at 08/15/06
    //	return false;

    return true;
}









/*-----------------------------------------------------------------------------------*/


void SCComponent::reset_edge_intersections()
{
    int i, j;
    Triangle *face;

    for(i = 0; i < nnodes; i++)
    {
        face = object->tlist.tris[i];

        for(j = 0; j < 3; j++)
        {
            face->edges[j]->reset_intersections();
        }
    }
}


/*------------------------------------------------------------------------------------*/
/*
*/
void PeriodicOrbit_Detector::detect_periodicorbit_morsedecomp()
{
    int i, j;

    ////Calculate the edge length here as the distance threshold
    //icVector3 temp;
    //temp.entry[0] = object->tlist.tris[0]->verts[0]->x -
    //	object->tlist.tris[0]->verts[1]->x;
    //temp.entry[1] = object->tlist.tris[0]->verts[0]->y -
    //	object->tlist.tris[0]->verts[1]->y;
    //temp.entry[2] = object->tlist.tris[0]->verts[0]->z -
    //	object->tlist.tris[0]->verts[1]->z;

    //ave_length = length(temp);
    /////////////////////////////////////////////////////////////

    ////File testing...
    g_start = clock();
    FILE *fp;

    //initialize the intersections on all edges
    Triangle *face;
    Edge *e;
    for(i = 0; i < object->tlist.ntris; i++)
    {
        face = object->tlist.tris[i];
        for(j = 0; j < 3; j++)
        {
            e = face->edges[j];
            e->num_intersections = 0;
            e->intersections[0].set(0, 0, 0);
            e->intersections[1].set(0, 0, 0);
        }
    }

    ndisplay_POs = 0;

    for(i = 0; i < morse_decomp->scclist->nsccs; i++)
    {

        if(!morse_decomp->scclist->scccomponents[i]->valid)
            continue;

        //_cprintf("start from Morse set %d.\n", i);

        //if (i==13822)
        //{
        //	int stop = 0;
        //}

        ////Find the possible attracting limit cycles in this SCC
        if(morse_decomp->scclist->scccomponents[i]->nattpts == 0)
        {
            ////use regular method to trace
            start_regular(i, 1);   //1 means attracting
        }

        else
        {
            start_seporatt_pts(i, 1);
        }

        ////Find the possible repelling limit cycles in this SCC
        if(morse_decomp->scclist->scccomponents[i]->nseppts == 0)
        {
            ////use regular method to trace
            start_regular(i, 0);   //0 means repelling
        }

        else
        {
            start_seporatt_pts(i, 0);
        }
    }

    //#ifdef WRITETIME
    //	finish = clock();
    //	double timespent = (double)(finish - start)/CLOCKS_PER_SEC;
    //	fp = fopen("detectprocess.txt", "a");
    //	fprintf(fp, "detection spent %f to detect %d limit cycles...\n", timespent, cur_limitcycle_index);
    //	fclose(fp);
    //#endif
}

void PeriodicOrbit_Detector::start_regular(int scc_index, int type)
{
}


/*
start tracing from attachment or separation points to detect periodic orbits
*/
void PeriodicOrbit_Detector::start_seporatt_pts(int scc_index, int type)
{
    //search all the edges inside the SCC[scc_index]
    //if we find an edge containing sep or att point, validate it
    //if it is a valid point, start tracing from it
    //1)for a sep point, trace backward
    //2)for a att point, trace forward

    int i, j;
    Triangle *face;
    Edge *cur_e;
    int flag = -1;
    chosen_edge = nullptr;
    int the_triangle;
    double t_s[3];


    for(i = 0; i < morse_decomp->scclist->scccomponents[scc_index]->nnodes; i++)
    {
        face = object->tlist.tris[morse_decomp->scclist->scccomponents[scc_index]->nodes[i]];

        //FILE *fp = fopen("detect_porbit.txt", "a");
        //fprintf(fp, "start %d Morse set at triangle %d.\n", scc_index,
        //	morse_decomp->scclist->scccomponents[scc_index]->nodes[i]);
        //fclose(fp);

        for(j = 0; j < 3; j++)
        {
            cur_e = face->edges[j];

            if(!cur_e->valid) //not valid edge
                continue;

            if(type == 0)
            {
                if(!cur_e->find_sep || cur_e->sep == nullptr)  //no separation point
                    continue;

                if(cur_e->sep_visit == 1) //this is a dead point
                    continue;

                if(!cur_e->is_valid_sep(face->index, scc_index)) //this is not a good separation point
                    continue;

                //reset the intersection information before starting tracing
                morse_decomp->scclist->scccomponents[scc_index]->reset_edge_intersections();
                flag = -1;

                //start tracing backward to locate repelling limit cycle if existed
                the_triangle = face->index;
                t_s[0] = cur_e->sep->x;
                t_s[1] = cur_e->sep->y;
                t_s[2] = cur_e->sep->z;

                //fp = fopen("detect_porbit.txt", "a");
                //fprintf(fp, "In %d Morse set, start tracing.\n", i);
                //fclose(fp);

                if(trace_to_find_recurrence(t_s, the_triangle, 1, scc_index, chosen_edge, flag))
                {
                    //we found a cycle, but we still need to return the intersection and the edge
                    //for calculating the closed streamline for the limit cycle

                    //fp = fopen("detect_porbit.txt", "a");
                    //fprintf(fp, "In %d Morse set, we find a recurrence.\n", i);
                    //fclose(fp);

                    //chosen_edge = g_chosenedge;
                    if(periodic_orbits->isFull())
                    {

                        if(!periodic_orbits->extend())
                        {
                            char rout[256], var[256];
                            sprintf(rout, "%s", "PeriodicOrbit_Detector::start_seporatt_pts");
                            sprintf(var, "%s", "periodic_orbits");

                            //write_mem_error(rout, var, 0);
                            return;
                        }


                        /*initialize the new elements*/
                        for(i = periodic_orbits->nporbits; i < periodic_orbits->curMaxNumPOrbits; i++)
                        {
                            periodic_orbits->polist[i] = new PeriodicOrbit();
                        }

                    }


                    if(cal_closedstreamline(chosen_edge, the_triangle, 1, scc_index))
                    {
                        //store the information for the periodic orbit
                        store_current_cellcycle(cellcycle->elems, cellcycle->nelems);

                        periodic_orbits->polist[periodic_orbits->nporbits]->index = periodic_orbits->nporbits;
                        periodic_orbits->polist[periodic_orbits->nporbits]->node_index = -1;

                        //store_cur_streamline(trajectories[cur_traj_index], num_linesegs_curtraj[_index]);

                        ////store the type of the limit cycle
                        periodic_orbits->polist[periodic_orbits->nporbits]->type = 0;

                        ////Initialize the connected list for the graph
                        periodic_orbits->polist[periodic_orbits->nporbits]->connected_POs = nullptr;
                        periodic_orbits->polist[periodic_orbits->nporbits]->num_connect_POs = 0;
                        periodic_orbits->polist[periodic_orbits->nporbits]->connected_fixedpts = nullptr;
                        periodic_orbits->polist[periodic_orbits->nporbits]->num_connect_fixedpts = 0;
                        //periodic_orbits->polist[periodic_orbits->nporbits]->connected_l = 0;
                        //periodic_orbits->polist[periodic_orbits->nporbits]->connected_r = 0;
                        periodic_orbits->nporbits++;

                        ndisplay_POs = periodic_orbits->nporbits-1;


                        ///*Testing code*/
                        //FILE *fp = fopen("nporbits.txt", "w");
                        //fprintf(fp, "current number of detected periodic orbits is :%d\n",
                        //	periodic_orbits->nporbits);
                        //fclose(fp);
                    }
                }

                ////probably we need to mark its neighboring separation points 07/25/06
                t_s[0] = cur_e->sep->x;
                t_s[1] = cur_e->sep->y;
                t_s[2] = cur_e->sep->z;
                mark_neighbor_pts(t_s, face->index, scc_index, 0);

                //set the point as dead
                cur_e->sep_visit = 1;
            }

            else
            {
                if(!cur_e->find_attp || cur_e->attp == nullptr)  //no attachment point
                    continue;

                if(cur_e->att_visit == 1) //this is a dead point
                    continue;

                if(!cur_e->is_valid_attp(face->index, scc_index)) //this is not a good attachment point
                    continue;


                morse_decomp->scclist->scccomponents[scc_index]->reset_edge_intersections();

                //start tracing forward to locate attracting limit cycle if existed
                the_triangle = face->index;
                t_s[0] = cur_e->attp->x;
                t_s[1] = cur_e->attp->y;
                t_s[2] = cur_e->attp->z;
                if(trace_to_find_recurrence(t_s, the_triangle, 0, scc_index, chosen_edge, flag))
                {
                    //we found a cycle, but we still need to return the intersection and the edge
                    //for calculating the closed streamline for the limit cycle
                    if(periodic_orbits->isFull())
                    {
                        //FILE *fp = fopen("extend_polist.txt", "w");
                        //fprintf(fp, "start calling the extend routine.\n");
                        //fclose(fp);

                        if(!periodic_orbits->extend())
                        {
                            char rout[256], var[256];
                            sprintf(rout, "%s", "PeriodicOrbit_Detector::start_seporatt_pts");
                            sprintf(var, "%s", "periodic_orbits");

                            //write_mem_error(rout, var, 0);
                            return;
                        }

                        //FILE *fp = fopen("extend_polist.txt", "w");
                        //fprintf(fp, "start allocating memeory for new periodic orbits\n");
                        //fprintf(fp, "Current maximum number of POs in the list is %d.\n",
                        //	periodic_orbits->nporbits);
                        //fclose(fp);

                        /*initialize the new elements*/
                        for(i = periodic_orbits->nporbits; i < periodic_orbits->curMaxNumPOrbits; i++)
                        {
                            periodic_orbits->polist[i] = new PeriodicOrbit();
                        }
                        //fp = fopen("extend_polist.txt", "w");
                        //fprintf(fp, "polist has been extended.\n");
                        //fprintf(fp, "Current maximum number of POs in the list is %d.\n",
                        //	periodic_orbits->nporbits);
                        //fclose(fp);

                    }


                    //chosen_edge = g_chosenedge;
                    if(cal_closedstreamline(chosen_edge, the_triangle, 0, scc_index))
                    {
                        //store the information for the periodic orbit
                        store_current_cellcycle(cellcycle->elems, cellcycle->nelems);

                        periodic_orbits->polist[periodic_orbits->nporbits]->index = periodic_orbits->nporbits;
                        periodic_orbits->polist[periodic_orbits->nporbits]->node_index = -1; ////10/16/05

                        //store_cur_streamline(trajectories[cur_traj_index], num_linesegs_curtraj[cur_traj_index]);

                        ////store the type of the limit cycle
                        periodic_orbits->polist[periodic_orbits->nporbits]->type = 1;

                        ////Initialize the connected list for the graph
                        periodic_orbits->polist[periodic_orbits->nporbits]->connected_POs = nullptr;
                        periodic_orbits->polist[periodic_orbits->nporbits]->num_connect_POs = 0;
                        periodic_orbits->polist[periodic_orbits->nporbits]->connected_fixedpts = nullptr;
                        periodic_orbits->polist[periodic_orbits->nporbits]->num_connect_fixedpts = 0;
                        //periodic_orbits->polist[periodic_orbits->nporbits]->connected_l = 0;
                        //periodic_orbits->polist[periodic_orbits->nporbits]->connected_r = 0;
                        periodic_orbits->nporbits++;

                        ndisplay_POs = periodic_orbits->nporbits-1;

                        ///*Testing code*/
                        //FILE *fp = fopen("nporbits.txt", "w");
                        //fprintf(fp, "current number of detected periodic orbits is :%d\n",
                        //	periodic_orbits->nporbits);
                        //fclose(fp);
                    }
                }

                ////probably we need to mark its neighboring separation points 07/25/06
                t_s[0] = cur_e->attp->x;
                t_s[1] = cur_e->attp->y;
                t_s[2] = cur_e->attp->z;
                mark_neighbor_pts(t_s, face->index, scc_index, 1);

                //set the point as dead
                cur_e->att_visit = 1;
            }
        }
    }
}


bool PeriodicOrbit_Detector::trace_to_find_recurrence(double g[3], int &triangle,
                                                      int type, int scc_index, Edge *chosen_edge,  int &flag)
{
    detect_traj = new Trajectory(-1);

    int i;
    flag = 0;
    double globalp[3];
    int pre_face, cur_face;

    pre_face = cur_face = triangle;

    globalp[0] = g[0];   globalp[1] = g[1]; globalp[2] = g[2];


    Edge *the_e = nullptr;
    icVector3 intersect;
    chosen_edge = nullptr;
    Edge *pre_e = nullptr;

    double seg_length = 0;

    cellcycle->nelems = 0;
    cellcycle->add_New(triangle);

    int NUMTRACETRIS = (int)sqrt((double)object->tlist.ntris);

    for(i = 0; i < 50*NUMTRACETRIS && i < object->tlist.ntris; i++)
    {

        if(cur_face == -1)
        {
            flag = 1;       //fail, reaches boundary of the mesh, no cycle exists on this path
            return false;
        }

        pre_face = cur_face;

        ////Here we need a sub routine that can return the intersection and the corresponding edge
        ////when it enters a new triangle

        pre_e = the_e;

        cur_face = trace_in_triangle_po_detect(globalp, cur_face, type, the_e, intersect, seg_length, flag);


        if(flag == 1 || pre_face == cur_face)
        {
            flag = 1;      //fail, reach singularity or other error, no cycle has been found on this path
            return false;
        }

        if(!is_repeated_elem(morse_decomp->scclist->scccomponents[scc_index]->nodes, cur_face,
                              morse_decomp->scclist->scccomponents[scc_index]->nnodes))
        {
            flag = 1;     //fail, goes out of current SCC
            return false;
        }

        /* we need to judge whether it forms a loop! 02/08/07 */
        if(!cellcycle->add_New(cur_face) && periodicorbit_detector->cellcycle->nelems > 1)
        {
            ////Move the repeated cycle judgement here 07/25/06
            int pre_limit_index;

            if(is_existing_porbit(cur_face, pre_limit_index, 1-type))
            {
                flag = 1;
                return false;
            }

            //the_e = g_theedge;
            the_e = chosen_edge = g_chosen_edge;

            if(the_e == nullptr)  //cross a vertex
                continue;

            if(the_e->num_intersections > 3)
            {
                //do not consider the tangent point now
                if(pre_e == the_e)
                    if(the_e->approx_tangent(seg_length, intersect))
                        continue;

                //validate the intersections on the edge
                the_e->intersections[1] = intersect;        //save the current intersection to the_e->intersections[1]
                if(the_e->is_good_intersections(pre_face, scc_index))  //this judgement may not be very good for highly curled field
                {
                    //g_chosenedge = chosen_edge = the_e;

                    chosen_edge = the_e;
                    chosen_edge->intersections[1] = intersect;
                    triangle = cur_face;
                    return true;                             //we find a cycle here 07/23/06
                }
            }

            the_e->intersections[0] = intersect;  //save the intersection with the edge
            the_e->num_intersections ++;
        }
    }

    return false;  //can not converge!

    delete detect_traj;
}

/*
Trace inside a triangle for periodic orbit detection
*/
int PeriodicOrbit_Detector::trace_in_triangle_po_detect(double globalp[3], int &face_id, int type,
                                                        Edge *cur_e, icVector3 &intersect, double &seg_length, int &flag)
{

    //similar to the regular tracing, except that you need to mark those special points too close to the cycle
    int i;
    double alpha[3];
    double cur_point[2], pre_point[2];
    double vert0[3];
    icVector3 VP, globalv;

    if(face_id < 0)
        return -1;

    Triangle *face = object->tlist.tris[face_id];
    Triangle *pre_f = face;


    icVector2 dis;
    seg_length = 0;

    ////initialize
    VP.entry[0] = globalp[0] - face->verts[0]->x;
    VP.entry[1] = globalp[1] - face->verts[0]->y;
    VP.entry[2] = globalp[2] - face->verts[0]->z;

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

            //if(detect_traj->cal_next_point_euler1(pre_point, cur_point, face_id, alpha, type))
            //if(detect_traj->cal_nextpt_2ndeuler(pre_point, cur_point, face_id, alpha, type))
            //if(detect_traj->cal_nextpt_RK4(pre_point, cur_point, face_id, alpha, type))
            if(detect_traj->get_next_pt(pre_point, cur_point, face_id, alpha, type, unsigned char(Integrator_opt)))
            {
                ////update the global point

                globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

                //Get global coordinates of the point
                globalp[0] = vert0[0] + globalv.entry[0];
                globalp[1] = vert0[1] + globalv.entry[1];
                globalp[2] = vert0[2] + globalv.entry[2];

            }

            else{  ////the curve reach a singularity
                flag = 1;
                return face_id;
            }
        }

        ////3. if the point is out of current triangle
        else{
            double t[2] = {0.};

            int PassVertornot = 0;

            int presave_face = face_id;

            detect_traj->cross_a_vertex(face_id, cur_point, pre_point, type, PassVertornot);

            //get_next_triangle(face_id, pre_point, cur_point, t, type, PassVertornot, alpha);

            ////update the globalpoint here (Modified on 01/30/07)
            if(PassVertornot > 0)
            {
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

                g_chosen_edge = chosen_edge = cur_e = nullptr;  //do not consider passing vertex cases now!!! 07/23/06
            }

            else{
                face_id = presave_face;

                int which_edge = -1;

                detect_traj->cross_boundary(pre_point, cur_point, face_id, alpha, which_edge, t);

                detect_traj->pass_edge(face_id, which_edge);

                //// transfer it to the global coordinates
                globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;

                intersect.entry[0] = globalp[0] = vert0[0] + globalv.entry[0];
                intersect.entry[1] = globalp[1] = vert0[1] + globalv.entry[1];
                intersect.entry[2] = globalp[2] = vert0[2] + globalv.entry[2];

                ////Get the corresponding edge
                for(int k = 0; k < 3; k++)
                {
                    int vertindex = face->verts[which_edge]->index;

                    cur_e = face->edges[k];
                    //g_theedge = cur_e;
                    g_chosen_edge = chosen_edge = cur_e;

                    if(cur_e->verts[0]->index != vertindex && cur_e->verts[1]->index != vertindex)
                        break;
                }
            }

            return face_id;
        }

    }

    return face_id;
}

/*
Calculate the closed streamline of the limit cycle according to the input edge and triangle
*/
bool PeriodicOrbit_Detector::cal_closedstreamline(Edge *the_e, int triangle, int type, int scc_index)
{
    int i;
    int flag = 0;
    double globalp[3];
    int pre_face, cur_face;

    icVector2 curp, dis;

    int stop = 1;

    ////Initialization part

    pre_face = cur_face = triangle;

    ////We use the latest intersection as the approximate fixed point to start the tracing
    globalp[0] = the_e->intersections[1].entry[0];
    globalp[1] = the_e->intersections[1].entry[1];
    globalp[2] = the_e->intersections[1].entry[2];

    cellcycle->nelems = 0;

    /*allocate space for the closed streamline of the periodic orbit*/
    periodic_orbits->polist[periodic_orbits->nporbits]->traj
        = new Trajectory(periodic_orbits->nporbits, 200);

    cellcycle->add_New(cur_face);


    int NUMTRACETRIS = (int)sqrt((double)object->tlist.ntris);

    for(i = 0; i < 20*NUMTRACETRIS && i < object->tlist.ntris; i++)
    {
        if(cur_face < 0)
        {
            delete periodic_orbits->polist[periodic_orbits->nporbits]->traj;
            return false;
        }

        pre_face = cur_face;
        cur_face = periodic_orbits->polist[periodic_orbits->nporbits]->traj->
                   trace_in_triangle(cur_face, globalp, type, flag);

        if(flag == 1 || flag == 2 || pre_face == cur_face )
        {
            delete periodic_orbits->polist[periodic_orbits->nporbits]->traj;
            return false;
        }

        if(!is_repeated_elem(cellcycle->elems, cur_face, cellcycle->nelems))
        {
            cellcycle->add_New(cur_face);
        }

        else{ ////form a cell cycle, probably we can stop tracing now
            if(!get_Cycle(cellcycle->elems, cellcycle->nelems, cur_face))
            {
                delete periodic_orbits->polist[periodic_orbits->nporbits]->traj;
                return false;
            }

            //we need to get the real cycle!
            if(stop >= 4)
                return true;
            stop ++;
        }

        /*
        The current tracing curve should not go outside current SCC !
        */
        if(!is_repeated_elem(morse_decomp->scclist->scccomponents[scc_index]->nodes, cur_face,
                              morse_decomp->scclist->scccomponents[scc_index]->nnodes))
        {
            flag = 1;     //fail, reach the boundary of the SCC
            delete periodic_orbits->polist[periodic_orbits->nporbits]->traj;
            return false;
        }

        /*
        We test whether it is previous detected cycle again here 07/25/06
        */
        int pre_limit_index = -1;
        if(is_existing_porbit(cur_face, pre_limit_index, 1-type))
        {
            //we need to build the connection later
            delete periodic_orbits->polist[periodic_orbits->nporbits]->traj;
            return false;
        }

        ////Mark those near by special points having the same type as the limit cycle as "dead" point
        mark_neighbor_pts(globalp, pre_face, scc_index, 1-type);

    }

    periodic_orbits->polist[periodic_orbits->nporbits]->traj;
    return false;
}


/*
During calculating the closed streamline of a detected limit cycle, we need to mark all its neighboring
special points as "dead" points
For attracting limit cycles, mark those close attachment points;
for repelling limit cycles, mark those close separation points.
*/
void PeriodicOrbit_Detector::mark_neighbor_pts(double p[3], int triangle, int scc_index, int type)
{
    int *NearbyTriangles = nullptr;
    int num_triangles = 0;
    int i, j;
    Triangle *face;
    Edge *cur_e;

    NearbyTriangles = get_disc_coarse(p, triangle, NearbyTriangles, num_triangles);
    for(i = 0; i < num_triangles; i++)
    {
        //Do not change the states of those special points not belonging to current SCC
        if(!is_repeated_elem(morse_decomp->scclist->scccomponents[scc_index]->nodes, NearbyTriangles[i],
                              morse_decomp->scclist->scccomponents[scc_index]->nnodes))
            continue;

        face = object->tlist.tris[NearbyTriangles[i]];

        for(j = 0; j < 3; j++)
        {
            cur_e = face->edges[j];

            ////No need to find out whether the edge contains a sep / att point or not
            if(type == 0)
                cur_e->sep_visit = 1;
            else
                cur_e->att_visit = 1;

        }
    }

    free(NearbyTriangles);

}

void  PeriodicOrbit_Detector::store_current_cellcycle(int *acellcycle, int num)
{
    //periodic_orbits->polist[periodic_orbits->nporbits]->cellcycle = (int*)malloc(sizeof(int) * num);
    periodic_orbits->polist[periodic_orbits->nporbits]->cellcycle = new int[num];

    for(int i = 0; i < num; i++)
    {
        periodic_orbits->polist[periodic_orbits->nporbits]->cellcycle[i] = acellcycle[i];
    }

    periodic_orbits->polist[periodic_orbits->nporbits]->ntris = num;

}


/*
As a special separatrices, the closed streamline of a periodic orbit
should be stored with the other separatrices
*/
void PeriodicOrbit_Detector::store_cur_streamline(LineSeg *streamline, int num)
{
    //limitcycles[cur_limitcycle_index].closed_streamline = (LineSeg*)malloc(sizeof(LineSeg)*num);

    //for(int i = 0; i < num; i++)
    //{
    //	limitcycles[cur_limitcycle_index].closed_streamline[i].gstart[0] = streamline[i].gstart[0];
    //	limitcycles[cur_limitcycle_index].closed_streamline[i].gstart[1] = streamline[i].gstart[1];
    //	limitcycles[cur_limitcycle_index].closed_streamline[i].gstart[2] = streamline[i].gstart[2];

    //	limitcycles[cur_limitcycle_index].closed_streamline[i].gend[0] = streamline[i].gend[0];
    //	limitcycles[cur_limitcycle_index].closed_streamline[i].gend[1] = streamline[i].gend[1];
    //	limitcycles[cur_limitcycle_index].closed_streamline[i].gend[2] = streamline[i].gend[2];

    //	limitcycles[cur_limitcycle_index].closed_streamline[i].start[0] = streamline[i].start[0];
    //	limitcycles[cur_limitcycle_index].closed_streamline[i].start[1] = streamline[i].start[1];
    //	limitcycles[cur_limitcycle_index].closed_streamline[i].start[2] = streamline[i].start[2];

    //	limitcycles[cur_limitcycle_index].closed_streamline[i].end[0] = streamline[i].end[0];
    //	limitcycles[cur_limitcycle_index].closed_streamline[i].end[1] = streamline[i].end[1];
    //	limitcycles[cur_limitcycle_index].closed_streamline[i].end[2] = streamline[i].end[2];

    //	limitcycles[cur_limitcycle_index].closed_streamline[i].Triangle_ID = streamline[i].Triangle_ID;
    //}

    //limitcycles[cur_limitcycle_index].num_linesegs = num;
}


bool PeriodicOrbit_Detector::is_existing_porbit(int triangle, int &pre_limit_index, unsigned char type)
{
    int i;

    for(i = 0; i < periodic_orbits->nporbits; i++)
    {

        if(periodic_orbits->polist[i]->type != type)
            continue;

        if(is_repeated_elem(periodic_orbits->polist[i]->cellcycle, triangle,
                             periodic_orbits->polist[i]->ntris))
        {
            pre_limit_index = i;  //return the limit cycle index
            return true;
        }
    }

    return false;
}

/*
Get the set of triangles fall inside the distance circle with the center being the input point p
In this routine, we use a coarse measure to find the disc.
We select the one ring neighborhood of the triangle containing the input point p
*/
int *PeriodicOrbit_Detector::get_disc_coarse(double p[3], int triangle,
                                             int *NearbyTriangles, int &num_triangles)
{
    Vertex *vert;
    Triangle *face = object->tlist.tris[triangle];
    Corner *c;

    int i, j;

    ////Find all the triangle
    num_triangles = 0;
    NearbyTriangles = extend_link(NearbyTriangles, num_triangles);
    NearbyTriangles[0] = triangle;
    num_triangles++;

    /* We here use the one ring neighborhood since it is different from the streamline tracing here */
    for(i = 0; i < 3; i++)
    {
        vert = face->verts[i];

        ////add the one ring neighborhood of the vertex to the list
        for(j = 0; j < vert->ncorners; j++)
        {
            c = vert->corners[j];

            if(is_repeated_elem(NearbyTriangles, c->t, num_triangles))
                continue;

            NearbyTriangles = extend_link(NearbyTriangles, num_triangles);
            NearbyTriangles[num_triangles] = c->t;
            num_triangles++;
        }
    }

    return NearbyTriangles;
}


/*
Repeat the routine in class Edge
*/
int *Edge::get_disc_coarse(double p[3], int triangle,
                           int *NearbyTriangles, int &num_triangles)
{
    Vertex *vert;
    Triangle *face = object->tlist.tris[triangle];
    Corner *c;

    int i, j;

    ////Find all the triangle
    num_triangles = 0;
    NearbyTriangles = extend_link(NearbyTriangles, num_triangles);
    NearbyTriangles[0] = triangle;
    num_triangles++;

    /* We here use the one ring neighborhood since it is different from the streamline tracing here */
    for(i = 0; i < 3; i++)
    {
        vert = face->verts[i];

        ////add the one ring neighborhood of the vertex to the list
        for(j = 0; j < vert->ncorners; j++)
        {
            c = vert->corners[j];

            if(is_repeated_elem(NearbyTriangles, c->t, num_triangles))
                continue;

            NearbyTriangles = extend_link(NearbyTriangles, num_triangles);
            NearbyTriangles[num_triangles] = c->t;
            num_triangles++;
        }
    }

    return NearbyTriangles;
}




