/*
This file contains routines of class MCG_Graph that we use to build the MCG graph based on
the obtain Morse decomposition

Created and Modified by Guoning Chen
copyright @2007
*/


#include "Analysis/MorseDecomp.h"
#include "Others/common_routines.h"
#include "VField.h"


extern Polyhedron *object;
extern MorseDecomp *morse_decomp;
extern int DebugOn;

MCG_Graph *mcg = NULL;
//local mcg
MCG_Graph *lmcg=NULL;//added 12/06/2009 by Qingqing

int *searchstack = NULL;
int nelem_in_stack;
int MaxStackElems;

int *repell_region = NULL;
int ntris_repell_region;
int curMaxNumTrisRepellRegion;


extern double used_tau ;
extern int Cal_Regions;//should the regions be computed?

bool EnConleyIndexComp = true;

bool RemRedundantMCGEdges = false;

void MCG_Graph::build_mcg()
{
    FILE *fp;
    //
    //fp=fopen("detect_porbit.txt", "a");
    //fprintf(fp, "start assigning the nodes in MCG\n");
    //fclose(fp);

    /*initialize all the triangle*/
    for(int i=0; i<object->tlist.ntris; i++)
    {
        object->tlist.tris[i]->visited = false;
        object->tlist.tris[i]->MS_type = 0;
    }

    assign_mcgnodes();


    //fp=fopen("detect_porbit.txt", "a");
    //fprintf(fp, "There are %d nodes in this MCG.\n", nlist->nmnodes);
    //fprintf(fp, "start computing the edges in MCG\n");
    //fclose(fp);

    build_mcg_edges_graph();
    layout_mcg();
    if(Cal_Regions)
    {
        cal_mcgedge_regions();
    }

}


/*
Assign nodes for the Morse sets in MCG
*/
void MCG_Graph::assign_mcgnodes()
{
    cur_mcgnode_index = 0;
    nlist->nmnodes = 0;

    r_counter = a_counter = s_counter = 0;

    int i;

    for(i = 0; i < morse_decomp->scclist->nsccs; i++)
    {
        //FILE *fp=fopen("detect_porbit.txt", "a");
        //fprintf(fp, "start assigning node %d.\n", i);
        //fclose(fp);

        if(morse_decomp->scclist->scccomponents[i]->nnodes <= 2
            && morse_decomp->scclist->scccomponents[i]->nfixedpoints < 1)
            continue;

        if(morse_decomp->scclist->scccomponents[i]->nnodes >= 2
            && morse_decomp->scclist->scccomponents[i]->valid == false)
            continue;

        //FILE *fp=fopen("detect_porbit.txt", "a");
        //fprintf(fp, "finish assigning node %d.\n", i);
        //fprintf(fp, "current number of nodes is %d.\n", cur_mcgnode_index);
        //fclose(fp);

        add_to_mcg_nodelist(i, cur_mcgnode_index);

        //fp=fopen("detect_porbit.txt", "a");
        //fprintf(fp, "finish adding to the list %d.\n", i);
        //fclose(fp);

        morse_decomp->scclist->scccomponents[i]->node_index = cur_mcgnode_index-1;

    }

}


/*lay out the nodes*/
//extern double ConleyGraphWin_x ;
//extern double ConleyGraphWin_y ;
//extern double types_interval_y ;
//double ConleyGraphWin_x ;
//double ConleyGraphWin_y ;
double ConleyGraphWin_x = 4;
double ConleyGraphWin_y = 2;
double types_interval_y = ConleyGraphWin_y / 4;

void MCG_Graph::layout_mcg()
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
    for(i = 0; i < cur_mcgnode_index; i++)
    {
        if(nlist->mnodes[i]->type == 0)
            num_repellers ++;
        else if(nlist->mnodes[i]->type == 1)
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
    for(i = 0; i < cur_mcgnode_index; i++)
    {
        if(nlist->mnodes[i]->type == 0)
        {
            nlist->mnodes[i]->pos.entry[0] = (cur_repeller_index+1)*repeller_interval;
            nlist->mnodes[i]->pos.entry[1] = repeller_y;
            cur_repeller_index ++;
        }
        else if(nlist->mnodes[i]->type == 1)
        {
            nlist->mnodes[i]->pos.entry[0] = (cur_attractor_index+1)*attractor_interval;
            nlist->mnodes[i]->pos.entry[1] = attractor_y;
            cur_attractor_index ++;
        }
        else{
            nlist->mnodes[i]->pos.entry[0] = (cur_saddle_index+1)*saddle_interval;
            nlist->mnodes[i]->pos.entry[1] = saddle_y;
            cur_saddle_index ++;
        }
    }

}


void MCG_Graph::build_mcg_edges_graph()
{
    cur_mcgedge_index = 0;

    /*Testing code*/
    FILE *fp;
    //fp = fopen("detect_porbit_cooling.txt", "a");
    //fp = fopen("detect_porbit_swirl.txt", "a");
    //fp = fopen("detect_porbit_tumble.txt", "a");
    if(DebugOn == 1){
    //fp = fopen("detect_porbit.txt", "a");
    //fprintf(fp, "current # of Morse sets is %d. \n", nlist->mnodes);
    //fprintf(fp, "start growing regions ... \n");
    //fclose(fp);
    }

    /*grow all nodes without knowing their types*/
    grow_all_mcgnodes();

    //fp = fopen("detect_porbit_cooling.txt", "a");
    //fp = fopen("detect_porbit_swirl.txt", "a");
    //fp = fopen("detect_porbit_tumble.txt", "a");
    if(DebugOn == 1){
    //fp = fopen("detect_porbit.txt", "a");
    //fprintf(fp, "finish growing edges ... \n");
    //fclose(fp);
    }

    /*remove any redundant edges*/
    if (RemRedundantMCGEdges )
        remove_redundant_edges();
}


unsigned char
get_region_type_2 (int *tris, int ntris)
{
    int i;

    int nsources=0, nsinks=0, nsaddles=0;

    unsigned char type = 0;

    for (i=0; i<ntris; i++)
    {
        Triangle * t = object->tlist.tris[tris[i]];

        if (t->singularityID>=0)
        {
            switch(object->slist.slist[t->singularityID]->type)
            {
            case SOURCE:
            case RFOCUS:
            case CWCENTER:
            case CCWCENTER:
                nsources++;
                break;

            case SINK:
            case AFOCUS:
                nsinks++;
                break;

            case SADDLE:
                nsaddles++;
                break;

            default:
                nsources++;
                nsinks++;
                break;
            }
        }
    }

    if (nsources+nsinks == nsaddles)
        type = 4;
    else
    {
        if (nsources+nsinks < nsaddles)
            type = SADDLE;

        else if (nsources > nsinks)
            type = SOURCE;
        else
            type = SINK;
    }

    return type;
}


void MCG_Graph::add_to_mcg_nodelist(int scc_index, int &index)
{
    if(index >= nlist->curMaxNumMNodes)
    {
        int oldnum = nlist->curMaxNumMNodes;

        /*extend the memory*/
        if(!nlist->extend(50))
        {
            /*report errors*/
            return;
        }

        for(int i = oldnum; i < nlist->curMaxNumMNodes; i++)
            nlist->mnodes[i] = new MCG_Node();
    }

    nlist->mnodes[index]->node_index = index;
    nlist->mnodes[index]->scc_index = scc_index;
    nlist->mnodes[index]->graph_edges = NULL;
    nlist->mnodes[index]->nedges = 0;
    nlist->mnodes[index]->visited = false;
    nlist->mnodes[index]->cancelled = false;
    nlist->mnodes[index]->boundaryN=0;
    nlist->mnodes[index]->boundary=NULL;

    /*  record the tau that we use to obtain this Morse set  */
    nlist->mnodes[index]->used_tau = used_tau;

    nlist->mnodes[index]->type = 3;  //non-determined MS type

    /*we also need to get the type of the node!!!*/

    if (!EnConleyIndexComp)
    {
        nlist->mnodes[index]->type = get_region_type(scc_index);


    /*we make the correction for the saddle region with triangles fewer than 5
    03/08/07*/

        if(nlist->mnodes[index]->type == 2
            && morse_decomp->scclist->scccomponents[scc_index]->nnodes < 5)
            nlist->mnodes[index]->type = make_correction_saddle(scc_index);

        if (	nlist->mnodes[index]->type < 0 ||
                nlist->mnodes[index]->type >= 3)
                nlist->mnodes[index]->type = 3;  //non-determined MS type




        //Assign the triangle labels based on the type of Morse set that each triangle belongs to
        for (int ii = 0; ii < morse_decomp->scclist->scccomponents[scc_index]->nnodes; ii++)
            object->tlist.tris[morse_decomp->scclist->scccomponents[scc_index]->nodes[ii]]->MS_type = nlist->mnodes[index]->type+1;
    }
    else
    {
        nlist->mnodes[index]->type =conley->GetMCGType(scc_index);
        /*
            correct the saddle triangles
        */
        nlist->mnodes[index]->priority=conley->GetMCGPriority();
        nlist->mnodes[index]->boundary=conley->GetMCGBoundary(nlist->mnodes[index]->boundaryN);
        conley->GetMCGConley(nlist->mnodes[index]->conley);

        if (morse_decomp->scclist->scccomponents[scc_index]->nnodes < 3 /*&&
            nlist->mnodes[index]->conley[0]==0 && nlist->mnodes[index]->conley[1]==0 && nlist->mnodes[index]->conley[2]==0*/)
        {
            /*
               get the type of the Morse set by counting the number of singularities
            */

            unsigned char correct_type = get_region_type_2 (morse_decomp->scclist->scccomponents[scc_index]->nodes,
                morse_decomp->scclist->scccomponents[scc_index]->nnodes);

            switch(correct_type)
            {
            case SOURCE:
                nlist->mnodes[index]->conley[2] = 1;
                nlist->mnodes[index]->type = 0;
                break;

            case SINK:
                nlist->mnodes[index]->conley[0] = 1;
                nlist->mnodes[index]->type = 1;
                break;

            case SADDLE:
                nlist->mnodes[index]->conley[1] = 1;
                nlist->mnodes[index]->type = 2;
                break;

            default:
                break;
            }

        }
        AddQueue(index);

        if (	nlist->mnodes[index]->type < 0 ||
                nlist->mnodes[index]->type >= 3)
                nlist->mnodes[index]->type = 3;  //non-determined MS type

        //Assign the triangle labels based on the type of Morse set that each triangle belongs to
        for (int ii = 0; ii < morse_decomp->scclist->scccomponents[scc_index]->nnodes; ii++)
            object->tlist.tris[morse_decomp->scclist->scccomponents[scc_index]->nodes[ii]]->MS_type = nlist->mnodes[index]->type+1;

    }

    /*if(nlist->mnodes[index]->boundaryE.size())nlist->mnodes[index]->boundaryE.clear();
    for(int i=0;i<conley->boundary_edgelist.size();i++)
        nlist->mnodes[index]->boundaryE.push_back(conley->boundary_edgelist[i]);*/






    /*------------------------------------------------------------------*/
    /*the inaccurate judgement could happen for larger Morse sets!!!!*/

    /*------------------------------------------------------------------*/

    if( nlist->mnodes[index]->type == 0)
    {
        nlist->mnodes[index]->labelindex = r_counter;
        r_counter++ ;
    }
    else if( nlist->mnodes[index]->type == 1)
    {
        nlist->mnodes[index]->labelindex = a_counter;
        a_counter++;
    }
    else
    {
        nlist->mnodes[index]->labelindex = s_counter;
        s_counter++;
    }

    nlist->nmnodes++;
    index++;
}


/*
Get the type of the region according to the incoming and outgoing edges of the nodes in
the SCC component
*/
int MCG_Graph::get_region_type(int scc_index)
{
    int i;
    int r_flag, a_flag;

    r_flag = a_flag = 0;


    if(get_boundary_nodes(scc_index)==1)
    {
        for(i = 0; i < morse_decomp->scclist->scccomponents[scc_index]->nnodes; i++)
        {
            if(get_triangle_type(morse_decomp->scclist->scccomponents[scc_index]->nodes[i], scc_index) == 0)
            {
                r_flag = 1;

                if(a_flag == 1)
                    return 2;
            }

            else if(get_triangle_type(morse_decomp->scclist->scccomponents[scc_index]->nodes[i], scc_index) == 1)
            {
                a_flag = 1;

                if(r_flag == 1)
                    return 2;
            }
            else
                return 2;
        }
    }

    else{
        for(i = 0; i < nboundarytris; i++)
        {
            if(get_triangle_type(boundary_nodes[i], scc_index) == 0)
            {
                r_flag = 1;

                if(a_flag == 1)
                    return 2;
            }

            else if(get_triangle_type(boundary_nodes[i], scc_index) == 1)
            {
                a_flag = 1;

                if(r_flag == 1)
                    return 2;
            }
            else
                return 2;
        }
    }

    if(r_flag == 1 && a_flag == 0)
        return 0;
    if(r_flag == 0 && a_flag == 1)
        return 1;
    return 3;
}

/*
The number of the triangles in the saddle-like regions we want to make correction to
is 4 now. 03/08/07
Probably this should happen during the judgement
*/
int MCG_Graph::make_correction_saddle(int scc_index)
{
    int j;
    int source, sink, saddle;

    source=sink=saddle = 0;

    for(j = 0; j < morse_decomp->scclist->scccomponents[scc_index]->nnodes; j++)
    {
        int t = morse_decomp->scclist->scccomponents[scc_index]->nodes[j];

        if( object->tlist.tris[t]->singularityID < 0)
            continue;

        if(object->slist.slist[object->tlist.tris[t]->singularityID]->type == SOURCE)
            source++;
        else if(object->slist.slist[object->tlist.tris[t]->singularityID]->type == SINK)
            sink++;
        else if(object->slist.slist[object->tlist.tris[t]->singularityID]->type == SADDLE)
            saddle++;
    }

    if(saddle == 0)
    {
        if(source > 0)
            return 0;
        else
            return 1;
    }
    else
        return 2;

    return 3;
}


/*
Grow all the nodes forward without considering its type
03/27/07
*/
void MCG_Graph::grow_all_mcgnodes()
{
    int i;

    clock_t l_start, l_finish;

    for(i = 0; i < nlist->nmnodes/*cur_mcgnode_index*/; i++)
    {
        //l_start = clock();

        if (nlist->mnodes[i]->cancelled) continue;

        grow_repeller_region_graph(nlist->mnodes[i]->scc_index, 0);

        //l_finish = clock();
        //FILE *fp = fopen("detect_porbit_cooling.txt", "a");
        //fprintf(fp, "the time for growing edges of node %d is %f seconds.\n",
        //	i, (double)(l_finish - l_start)/CLOCKS_PER_SEC);
        //fclose(fp);
    }

    /*reverse the edges*/
 //   ReverseEdges();

    //
    ///*grow based on reversed graph*/
    //for(i = 0; i < cur_mcgnode_index; i++)
    //{
    //	grow_repeller_region_graph_2(mcgnodes[i].scc_index, 1);
    //}

    //ReverseEdges();
}


/*
remove redundant edges in MCG and make it minimum
*/
void MCG_Graph::remove_redundant_edges()
{
    int i, j;
    int *expandnodes = new int[cur_mcgnode_index];
    int num_expanded_nodes;
    int cur_node, other_node;

    for(i = 0; i < cur_mcgnode_index; i++)
    {
        for(j = 0; j < nlist->mnodes[i]->nedges; j++)
        {
            other_node = elist->edges[nlist->mnodes[i]->graph_edges[j]]->node_index2;
            if(other_node == i)
                continue;  /*consider outgoing edges only*/

            if(has_longer_path_in_mcg(i, other_node, nlist->mnodes[i]->graph_edges[j]))
            {
                elist->edges[nlist->mnodes[i]->graph_edges[j]]->cancel = true;
            }
        }
    }

    /*update the edge list in the nodes*/
    for(i = 0; i < cur_mcgedge_index; i++)
    {
        if(elist->edges[i]->cancel)
        {
            del_one_edge_from(nlist->mnodes[elist->edges[i]->node_index1]->graph_edges,
                nlist->mnodes[elist->edges[i]->node_index1]->nedges, i);
            del_one_edge_from(nlist->mnodes[elist->edges[i]->node_index2]->graph_edges,
                nlist->mnodes[elist->edges[i]->node_index2]->nedges, i);
        }
    }
    delete [] expandnodes;
}


/*
Get the boundary nodes of the specified SCC
*/
int MCG_Graph::get_boundary_nodes(int scc_index)
{
    int MaxTris;
    if(morse_decomp->scclist->scccomponents[scc_index]->nnodes < 5)
    {
        if(boundary_nodes != NULL)
            free(boundary_nodes);
        MaxTris = morse_decomp->scclist->scccomponents[scc_index]->nnodes;
        boundary_nodes = (int*)malloc(sizeof(int)*MaxTris);
        nboundarytris = 0;
        int i;
        for(i = 0; i < morse_decomp->scclist->scccomponents[scc_index]->nnodes; i++)
        {
            boundary_nodes[nboundarytris]=morse_decomp->scclist->scccomponents[scc_index]->nodes[i];
            nboundarytris++;
        }
        return 1;
    }


    MaxTris = morse_decomp->scclist->scccomponents[scc_index]->nnodes;
    boundary_nodes = (int*)malloc(sizeof(int)*MaxTris);
    nboundarytris = 0;

    int i, j;
    Triangle *face;
    Edge *e;

    for(i = 0; i < morse_decomp->scclist->scccomponents[scc_index]->nnodes; i++)
    {
        face = object->tlist.tris[morse_decomp->scclist->scccomponents[scc_index]->nodes[i]];
        face->visited=true;
    }

    for(i = 0; i < morse_decomp->scclist->scccomponents[scc_index]->nnodes; i++)
    {
        face = object->tlist.tris[morse_decomp->scclist->scccomponents[scc_index]->nodes[i]];
        for(j = 0; j < 3; j++)
        {
            e = face->edges[j];

            if(e->tris[0]==NULL||e->tris[1]==NULL)/*this is a boundary face*/
            {
                //boundary_nodes[nboundarytris] = face->index;
                //nboundarytris++;
                continue;
            }

            /*this will be really slow*/
            //if(is_repeated_elem(morse_decomp->scclist->scccomponents[scc_index]->nodes, e->tris[0]->index,
            //	morse_decomp->scclist->scccomponents[scc_index]->nnodes)
            //	&& is_repeated_elem(morse_decomp->scclist->scccomponents[scc_index]->nodes, e->tris[1]->index,
            //	morse_decomp->scclist->scccomponents[scc_index]->nnodes))
            //	continue;  //this is an inner edge
            if(object->tlist.tris[e->tris[0]->index]->visited
                && object->tlist.tris[e->tris[1]->index]->visited)
                continue;  //this is an inner edge

            /*this is a boundary triangle*/
            boundary_nodes[nboundarytris] = face->index;
            nboundarytris++;
            break;
        }
    }

    /*reset the nodes in this SCC*/
    for(i = 0; i < morse_decomp->scclist->scccomponents[scc_index]->nnodes; i++)
    {
        face = object->tlist.tris[morse_decomp->scclist->scccomponents[scc_index]->nodes[i]];
        face->visited=false;
    }

    return nboundarytris;
}


/*
the routine return the type of the specified triangle "tri"
0-- repelling, means that all its incoming edges are inside the
same SCC
1-- attracting,
2-- saddle
The input is the directed graph obtained during Morse decomposition
The problem is the triangle having tangential edge may be treated
as saddle triangle
*/
int MCG_Graph::get_triangle_type(int tri, int scc_index)
{
    int other_node;
    int i;
    int r_flag, a_flag;
    r_flag = a_flag = 0;
    DirGraph_Node **sccnodes = morse_decomp->dg->nlist->dirnodes;
    Graph_Edge **sccedges = morse_decomp->dg->elist->edges;

    for(i = 0; i < sccnodes[tri]->nedges; i++)
    {
        /*we do not consider the edge pointing to itself*/
        if(sccedges[sccnodes[tri]->edges[i]]->node_index1 == tri
            && sccedges[sccnodes[tri]->edges[i]]->node_index2 == tri)
            continue;

        other_node = sccedges[sccnodes[tri]->edges[i]]->node_index2;

        /*this is an outward edge*/
        if(other_node != tri)
        {
            if(sccnodes[other_node]->sscomp_index != scc_index)
            {
                r_flag ++;
            }
            continue;
        }

        /*this is an incoming edge*/
        other_node = sccedges[sccnodes[tri]->edges[i]]->node_index1;

        /*this edge is coming from the other node in the same SCC*/
        if(sccnodes[other_node]->sscomp_index == scc_index)
            continue;
        else{
            a_flag++;
        }
    }

    if(r_flag > 0 && a_flag == 0)
        return 0;
    else if(r_flag == 0 && a_flag > 0)
        return 1;
    else
        return 2;
}


/*
grow the repelling region using the obtained directed graph
New method but still based BFS
03/27/07
*/
void MCG_Graph::grow_repeller_region_graph(int scc_index, int inverse)
{

    int i, j, k, cur_node, other_node;
    DirGraph_Node **sccnodes = morse_decomp->dg->nlist->dirnodes;
    Graph_Edge **sccedges = morse_decomp->dg->elist->edges;

    /*reset all the flags of the visited nodes*/
    for(i = 0; i < ntris_repell_region; i++)
    {
        sccnodes[repell_region[i]]->visited = 0;
    }

    if(repell_region != NULL)
        free(repell_region);

    curMaxNumTrisRepellRegion = object->tlist.ntris;
    repell_region = (int*)malloc(sizeof(int)*curMaxNumTrisRepellRegion);
    ntris_repell_region = 0;

    /*
    use all the triangles in the region
    */
    copy_array_Int(morse_decomp->scclist->scccomponents[scc_index]->nodes, repell_region,
        morse_decomp->scclist->scccomponents[scc_index]->nnodes);
    ntris_repell_region = morse_decomp->scclist->scccomponents[scc_index]->nnodes;


    //reset_sccnodes_flags();

    for(i = 0; i < ntris_repell_region; i++)
    {
        /*expand current node along its outgoing edges*/
        cur_node = repell_region[i];

        if(sccnodes[cur_node]->visited == 2)
            continue;

        sccnodes[cur_node]->visited = 2;

        for(j = 0; j < sccnodes[cur_node]->nedges; j++)
        {
            if(sccedges[sccnodes[cur_node]->edges[j]]->node_index2 == cur_node)
                continue;

            other_node = sccedges[sccnodes[cur_node]->edges[j]]->node_index2;



            if(sccnodes[other_node]->sscomp_index != scc_index
                && sccnodes[other_node]->visited == 0
                /*&& !is_repeated_elem(repell_region, other_node, ntris_repell_region)*/
                )
            {
                /*judge whether it corresponds to other valid Morse set*/
                /*if((morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->nnodes <= 2
                    && morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->nfixedpoints > 0)
                    ||(morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->nnodes > 2
                    && morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->valid == true))*/
                if(IsMCGN(sccnodes[other_node]->sscomp_index)
                    && morse_decomp->scclist->scccomponents[scc_index]->node_index>=0
                && morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index>=0)
                {
                    if(!is_connected_mcg(morse_decomp->scclist->scccomponents[scc_index]->node_index,
                        morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index)
                        //&& !has_interval_mcg(scclist.scccomponents[scc_index].node_index,
                        //scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index)

                        && !nlist->mnodes[morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index]->cancelled
                        )
                    {
                        if(inverse == 0)
                        {
                            if(is_valid_link(morse_decomp->scclist->scccomponents[scc_index]->node_index,
                                morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index))
                            {
                                add_to_mcg_edgelist(morse_decomp->scclist->scccomponents[scc_index]->node_index,
                                    morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index);
                                /*add the edge to the corresponding nodes*/
                                add_edge_to_node(morse_decomp->scclist->scccomponents[scc_index]->node_index, cur_mcgedge_index-1);
                                add_edge_to_node(morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index,
                                    cur_mcgedge_index-1);
                            }
                        }
                        else{
                            if(is_valid_link(morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index,
                                morse_decomp->scclist->scccomponents[scc_index]->node_index))
                            {
                                add_to_mcg_edgelist(morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index,
                                    morse_decomp->scclist->scccomponents[scc_index]->node_index);
                                /*add the edge to the corresponding nodes*/
                                add_edge_to_node(morse_decomp->scclist->scccomponents[scc_index]->node_index, cur_mcgedge_index-1);
                                add_edge_to_node(morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index,
                                    cur_mcgedge_index-1);
                            }
                        }

                        /*add all the nodes in this Morse sets into the array*/
                        for(k = 0; k < morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->nnodes; k++)
                        {
                            if(sccnodes[morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->nodes[k]]->visited
                                == 1)
                                continue;
                            repell_region[ntris_repell_region] =
                                morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->nodes[k];
                            /*set the flag*/
                            sccnodes[repell_region[ntris_repell_region]]->visited = 1;
                            ntris_repell_region++;
                        }

                    }
                }
                else{
                    /*add the node into the array*/
                    repell_region[ntris_repell_region] = other_node;
                            sccnodes[repell_region[ntris_repell_region]]->visited = 1;
                    ntris_repell_region++;
                }
            }
        }
    }
}



bool MCG_Graph::has_longer_path_in_mcg(int n1, int n2, int edgeindex)
{
    //the maximum searching number is the number of the edges in the graph

    int i;
    //reset the 'visited' flags of all the edges
    for(i = 0; i < cur_mcgedge_index; i++)
    {
        elist->edges[i]->visited = false;
    }

    /*set the direct edge from 'from' to 'to' to be visited
    */
    elist->edges[edgeindex]->visited = true;

    //the maximum searching number is the number of the edges in the graph
    int count = 0;
    int cur_node = n1;
    int cur_node_id = 0;

    if(searchstack != NULL)
        free(searchstack);

    searchstack = (int *)malloc(sizeof(int) * (cur_mcgnode_index*2));
    MaxStackElems = cur_mcgnode_index*2;

    searchstack[0] = n1;
    nelem_in_stack = 1;

    while(cur_node != n2 && count < cur_mcgedge_index)
    {
        /* pick the last element of the stack */
        cur_node = searchstack[nelem_in_stack-1];

        nelem_in_stack --;

        if(cur_node == n2)
            return true;

        /***********************************************************/

        //pick one unvisited edges
        for(i = 0; i < nlist->mnodes[cur_node]->nedges; i++)
        {
            if(elist->edges[nlist->mnodes[cur_node]->graph_edges[i]]->visited
                ||elist->edges[nlist->mnodes[cur_node]->graph_edges[i]]->cancel)
                continue;

            /*consider outgoing edges only*/
            if(elist->edges[nlist->mnodes[cur_node]->graph_edges[i]]->node_index2 == cur_node)
                continue;

            //set the node_index2 of the edge as the next cur_node

            if(nelem_in_stack >= MaxStackElems) /*extend it*/
            {
                searchstack = (int *)realloc(searchstack, sizeof(int)*(MaxStackElems+50));

                if(searchstack == NULL)
                    exit(-1);

                MaxStackElems += 50;
            }

            searchstack[nelem_in_stack] = elist->edges[nlist->mnodes[cur_node]->graph_edges[i]]->node_index2;
            nelem_in_stack++;

            elist->edges[nlist->mnodes[cur_node]->graph_edges[i]]->visited = true;
            count++;
        }

        if(nelem_in_stack == 0)
            return false;

    }

    if(count < cur_mcgedge_index)
    {
        return true;
    }
    return false;
}


/*
Add a new edge into the edge list
*/
void MCG_Graph::add_to_mcg_edgelist(int node1, int node2)
{
    if(elist->isFull())
    {
        int oldnum = elist->curMaxNumGedges;
        if(!elist->extend())
        {
            /*report errors*/
            char rout[256], var[256];
            sprintf(rout, "%s", "MCG_Graph::add_to_mcg_edgelist");
            sprintf(var, "%s", "elist");

            //write_mem_error(rout, var, 1);
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

    elist->edges[cur_mcgedge_index]->edge_index = cur_mcgedge_index;
    elist->edges[cur_mcgedge_index]->node_index1 = node1;
    elist->edges[cur_mcgedge_index]->node_index2 = node2;
    elist->edges[cur_mcgedge_index]->cancel = false;
    elist->edges[cur_mcgedge_index]->visited = false;
    elist->edges[cur_mcgedge_index]->ntris=0;
    elist->edges[cur_mcgedge_index]->triangles=NULL;

    elist->nedges++;
    cur_mcgedge_index++;
}


void MCG_Graph::add_edge_to_node(int node, int edgeindex)
{
    if(node < 0)
        return;

    nlist->mnodes[node]->graph_edges = extend_link(nlist->mnodes[node]->graph_edges,
        nlist->mnodes[node]->nedges);

    nlist->mnodes[node]->graph_edges[nlist->mnodes[node]->nedges] = edgeindex;
    nlist->mnodes[node]->nedges++;

}


void MCG_Graph::reset_sccnodes_flags()
{
    int i;
    DirGraph_Node **sccnodes = morse_decomp->dg->nlist->dirnodes;

    for(i = 0; i < morse_decomp->dg->nlist->ndirnodes; i++)
    {
        sccnodes[i]->visited = 0;
        sccnodes[i]->rev_visited = false;
    }
}


/*Judge whether there has already one edge between them
Note: not consider double connection now 03/08/07
*/
bool MCG_Graph::is_connected_mcg(int node1, int node2)
{
    int i, cur_edge;

    for(i = 0; i < nlist->mnodes[node1]->nedges; i++)
    {
        cur_edge = nlist->mnodes[node1]->graph_edges[i];

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

/*
judge whether a new added link is proper or not.
The criterion for "proper" now is
source->saddle/sink
saddle->saddle/sink
The inputs are the node indices of the MCG
*/

bool MCG_Graph::is_valid_link(int from, int to)
{
    /*
       This is only a temporary solution
    */
    if (from < 0 || to < 0) return false;
    if(nlist->mnodes[from]->type == 1 || nlist->mnodes[to]->type == 0)
        return false;

    return true;
}
bool MCG_Graph::IsMCGN(int scc)
{
    for(int i=0;i<nlist->nmnodes;i++)
    {
        if(nlist->mnodes[i]->scc_index==scc)
            return true;
    }
    return false;
}

void MCG_Graph::cal_mcgedge_regions(void)
{
        /* The following provides the description of the algorithm:

    for each saddle like Morse set in MCG
        if it has edges associated with it
        compute the regions of it using forward and backward graph of MFIG
          consider all its edges
              case 1: if it is saddle to sink
                      grow the sink region using backward graph of MFIG
                      intersect the obtained region with the forward region of the saddle
                      save the intersect regions to the member variable *triangles of MCGEdge
              case 2: if it is source to saddle
                      grow the source region using forward graph of MFIG
                      intersect the obtained region with the backward region of the saddle
                      save the intersect regions to the member variable *triangles of MCGEdge
              case 3: if it is saddle to saddle
                      according to the direction of connection
                      use the forward region of the starting saddle to intersect
                      with the backward region of the ending saddle
                      save the intersect regions to the member variable *triangles of MCGEdge
    */

    /*declare variables*/

    int *forward_region = NULL;
    int curMaxNumForward;
    int nforward_tris = 0;
    int *backward_region = NULL;
    int nbackward_tris = 0;
    int curMaxNumBackward;
    int *other_region = NULL;
    int curMaxNumOther;
    int nother_tris = 0;
    int *intersect_region = NULL;
    int curMaxNumIntersect;
    int nintersect_tris = 0;

    /*Allocate the initial space for those variables*/
    curMaxNumForward = object->tlist.ntris;
    forward_region = (int*)malloc(sizeof(int)*curMaxNumForward);

    curMaxNumBackward = object->tlist.ntris;
    backward_region = (int*)malloc(sizeof(int)*curMaxNumBackward);

    curMaxNumOther = object->tlist.ntris;
    other_region = (int*)malloc(sizeof(int)*curMaxNumOther);

    curMaxNumIntersect = object->tlist.ntris;
    intersect_region = (int*)malloc(sizeof(int)*curMaxNumIntersect);

    /*initialize the mcg edges*/
    //int i, j;

    for(int i=0; i<cur_mcgedge_index; i++)
    {
        elist->edges[i]->visited = false;
        /*mcgedges[i].ntris = 0;
        mcgedges[i].triangles = NULL;*/

        if (elist->edges[i]->cancel)
            elist->edges[i]->visited = true;
    }

    /**/

    bool forward_backward = false;

    int othernode;

    /*  reset the flags of the directed nodes for region growing  */
    reset_sccnodes_flags();


    for(int i=0; i<cur_mcgnode_index; i++)
    {
        //compute only for the new edges
        //if(elist->edges[i]->ntris!=0)
        //	continue;
        if(nlist->mnodes[i]->type != 2)
            continue;

        if (nlist->mnodes[i]->cancelled)
            continue;

        /**/
        if(nlist->mnodes[i]->nedges == 0)
            continue;

        nforward_tris = nbackward_tris = 0;

        /*grow two regions for the saddle*/
        get_saddle_two_regions(nlist->mnodes[i]->scc_index, forward_region, nforward_tris, curMaxNumForward,
            backward_region, nbackward_tris, curMaxNumBackward);

        for(int j=0; j<nlist->mnodes[i]->nedges; j++)
        {
            if(elist->edges[nlist->mnodes[i]->graph_edges[j]]->visited)
                continue;

            elist->edges[nlist->mnodes[i]->graph_edges[j]]->visited = true;

            othernode = elist->edges[nlist->mnodes[i]->graph_edges[j]]->node_index2;
            if(othernode == i) /*if it is current saddle node*/
                othernode =elist->edges[nlist->mnodes[i]->graph_edges[j]]->node_index1;

            if (nlist->mnodes[othernode]->cancelled)
                continue;

            nintersect_tris = nother_tris = 0;

            /*case 1: the other node is a sink*/
            if(nlist->mnodes[othernode]->type == 1)
            {
                /*grow backward*/

                cal_forward_region_MFIG(nlist->mnodes[othernode]->scc_index, 1,
                    other_region, nother_tris, curMaxNumOther);

                forward_backward = true;

            }

            /*case 2: the other node is a source*/
            else if(nlist->mnodes[othernode]->type == 0)
            {
                /*grow forward*/
                cal_forward_region_MFIG(nlist->mnodes[othernode]->scc_index, 0,
                    other_region, nother_tris, curMaxNumOther);

                forward_backward = false;
            }

            /*case 3: the other node is a saddle*/
            else
            {
                /*use the direction to decide*/
                if(othernode == elist->edges[nlist->mnodes[i]->graph_edges[j]]->node_index1)
                {
                    /*othernode is a source*/
                    /*grow forward from it*/
                    cal_forward_region_MFIG(nlist->mnodes[i]->scc_index, 0,
                        other_region, nother_tris, curMaxNumOther);

                    /*compute the intersection*/

                    forward_backward = false;
                }
                else
                {
                    /*othernode is a sink*/
                    /*grow backward from it*/
                    cal_forward_region_MFIG(nlist->mnodes[othernode]->scc_index, 1,
                        other_region, nother_tris, curMaxNumOther);

                    /*compute the intersection*/
                    forward_backward = true;
                }
            }

            /*compute the intersection*/

            if(!forward_backward) /*forward, should use the backward region of the saddle*/
                    intersect_twoRegions(other_region, nother_tris,
                        backward_region, nbackward_tris,
                        intersect_region, nintersect_tris);
                    //intersect_twoRegions(other_region, nother_tris,
                    //	forward_region, nforward_tris,
                    //	intersect_region, nintersect_tris);
            else
                    intersect_twoRegions(other_region, nother_tris,
                        forward_region, nforward_tris,
                        intersect_region, nintersect_tris);
                    //intersect_twoRegions(other_region, nother_tris,
                    //	backward_region, nbackward_tris,
                    //	intersect_region, nintersect_tris);


            /*remove the redundant triangles*/
            copy_array_Int(morse_decomp->scclist->scccomponents[nlist->mnodes[othernode]->scc_index]->nodes, other_region,
                morse_decomp->scclist->scccomponents[nlist->mnodes[othernode]->scc_index]->nnodes);
            nother_tris = morse_decomp->scclist->scccomponents[nlist->mnodes[othernode]->scc_index]->nnodes;

            remove_redundant(intersect_region, nintersect_tris,
                other_region, nother_tris);

            copy_array_Int(morse_decomp->scclist->scccomponents[nlist->mnodes[i]->scc_index]->nodes, other_region,
                morse_decomp->scclist->scccomponents[nlist->mnodes[i]->scc_index]->nnodes);
            nother_tris = morse_decomp->scclist->scccomponents[nlist->mnodes[i]->scc_index]->nnodes;

            remove_redundant(intersect_region, nintersect_tris,
                other_region, nother_tris);

            /*copy the intersection region to the edge*/

            elist->edges[nlist->mnodes[i]->graph_edges[j]]->triangles =
                (int *)malloc(sizeof(int)*nintersect_tris);
            copy_array_Int(intersect_region, elist->edges[nlist->mnodes[i]->graph_edges[j]]->triangles, nintersect_tris);
            elist->edges[nlist->mnodes[i]->graph_edges[j]]->ntris = nintersect_tris;
        }
    }

    free(forward_region);
    free(backward_region);
    free(intersect_region);
    free(other_region);
}

void  MCG_Graph::get_saddle_two_regions(int scc_index,
    int *forward_region, int &nforward_tris, int &curMaxNumForward,
    int *backward_region, int &nbackward_tris, 	int &curMaxNumBackward)
{
    /*grow the forward region*/
    cal_forward_region_MFIG(scc_index, 0, forward_region, nforward_tris, curMaxNumForward);

    /*grow the backward region*/
    /*reverse the edges of MFIG*/

    cal_forward_region_MFIG(scc_index, 1, backward_region, nbackward_tris, curMaxNumForward);

    /*reverse the edges of MFIG*/
}

void  MCG_Graph::cal_forward_region_MFIG(int scc_index, int inverse,
                             int *repell_region, int &ntris_repell_region, int &curMaxNumTrisRepellRegion)
{

    //curMaxNumTrisRepellRegion = Object.nfaces;
    //if(repell_region != NULL)
    //	free(repell_region);

    //repell_region = (int*)malloc(sizeof(int)*curMaxNumTrisRepellRegion);


    if(inverse == 1)
        ReverseEdges();

    /*  we DO NOT need to set the flag for the whole mesh again and again!  */
    //reset_sccnodes_flags();
    int i, j, /*k,*/ cur_node, other_node;


    /*
    use all the triangles in the region
    */
    ntris_repell_region = 0;
    copy_array_Int(morse_decomp->scclist->scccomponents[scc_index]->nodes, repell_region,
        morse_decomp->scclist->scccomponents[scc_index]->nnodes);
    ntris_repell_region = morse_decomp->scclist->scccomponents[scc_index]->nnodes;




    /*  we may grow too much if we can not stop earlier  */

    for(i = 0; i < ntris_repell_region; i++)
    {
        /*expand current node along its outgoing edges*/
        cur_node = repell_region[i];

        if(morse_decomp->dg->nlist->dirnodes[cur_node]->visited == 1)
            continue;

        morse_decomp->dg->nlist->dirnodes[cur_node]->visited = 1;

        for(j = 0; j < morse_decomp->dg->nlist->dirnodes[cur_node]->nedges; j++)
        {
            if(morse_decomp->dg->elist->edges[morse_decomp->dg->nlist->dirnodes[cur_node]->edges[j]]->node_index2 == cur_node)
                continue;

            other_node = morse_decomp->dg->elist->edges[morse_decomp->dg->nlist->dirnodes[cur_node]->edges[j]]->node_index2;

            if(morse_decomp->dg->nlist->dirnodes[other_node]->sscomp_index != scc_index
                //&& morse_decomp->scclist->scccomponents[morse_decomp->dg->nlist->dirnodes[other_node]->sscomp_index]->node_index < 0 // Do not cross other Morse sets
                && morse_decomp->dg->nlist->dirnodes[other_node]->visited == 0
                /*&& !IsRepeated(repell_region, other_node, ntris_repell_region)*/
                && !morse_decomp->dg->nlist->dirnodes[other_node]->rev_visited)
            {

                    /*add the node into the array*/
                    repell_region[ntris_repell_region] = other_node;
                    ntris_repell_region++;

                    morse_decomp->dg->nlist->dirnodes[other_node]->rev_visited = true;
            }

            //else if (morse_decomp->dg->nlist->dirnodes[other_node]->sscomp_index == target_scc_index
            //	&& !morse_decomp->dg->nlist->dirnodes[other_node]->rev_visited)
            //{
            //		/*add the node into the array*/
            //		repell_region[ntris_repell_region] = other_node;
            //		ntris_repell_region++;

            //		morse_decomp->dg->nlist->dirnodes[other_node]->rev_visited = true;
            //}
        }
    }

    if(inverse == 1)
        ReverseEdges();

    /*  reset all the flag  */
    DirGraph_Node **dg_nlist = morse_decomp->dg->nlist->dirnodes;
    for(i=0; i<ntris_repell_region; i++)
    {
        dg_nlist[repell_region[i]]->visited = 0;
        dg_nlist[repell_region[i]]->rev_visited = false;
    }
}


void  MCG_Graph::intersect_twoRegions(int *source1, int num1,
                      int *source2, int num2,
                      int *dest, int &num_dest)
{
    int i, j;
    int cur_index = 0;
    int cur_triangle;


    for(i = 0; i < num1; i++)
    {
        cur_triangle = source1[i];

        for(j = 0; j < num2; j++)
        {
            if(cur_triangle == source2[j])
            {
                dest[cur_index] = cur_triangle;
                cur_index++;
            }
        }

    }

    num_dest = cur_index;
}

void  MCG_Graph::remove_redundant(int *source1, int &num1, int *source2, int num2)
{
    int i, j;
    int *temp;
    int cur_elems = 0;
    if(num1 > num2) temp=(int*)malloc(sizeof(int)*num1);
    else temp = (int*)malloc(sizeof(int)*num2);
    for(i=0; i<num1; i++)
    {
        for(j=0; j<num2; j++)
        {
            if(source1[i]==source2[j])
                break;
        }

        if(j>=num2)
        {
            temp[cur_elems] = source1[i];
            cur_elems ++;
        }
    }

    for(i=0; i<cur_elems; i++)
        source1[i] = temp[i];
    num1 = cur_elems;
}

void  MCG_Graph::ReverseEdges()
{
    int i, temp;

    for(i = 0; i < morse_decomp->dg->elist->nedges; i++)
    {
        temp =  morse_decomp->dg->elist->edges[i]->node_index1;
        morse_decomp->dg->elist->edges[i]->node_index1 = morse_decomp->dg->elist->edges[i]->node_index2;
        morse_decomp->dg->elist->edges[i]->node_index2 = temp;
    }
}
bool MCG_Graph::IsRepeated(int *a, int b, int num)
{
    for(int i = 0; i < num; i++)
    {
        if(a[i] == b)
            return true;
    }

    return false;
}

int  MCG_Graph::select_morse(double min_pri)
{
    int morse=-1;

    //no pq

    //now the variable "morse" is useless.
    if(NeedRefined.size())NeedRefined.clear();

    //Now add morse set into this queue
    //now we need to scan all the morse sets
    for(int i=0;i<nlist->nmnodes;i++)
    {
        //if discard, don't add
        if(NeedRefinement(i))
            NeedRefined.push_back(i);
    }



    //double max_pri=0/*min_pri*/;
    double max_pri=min_pri;
    for(int i=0;i<NeedRefined.size();i++)
    {
        if(nlist->mnodes[NeedRefined[i]]->priority>max_pri)
        {
            max_pri=nlist->mnodes[NeedRefined[i]]->priority;
            morse=NeedRefined[i];
        }
    }
    return morse;

}

void  MCG_Graph::AddQueue(int morse)
{
    //now the variable "morse" is useless.
    if(NeedRefined.size())NeedRefined.clear();

    //Now add morse set into this queue
    //now we need to scan all the morse sets
    for(int i=0;i<nlist->nmnodes;i++)
    {
        //if discard, don't add
        //if(NeedRefinement(i))
            NeedRefined.push_back(i);
    }
}


bool MCG_Graph::NeedRefinement(int morse)
{
    //if exclude
    if (nlist->mnodes[morse]->cancelled) return false;

    int cur_scc=nlist->mnodes[morse]->scc_index;
    for(int i=0;i<morse_decomp->scclist->scccomponents[cur_scc]->nnodes;i++)
    {
        if(object->tlist.tris[morse_decomp->scclist->scccomponents[cur_scc]->nodes[i]]->exclude)
            return false;
    }
    if(morse_decomp->scclist->scccomponents[cur_scc]->nnodes<3)
        return false;
    return true;
}


void
MCG_Graph::grow_saddle_mcgnodes()
{
    int i;

    clock_t l_start, l_finish;

    for(i = 0; i < nlist->nmnodes/*cur_mcgnode_index*/; i++)
    {
        if (nlist->mnodes[i]->type == 2)
        {
            grow_from_MorseSet (i, true,  morse_decomp);  // grow forward
            grow_from_MorseSet (i, false, morse_decomp);  // grow backward
        }

    }
}



void
MCG_Graph::grow_other_mcgnodes()
{
    int i;

    for (i=0; i<nlist->nmnodes; i++)
    {
        if (nlist->mnodes[i]->type == 2)
            continue;

        else if (nlist->mnodes[i]->type == 0) // source like
        {
            if (nlist->mnodes[i]->nedges == 0)
                grow_from_MorseSet (i, true, morse_decomp);
        }

        else if (nlist->mnodes[i]->type == 1) // sink like
        {
            if (nlist->mnodes[i]->nedges == 0)
                grow_from_MorseSet (i, false, morse_decomp);
        }
    }
}

/*
   Remove indirect links in MCG 02/21/2010
*/
//void
//MCG_Graph::remove_indir_mcgedges()
//{
//	int i, j;
//	int *expandnodes = new int[cur_mcgnode_index];
//	int num_expanded_nodes;
//	int cur_node, other_node;
//
//	for(i = 0; i < cur_mcgnode_index; i++)
//	{
//		//if (nlist->mnodes[i]->type != 2)
//		//	continue;
//
//		for(j = 0; j < nlist->mnodes[i]->nedges; j++)
//		{
//			other_node = elist->edges[nlist->mnodes[i]->graph_edges[j]]->node_index2;
//			if(other_node == i)
//				continue;  /*consider outgoing edges only*/
//
//			if (nlist->mnodes[i]->type != 2 && nlist->mnodes[other_node]->type != 2)
//			{
//				if(has_longer_path_in_mcg(i, other_node, nlist->mnodes[i]->graph_edges[j])) // we need to change this !!
//					elist->edges[nlist->mnodes[i]->graph_edges[j]]->cancel = true;
//			}
//
//			// if there is a direct path using dfs (not passing other mcg nodes) between two nodes, leave it; otherwise, cancel it
//
//			else if (!has_dir_path_between_MorseSets(i, other_node, morse_decomp))
//				elist->edges[nlist->mnodes[i]->graph_edges[j]]->cancel = true;
//		}
//	}
//
//	/*update the edge list in the nodes*/
//	for(i = 0; i < cur_mcgedge_index; i++)
//	{
//		if(elist->edges[i]->cancel)
//		{
//			del_one_edge_from(nlist->mnodes[elist->edges[i]->node_index1]->graph_edges,
//				nlist->mnodes[elist->edges[i]->node_index1]->nedges, i);
//			del_one_edge_from(nlist->mnodes[elist->edges[i]->node_index2]->graph_edges,
//				nlist->mnodes[elist->edges[i]->node_index2]->nedges, i);
//		}
//	}
//	delete [] expandnodes;
//}

/*
Note that this assumes that the strongly connected components are from the directed graph stored in the input Morse decomposition variable
Therefore, we need to be careful if the input morse_decomp is a local result
*/
void
MCG_Graph::grow_from_MorseSet(int MorseID, bool back_forth, MorseDecomp *morse_decomp)
{
    /*
        obtain the corresponding triangle region
    */

    if (MorseID < 0 || MorseID >= nlist->nmnodes) return;

    int scc_id = nlist->mnodes[MorseID]->scc_index;

    unsigned char start_type = nlist->mnodes[MorseID]->type;

    DirGraph_NodeList *dg_nlist = morse_decomp->dg->nlist;
    Graph_EdgeList *dg_elist = morse_decomp->dg->elist;

    SCComponent *start_scc = morse_decomp->scclist->scccomponents[scc_id];

    int i, j;

    for (i=0; i<dg_nlist->ndirnodes; i++)
        dg_nlist->dirnodes[i]->visited = 0;

    /*
        copy the triangles of the strongly connected component corresponding to the given Morse set
    */

    DynList_Int *reg1 = new DynList_Int(start_scc->nnodes*3);

    for (i=0; i<start_scc->nnodes; i++)
        reg1->add_New_2 (start_scc->nodes[i]);

    for (i=0; i<start_scc->nnodes; i++)
        dg_nlist->dirnodes[start_scc->nodes[i]]->visited = 1;

    /*
       start growing from the initial region
    */

    int conn_scc_id;

    for (i=0; i<reg1->nelems; i++)
    {
        //if (dg_nlist->dirnodes[reg1->elems[i]]->visited > 0)
        //	continue;

        /*  get the current element  */
        DirGraph_Node *cur_n = dg_nlist->dirnodes[reg1->elems[i]];

        for (j=0; j<cur_n->nedges; j++)
        {
            Graph_Edge *cur_adj_e = dg_elist->edges[cur_n->edges[j]];

            int adj_n;

            if (back_forth)
                adj_n = cur_adj_e->node_index2;  // consider outgoing direction

            else
                adj_n = cur_adj_e->node_index1;  // consider the incoming direction

            if (adj_n == cur_n->node_index || dg_nlist->dirnodes[adj_n]->visited > 0)
                continue;

            /*
                add to the list of the growing region
            */
            reg1->add_New_2(adj_n);

            dg_nlist->dirnodes[adj_n]->visited = 1;  // to mark that the node is already in the region

            /*
                check and see whether we reach another Morse set
            */

            conn_scc_id = dg_nlist->dirnodes[adj_n]->sscomp_index;

            if (conn_scc_id != scc_id && conn_scc_id >= 0
                && morse_decomp->scclist->scccomponents[conn_scc_id]->node_index>=0 /* need to verify it is a valid Morse set*/)
            {
                /*
                    the growing region reaches another Morse set
                */
                /*
                    In order to get a better connection region, we need to grow the region further until it contains the second Morse set
                    When should we stop?
                    Let us not consider that right now and return to this issue later.
                */

                /*
                    we need to see whether one of these two Morse sets is saddle-like
                */

                SCComponent *con_scc = morse_decomp->scclist->scccomponents[conn_scc_id];
                DynList_Int *con_reg = new DynList_Int();
                bool add_con_reg = false;

                if (start_type == 2 || nlist->mnodes[con_scc->node_index]->type == 2)
                {
                    /*
                      if one of the two Morse sets is saddle-like, we need to compute the connection region between them
                    */
                    add_con_reg = true;
                    /*
                        We actually grow region backward until we reach the start Morse set
                    */
                    DynList_Int *rev_reg = new DynList_Int();
                    grow_reg_from_to (conn_scc_id, scc_id, !back_forth, morse_decomp, rev_reg);


                    /*
                        we then compute the intersection region as the desired connection region
                        NOTE: we need to remove the triangles of the two Morse sets for better visualization! 02/21/2010
                    */
                    con_reg = get_intersect_reg(rev_reg, reg1, morse_decomp->dg, con_reg);

                    delete rev_reg;
                }

                /*
                    construct a new edge if there is no one connecting them
                    if it is a saddle-like connection, we need to store the connection region as well
                */

                /*
                    please make sure there is no redundant edges!
                */


                if (add_con_reg)
                {
                    if (back_forth)
                        if(is_valid_link(MorseID, con_scc->node_index) && !is_connected_mcg(MorseID, con_scc->node_index))
                        {
                            add_to_mcg_edge_list(MorseID, con_scc->node_index, con_reg);
                            add_edge_to_node(MorseID, cur_mcgedge_index-1);
                            add_edge_to_node(con_scc->node_index, cur_mcgedge_index-1);
                        }
                    else
                        if(is_valid_link(con_scc->node_index, MorseID) && !is_connected_mcg(con_scc->node_index, MorseID))
                        {
                            add_to_mcg_edge_list(con_scc->node_index, MorseID, con_reg);
                            add_edge_to_node(MorseID, cur_mcgedge_index-1);
                            add_edge_to_node(con_scc->node_index, cur_mcgedge_index-1);
                        }
                }
                else
                {
                    if (back_forth)
                        if(is_valid_link(MorseID, con_scc->node_index) && !is_connected_mcg(MorseID, con_scc->node_index))
                        {
                            add_to_mcg_edge_list(MorseID, con_scc->node_index);
                            add_edge_to_node(MorseID, cur_mcgedge_index-1);
                            add_edge_to_node(con_scc->node_index, cur_mcgedge_index-1);
                        }
                    else
                        if(is_valid_link(con_scc->node_index, MorseID) && !is_connected_mcg(con_scc->node_index, MorseID))
                        {
                            add_to_mcg_edge_list(con_scc->node_index, MorseID);
                            add_edge_to_node(MorseID, cur_mcgedge_index-1);
                            add_edge_to_node(con_scc->node_index, cur_mcgedge_index-1);
                        }
                }

                /*
                    after copying the connection region to the local variable of the MCG edge, we can delete it
                */
                delete con_reg;

                /*
                    we need to copy the whole region of the reaching Morse set to the region variable reg1
                */

                for (int k=0; k<con_scc->nnodes; k++)
                {
                    if (dg_nlist->dirnodes[con_scc->nodes[k]]->visited > 0)
                        continue;

                    reg1->add_New_2 (con_scc->nodes[k]);
                    dg_nlist->dirnodes[con_scc->nodes[k]]->visited = 1;
                }
            }
        }
    }

    delete reg1;
}


DynList_Int
*MCG_Graph::get_intersect_reg (DynList_Int *reg1, DynList_Int *reg2, DirGraph *dg, DynList_Int *intersect_reg)
{
    int i;

    for (i=0; i<reg1->nelems; i++)
        dg->nlist->dirnodes[reg1->elems[i]]->parent = 0;

    for (i=0; i<reg1->nelems; i++)
        dg->nlist->dirnodes[reg1->elems[i]]->parent ++;

    for (i=0; i<reg2->nelems; i++)
        dg->nlist->dirnodes[reg2->elems[i]]->parent ++;

    for (i=0; i<reg1->nelems; i++)
    {
        if (dg->nlist->dirnodes[reg1->elems[i]]->parent == 2)
            intersect_reg->add_New_2(reg1->elems[i]);
    }

    return intersect_reg;
}



DynList_Int
*MCG_Graph::grow_reg_from_to (int scc1, int scc2, bool back_forth,MorseDecomp *morse_decomp, DynList_Int *rev_reg)
{
    int i, j;
    /*
        First, initialize the region
    */
    DirGraph_NodeList *dg_nlist = morse_decomp->dg->nlist;
    Graph_EdgeList *dg_elist = morse_decomp->dg->elist;

    SCComponent *scc_tris1 = morse_decomp->scclist->scccomponents[scc1];
    SCComponent *scc_tris2 = morse_decomp->scclist->scccomponents[scc2];

    for (i=0; i<dg_nlist->ndirnodes; i++)
        dg_nlist->dirnodes[i]->rev_visited = false;

    /*
        Now growing the region according to the back_forth flag
    */
    for (i=0; i<scc_tris1->nnodes; i++)
    {
        rev_reg->add_New_2 (scc_tris1->nodes[i]);
        dg_nlist->dirnodes[scc_tris1->nodes[i]]->rev_visited = true;
    }

    int reach_ntris_scc2 = 0;

    for (i=0; i<rev_reg->nelems; i++)
    {
        DirGraph_Node *cur_n = dg_nlist->dirnodes[rev_reg->elems[i]];

        //if (cur_n->rev_visited)
        //	continue;

        for (j=0; j<cur_n->nedges; j++)
        {
            Graph_Edge *cur_adj_e = dg_elist->edges[cur_n->edges[j]];

            int adj_n;
            if (back_forth)
            {
                adj_n = cur_adj_e->node_index2;
            }
            else
                adj_n = cur_adj_e->node_index1;

            if (adj_n == cur_n->node_index || dg_nlist->dirnodes[adj_n]->rev_visited)
                continue;

            rev_reg->add_New_2(adj_n);

            dg_nlist->dirnodes[adj_n]->rev_visited = true;

            if (dg_nlist->dirnodes[adj_n]->sscomp_index == scc2)
                reach_ntris_scc2 ++;

            if (reach_ntris_scc2 >= scc_tris2->nnodes/2+1)
                return rev_reg;
        }
    }

    return rev_reg;
}

void
MCG_Graph::add_to_mcg_edge_list (int from, int to, DynList_Int *con_reg )
{
    if(elist->isFull())
    {
        int oldnum = elist->curMaxNumGedges;
        if(!elist->extend())
        {
            /*report errors*/
            char rout[256], var[256];
            sprintf(rout, "%s", "MCG_Graph::add_to_mcg_edgelist");
            sprintf(var, "%s", "elist");

            //write_mem_error(rout, var, 1);
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

    elist->edges[cur_mcgedge_index]->edge_index = cur_mcgedge_index;
    elist->edges[cur_mcgedge_index]->node_index1 = from;
    elist->edges[cur_mcgedge_index]->node_index2 = to;
    elist->edges[cur_mcgedge_index]->cancel = false;
    elist->edges[cur_mcgedge_index]->visited = false;

    if (con_reg == NULL)
    {
        elist->edges[cur_mcgedge_index]->ntris=0;
        elist->edges[cur_mcgedge_index]->triangles=NULL;
    }

    else
    {
        int i;

        elist->edges[cur_mcgedge_index]->triangles = (int*)malloc(sizeof(int)*con_reg->nelems);

        for (i=0; i<con_reg->nelems; i++)
        {
            elist->edges[cur_mcgedge_index]->triangles[i] = con_reg->elems[i];
        }

        elist->edges[cur_mcgedge_index]->ntris = con_reg->nelems;
    }

    elist->nedges++;
    cur_mcgedge_index++;
}




/*
    construct local MCG within a given triangular region
*/

void
MCG_Graph::grow_all_mcgnodes_local(MorseDecomp *local_decomp)
{
    int i;

    clock_t l_start, l_finish;

    for(i = 0; i < nlist->nmnodes/*cur_mcgnode_index*/; i++)
    {
        //l_start = clock();

        if (nlist->mnodes[i]->cancelled) continue;

        grow_repeller_region_graph_local(local_decomp, nlist->mnodes[i]->scc_index, 0);

    }
}

void
MCG_Graph::grow_repeller_region_graph_local(MorseDecomp *morse_decomp, int scc_index, int inverse)
{
    int i, j, k, cur_node, other_node;
    DirGraph_Node **sccnodes = morse_decomp->dg->nlist->dirnodes;
    Graph_Edge **sccedges = morse_decomp->dg->elist->edges;

    /*reset all the flags of the visited nodes*/
    for(i = 0; i < ntris_repell_region; i++)
    {
        sccnodes[repell_region[i]]->visited = 0;
    }

    if(repell_region != NULL)
        free(repell_region);

    curMaxNumTrisRepellRegion = object->tlist.ntris;
    repell_region = (int*)malloc(sizeof(int)*curMaxNumTrisRepellRegion);
    ntris_repell_region = 0;

    /*
    use all the triangles in the region
    */
    copy_array_Int(morse_decomp->scclist->scccomponents[scc_index]->nodes, repell_region,
        morse_decomp->scclist->scccomponents[scc_index]->nnodes);
    ntris_repell_region = morse_decomp->scclist->scccomponents[scc_index]->nnodes;


    //reset_sccnodes_flags();

    for(i = 0; i < ntris_repell_region; i++)
    {
        /*expand current node along its outgoing edges*/
        cur_node = repell_region[i];

        if(sccnodes[cur_node]->visited == 2)
            continue;

        sccnodes[cur_node]->visited = 2;

        for(j = 0; j < sccnodes[cur_node]->nedges; j++)
        {
            if(sccedges[sccnodes[cur_node]->edges[j]]->node_index2 == cur_node)
                continue;

            other_node = sccedges[sccnodes[cur_node]->edges[j]]->node_index2;



            if(sccnodes[other_node]->sscomp_index != scc_index
                && sccnodes[other_node]->visited == 0
                /*&& !is_repeated_elem(repell_region, other_node, ntris_repell_region)*/
                )
            {
                /*judge whether it corresponds to other valid Morse set*/
                /*if((morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->nnodes <= 2
                    && morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->nfixedpoints > 0)
                    ||(morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->nnodes > 2
                    && morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->valid == true))*/
                if(IsMCGN(sccnodes[other_node]->sscomp_index)
                    && morse_decomp->scclist->scccomponents[scc_index]->node_index>=0
                && morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index>=0)
                {
                    if(!is_connected_mcg(morse_decomp->scclist->scccomponents[scc_index]->node_index,
                        morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index)
                        //&& !has_interval_mcg(scclist.scccomponents[scc_index].node_index,
                        //scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index)

                        && !nlist->mnodes[morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index]->cancelled
                        )
                    {
                        if(inverse == 0)
                        {
                            if(is_valid_link(morse_decomp->scclist->scccomponents[scc_index]->node_index,
                                morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index))
                            {
                                add_to_mcg_edgelist(morse_decomp->scclist->scccomponents[scc_index]->node_index,
                                    morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index);
                                /*add the edge to the corresponding nodes*/
                                add_edge_to_node(morse_decomp->scclist->scccomponents[scc_index]->node_index, cur_mcgedge_index-1);
                                add_edge_to_node(morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index,
                                    cur_mcgedge_index-1);
                            }
                        }
                        else{
                            if(is_valid_link(morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index,
                                morse_decomp->scclist->scccomponents[scc_index]->node_index))
                            {
                                add_to_mcg_edgelist(morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index,
                                    morse_decomp->scclist->scccomponents[scc_index]->node_index);
                                /*add the edge to the corresponding nodes*/
                                add_edge_to_node(morse_decomp->scclist->scccomponents[scc_index]->node_index, cur_mcgedge_index-1);
                                add_edge_to_node(morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index,
                                    cur_mcgedge_index-1);
                            }
                        }

                        /*add all the nodes in this Morse sets into the array*/
                        for(k = 0; k < morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->nnodes; k++)
                        {
                            if(sccnodes[morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->nodes[k]]->visited
                                == 1)
                                continue;
                            repell_region[ntris_repell_region] =
                                morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->nodes[k];
                            /*set the flag*/
                            sccnodes[repell_region[ntris_repell_region]]->visited = 1;
                            ntris_repell_region++;
                        }

                    }
                }
                else{
                    /*add the node into the array*/
                    repell_region[ntris_repell_region] = other_node;
                            sccnodes[repell_region[ntris_repell_region]]->visited = 1;
                    ntris_repell_region++;
                }
            }
        }
    }
}

void
MCG_Graph::build_mcg_edges_graph_local(MorseDecomp *local_decomp)
{
    cur_mcgedge_index = 0;

    ///*
    //    set the triangles inside the given triangle region
    //*/
    //int i;
    //for (i=0; i<ntris; i++)
    //{
    //	//morse_decomp
    //}

    /*grow all nodes without knowing their types*/
    grow_all_mcgnodes_local(local_decomp);

    /*remove any redundant edges*/
    remove_redundant_edges();
}



/**********************************************************************
  02/28/2010
  save the obtained global MCG
  Format:
  #FG nodes
  #mnodes
  #medges
  // node list
  # triangles
  triangle list
  ...
  # triangles
  triangle list
  ...

  ...
  // edge list
  node1 node2
  node1 node2
  node1 node2
  ...
*/
void
MCG_Graph::save_MCG(char *filename)
{
}

void
MCG_Graph::load_MCG(char *filename)
{
}



/*
    Added by Guoning to output the Morse sets 04/04/2011
*/
void output_MS(char *filename, MCG_Graph *mcg, MorseDecomp *md)
{
    FILE *fp = fopen (filename, "w");

    fprintf (fp, "#Morse_sets: %d\n", mcg->nlist->nmnodes);

    int i, j;

    for (i=0; i<mcg->nlist->nmnodes; i++)
    {
        MCG_Node *mn=mcg->nlist->mnodes[i];

        fprintf(fp, "%d\n", i);

        SCComponent *sc = md->scclist->scccomponents[mn->scc_index];

        //fprintf(fp, "%d %d %d\n", sc->Conley0, sc->Conley1, sc->Conley2);

        //if (sc->nnodes<=2 && mn->conley[0]==0 && mn->conley[1]==0 && mn->conley[2]==0)
        //{
        //	fprintf(fp, "%d %d %d\n", 0, 1, 0);
        //}
        //else
            fprintf(fp, "%d %d %d\n", mn->conley[0], mn->conley[1], mn->conley[2]);

        fprintf(fp, "%d\n", sc->nnodes);

        for (j=0; j<sc->nnodes; j++)
        {
            fprintf (fp, "%d ", sc->nodes[j]);
        }

        fprintf(fp, "\n");
    }

    fclose(fp);
}

