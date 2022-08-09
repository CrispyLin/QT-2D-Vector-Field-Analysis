#include "RegionTauMap.h"
#include "Others/common_routines.h"
#include <vector>

using namespace std;

#define COMP_LOCAL_MCG

extern MorseDecomp *morse_decomp;
extern MCG_Graph *mcg;
extern Polyhedron *object;
extern MorseDecomp *local_decomp;
int refined_morse;//added 12/04/2009 by Qingqing
extern SCCList *tscclist;
int newSCCN;
int orEndSCC;
int or_scc;
int orEndDirE;//the id of dirnodes after delete the edges inside the morse set
int orEndMCGN;//the id of mcgnodes before refinement
int orEndMCGE;//the id of mcgedges before refinement

extern MCG_Graph *lmcg;

extern double used_tau ;
extern int Cal_Regions;

bool enFastRefinement = false;
extern double edge_sample_error ;


RegionTauMap::RegionTauMap(void)
    : /*conley(nullptr)
, */tedges(nullptr)
    , curMaxNumTempEdges(0)
    , num_tedges(0)
    , cur_end_tedge(0)
{
    //conley= new ConleyIndex();
}

RegionTauMap::~RegionTauMap(void)
{
    //delete conley;
}

void RegionTauMap::RefineRegion(int id, double tau)
{

    //conley->build_morse_nlist(id);
    int scc=mcg->nlist->mnodes[id]->scc_index;
    //conley->build_morse_nlist(scc,0);
    build_morse_nlist(scc);
    //output nodes done
    refined_morse=id;
    //delete the edges in the directed graph
    CleanRegion(id);

    //Flow Combinatorialization

    local_flow_comb(tau);


    /*
        Guoning's comments 02/18/2010
        Add the checking of connectivity here
        call the function local_decomp->obtain_connected_imgs_reg (...)
    */

    //build local object
    BuildLocalRegion();

    //Morse Set Identification

    local_decomp->dg->find_SCCS();
    //object->capture_Singularities();  // added by Guoning 02/16/2010
    local_decomp->build_SCCElemList/*build_SCCElemList_local*/();
    //local_decomp->mark_all_valid_SCCS();  // when we verify the validity of the SCCs, we need to consider its original triangle ids !!

    /*   make a triangle list   */
    int *tris = new int[morse_nlist.size()];
    int ntris = morse_nlist.size();

    for (int i=0; i<morse_nlist.size(); i++)
        tris[i] = morse_nlist[i];

    local_decomp->mark_all_valid_SCCS(tris, ntris);

    delete [] tris;


    //////////////////////////////////////////////
    //Print the local scc
    //////////////////////////////////////////////
    //FILE *local=fopen("local.txt","w");
    //fprintf(local,"scc triangles:\n");
    //for(int i=0;i<local_decomp->scclist->nsccs;i++)
    //{
    //	fprintf(local,"\nscc%d:\n",i);
    //	for(int j=0;j<local_decomp->scclist->scccomponents[i]->nnodes;j++)
    //	{
    //		fprintf(local,"%d\t", local_decomp->scclist->scccomponents[i]->nodes[j]);
    //	}
    //}
    //fclose(local);
    //////////////////////////////////////////////
    //Print the local scc
    //////////////////////////////////////////////

    //copy scclist to global scclist
    LocalToGlobalSCC();


#ifdef COMP_LOCAL_MCG
    //Morse Set Classification
    //MCG Construction
    BuildLocalMCG();

    ////////////////////////////////
    //Print the local mcg
    ////////////////////////////////
    //FILE *locmcg=fopen("local_mcg.txt","w");
    //	fprintf(locmcg,"local_mcg %d\n", lmcg->nlist->nmnodes);
    //	for(int i=0;i<lmcg->nlist->nmnodes;i++)
    //	{
    //		fprintf(locmcg,"%d\t%d\n",i,lmcg->nlist->mnodes[i]->scc_index);
    //	}
    //	fclose(locmcg);
    ////////////////////////////////
    //Print the local mcg
    ////////////////////////////////

    MergeMCG();


//FILE* mcgN=fopen("local_mcg.txt","w");
//fprintf(mcgN,"source#\tsink#\tsaddle#\n");
//fprintf(mcgN,"%d\t%d\t%d\n\n",lmcg->r_counter,lmcg->a_counter,lmcg->s_counter);

//fprintf(mcgN,"mcgNid\tSCCid\ttype\tmcg_labelindex\n");
//for(int i=0;i<lmcg->nlist->nmnodes;i++)
//	fprintf(mcgN,"%d\t%d\t%u\t%d\n",i,lmcg->nlist->mnodes[i]->scc_index,lmcg->nlist->mnodes[i]->type,lmcg->nlist->mnodes[i]->labelindex);


//fprintf(mcgN,"\nmcgEid\tNodeId1\tNodeId2\n");
//for(int i=0;i<lmcg->elist->nedges;i++)
//	fprintf(mcgN,"%d\t%d\t%d\n",i,lmcg->nlist->mnodes[lmcg->elist->edges[i]->node_index1]->global_index,lmcg->nlist->mnodes[lmcg->elist->edges[i]->node_index2]->global_index);
//fclose(mcgN);


/////////////////////////
//global mcg computing
/////////////////////////
//if(mcg != nullptr)
//	delete mcg;
//mcg = new MCG_Graph();
//mcg->init_MCG();
//mcg->build_mcg();
/////////////////////////
//global mcg computing
/////////////////////////

// construct the connection region for planar example? 03/07/2010
//mcg->cal_mcgedge_regions();
#endif

    free (tedges);
    tedges = nullptr;

    //FILE *local=fopen("localmcginf.txt","w");
    //fprintf(local,"Finish merging MCG: %d\n", refined_morse);
    //fclose(local);
}


bool RegionTauMap::RefineRegion_forAuto(int id, double tau)
{

    //conley->build_morse_nlist(id);
    int scc=mcg->nlist->mnodes[id]->scc_index;
    //conley->build_morse_nlist(scc,0);
    build_morse_nlist(scc);
    //output nodes done
    refined_morse=id;
    //delete the edges in the directed graph
    CleanRegion(id);

    //Flow Combinatorialization


    //  for debugging purpose 07/05/2010
    clock_t start, end;
    FILE *fp=nullptr;

    //#ifdef OUTPUT_PERFORMANCE
    //		fp = fopen("refining_morseSet.txt", "a");
    //		fprintf(fp, "*** start computing local flow combinatorialization.\n");
    //		fclose(fp);
    //		start = clock();
    //	//}
    //#endif

    //local_flow_comb(tau);
    num_tedges = 0;
    build_local_connected_graph(tau);
    //build_local_connected_graph_2(tau);

    //#ifdef OUTPUT_PERFORMANCE
    //		end = clock();
    //		fp = fopen("refining_morseSet.txt", "a");
    //		fprintf(fp, "finish computing local flow combinatorialization.\n");
    //		fprintf(fp, "it took %f seconds.\n",(double)(end - start)/CLOCKS_PER_SEC);
    //		fclose(fp);
    //#endif

    /*
        Guoning's comments 02/18/2010
        Add the checking of connectivity here
        call the function local_decomp->obtain_connected_imgs_reg (...)
    */

    //#ifdef OUTPUT_PERFORMANCE
    //		fp = fopen("refining_morseSet.txt", "a");
    //		fprintf(fp, "*** start building a local directed graph.\n");
    //		fclose(fp);
    //		start = clock();
    //#endif

    //build local object
    //BuildLocalRegion();
    BuildLocalRegion_2();

    //#ifdef OUTPUT_PERFORMANCE
    //		end = clock();
    //		fp = fopen("refining_morseSet.txt", "a");
    //		fprintf(fp, "finish building the local directed graph.\n");
    //		fprintf(fp, "it took %f seconds.\n",(double)(end - start)/CLOCKS_PER_SEC);
    //		fclose(fp);
    //#endif

    //Morse Set Identification

    //#ifdef OUTPUT_PERFORMANCE
    //		fp = fopen("refining_morseSet.txt", "a");
    //		fprintf(fp, "*** start extracting SCCs.\n");
    //		fclose(fp);
    //		start = clock();
    //#endif

    local_decomp->dg->find_SCCS();

    //#ifdef OUTPUT_PERFORMANCE
    //		end = clock();
    //		fp = fopen("refining_morseSet.txt", "a");
    //		fprintf(fp, "finish extracting SCCs.\n");
    //		fprintf(fp, "it took %f seconds.\n",(double)(end - start)/CLOCKS_PER_SEC);
    //		fclose(fp);
    //#endif

    //object->capture_Singularities();  // added by Guoning 02/16/2010
    //#ifdef OUTPUT_PERFORMANCE
    //		fp = fopen("refining_morseSet.txt", "a");
    //		fprintf(fp, "*** start building local SCC list.\n");
    //		fclose(fp);
    //		start = clock();
    //#endif

    local_decomp->build_SCCElemList_local();

    //#ifdef OUTPUT_PERFORMANCE
    //		end = clock();
    //		fp = fopen("refining_morseSet.txt", "a");
    //		fprintf(fp, "finish building local SCC list.\n");
    //		fprintf(fp, "it took %f seconds.\n",(double)(end - start)/CLOCKS_PER_SEC);
    //		fclose(fp);
    //#endif

    /*   make a triangle list   */
    int *tris = new int[morse_nlist.size()];
    int ntris = morse_nlist.size();

    for (int i=0; i<morse_nlist.size(); i++)
        tris[i] = morse_nlist[i];

    if (!local_decomp->mark_all_valid_SCCS(tris, ntris))
    {
        delete [] tris;
        //free (tedges);
        //tedges = nullptr;
        return false;
    }

    delete [] tris;

    /// if there is a SCC which is supposed to be a Morse set formed by a disconnected region, we return fail

    //#ifdef OUTPUT_PERFORMANCE
    //		fp = fopen("refining_morseSet.txt", "a");
    //		fprintf(fp, "*** start merging SCCs.\n");
    //		fclose(fp);
    //		start = clock();
    //#endif

    //copy scclist to global scclist
    LocalToGlobalSCC();

//#ifdef OUTPUT_PERFORMANCE
//		end = clock();
//		fp = fopen("refining_morseSet.txt", "a");
//		fprintf(fp, "finish merging SCCs.\n");
//		fprintf(fp, "it took %f seconds.\n",(double)(end - start)/CLOCKS_PER_SEC);
//		fclose(fp);
//#endif

#ifdef COMP_LOCAL_MCG
    //Morse Set Classification
    //MCG Construction
    //#ifdef OUTPUT_PERFORMANCE
    //		fp = fopen("refining_morseSet.txt", "a");
    //		fprintf(fp, "*** start merging MCG.\n");
    //		fclose(fp);
    //		start = clock();
    //#endif

    BuildLocalMCG();

    //#ifdef OUTPUT_PERFORMANCE
    //		end = clock();
    //		fp = fopen("refining_morseSet.txt", "a");
    //		fprintf(fp, "finish building local MCG and start merging.\n");
    //		fprintf(fp, "there are %d nodes and %d edges.\n", lmcg->nlist->nmnodes, lmcg->elist->nedges);
    //		fprintf(fp, "the refined Morse set has %d edges.\n", mcg->nlist->mnodes[refined_morse]->nedges);
    //		fprintf(fp, "the total number of edges in the original MCG is %d\n", mcg->elist->nedges);
    //		fprintf(fp, "it took %f seconds.\n",(double)(end - start)/CLOCKS_PER_SEC);
    //		fclose(fp);
    //#endif

    MergeMCG();

//#ifdef OUTPUT_PERFORMANCE
//		end = clock();
//		fp = fopen("refining_morseSet.txt", "a");
//		fprintf(fp, "finish merging MCG.\n");
//		fprintf(fp, "the total number of edges in current MCG is %d\n", mcg->elist->nedges);
//		fprintf(fp, "it took %f seconds.\n",(double)(end - start)/CLOCKS_PER_SEC);
//		fclose(fp);
//#endif

#endif

    /////////////////////////

    // construct the connection region for planar example? 03/07/2010

    //mcg->cal_mcgedge_regions();  // comment out by Guoning for a quick test 06/29/2010

    //free (tedges);
    //tedges = nullptr;

    return true;
}


void RegionTauMap::CleanRegion(int id)
{
    DirGraph_Node **sccnodes = morse_decomp->dg->nlist->dirnodes;
    Graph_Edge **sccedges = morse_decomp->dg->elist->edges;
    int num_sccedges=morse_decomp->dg->elist->nedges;

    cur_SCC_id = mcg->nlist->mnodes[id]->scc_index;

    for(int nodeID=0;nodeID</*conley->*/morse_nlist.size();nodeID++)
    {

        int cur_node=/*conley->*/morse_nlist[nodeID];
        for(int edgeID=0;edgeID<sccnodes[cur_node]->nedges;edgeID++)
        {
            num_sccedges = morse_decomp->dg->elist->nedges;  //modified by Guoning 09/24/09

            int cur_edge=sccnodes[cur_node]->edges[edgeID];
            //get two ends of the edge, these two variable should not be changed!
            int node1ID=sccedges[cur_edge]->node_index1;
            int node2ID=sccedges[cur_edge]->node_index2;

            //get the scc_id of the two ends
            int n1sccID=sccnodes[node1ID]->sscomp_index;
            int n2sccID=sccnodes[node2ID]->sscomp_index;

            if(n1sccID==n2sccID)
            {
                //delete edge from sccedges
                sccedges[cur_edge]->cancel=sccedges[num_sccedges-1]->cancel;
                sccedges[cur_edge]->node_index1=sccedges[num_sccedges-1]->node_index1;
                sccedges[cur_edge]->node_index2=sccedges[num_sccedges-1]->node_index2;
                sccedges[cur_edge]->visited=sccedges[num_sccedges-1]->visited;

                //change the node.edge
                int n1c=sccedges[cur_edge]->node_index1;
                int n2c=sccedges[cur_edge]->node_index2;

                for(int edgeN1=0;edgeN1<sccnodes[n1c]->nedges;edgeN1++)
                {
                    if(sccnodes[n1c]->edges[edgeN1]==num_sccedges-1)
                    {
                        sccnodes[n1c]->edges[edgeN1]=cur_edge;
                        break;
                    }
                }

                for(int edgeN2=0;edgeN2<sccnodes[n2c]->nedges;edgeN2++) // corrected by Guoning 09/24/09
                {
                    if(sccnodes[n2c]->edges[edgeN2]==num_sccedges-1)
                    {
                        sccnodes[n2c]->edges[edgeN2]=cur_edge;
                        break;
                    }
                }
                morse_decomp->dg->elist->nedges--;

                //delete cur_edge from its original two end nodes
                for(int j=0;j<sccnodes[node1ID]->nedges;j++)
                {

                    if(sccnodes[node1ID]->edges[j]==cur_edge)
                    {
                        if(j!=sccnodes[node1ID]->nedges-1)
                            sccnodes[node1ID]->edges[j]=sccnodes[node1ID]->edges[sccnodes[node1ID]->nedges-1];
                        sccnodes[node1ID]->nedges--;

                        break;
                    }
                }

                for(int j=0;j<sccnodes[node2ID]->nedges;j++)
                {

                    if(sccnodes[node2ID]->edges[j]==cur_edge)
                    {
                        if(j!=sccnodes[node2ID]->nedges-1)
                            sccnodes[node2ID]->edges[j]=sccnodes[node2ID]->edges[sccnodes[node2ID]->nedges-1];
                        sccnodes[node2ID]->nedges--;

                        break;
                    }
                }
                edgeID--;
            }
        }
    }
}


void RegionTauMap::local_flow_comb(double tau)
{
    init_Region();
    build_local_connected_graph(tau);

    //output temp edgelist
    /*ofstream tempelist("tempelist.txt");
    for(int i=0;i<num_tedges;i++)
    {
        tempelist<<i<<"\t"<<"n1:"<<"\t"<<tedges[i].node_index1<<"\t"<<"n2:"<<"\t"<<tedges[i].node_index2<<endl;
    }
    tempelist.close();*/
    //output temp edgelist


    /*add_temp_to_edgelist();*/

    //output directed graph edgelist
    /*ofstream dgelist("dgelist.txt");
    for(int i=0;i<morse_decomp->dg->elist->nedges;i++)
    {
        dgelist<<i<<"\t"<<"n1:"<<"\t"<<morse_decomp->dg->elist->edges[i]->node_index1<<"\t"<<"n2:"<<"\t"<<morse_decomp->dg->elist->edges[i]->node_index2<<endl;
    }
    dgelist.close();*/
    //output directed graph edgelist
}


void RegionTauMap::init_Region(void)
{
    tedges = (TempEdge*) malloc(sizeof(TempEdge)* object->elist.nedges*2); //it requires the construction of the edge list
    num_tedges = 0;
    curMaxNumTempEdges = object->elist.nedges*2;
}

void RegionTauMap::build_local_connected_graph(double tau)
{

    ///***1. first trace backward*/
    trace_local_Verts_b2(tau);//rewrite
    //trace_local_Verts_b(tau);//rewrite

    ///*1.2. build the edges according to the result*/
    //build_edges_local_Vers(1);
    build_edges_local_Vers_b();

    trace_all_centers_tris_build_edges(tau, 1);
    trace_all_edges_build_di_edges_adp(tau, 1);
    //trace_all_edges_build_di_edges_adp_b(tau);

    ///***2. perform forward tracing*/
    trace_local_Verts_f2(tau);

    ///*2.2. build the edges according to the result*/
    //build_edges_local_Vers(0);
    build_edges_local_Vers_f();

    trace_all_centers_tris_build_edges(tau, 0);
    trace_all_edges_build_di_edges_adp(tau, 0);
    //trace_all_edges_build_di_edges_adp_f(tau);
}

void RegionTauMap::add_temp_to_edgelist(void)
{
    int node1ID;
    int node2ID;
    int n1sccID;
    int n2sccID;

    for(int i=0;i<num_tedges;i++)
    {
        //get two ends of the edge
        node1ID=tedges[i].node_index1;
        node2ID=tedges[i].node_index2;

        /*if(node1ID<0 || node2ID<0)
        {
            less0<<"node1ID<0 || node2ID<0"<<endl;
            continue;
        }*/

        //get the scc_id of the two ends
        n1sccID=morse_decomp->dg->nlist->dirnodes[node1ID]->sscomp_index;
        n2sccID=morse_decomp->dg->nlist->dirnodes[node2ID]->sscomp_index;

        /*if(n1sccID<0 || n1sccID<0)
        {
            less0<<"n1sccID<0 || n1sccID<0"<<endl;
            continue;
        }*/

        //if(n1sccID==n2sccID && cur_SCC_id == n1sccID && node1ID !=node2ID)  // comment out by Guoning 03/08/2010

        if (node1ID != node2ID)
        {
            //??????????????????????
            //if(!has_Edge_From_To(node1ID, node2ID))
            if(!morse_decomp->dg->is_repeated_edge_2(node1ID, node2ID))
            {
                //????????????????????
                morse_decomp->dg->add_to_edgelist(node1ID, node2ID, morse_decomp->dg->elist->nedges);

                /*SCC_AddToEdge(node1ID, node2ID, num_sccedges);*/

                /*add the edge into the node's edge list*/
                //?????????????????????
                morse_decomp->dg->add_edge_to_node(node1ID, morse_decomp->dg->elist->nedges-1);
                morse_decomp->dg->add_edge_to_node(node2ID, morse_decomp->dg->elist->nedges-1);
                /*SCC_AddEdgeToNode(node1ID, num_sccedges-1);
                SCC_AddEdgeToNode(node2ID, num_sccedges-1);*/
            }
        }
    }

    //free(tedges);

}

void RegionTauMap::trace_local_Verts(double tau, int backward)
{
    int i, k;
    int tri_id, end_tris;
    Vertex *v;
    icVector3 stP, newP;

    double lp[2];
    double alpha[3];
    Triangle *face;
    Triangle *t;

    Trajectory *temp= new Trajectory(-1, 1);

    /////////////////////////////////////////////
    //for(i = 0; i < object->vlist.nverts; i++)
    //	object->vlist.verts[i]->visited=false;
    for(int j=0;j</*conley->*/morse_nlist.size();j++)
    {
        t=object->tlist.tris[/*conley->*/morse_nlist[j]];
        for(k=0;k<t->nverts;k++)
            t->verts[k]->visited = false;
    }
    /////////////////////////////////////////////

    //for(i = 0; i < object->vlist.nverts; i++)

    /////////////////////////////////////////////
    for(int j=0;j</*conley->*/morse_nlist.size();j++)
    /////////////////////////////////////////////
    {
        /////////////////////////////////////////
        t=object->tlist.tris[/*conley->*/morse_nlist[j]];
        for(int n=0;n<t->nverts;n++)
        {
            v=t->verts[n];
            if(v->visited==false)
            {
                /////////////////////////////////////////
                v->visited=true;

                //v = object->vlist.verts[i];



                //if(length(v->t_vec) < 1e-10) continue;   //probably no vector value on it
                /* we can choose the triangle that the flow will lead the vertex go into */

                //temp->pass_vertex(v->index, tri_id, backward);

                //if(tri_id < 0)
                //{
                //	tri_id = v->corners[0]->t;
                //}

                if (backward==0)
                    tri_id = v->tauPt_f.which_tri;
                else
                    tri_id = v->tauPt_b.which_tri;

                stP.entry[0] = v->x;
                stP.entry[1] = v->y;
                stP.entry[2] = v->z;


                morse_decomp->trace_Ver(tri_id, stP.entry, newP.entry, end_tris, tau, backward);

                v->img_tau[0] = newP.entry[0];
                v->img_tau[1] = newP.entry[1];
                v->img_tau[2] = newP.entry[2];

                v->imgtri = end_tris;
                ///////////////
                v->visited=true;
            }
        }
        ///////////////
    }
    delete temp;

}


void RegionTauMap::trace_local_Verts_f2(double tau)
{
    int i, k;
    int tri_id, end_tris;
    Vertex *v;
    icVector3 stP, newP;

    double lp[2];
    double alpha[3];
    Triangle *face;
    Triangle *t;

    Trajectory *temp= new Trajectory(-1, 1);

    /////////////////////////////////////////////
    //for(i = 0; i < object->vlist.nverts; i++)
    //	object->vlist.verts[i]->visited=false;
    for(int j=0;j</*conley->*/morse_nlist.size();j++)
    {
        t=object->tlist.tris[/*conley->*/morse_nlist[j]];
        for(k=0;k<t->nverts;k++)
            t->verts[k]->visited = false;
    }
    /////////////////////////////////////////////

    //for(i = 0; i < object->vlist.nverts; i++)

    /////////////////////////////////////////////
    for(int j=0;j</*conley->*/morse_nlist.size();j++)
    /////////////////////////////////////////////
    {
        /////////////////////////////////////////
        t=object->tlist.tris[/*conley->*/morse_nlist[j]];
        for(int n=0;n<t->nverts;n++)
        {
            v=t->verts[n];
            if(v->visited==false)
            {
                /////////////////////////////////////////
                v->visited=true;


                tri_id = v->tauPt_f.which_tri;

                stP.entry[0] = v->x;
                stP.entry[1] = v->y;
                stP.entry[2] = v->z;


                morse_decomp->trace_Ver_f(tri_id, stP.entry, newP.entry, end_tris, tau);

                v->img_tau[0] = newP.entry[0];
                v->img_tau[1] = newP.entry[1];
                v->img_tau[2] = newP.entry[2];

                v->imgtri = end_tris;
                ///////////////
                v->visited=true;
            }
        }
        ///////////////
    }
    delete temp;
}

void RegionTauMap::trace_local_Verts_b2(double tau)
{
    int i, k;
    int tri_id, end_tris;
    Vertex *v;
    icVector3 stP, newP;

    double lp[2];
    double alpha[3];
    Triangle *face;
    Triangle *t;

    Trajectory *temp= new Trajectory(-1, 1);

    /////////////////////////////////////////////
    //for(i = 0; i < object->vlist.nverts; i++)
    //	object->vlist.verts[i]->visited=false;
    for(int j=0;j</*conley->*/morse_nlist.size();j++)
    {
        t=object->tlist.tris[/*conley->*/morse_nlist[j]];
        for(k=0;k<t->nverts;k++)
            t->verts[k]->visited = false;
    }
    /////////////////////////////////////////////

    //for(i = 0; i < object->vlist.nverts; i++)

    /////////////////////////////////////////////
    for(int j=0;j</*conley->*/morse_nlist.size();j++)
    /////////////////////////////////////////////
    {
        /////////////////////////////////////////
        t=object->tlist.tris[/*conley->*/morse_nlist[j]];
        for(int n=0;n<t->nverts;n++)
        {
            v=t->verts[n];
            if(v->visited==false)
            {
                /////////////////////////////////////////
                v->visited=true;


                tri_id = v->tauPt_b.which_tri;

                stP.entry[0] = v->x;
                stP.entry[1] = v->y;
                stP.entry[2] = v->z;


                morse_decomp->trace_Ver_b(tri_id, stP.entry, newP.entry, end_tris, tau);

                v->img_tau[0] = newP.entry[0];
                v->img_tau[1] = newP.entry[1];
                v->img_tau[2] = newP.entry[2];

                v->imgtri = end_tris;
                ///////////////
                v->visited=true;
            }
        }
        ///////////////
    }
    delete temp;

}

void RegionTauMap::build_edges_local_Vers(int backward)
{
    int i;
    Vertex *v;

    Triangle *t;

    //for(i = 0; i < object->vlist.nverts; i++)
    //	object->vlist.verts[i]->visited=false;

    for(int j=0;j</*conley->*/morse_nlist.size();j++)
    {
        t=object->tlist.tris[/*conley->*/morse_nlist[j]];
        for(int k=0;k<t->nverts;k++)
            t->verts[k]->visited = false;
    }

    for(int j=0;j</*conley->*/morse_nlist.size();j++)
    {
        t=object->tlist.tris[/*conley->*/morse_nlist[j]];

        for(int k=0;k<t->nverts;k++)
        {
            v=t->verts[k];
            if(v->visited==false)
            {
                cur_end_tedge = 0;
                build_edges_Ver(t->verts[k]->index, backward);//rewrite
                v->visited=true;
            }
        }
    }
}


void
RegionTauMap::build_edges_local_Vers_f()
{
    int i;
    Vertex *v;

    Triangle *t;

    //for(i = 0; i < object->vlist.nverts; i++)
    //	object->vlist.verts[i]->visited=false;

    for(int j=0;j</*conley->*/morse_nlist.size();j++)
    {
        t=object->tlist.tris[/*conley->*/morse_nlist[j]];
        for(int k=0;k<t->nverts;k++)
            t->verts[k]->visited = false;
    }

    for(int j=0;j</*conley->*/morse_nlist.size();j++)
    {
        t=object->tlist.tris[/*conley->*/morse_nlist[j]];

        for(int k=0;k<t->nverts;k++)
        {
            v=t->verts[k];
            if(v->visited==false)
            {
                cur_end_tedge = 0;
                build_edges_Ver_f(t->verts[k]->index);//rewrite
                v->visited=true;
            }
        }
    }
}


void
RegionTauMap::build_edges_local_Vers_b()
{
    int i;
    Vertex *v;

    Triangle *t;

    //for(i = 0; i < object->vlist.nverts; i++)
    //	object->vlist.verts[i]->visited=false;

    for(int j=0;j</*conley->*/morse_nlist.size();j++)
    {
        t=object->tlist.tris[/*conley->*/morse_nlist[j]];
        for(int k=0;k<t->nverts;k++)
            t->verts[k]->visited = false;
    }

    for(int j=0;j</*conley->*/morse_nlist.size();j++)
    {
        t=object->tlist.tris[/*conley->*/morse_nlist[j]];

        for(int k=0;k<t->nverts;k++)
        {
            v=t->verts[k];
            if(v->visited==false)
            {
                cur_end_tedge = 0;
                build_edges_Ver_b(t->verts[k]->index);//rewrite
                v->visited=true;
            }
        }
    }
}


void RegionTauMap::trace_all_centers_tris_build_edges(double tau, int backward)
{


    for(int j=0;j</*conley->*/morse_nlist.size();j++)
    {
        trace_center_tris_build_edge(/*conley->*/morse_nlist[j], tau, backward);

    }
}

void RegionTauMap::trace_all_edges_build_di_edges_adp(double tau, int backward)
{
    int i, j;
    Triangle *t;
    Edge *e;
    Vertex *v1, *v2;
    double st1[3], st2[3];
    int endtri;

    init_local_Edges();

    for(i=0;i</*conley->*/morse_nlist.size();i++)
    {
        t=object->tlist.tris[/*conley->*/morse_nlist[i]];

        for(j = 0; j < 3; j++)
        {
            e = t->edges[j];
            if(e->visited )/*we process each edge only once*/
                continue;

            e->visited = true;

            v1 = e->verts[0];
            v2 = e->verts[1];

            /*currently, we don't deal with the boundary edges!!!!*/
            if(v1->imgtri < 0 || v2->imgtri < 0)
                continue;

            /*if the images of the two vertices are already continuous, need not
            consider this edge any more*/
            if(v1->imgtri == v2->imgtri
                || morse_decomp->are_close_neighbors(v1->imgtri, v2->imgtri))
                continue;

            st1[0] = v1->x;
            st1[1] = v1->y;
            st1[2] = v1->z;
            st2[0] = v2->x;
            st2[1] = v2->y;
            st2[2] = v2->z;

            Triangle *nei_tri = e->tris[0];

            if (nei_tri == t)
                nei_tri = e->tris[1];

            if (nei_tri == nullptr)
                continue;

            int neighbor_tri = nei_tri->index;


            //if(e->tris[0]!=nullptr)
            //{
            //	neighbor_tri = e->tris[0]->index;
            //	if(neighbor_tri == i && e->tris[1] == nullptr)
            //		continue;
            //	else if(neighbor_tri == i && e->tris[1] !=nullptr)
            //		neighbor_tri = e->tris[1]->index;
            //}
            //else if(e->tris[1] !=nullptr && e->tris[1]->index != i)
            //	neighbor_tri=e->tris[1]->index;
            //else
            //	continue;

            //trace_an_edge_build_di_edges_adp(st1, st2, v1->imgtri, v2->imgtri, /*conley->*/morse_nlist[i], neighbor_tri,tau, backward);

            if (enFastRefinement)
                trace_an_edge_build_di_edges_adp(st1, st2, v1->imgtri, v2->imgtri, /*conley->*/morse_nlist[i], neighbor_tri,tau, backward);
            else
                adp_edge_sampling(st1, st2, v1->imgtri, v2->imgtri, morse_nlist[i], neighbor_tri,
                                  tau, backward);
            //adp_edge_sampling_2(st1, st2, v1->imgtri, v2->imgtri, morse_nlist[i], neighbor_tri,
            //	tau, backward);

            /*
                modified by Guoning 07/22/2010
            */
            //if (enFastRefinement)
            //	morse_decomp->trace_an_edge_build_di_edges_adp(st1, st2, v1->imgtri, v2->imgtri, i, neighbor_tri,
            //		tau, backward);

            //else
            //	morse_decomp->adp_edge_sampling(st1, st2, v1->imgtri, v2->imgtri, i, neighbor_tri,
            //		tau, backward);

        }
    }
}



void
RegionTauMap::trace_all_edges_build_di_edges_adp_f(double tau)
{
    int i, j;
    Triangle *t;
    Edge *e;
    Vertex *v1, *v2;
    double st1[3], st2[3];
    int endtri;

    init_local_Edges();

    for(i=0;i</*conley->*/morse_nlist.size();i++)
    {
        t=object->tlist.tris[/*conley->*/morse_nlist[i]];

        for(j = 0; j < 3; j++)
        {
            e = t->edges[j];
            if(e->visited )/*we process each edge only once*/
                continue;

            e->visited = true;

            v1 = e->verts[0];
            v2 = e->verts[1];

            /*currently, we don't deal with the boundary edges!!!!*/
            if(v1->imgtri < 0 || v2->imgtri < 0)
                continue;

            /*if the images of the two vertices are already continuous, need not
            consider this edge any more*/
            if(v1->imgtri == v2->imgtri
                || morse_decomp->are_close_neighbors(v1->imgtri, v2->imgtri))
                continue;

            st1[0] = v1->x;
            st1[1] = v1->y;
            st1[2] = v1->z;
            st2[0] = v2->x;
            st2[1] = v2->y;
            st2[2] = v2->z;

            Triangle *nei_tri = e->tris[0];

            if (nei_tri == t)
                nei_tri = e->tris[1];

            if (nei_tri == nullptr)
                continue;

            int neighbor_tri = nei_tri->index;


            if (enFastRefinement)
                trace_an_edge_build_di_edges_adp_f(st1, st2, v1->imgtri, v2->imgtri, /*conley->*/morse_nlist[i], neighbor_tri,tau);
            else
                adp_edge_sampling(st1, st2, v1->imgtri, v2->imgtri, morse_nlist[i], neighbor_tri,
                                  tau, 0);
            //adp_edge_sampling_2(st1, st2, v1->imgtri, v2->imgtri, morse_nlist[i], neighbor_tri,
            //	tau, 0);

        }
    }
}


void
RegionTauMap::trace_all_edges_build_di_edges_adp_b(double tau)
{
    int i, j;
    Triangle *t;
    Edge *e;
    Vertex *v1, *v2;
    double st1[3], st2[3];
    int endtri;

    init_local_Edges();

    for(i=0;i</*conley->*/morse_nlist.size();i++)
    {
        t=object->tlist.tris[/*conley->*/morse_nlist[i]];

        for(j = 0; j < 3; j++)
        {
            e = t->edges[j];
            if(e->visited )/*we process each edge only once*/
                continue;

            e->visited = true;

            v1 = e->verts[0];
            v2 = e->verts[1];

            /*currently, we don't deal with the boundary edges!!!!*/
            if(v1->imgtri < 0 || v2->imgtri < 0)
                continue;

            /*if the images of the two vertices are already continuous, need not
            consider this edge any more*/
            if(v1->imgtri == v2->imgtri
                || morse_decomp->are_close_neighbors(v1->imgtri, v2->imgtri))
                continue;

            st1[0] = v1->x;
            st1[1] = v1->y;
            st1[2] = v1->z;
            st2[0] = v2->x;
            st2[1] = v2->y;
            st2[2] = v2->z;

            Triangle *nei_tri = e->tris[0];

            if (nei_tri == t)
                nei_tri = e->tris[1];

            if (nei_tri == nullptr)
                continue;

            int neighbor_tri = nei_tri->index;



            if (enFastRefinement)
                trace_an_edge_build_di_edges_adp_b(st1, st2, v1->imgtri, v2->imgtri, morse_nlist[i], neighbor_tri,tau);
            else
                adp_edge_sampling(st1, st2, v1->imgtri, v2->imgtri, morse_nlist[i], neighbor_tri,
                                  tau, 1);
            //adp_edge_sampling_2(st1, st2, v1->imgtri, v2->imgtri, morse_nlist[i], neighbor_tri,
            //	tau, 1);

        }
    }
}


void RegionTauMap::build_edges_Ver(int vertid, int backward)
{
    Vertex *v = object->vlist.verts[vertid];

    int i;

    for(i = 0; i < v->ncorners; i++)
    {
        if(v->corners[i]->t < 0 || v->imgtri < 0)
            continue;

        if(backward == 0) /*forward tracing*/
        {
            //if(morse_decomp->dg->is_repeated_edge_2(v->corners[i]->t, v->imgtri))
            //	continue;

            build_one_edge(v->corners[i]->t, v->imgtri);
        }

        else /*backward tracing*/
        {
            //if(morse_decomp->dg->is_repeated_edge_2(v->imgtri, v->corners[i]->t))
            //	continue;

            build_one_edge(v->imgtri, v->corners[i]->t);
        }
    }
}


void
RegionTauMap::build_edges_Ver_f(int vertid)
{
    Vertex *v = object->vlist.verts[vertid];

    int i;

    for(i = 0; i < v->ncorners; i++)
    {
        if(v->corners[i]->t < 0 || v->imgtri < 0)
            continue;

        build_one_edge(v->corners[i]->t, v->imgtri);
    }
}


void
RegionTauMap::build_edges_Ver_b(int vertid)
{
    Vertex *v = object->vlist.verts[vertid];

    int i;

    for(i = 0; i < v->ncorners; i++)
    {
        if(v->corners[i]->t < 0 || v->imgtri < 0)
            continue;


        build_one_edge(v->imgtri, v->corners[i]->t);
    }
}


void RegionTauMap::build_one_edge(int node_from, int node_to)
{
    Region_AddToEdge(node_from, node_to, num_tedges);
}

void RegionTauMap::Region_AddToEdge(int node1, int node2, int& cur_index)
{
    // Try to modify these using a dynamic structure
    if(cur_index >= curMaxNumTempEdges)
    {
        /*
           avoid using realloc !!
        */

        TempEdge *temp = tedges;
        //tedges = (TempEdge *)realloc(tedges, sizeof(TempEdge)*(curMaxNumTempEdges+object->elist.nedges/2));
        //tedges = (TempEdge *)malloc(sizeof(TempEdge)*(curMaxNumTempEdges*2/*+object->elist.nedges/2*/));
        //tedges = (TempEdge *)malloc(sizeof(TempEdge)*(curMaxNumTempEdges+object->elist.nedges/2));
        tedges = (TempEdge *)malloc(sizeof(TempEdge)*(curMaxNumTempEdges+curMaxNumTempEdges/2));
        if(tedges == nullptr)
        {
            //MessageBox(nullptr, "failed to reallocate memory!", "Error", MB_OK);
            tedges = temp;
            exit(-1);
        }

        for (int i=0; i<curMaxNumTempEdges; i++)
        {
            tedges[i].edge_index = temp[i].edge_index;
            tedges[i].node_index1 = temp[i].node_index1;
            tedges[i].node_index2 = temp[i].node_index2;
            tedges[i].cancelled = temp[i].cancelled;
            tedges[i].visited = temp[i].visited;
        }
        //curMaxNumTempEdges += Object.nedges/2;
        //curMaxNumTempEdges += object->elist.nedges/2;
        //curMaxNumTempEdges += curMaxNumTempEdges;
        //
        curMaxNumTempEdges += (curMaxNumTempEdges/2);

        free (temp);
    }

    tedges[cur_index].node_index1 = node1;
    tedges[cur_index].node_index2 = node2;
    //tedges[cur_index].node_index1 = node2;
    //tedges[cur_index].node_index2 = node1;
    tedges[cur_index].edge_index = cur_index;
    tedges[cur_index].cancelled = false;
    cur_index++;

    cur_end_tedge = cur_index;
}

void RegionTauMap::trace_center_tris_build_edge(int tri, double tau, int backward)
{
    double center[3] = {0., 0., 0.};

    Triangle *t=object->tlist.tris[tri];

    int i, endtri;

    //for(i = 0; i < 3; i++)
    //{

    //	center[0] += object->vlist.verts[t->verts[i]->index]->x;
    //	center[1] += object->vlist.verts[t->verts[i]->index]->y;
    //	center[2] += object->vlist.verts[t->verts[i]->index]->z;
    //}

    //center[0] /= 3.;
    //center[1] /= 3.;
    //center[2] /= 3.;

    center[0] = t->tauPt_f.x;
    center[1] = t->tauPt_f.y;
    center[2] = t->tauPt_f.z;

    /*start tracing here*/
    morse_decomp->trace_Ver(tri, center, center, endtri, tau, backward);

    if(endtri < 0)
        return;

    if(backward == 0)  /*forward tracing*/
    {
        build_one_edge(tri, endtri);
    }

    else /*backward tracing*/
    {
        build_one_edge(endtri, tri);
    }
}

void RegionTauMap::init_local_Edges(void)
{
    int i, j;
    Triangle *t;
    Edge *cur_edge;
    int num_region_node;
    num_region_node=/*conley->*/morse_nlist.size();

    for(i=0;i<num_region_node;i++)
    {
        t=object->tlist.tris[/*conley->*/morse_nlist[i]];
        for(j = 0; j < 3; j++)
        {
            cur_edge = t->edges[j];
            cur_edge->OnBoundary = false;
            cur_edge->visited = false;
        }
    }
}

void RegionTauMap::trace_an_edge_build_di_edges_adp(double st1[3], double st2[3], int t1, int t2, int tri, int neighbor_tri, double tau, int backward)
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



void
RegionTauMap::trace_an_edge_build_di_edges_adp_f(double st1[3],
                                                 double st2[3],
                                                 int t1,
                                                 int t2,
                                                 int tri,
                                                 int neighbor_tri,
                                                 double tau)
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
    trace_recursive_an_edge_f(st1, st2, t1, tri, neighbor_tri, tau, level);

    /*call the recursive adaptive edge sampling for the second half of the edge*/
    trace_recursive_an_edge_f(stack_st, middle_p, t2, tri, neighbor_tri, tau, level);
}

void
RegionTauMap::trace_an_edge_build_di_edges_adp_b(double st1[3],
                                                 double st2[3],
                                                 int t1,
                                                 int t2,
                                                 int tri,
                                                 int neighbor_tri,
                                                 double tau)
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
    trace_recursive_an_edge_b(st1, st2, t1, tri, neighbor_tri, tau, level);

    /*call the recursive adaptive edge sampling for the second half of the edge*/
    trace_recursive_an_edge_b(stack_st, middle_p, t2, tri, neighbor_tri, tau, level);
}



void RegionTauMap::trace_recursive_an_edge(double v1[3],
                                           double v2[3],
                                           int& t1,
                                           int tri,
                                           int neighbor_tri,
                                           double tau,
                                           int backward,
                                           int& level)
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

    if (backward==0)
        morse_decomp->trace_Ver_f(tri, v2, v2_end, t2, tau);
    else
        morse_decomp->trace_Ver_b(tri, v2, v2_end, t2, tau);

    //morse_decomp->trace_Ver(tri, v2, v2_end, t2, tau, backward);

    /*boundary!!*/
    if(t2 < 0)
        return;

    icVector3 dis;
    dis.entry[0] = v1[0]-v2[0];
    dis.entry[1] = v1[1]-v2[1];
    dis.entry[2] = v1[2]-v2[2];

    if(length(dis) < edge_sample_error/*1e-8*/)
    {
        v1[0] = v2[0];
        v1[1] = v2[1];
        v1[2] = v2[2];
        t1 = t2;
        /*need to build edges here*/
        if(backward == 0) /*forward tracing*/
        {
            //if(!morse_decomp->dg->is_repeated_edge_2(tri, t1))
            build_one_edge(tri, t1);
            if(neighbor_tri>=0/*&&!morse_decomp->dg->is_repeated_edge_2(neighbor_tri, t1)*/)
                build_one_edge(neighbor_tri, t1);
        }
        else /*backward tracing*/
        {

            //if(!morse_decomp->dg->is_repeated_edge_2(t1,tri))
            build_one_edge(t1, tri);
            if(neighbor_tri>=0/*&&!morse_decomp->dg->is_repeated_edge_2(t1,neighbor_tri)*/)
                build_one_edge(t1, neighbor_tri);
        }
        return;
    }


    /*if the two points are too close, stop!*/
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

            //if(!morse_decomp->dg->is_repeated_edge_2(tri, t1))
            build_one_edge(tri, t1);
            if(neighbor_tri>=0/*&&!morse_decomp->dg->is_repeated_edge_2(neighbor_tri, t1)*/)
                build_one_edge(neighbor_tri, t1);
        }
        else /*backward tracing*/
        {

            //if(!morse_decomp->dg->is_repeated_edge_2(t1, tri))
            build_one_edge(t1, tri);
            if(neighbor_tri>=0 /*&&!morse_decomp->dg->is_repeated_edge_2(t1,neighbor_tri)*/)
                build_one_edge(t1, neighbor_tri);
        }
        return;
    }

    if(t1 == t2 || morse_decomp->are_close_neighbors(t1, t2))
    {
        v1[0] = v2[0];
        v1[1] = v2[1];
        v1[2] = v2[2];
        t1 = t2;

        /*need to build edges here*/
        if(backward == 0) /*forward tracing*/
        {
            //if(!morse_decomp->dg->is_repeated_edge_2(tri, t1))
            build_one_edge(tri, t1);
            if(neighbor_tri>=0/*&&!morse_decomp->dg->is_repeated_edge_2(neighbor_tri, t1)*/)
                build_one_edge(neighbor_tri, t1);
        }
        else /*backward tracing*/
        {
            //if(!morse_decomp->dg->is_repeated_edge_2(t1, tri))
            build_one_edge(t1, tri);
            if(neighbor_tri>=0/*&&!morse_decomp->dg->is_repeated_edge_2(t1, neighbor_tri)*/)
                build_one_edge(t1, neighbor_tri);
        }
        return;
    }



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


void RegionTauMap::trace_recursive_an_edge_f(double v1[3],
                                             double v2[3],
                                             int& t1,
                                             int tri,
                                             int neighbor_tri,
                                             double tau,
                                             int& level)
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

    morse_decomp->trace_Ver_f(tri, v2, v2_end, t2, tau);

    //morse_decomp->trace_Ver(tri, v2, v2_end, t2, tau, backward);

    /*boundary!!*/
    if(t2 < 0)
        return;

    icVector3 dis;
    dis.entry[0] = v1[0]-v2[0];
    dis.entry[1] = v1[1]-v2[1];
    dis.entry[2] = v1[2]-v2[2];

    if(length(dis) < edge_sample_error/*1e-8*/)
    {
        v1[0] = v2[0];
        v1[1] = v2[1];
        v1[2] = v2[2];
        t1 = t2;
        /*need to build edges here*/

        build_one_edge(tri, t1);
        if(neighbor_tri>=0)
            build_one_edge(neighbor_tri, t1);
        return;
    }


    /*if the two points are too close, stop!*/
    if(level > 10)
    {
        /*trace v2, get t2 */
        v1[0] = v2[0];
        v1[1] = v2[1];
        v1[2] = v2[2];
        t1 = t2;
        /*need to build edges here*/

        build_one_edge(tri, t1);
        if(neighbor_tri>=0)
            build_one_edge(neighbor_tri, t1);
        return;
    }

    if(t1 == t2 || morse_decomp->are_close_neighbors(t1, t2))
    {
        v1[0] = v2[0];
        v1[1] = v2[1];
        v1[2] = v2[2];
        t1 = t2;

        /*need to build edges here*/
        //if(!morse_decomp->dg->is_repeated_edge_2(tri, t1))
        build_one_edge(tri, t1);
        if(neighbor_tri>=0/*&&!morse_decomp->dg->is_repeated_edge_2(neighbor_tri, t1)*/)
            build_one_edge(neighbor_tri, t1);
        return;
    }



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
        trace_recursive_an_edge_f(v1, v2, t1, tri, neighbor_tri, tau, level);
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
        trace_recursive_an_edge_f(v1, v2, t1, tri, neighbor_tri, tau, level);
        level--;
    }
}







void RegionTauMap::trace_recursive_an_edge_b(double v1[3],
                                             double v2[3],
                                             int& t1,
                                             int tri,
                                             int neighbor_tri,
                                             double tau,
                                             int& level)
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

    morse_decomp->trace_Ver_b(tri, v2, v2_end, t2, tau);

    //morse_decomp->trace_Ver(tri, v2, v2_end, t2, tau, backward);

    /*boundary!!*/
    if(t2 < 0)
        return;

    icVector3 dis;
    dis.entry[0] = v1[0]-v2[0];
    dis.entry[1] = v1[1]-v2[1];
    dis.entry[2] = v1[2]-v2[2];

    if(length(dis) < edge_sample_error/*1e-8*/)
    {
        v1[0] = v2[0];
        v1[1] = v2[1];
        v1[2] = v2[2];
        t1 = t2;
        /*need to build edges here*/

        build_one_edge(t1, tri);
        if(neighbor_tri>=0)
            build_one_edge(t1, neighbor_tri);
        return;
    }


    /*if the two points are too close, stop!*/
    if(level > 10)
    {
        /*trace v2, get t2 */
        v1[0] = v2[0];
        v1[1] = v2[1];
        v1[2] = v2[2];
        t1 = t2;
        /*need to build edges here*/

        build_one_edge(t1, tri);
        if(neighbor_tri>=0)
            build_one_edge(t1, neighbor_tri);
        return;
    }

    if(t1 == t2 || morse_decomp->are_close_neighbors(t1, t2))
    {
        v1[0] = v2[0];
        v1[1] = v2[1];
        v1[2] = v2[2];
        t1 = t2;

        /*need to build edges here*/
        //if(!morse_decomp->dg->is_repeated_edge_2(t1, tri))
        build_one_edge(t1, tri);
        if(neighbor_tri>=0/*&&!morse_decomp->dg->is_repeated_edge_2(t1, neighbor_tri)*/)
            build_one_edge(t1, neighbor_tri);
        return;
    }



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
        trace_recursive_an_edge_b(v1, v2, t1, tri, neighbor_tri, tau, level);
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
        trace_recursive_an_edge_b(v1, v2, t1, tri, neighbor_tri, tau, level);
        level--;
    }
}



/*
 non-recursive adaptive edge sampling
*/
void
RegionTauMap::adp_edge_sampling(double st1[3],
                                double st2[3],
                                int t1, int t2,
                                int tri,
                                int neighbor_tri,
                                double tau,
                                int backward)
{
    /*do one more recursive here*/
    int level = 1;
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
    //morse_decomp->trace_Ver(tri, st1, v1_end, t1, tau, backward);

    if (backward == 0)
        morse_decomp->trace_Ver_f(tri, st1, v1_end, t1, tau);
    else
        morse_decomp->trace_Ver_b(tri, st1, v1_end, t1, tau);

    while(!stack_pts.empty())
    {
        point3 top_p = stack_pts.top();
        stack_pts.pop();

        st2[0] = top_p.p[0];
        st2[1] = top_p.p[1];
        st2[2] = top_p.p[2];

        /*   conduct a tracing from this point  */
        double v2_end[3];
        //morse_decomp->trace_Ver(tri, st2, v2_end, t2, tau, backward);

        if (backward == 0)
            morse_decomp->trace_Ver_f(tri, st2, v2_end, t2, tau);
        else
            morse_decomp->trace_Ver_b(tri, st2, v2_end, t2, tau);

        //if (t2 < 0 || t1 < 0)
        //	continue;

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

            if (t2 >= 0)
            {
                /*need to build edges here*/
                if(backward == 0) /*forward tracing*/
                {
                    //if(!morse_decomp->dg->is_repeated_edge_2(tri, t1))
                    build_one_edge(tri, t1);
                    //if(neighbor_tri>=0&&!dg->is_repeated_edge(neighbor_tri, t1))
                    if(neighbor_tri>=0/*&&!morse_decomp->dg->is_repeated_edge_2(neighbor_tri, t1)*/)
                        build_one_edge(neighbor_tri, t1);
                }
                else /*backward tracing*/
                {
                    //if(!morse_decomp->dg->is_repeated_edge_2(t1, tri))
                    build_one_edge(t1, tri);
                    if(neighbor_tri>=0/*&&!morse_decomp->dg->is_repeated_edge_2(t1, neighbor_tri)*/)
                        build_one_edge(t1, neighbor_tri);
                }
            }

            continue;
        }

        /*   check connectivity with the image of the previous point  */

        /*   if they intersect, */
        if (t1 == t2 /*&& t1 > 0 && t2 > 0*/ || morse_decomp->are_close_neighbors(t1, t2))
        {
            /*need to build edges here*/
            if (t1 > 0 && t2 > 0)
            {
                if(backward == 0) /*forward tracing*/
                {
                    //if(!morse_decomp->dg->is_repeated_edge_2(tri, t1))
                    build_one_edge(tri, t1);
                    if(neighbor_tri>=0/*&&!morse_decomp->dg->is_repeated_edge_2(neighbor_tri, t1)*/)
                        build_one_edge(neighbor_tri, t1);
                }
                else /*backward tracing*/
                {
                    //if(!morse_decomp->dg->is_repeated_edge_2(t1, tri))
                    build_one_edge(t1, tri);
                    if(neighbor_tri>=0/*&&!morse_decomp->dg->is_repeated_edge_2(t1, neighbor_tri)*/)
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



void
RegionTauMap::adp_edge_sampling_2(double st1[3],
                                  double st2[3],
                                  int t1, int t2,
                                  int tri,
                                  int neighbor_tri,
                                  double tau,
                                  int backward)
{
    /*do one more recursive here*/
    int level = 1;
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
    //morse_decomp->trace_Ver(tri, st1, v1_end, t1, tau, backward);
    morse_decomp->trace_Ver_local(tri, st1, v1_end, t1, tau, backward, cur_SCC_id);

    while(!stack_pts.empty())
    {
        point3 top_p = stack_pts.top();
        stack_pts.pop();

        st2[0] = top_p.p[0];
        st2[1] = top_p.p[1];
        st2[2] = top_p.p[2];

        /*   conduct a tracing from this point  */
        double v2_end[3];
        //morse_decomp->trace_Ver(tri, st2, v2_end, t2, tau, backward);
        morse_decomp->trace_Ver_local(tri, st2, v2_end, t2, tau, backward, cur_SCC_id);

        //if (t2 < 0 || t1 < 0)
        //	continue;

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

            //if (t2 >= 0)
            //{
            /*need to build edges here*/
            if(backward == 0) /*forward tracing*/
            {
                //if(!dg->is_repeated_edge(tri, t1))
                //if(!morse_decomp->dg->is_repeated_edge_2(tri, t1))
                build_one_edge(tri, t1);
                if(neighbor_tri>=0/*&&!morse_decomp->dg->is_repeated_edge_2(neighbor_tri, t1)*/)
                    build_one_edge(neighbor_tri, t1);
            }
            else /*backward tracing*/
            {
                //if(!morse_decomp->dg->is_repeated_edge_2(t1, tri))
                build_one_edge(t1, tri);
                if(neighbor_tri>=0/*&&!morse_decomp->dg->is_repeated_edge_2(t1, neighbor_tri)*/)
                    build_one_edge(t1, neighbor_tri);
            }
            //}
            continue;
        }

        /*   check connectivity with the image of the previous point  */

        /*   if they intersect, */
        if (t1 == t2 /*&& t1 > 0 && t2 > 0*/ || morse_decomp->are_close_neighbors(t1, t2))
        {
            /*need to build edges here*/
            //if (t1 > 0 && t2 > 0)
            //{
            if(backward == 0) /*forward tracing*/
            {
                //if(!morse_decomp->dg->is_repeated_edge_2(tri, t1))
                build_one_edge(tri, t1);
                if(neighbor_tri>=0/*&&!morse_decomp->dg->is_repeated_edge_2(neighbor_tri, t1)*/)
                    build_one_edge(neighbor_tri, t1);
            }
            else /*backward tracing*/
            {
                //if(!morse_decomp->dg->is_repeated_edge_2(t1, tri))
                build_one_edge(t1, tri);
                if(neighbor_tri>=0/*&&!morse_decomp->dg->is_repeated_edge_2(t1, neighbor_tri)*/)
                    build_one_edge(t1, neighbor_tri);
            }
            //}

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


void RegionTauMap::Auto_Refine(double tau_max,double ini_tau,double min_pri, int max_iter)
{
    //initializetau
    vector<int> old_nodelist;//to see if the morse set has already been refined
    int id;
    double tau;

    //initialize the mesh
    int i;
    for(i=0;i<object->tlist.ntris;i++)
        object->tlist.tris[i]->exclude=false;//first no face is discard

    for (i=0; i<object->tlist.ntris; i++)
    {
        object->tlist.tris[i]->in_img = false;
        object->tlist.tris[i]->visited = false;
        object->tlist.tris[i]->used_tau = 0.;
    }


    FILE *morse_id = nullptr;

#ifdef OUTPUT_PERFORMANCE
    morse_id = fopen ("refining_morseSet.txt", "w");
    fprintf(morse_id,"So it starts ...\n");
    fclose(morse_id);
#endif

    init_Region();

    //auto refine pipeline
    clock_t start = clock();

    for(int iter=0;iter<max_iter;iter++)
    {
        //get current morse_id & local_tau
        id=mcg->select_morse(min_pri);

#ifdef OUTPUT_PERFORMANCE
        morse_id = fopen ("refining_morseSet.txt", "a");
        fprintf(morse_id,"\n ------------------------------------------------------\n");
        fprintf (morse_id, "%d: The current extraced Morse set to refine is %d\n", iter+1, id);
        fprintf (morse_id, "The number of edges in current MCG is %d and nodes is %d.\n", mcg->elist->nedges, mcg->nlist->nmnodes);
        fclose(morse_id);
#endif \
    //if no morse set in the quene, over
        if(id==-1) break/*continue*/;

        //tau=conley->get_tau(0);//should the morse set be 0 forever????
        //tau=min(ini_tau, tau/100.);

        tau=ini_tau;

        int scc_id = mcg->nlist->mnodes[id]->scc_index;

        //each time refine one morse set
        bool done=false;


        //get the morse set nodelist
        old_nodelist.clear();
        for(int nodeID=0;nodeID<morse_decomp->scclist->scccomponents[scc_id]->nnodes;nodeID++)
        {
            int cur_node=morse_decomp->scclist->scccomponents[scc_id]->nodes[nodeID];
            old_nodelist.push_back(cur_node);
            //object->tlist.tris[old_nodelist[nodeID]]->exclude=true;
        }

        /*
            obtain the previously used tau from the current Morse set
        */
        double pre_tau = object->tlist.tris[old_nodelist[0]]->used_tau;

        if (/*mcg->nlist->mnodes[id]->used_tau*/pre_tau >= tau_max)
        {
            set_exclude_MorseSet (id);  // we exclude this Morse set from the refinement list
            /*  In the improved implementation, this Morse set can be simply dequeue from the priority queue  */
            continue;
        }

        else
            tau = max (ini_tau, 2*pre_tau/*mcg->nlist->mnodes[id]->used_tau*/);

        //for visualization
        //cur_morse=id;

        //refine the morse set
        while(tau<tau_max)
        {
            used_tau = tau;

            ////refine morset first
            //RefineRegion(id,tau);

            //if (id == 85) // for debugging purpose
            //{
            //	int test = 0;
            //}

            // Modified by Guoning on 06/29/2010
            // Here we check whether the computation creates some disconnected SCCs or not
            if (RefineRegion_forAuto(id, tau))
            {

                //RefineRegion(id, tau);

                //if (mcg->nlist->mnodes[id]->cancelled) // removed by Guoning on 07/07/2010
                //	break;

                //check if the morse set refined
                //get the first morse id in old_nodelist
                //int morse1=morse_decomp->scclist->scccomponents[morse_decomp->dg->nlist->dirnodes[old_nodelist[0]]->sscomp_index]->node_index;

                // we need to test whether the decomposed Morse sets are still connected or not  06/29/2010

                //compare the rest of the old_nodelist
                //int pre_id = -1;


                /*
                    Since we already compute the local MCG, we can simply use it to judge whether the Morse set
                    has been refined or dis-regarded.
                    Modified by Guoning on 07/06/2010
                */
                //for(int nodeID=0;nodeID<old_nodelist.size();nodeID++)
                //{
                //	int cur_id=morse_decomp->scclist->scccomponents[morse_decomp->dg->nlist->dirnodes[old_nodelist[nodeID]]->sscomp_index]->node_index;

                //	if (cur_id < 0)
                //		continue;

                //	else if (pre_id < 0)
                //		pre_id = cur_id;

                //	//refined already, need to jump out
                //	else if(cur_id!=pre_id)
                //	{
                //		done=true;
                //		break;
                //	}
                //}

                //if (pre_id < 0) // this is not a Morse set any more, added by Guoning 07/05/2010
                //{
                //	tau = tau_max+1.;
                //	break;
                //}

                if (lmcg->nlist->nmnodes>1)
                {
                    done=true;
                    /*
                       We need to record the \tau value for each triangle inside this region
                    */
                    for(int nodeID=0;nodeID<old_nodelist.size();nodeID++)
                    {
                        object->tlist.tris[old_nodelist[nodeID]]->used_tau = tau;
                    }
                    break;
                }

                else if (lmcg->nlist->nmnodes==0) // This is not a Morse set any more
                {
                    tau=tau_max+1.;
                    break;
                }

                //if refined break
                //if(done)break;
            }


            //get the new tau
            tau=tau*2;
        }

        /*   We exclude the Morse set from the refinement list  */
        if (/*!done && */tau>= tau_max)  // for this kind of refinement, if the being used tau value is larger than the tau_max, we can safely remove this Morse set
        {                                // from the priority queue
            set_exclude_MorseSet (id);
        }

#ifdef OUTPUT_PERFORMANCE
        morse_id = fopen ("refining_morseSet.txt", "a");
        fprintf (morse_id, "Finish refining Morse set %d (containing %d triangles)\n", id, old_nodelist.size());
        fclose(morse_id);
#endif
    }

    //clock_t end = clock();
    //	morse_id = fopen ("refining_morseSet.txt", "w");
    //	fprintf (morse_id, "\n\n The time spent on refining Morse sets is %f\n", (double)(end - start)/CLOCKS_PER_SEC);
    //	fclose(morse_id);

    /*
       Guoning: moved the release of the local edge list to the end of refinement 07/05/2010
    */
    free (tedges);
    tedges = nullptr;

    morse_decomp->save_dirGraph_FG((char*)"temp_dirGraph.dg");

    Cal_Regions=1;

    /*initialize all the triangle*/

    /*  The recomputation of the MCG is conducted after the refinement of Morse sets
    (03/02/2010)
        This should be changed when the bug of merging MCG is found.
        This needs to be modified after submission.
    */

    /*  first, remove all the edges and the edge lists in the Morse sets  */
    for (int i=0; i<mcg->nlist->nmnodes; i++)
    {
        if (mcg->nlist->mnodes[i]->graph_edges != nullptr)
        {
            delete mcg->nlist->mnodes[i]->graph_edges;
            mcg->nlist->mnodes[i]->graph_edges = nullptr;
        }
        mcg->nlist->mnodes[i]->nedges = 0;
    }

    for (int i=0; i<mcg->elist->curMaxNumGedges; i++)
    {
        if (mcg->elist->edges[i] != nullptr)
            delete mcg->elist->edges[i];
    }

    mcg->init_MCG();
    mcg->build_mcg_edges_graph();
    //mcg->cal_mcgedge_regions(); // the computation of the connection regions has some problems right now for cooling jacket dataset

    //RefineConnRegion(0.1);
}





void RegionTauMap::Auto_Refine_2(double tau_max)
{
    //initializetau
    vector<int> old_nodelist;//to see if the morse set has already been refined
    int id;
    double tau;

    //initialize the mesh
    int i;
    for(i=0;i<object->tlist.ntris;i++)
        object->tlist.tris[i]->exclude=false;//first no face is discard

    for (i=0; i<object->tlist.ntris; i++)
    {
        object->tlist.tris[i]->in_img = false;
        object->tlist.tris[i]->visited = false;
        object->tlist.tris[i]->used_tau = 0.;
    }


    FILE *morse_id = nullptr;

#ifdef OUTPUT_PERFORMANCE
    morse_id = fopen ("refining_morseSet.txt", "w");
    fprintf(morse_id,"So it starts ...\n");
    fclose(morse_id);
#endif

    init_Region();

    //auto refine pipeline
    clock_t start = clock();

    int nMorse = mcg->nlist->nmnodes;

    for(int iter=0;iter<nMorse;iter++)
    {
        //get current morse_id & local_tau
        //id=mcg->select_morse(min_pri);

        id = iter;

        if (!mcg->NeedRefinement(id))
            continue;

#ifdef OUTPUT_PERFORMANCE
        morse_id = fopen ("refining_morseSet.txt", "a");
        fprintf(morse_id,"\n ------------------------------------------------------\n");
        fprintf (morse_id, "%d: The current extraced Morse set to refine is %d\n", iter+1, id);
        fprintf (morse_id, "The number of edges in current MCG is %d and nodes is %d.\n", mcg->elist->nedges, mcg->nlist->nmnodes);
        fclose(morse_id);
#endif

        //tau=conley->get_tau(0);//should the morse set be 0 forever????
        //tau=min(ini_tau, tau/100.);

        tau=tau_max;

        int scc_id = mcg->nlist->mnodes[id]->scc_index;

        //each time refine one morse set
        bool done=false;


        //get the morse set nodelist
        old_nodelist.clear();
        for(int nodeID=0;nodeID<morse_decomp->scclist->scccomponents[scc_id]->nnodes;nodeID++)
        {
            int cur_node=morse_decomp->scclist->scccomponents[scc_id]->nodes[nodeID];
            old_nodelist.push_back(cur_node);
            //object->tlist.tris[old_nodelist[nodeID]]->exclude=true;
        }

        if (RefineRegion_forAuto(id, tau))
        {

            if (lmcg->nlist->nmnodes>1)
            {
                done=true;
                /*
                       We need to record the \tau value for each triangle inside this region
                    */
                for(int nodeID=0;nodeID<old_nodelist.size();nodeID++)
                {
                    object->tlist.tris[old_nodelist[nodeID]]->used_tau = tau;
                }
            }

            //if refined break
            //if(done)break;
        }


#ifdef OUTPUT_PERFORMANCE
        morse_id = fopen ("refining_morseSet.txt", "a");
        fprintf (morse_id, "Finish refining Morse set %d (containing %d triangles)\n", id, old_nodelist.size());
        fclose(morse_id);
#endif
    }

    //clock_t end = clock();
    //	morse_id = fopen ("refining_morseSet.txt", "w");
    //	fprintf (morse_id, "\n\n The time spent on refining Morse sets is %f\n", (double)(end - start)/CLOCKS_PER_SEC);
    //	fclose(morse_id);

    /*
       Guoning: moved the release of the local edge list to the end of refinement 07/05/2010
    */
    free (tedges);
    tedges = nullptr;

    morse_decomp->save_dirGraph_FG((char*)"temp_dirGraph.dg");

    Cal_Regions=1;

    /*initialize all the triangle*/

    /*  The recomputation of the MCG is conducted after the refinement of Morse sets
    (03/02/2010)
        This should be changed when the bug of merging MCG is found.
        This needs to be modified after submission.
    */

    /*  first, remove all the edges and the edge lists in the Morse sets  */
    for (int i=0; i<mcg->nlist->nmnodes; i++)
    {
        if (mcg->nlist->mnodes[i]->graph_edges != nullptr)
        {
            delete mcg->nlist->mnodes[i]->graph_edges;
            mcg->nlist->mnodes[i]->graph_edges = nullptr;
        }
        mcg->nlist->mnodes[i]->nedges = 0;
    }

    for (int i=0; i<mcg->elist->curMaxNumGedges; i++)
    {
        if (mcg->elist->edges[i] != nullptr)
            delete mcg->elist->edges[i];
    }

    mcg->init_MCG();
    mcg->build_mcg_edges_graph();
    //mcg->cal_mcgedge_regions(); // the computation of the connection regions has some problems right now for cooling jacket dataset

    //RefineConnRegion(0.1);
}





/*
    No limit for tau_max and maximal iterations
*/

void
RegionTauMap::Auto_Refine(double ini_tau,double min_pri, int ntrials)
{
    //initializetau
    vector<int> old_nodelist;//to see if the morse set has already been refined
    int id;
    double tau;

    //initialize the mesh
    for(int i=0;i<object->tlist.ntris;i++)
        object->tlist.tris[i]->exclude=false;//first no face is discard

    for (int i=0; i<object->tlist.ntris; i++)
    {
        object->tlist.tris[i]->in_img = false;
        object->tlist.tris[i]->visited = false;
        object->tlist.tris[i]->used_tau = 0.;
    }




    FILE *morse_id = nullptr;

    morse_id = fopen ("refining_morseSet.txt", "w");
    fprintf(morse_id,"");
    fclose(morse_id);

    init_Region();
    clock_t start = clock();

    //auto refine pipeline
    while(1)
    {
        //get current morse_id & local_tau
        id=mcg->select_morse(min_pri);

        morse_id = fopen ("refining_morseSet.txt", "a");
        fprintf (morse_id, "The current extraced Morse set to refine is %d\n", id);
        fclose(morse_id);

        //if no morse set in the quene, over
        if(id==-1) break;

        //tau=ini_tau;

        //if (mcg->nlist->mnodes[id]->used_tau >= tau_max)
        //{
        //	set_exclude_MorseSet (id);  // we exclude this Morse set from the refinement list
        //	/*  In the improved implementation, this Morse set can be simply dequeue from the priority queue  */
        //	continue;
        //}

        //else
        tau = max (ini_tau, mcg->nlist->mnodes[id]->used_tau);

        int scc_id = mcg->nlist->mnodes[id]->scc_index;

        //each time refine one morse set
        bool done=false;


        //get the morse set nodelist
        old_nodelist.clear();
        for(int nodeID=0;nodeID<morse_decomp->scclist->scccomponents[scc_id]->nnodes;nodeID++)
        {
            int cur_node=morse_decomp->scclist->scccomponents[scc_id]->nodes[nodeID];
            old_nodelist.push_back(cur_node);
            //object->tlist.tris[old_nodelist[nodeID]]->exclude=true;
        }

        /*
            obtain the previously used tau from the current Morse set
        */
        double pre_tau = object->tlist.tris[old_nodelist[0]]->used_tau;

        tau = max (ini_tau, 2*pre_tau/*mcg->nlist->mnodes[id]->used_tau*/);

        //each time refine one morse set
        //bool done=false;

        //int scc_id = mcg->nlist->mnodes[id]->scc_index;

        //get the morse set nodelist
        //if(old_nodelist.size())old_nodelist.clear();
        //for(int nodeID=0;nodeID<morse_decomp->scclist->scccomponents[scc_id]->nnodes;nodeID++)
        //{
        //	int cur_node=morse_decomp->scclist->scccomponents[scc_id]->nodes[nodeID];
        //	old_nodelist.push_back(cur_node);
        //}

        //for visualization
        //cur_morse=id;

        //refine the morse set
        int count = 0;
        while(count < ntrials)
        {
            used_tau = tau;

            //refine morset first
            //RefineRegion(id,tau);
            RefineRegion_forAuto(id,tau);

            //if (mcg->nlist->mnodes[id]->cancelled)
            //	break;

            //check if the morse set refined
            //get the first morse id in old_nodelist
            //int morse1=morse_decomp->scclist->scccomponents[morse_decomp->dg->nlist->dirnodes[old_nodelist[0]]->sscomp_index]->node_index;

            //compare the rest of the old_nodelist
            //int pre_id = -1;
            //for(int nodeID=0;nodeID<old_nodelist.size();nodeID++)
            //{
            //	int cur_id=morse_decomp->scclist->scccomponents[morse_decomp->dg->nlist->dirnodes[old_nodelist[nodeID]]->sscomp_index]->node_index;

            //	if (cur_id < 0)
            //		continue;

            //	else if (pre_id < 0)
            //		pre_id = cur_id;

            //	//refined already, need to jump out
            //	else if(cur_id!=pre_id)
            //	{
            //		done=true;
            //		break;
            //	}
            //}

            //if refined break
            //if(done)break;

            /*
               We need to record the \tau value for each triangle inside this region
            */
            //for(int nodeID=0;nodeID<morse_nlist.size();nodeID++)
            //{
            //	object->tlist.tris[morse_nlist[nodeID]]->used_tau = tau;
            //}

            if (lmcg->nlist->nmnodes>1)
            {
                done=true;
                /*
                       We need to record the \tau value for each triangle inside this region
                    */
                for(int nodeID=0;nodeID<old_nodelist.size();nodeID++)
                {
                    object->tlist.tris[old_nodelist[nodeID]]->used_tau = tau;
                }
                break;
            }

            else if (lmcg->nlist->nmnodes==0) // This is not a Morse set any more
            {
                //tau=tau_max+1.;
                set_exclude_MorseSet(id);
                break;
            }

            /*
               update the region triangle list
            */
            DynList_Int *temp_tris = new DynList_Int(old_nodelist.size());
            for (int nodeID=0; nodeID<old_nodelist.size(); nodeID++)
            {
                int cur_id=morse_decomp->scclist->scccomponents[morse_decomp->dg->nlist->dirnodes[old_nodelist[nodeID]]->sscomp_index]->node_index;

                if (cur_id < 0)
                    continue;
                temp_tris->add_New(old_nodelist[nodeID]);
            }

            /**/
            old_nodelist.clear();
            for (int nodeID=0; nodeID<temp_tris->nelems; nodeID++)
                old_nodelist.push_back(temp_tris->elems[nodeID]);

            delete temp_tris;


            //get the new tau
            tau=tau*2;

            count ++;
        }

        /*   We exclude the Morse set from the refinement list  */
        if (!done && count>= ntrials)
        {
            set_exclude_MorseSet (id);
        }

        morse_id = fopen ("refining_morseSet.txt", "a");
        fprintf (morse_id, "Finish refining Morse set %d\n", id);
        fclose(morse_id);

    }

    clock_t end = clock();
    morse_id = fopen ("refining_morseSet.txt", "a");
    fprintf (morse_id, "\n\n The time spent on refining Morse sets is %f\n", (double)(end - start)/CLOCKS_PER_SEC);
    fclose(morse_id);

    /*
        Guoning: we need to output the directed graph into a temporary file in case the MCG computatio can't finish
    */
    morse_decomp->save_dirGraph_FG((char*)"temp_dirGraph.dg");

    /*
       Guoning: moved the release of the local edge list to the end of refinement 07/05/2010
    */
    free (tedges);
    tedges = nullptr;


    Cal_Regions=1;

    /*initialize all the triangle*/

    /*  The recomputation of the MCG is conducted after the refinement of Morse sets
    (03/02/2010)
        This should be changed when the bug of merging MCG is found.
        This needs to be modified after submission.
    */

    ///////////////////////////////////////////////////////////////
    /*
         comment out by Guoning 03/06/2010
    */

    ///*  first, remove all the edges and the edge lists in the Morse sets  */
    //for (int i=0; i<mcg->nlist->nmnodes; i++)
    //{
    //	if (mcg->nlist->mnodes[i]->graph_edges != nullptr)
    //		delete mcg->nlist->mnodes[i]->graph_edges;
    //	mcg->nlist->mnodes[i]->nedges = 0;
    //}

    //for (int i=0; i<mcg->elist->curMaxNumGedges; i++)
    //{
    //	if (mcg->elist->edges[i] != nullptr)
    //		delete mcg->elist->edges[i];
    //}
    //mcg->init_MCG();
    //mcg->build_mcg_edges_graph();
    //mcg->cal_mcgedge_regions(); // the computation of the connection regions has some problems right now for cooling jacket dataset

    //RefineConnRegion(0.1);
    /////////////////////////////////////////////////////////////

    // For this computation, the labeling of the Morse sets is not that crucial
    //if(mcg != nullptr)
    //	delete mcg;
    //
    ////fp=fopen("detect_porbit.txt","w");
    ////fprintf(fp, "start initializing the MCG\n");
    ////fclose(fp);

    //mcg = new MCG_Graph();
    //mcg->init_MCG();
    //mcg->build_mcg();

    /*  first, remove all the edges and the edge lists in the Morse sets  */

    // Commented out by Guoning for debugging 07/15/2010
    //for (int i=0; i<mcg->nlist->nmnodes; i++)
    //{
    //	if (mcg->nlist->mnodes[i]->graph_edges != nullptr)
    //	{
    //		delete mcg->nlist->mnodes[i]->graph_edges;
    //		mcg->nlist->mnodes[i]->graph_edges = nullptr;
    //	}
    //	mcg->nlist->mnodes[i]->nedges = 0;
    //}

    //for (int i=0; i<mcg->elist->curMaxNumGedges; i++)
    //{
    //	if (mcg->elist->edges[i] != nullptr)
    //		delete mcg->elist->edges[i];
    //}

    //mcg->init_MCG();
    //mcg->build_mcg_edges_graph();
}



/*
    The Morse set and SCC information is stored in the global Morse decomposition variable "morse_decomp"
*/
void
RegionTauMap::set_exclude_MorseSet(int id)
{
    int scc_id = mcg->nlist->mnodes[id]->scc_index;
    SCComponent *scc = morse_decomp->scclist->scccomponents[scc_id];

    int i;
    for (i=0; i<scc->nnodes; i++)
    {
        object->tlist.tris[scc->nodes[i]]->exclude = true;
    }
}

void
RegionTauMap::set_MorseSet_id_tri(int id)
{
    int scc_id = mcg->nlist->mnodes[id]->scc_index;
    SCComponent *scc = morse_decomp->scclist->scccomponents[scc_id];

    int i;
    for (i=0; i<scc->nnodes; i++)
    {
        object->tlist.tris[scc->nodes[i]]->local_index = scc_id;
    }
}

void RegionTauMap::BuildLocalRegion(void)
{
    init_local_graph();
    AddLocalNodes();
    AddLocalEdges();
}

void RegionTauMap::init_local_graph(void)
{
    //build new local_decomp
    if(local_decomp != nullptr)
        delete local_decomp;
    local_decomp = new MorseDecomp(); /*initialize the Morse Decomposition component*/

    //build new diagraph
    //if(local_decomp->dg != nullptr)
    //	delete local_decomp->dg;

    //initialize dg
    int n=/*conley->*/morse_nlist.size();
    local_decomp->dg = new DirGraph(n,num_tedges/*300*/);

    /*allocate the real space for the node list and edge list*/
    for(int i = 0; i < local_decomp->dg->nlist->curMaxNumENodes; i++)
    {
        local_decomp->dg->nlist->dirnodes[i] = new DirGraph_Node();
        //local_decomp->dg->nlist->dirnodes[i]->edges = nullptr;
        //local_decomp->dg->nlist->dirnodes[i]->nedges = 0;
    }

    for(int i = 0; i < local_decomp->dg->elist->curMaxNumGedges; i++)
    {
        local_decomp->dg->elist->edges[i] = new Graph_Edge();
    }

}

void RegionTauMap::AddLocalEdges(void)
{
    //add edges to global morse object as well. Make sure the edges add to both global and local
    //do not add repeated edge, or the two ends are in the same triangle, or the edge across the boundary of the morse set
    int node1ID;
    int node2ID;
    int n1sccID;
    int n2sccID;
    int localN1;
    int localN2;

    //record the original end of elist
    orEndDirE=morse_decomp->dg->elist->nedges-1;
    //add edge to local morse object
    for(int i=0;i<num_tedges;i++)
    {
        //get two ends of the edge
        node1ID=tedges[i].node_index1;
        node2ID=tedges[i].node_index2;

        //get two ends indices in the local morse set
        localN1=object->tlist.tris[node1ID]->local_index;
        localN2=object->tlist.tris[node2ID]->local_index;

        //get the scc_id of the two ends
        n1sccID=morse_decomp->dg->nlist->dirnodes[node1ID]->sscomp_index;
        n2sccID=morse_decomp->dg->nlist->dirnodes[node2ID]->sscomp_index;

        if(n1sccID==n2sccID && cur_SCC_id == n1sccID && node1ID !=node2ID)
        {
            if(!morse_decomp->dg->is_repeated_edge_2(node1ID, node2ID))
            {
                //global
                morse_decomp->dg->add_to_edgelist(node1ID, node2ID, morse_decomp->dg->elist->nedges);
                morse_decomp->dg->add_edge_to_node(node1ID, morse_decomp->dg->elist->nedges-1);
                morse_decomp->dg->add_edge_to_node(node2ID, morse_decomp->dg->elist->nedges-1);
                //local
                local_decomp->dg->add_to_edgelist(localN1,localN2,local_decomp->dg->elist->nedges);
                local_decomp->dg->add_edge_to_node(localN1,local_decomp->dg->elist->nedges-1);
                local_decomp->dg->add_edge_to_node(localN2,local_decomp->dg->elist->nedges-1);
            }
        }

        /*
            Modified by Guoning on 07/02/2010
        */
        //if(/*n1sccID==n2sccID && cur_SCC_id == n1sccID &&*/ node1ID !=node2ID)
        //{
        //	if(!morse_decomp->dg->is_repeated_edge_2(node1ID, node2ID))
        //	{
        //		//global
        //		morse_decomp->dg->add_to_edgelist(node1ID, node2ID, morse_decomp->dg->elist->nedges);
        //		morse_decomp->dg->add_edge_to_node(node1ID, morse_decomp->dg->elist->nedges-1);
        //		morse_decomp->dg->add_edge_to_node(node2ID, morse_decomp->dg->elist->nedges-1);
        //		//local

        //		if (localN1>=0&&localN1<morse_nlist.size()&&localN2>=0&&localN2<morse_nlist.size())
        //		{
        //		local_decomp->dg->add_to_edgelist(localN1,localN2,local_decomp->dg->elist->nedges);
        //		local_decomp->dg->add_edge_to_node(localN1,local_decomp->dg->elist->nedges-1);
        //		local_decomp->dg->add_edge_to_node(localN2,local_decomp->dg->elist->nedges-1);
        //		}
        //	}
        //}

    }
}

void RegionTauMap::AddLocalNodes(void)
{
    //record the local index of vertex to the global index
    for(int i=0;i<object->tlist.ntris;i++)
        object->tlist.tris[i]->cur_morse=false;

    for(int i=0;i</*conley->*/morse_nlist.size();i++)
    {
        //record the local index of vertex to the global index
        Triangle* t= object->tlist.tris[/*conley->*/morse_nlist[i]];
        t->local_index=i;
        t->cur_morse=true;
        //Add Local Nodes to the new morse object
        local_decomp->dg->nlist->dirnodes[i]->node_index=i;
        local_decomp->dg->nlist->dirnodes[i]->global_index=/*conley->*/morse_nlist[i];
    }
}

/////////////////////////////////////////////////////////////////////////////////////
/*
    Guoning's version
*/

void RegionTauMap::BuildLocalRegion_2(void)
{
    reset_local_region_2();
    AddLocalNodes();
    AddLocalEdges();
}

void RegionTauMap::init_local_graph_2(void)
{
    //build new local_decomp
    if(local_decomp != nullptr)
        delete local_decomp;
    local_decomp = new MorseDecomp(); /*initialize the Morse Decomposition component*/

    //initialize dg
    local_decomp->dg = new DirGraph(object->tlist.ntris,object->elist.nedges);

    /*allocate the real space for the node list and edge list*/
    for(int i = 0; i < local_decomp->dg->nlist->curMaxNumENodes; i++)
    {
        local_decomp->dg->nlist->dirnodes[i] = new DirGraph_Node();
    }

    for(int i = 0; i < local_decomp->dg->elist->curMaxNumGedges; i++)
    {
        local_decomp->dg->elist->edges[i] = new Graph_Edge();
    }

}

void RegionTauMap::reset_local_region_2(void)
{
    unsigned int i;
    /*
        we need to clear the edge list any how
    */
    for (i=0; i<local_decomp->dg->nlist->ndirnodes; i++)
    {
        if (local_decomp->dg->nlist->dirnodes[i]->edges != nullptr)
        {
            free(local_decomp->dg->nlist->dirnodes[i]->edges);
            local_decomp->dg->nlist->dirnodes[i]->edges = nullptr;
            local_decomp->dg->nlist->dirnodes[i]->nedges = 0;
        }
    }
    local_decomp->dg->nlist->ndirnodes = morse_nlist.size();
    local_decomp->dg->elist->nedges = 0;
}



/////////////////////////////////////////////////////////////////////////////////////


void RegionTauMap::LocalToGlobalSCC(void)
{
    //copy the first scc to the original index
    or_scc=mcg->nlist->mnodes[refined_morse]->scc_index;
    CopySCC(0,or_scc);

    //copy the rest of the SCCs to the global
    //create new elements in the morse_decomp, and get the new dg->nsscs, scclist->nsccs
    //dg->

    newSCCN=local_decomp->scclist->nsccs-1;
    orEndSCC=morse_decomp->scclist->nsccs-1;
    morse_decomp->dg->num_sccomps+=newSCCN;

    //scclist->nsccs
    morse_decomp->scclist->nsccs=morse_decomp->dg->num_sccomps;

    //copy the scclist to a temp scclist
    //create temp scclist
    //if(tscclist!=nullptr)delete tscclist;

    tscclist = new SCCList(morse_decomp->dg->num_sccomps);
    tscclist->nsccs = morse_decomp->scclist->nsccs;
    tscclist->curMaxNumSCCS=tscclist->nsccs;
    for(int i = 0; i < tscclist->nsccs; i++)
        tscclist->scccomponents[i] = new SCComponent();

    //morse_decomp->dg->num_sccomps=morse_decomp->scclist->nsccs=orEndSCC+1;

    //copy the morse_decomp part
    for(int i=0;i<orEndSCC+1;i++)
        CopySCCTemp(i);
    //copy temp to morse_decomp
    delete morse_decomp->scclist;
    morse_decomp->scclist=tscclist;






    //reallocate the space for the new SCCs
    //morse_decomp->scclist->scccomponents=(SCComponent **)realloc(morse_decomp->scclist->scccomponents,sizeof(SCComponent *)*newSCCN);

    for(int i=1;i<local_decomp->scclist->nsccs;i++)
    {
        //create new elements
        int g=orEndSCC+i;
        //morse_decomp->scclist->scccomponents[g]=new SCComponent();
        //copy
        CopySCC(i,g);
    }


    ///*FILE* s=fopen("scclist.txt","w");
    //for(int i=0;i<morse_decomp->scclist->nsccs;i++)
    //	fprintf(s,"%d\t%d\t%d\n",i,morse_decomp->scclist->scccomponents[i]->nnodes,morse_decomp->scclist->scccomponents[i]->nfixedpoints);

    //fclose(s);*/

}
//local to global
void RegionTauMap::CopySCC(int l, int g)
{

    SCComponent *G=morse_decomp->scclist->scccomponents[g];
    SCComponent *L=local_decomp->scclist->scccomponents[l];


    //int nnodes;int *nodes; //the list of triangle lists belonging to this component
    G->nnodes=L->nnodes;
    if(G->nodes!=nullptr)
        free(G->nodes);
    if(G->nnodes)
    {
        G->nodes=(int*) malloc(sizeof(int)*G->nnodes);


        for(int i=0;i<G->nnodes;i++)
        {
            G->nodes[i]=local_decomp->dg->nlist->dirnodes[L->nodes[i]]->global_index;
            morse_decomp->dg->nlist->dirnodes[G->nodes[i]]->sscomp_index=g;
        }
    }
    else
        G->nodes=nullptr;

    //other
    G->node_index=L->node_index;
    G->valid=L->valid;
    G->nseppts=L->nseppts;
    G->nattpts=L->nattpts;
    G->nboundaries=L->nboundaries;

    //int nfixedpoints; int *singular_tri;    //the list of triangles containing singularity
    G->nfixedpoints=L->nfixedpoints;
    if(G->singular_tri!=nullptr)
        free(G->singular_tri);
    if(G->nfixedpoints)
    {
        G->singular_tri=(int*)malloc(sizeof(int)*G->nfixedpoints);
        for(int i=0;i<G->nfixedpoints;i++)
            G->singular_tri[i]=L->singular_tri[i];
    }
    else
        G->singular_tri=nullptr;

    //didn't use, so first comment it out
    //int *periodicorbits;	int nperiodicorbits;
    /*G->nperiodicorbits=L->nperiodicorbits;
    if(G->periodicorbits!=nullptr)
        free(G->periodicorbits);
    G->periodicorbits=(int*)malloc(sizeof(int)*G->nperiodicorbits);
    for(int i=0;i<G->nperiodicorbits;i++)
        G->periodicorbits[i]=L->periodicorbits[i];*/

    //conley index
    G->XM=L->XM;
    G->XL=L->XL;
    G->Conley0=L->Conley0;
    G->Conley1=L->Conley1;
    G->Conley2=L->Conley2;
    G->classification=L->classification;
    G->priority=L->priority;
    G->variance_vector=L->variance_vector;

    //record the local morse set#
    G->local_SCC=l;
    L->global_SCC=g;
}

void RegionTauMap::CopySCCTemp(int i)
{
    SCComponent *O=morse_decomp->scclist->scccomponents[i];
    SCComponent *T=tscclist->scccomponents[i];


    //int nnodes;int *nodes; //the list of triangle lists belonging to this component
    T->nnodes=O->nnodes;
    if(T->nodes!=nullptr)
        free(T->nodes);
    if(T->nnodes)
    {
        T->nodes=(int*) malloc(sizeof(int)*T->nnodes);


        for(int i=0;i<T->nnodes;i++)
        {
            T->nodes[i]=O->nodes[i];
            //morse_decomp->dg->nlist->dirnodes[T->nodes[i]]->sscomp_index=g;
        }
    }
    else
        T->nodes=nullptr;

    //other
    T->node_index=O->node_index;  // Do we copy the right Morse set ID for this SCC??!! (03/01/2010)
    T->valid=O->valid;
    T->nseppts=O->nseppts;
    T->nattpts=O->nattpts;
    T->nboundaries=O->nboundaries;

    //int nfixedpoints; int *singular_tri;    //the list of triangles containing singularity
    T->nfixedpoints=O->nfixedpoints;
    if(T->singular_tri!=nullptr)
        free(T->singular_tri);
    if(T->nfixedpoints)
    {
        T->singular_tri=(int*)malloc(sizeof(int)*T->nfixedpoints);
        for(int i=0;i<T->nfixedpoints;i++)
            T->singular_tri[i]=O->singular_tri[i];
    }
    else
        T->singular_tri=nullptr;

    //didn't use, so first comment it out
    //int *periodicorbits;	int nperiodicorbits;
    /*G->nperiodicorbits=L->nperiodicorbits;
    if(G->periodicorbits!=nullptr)
        free(G->periodicorbits);
    G->periodicorbits=(int*)malloc(sizeof(int)*G->nperiodicorbits);
    for(int i=0;i<G->nperiodicorbits;i++)
        G->periodicorbits[i]=L->periodicorbits[i];*/

    //conley index
    T->XM=O->XM;
    T->XL=O->XL;
    T->Conley0=O->Conley0;
    T->Conley1=O->Conley1;
    T->Conley2=O->Conley2;
    T->classification=O->classification;
    T->priority=O->priority;
    T->variance_vector=O->variance_vector;

}
void RegionTauMap::build_morse_nlist(int scc)
{
    if(morse_nlist.size())morse_nlist.clear();


    for(int nodeID=0;nodeID<morse_decomp->scclist->scccomponents[scc]->nnodes;nodeID++)
        morse_nlist.push_back(morse_decomp->scclist->scccomponents[scc]->nodes[nodeID]);


    //printf the list  // remove the print out debugging 06/29/2010
    //FILE *f=fopen("morse_nlist.txt","w");
    //fprintf(f,"%d\n",morse_nlist.size());
    //for(int i=0;i<morse_nlist.size();i++)
    //	fprintf(f,"%d\t%d\n",i,morse_nlist[i]);
    //fclose(f);

}


///////////////////////////////////////////////////////////////////////////////////////
//BEGIN MCG
///////////////////////////////////////////////////////////////////////////////////////
//All
void RegionTauMap::BuildLocalMCG(void)
{
    if(lmcg!=nullptr)
        delete lmcg;
    lmcg=new MCG_Graph();
    lmcg->init_MCG();
    build_mcg_local();
}

void RegionTauMap::MergeMCG(void)
{

    //make sure that there are some mcg nodes in the local MCG
    if(lmcg->nlist->nmnodes>0)
    {
        MergeMCGNodes();

        /////////////////////////////
        //printf out local morse set.
        /////////////////////////////
        //FILE *local=fopen("localmcginf.txt","w");
        //fprintf(local,"Nodes:\n");
        //for(int i=0;i<lmcg->nlist->nmnodes;i++)
        //{
        //	fprintf(local,"%d\t", lmcg->nlist->mnodes[i]->global_index);
        //}
        //fprintf(local,"\nEdges:\n");
        //for(int i=0;i<lmcg->elist->nedges;i++)
        //{
        //	fprintf(local,"%d\t%d\t%d\n",i,lmcg->nlist->mnodes[lmcg->elist->edges[i]->node_index1]->global_index,lmcg->nlist->mnodes[lmcg->elist->edges[i]->node_index2]->global_index);
        //}
        //fclose(local);

        /////////////////////////////
        //printf out local morse set.
        /////////////////////////////

        //MergeMCGEdges();  // remove the merge of the MCG edges (there was a bug here) 07/01/2010
    }
    //If there is no MCG nodes left,
    else
    {
        //FILE *local=fopen("localmcginf.txt","w");
        //fprintf(local,"Remove Morse set: %d\n", refined_morse);
        //fclose(local);

        //cancel the mcg node and its edges
        MCGNodeVanish();

        //local=fopen("localmcginf.txt","a");
        //fprintf(local,"Finish removing Morse set: %d\n", refined_morse);
        //fclose(local);
    }



    //mcg->build_mcg_edges_graph();

    mcg->layout_mcg();

    //FILE *m=fopen("mcg.txt","w");
    /*for(int i=0;i<mcg->nlist->nmnodes;i++)
        fprintf(m,"%d\t%d\t%u\n",i,mcg->nlist->mnodes[i]->scc_index,mcg->nlist->mnodes[i]->type);
    fclose(m);*/
}

void RegionTauMap::build_mcg_local(void)
{
    assign_mcgnodes_local();

    //lmcg->build_mcg_edges_graph();  // remove all the edges computation for local MCG. These will be updated after all the refinement

    //output
    /*FILE* mcgN=fopen("local_mcg.txt","w");
    fprintf(mcgN,"source#\tsink#\tsaddle#\n");
    fprintf(mcgN,"%d\t%d\t%d\n\n",lmcg->r_counter,lmcg->a_counter,lmcg->s_counter);

    fprintf(mcgN,"mcgNid\tSCCid\ttype\tmcg_labelindex\n");
    for(int i=0;i<lmcg->nlist->nmnodes;i++)
        fprintf(mcgN,"%d\t%d\t%u\t%d\n",i,lmcg->nlist->mnodes[i]->scc_index,lmcg->nlist->mnodes[i]->type,lmcg->nlist->mnodes[i]->labelindex);


    fprintf(mcgN,"\nmcgEid\tNodeId1\tNodeId2\n");
    for(int i=0;i<lmcg->elist->nedges;i++)
        fprintf(mcgN,"%d\t%d\t%d\n",i,lmcg->elist->edges[i]->node_index1,lmcg->elist->edges[i]->node_index2);
    fclose(mcgN);*/
}

//Node
void RegionTauMap::assign_mcgnodes_local(void)
{
    lmcg->cur_mcgedge_index=0;
    lmcg->cur_mcgnode_index=0;
    lmcg->nlist->nmnodes=0;
    lmcg->r_counter=lmcg->a_counter=lmcg->s_counter=0;

    scc_mcg_nodes(or_scc);//the original SCC refined in global

    for(int i=orEndSCC+1;i<morse_decomp->scclist->nsccs;i++)
        scc_mcg_nodes(i);
    //int a=0;
}

void RegionTauMap::scc_mcg_nodes(int scc_index)
{
    int i=scc_index;

    if(morse_decomp->scclist->scccomponents[i]->nnodes <= 2
        && morse_decomp->scclist->scccomponents[i]->nfixedpoints < 1)
    {;}
    else
    {
        if(morse_decomp->scclist->scccomponents[i]->nnodes >= 2
            && morse_decomp->scclist->scccomponents[i]->valid == false)
        {;}
        else{

            lmcg->add_to_mcg_nodelist(i, lmcg->cur_mcgnode_index);
            morse_decomp->scclist->scccomponents[i]->node_index=lmcg->cur_mcgnode_index-1;
        }

    }

}

void RegionTauMap::MergeMCGNodes(void)
{
    //1. get the type of the refined morse set
    unsigned char cur_type=mcg->nlist->mnodes[refined_morse]->type;
    orEndMCGN=mcg->nlist->nmnodes-1;
    orEndMCGE=mcg->elist->nedges-1;


    //2. get the first new morse set with the same type of the refined one.
    int FillNode=-1;
    int i;
    for( i=0;i<lmcg->nlist->nmnodes;i++)
    {
        if(lmcg->nlist->mnodes[i]->type==cur_type)
        {
            FillNode=i;
            break;
        }
    }

    //if there's no same node with the refined morse set
    if(i==lmcg->nlist->nmnodes)
    {
        MCGNoSameType();
        FillNode=0;//copy the first node to the refined morse set, then copy the other nodes.
    }

    //if there is at least one node has the same type as the refined morse set
    else
        //3.4.
        MCGRefinedMorseCopy(FillNode);

    //5.6.copy the rest of the local mcg nodes to the global mcg node list
    MCGNtoN(FillNode);
}






void RegionTauMap::MCGNtoN(/*int FillSCC,*/int FillNode)
{




    //5. if(( curMaxNumberOfMCGNodes-mcg->nlist->nmnodes)<lmcg->nlist->nmnodes) reallocate mcg nodes space, size=lmcg->nlist->nmnodes
    if((mcg->nlist->curMaxNumMNodes-mcg->nlist->nmnodes)<lmcg->nlist->nmnodes)
    {
        int oldnum = mcg->nlist->curMaxNumMNodes;

        /*extend the memory*/
        if(!mcg->nlist->extend(lmcg->nlist->nmnodes))
        {
            /*report errors*/
            return;
        }

        for(int i = oldnum; i < mcg->nlist->curMaxNumMNodes; i++)
            mcg->nlist->mnodes[i] = new MCG_Node();
    }


    //6. for each remaining  mcg nodes

    for(int i=0;i<lmcg->nlist->nmnodes;i++)
    {
        if(i==FillNode)continue;
        //(a)copy to the global mcg nodes list
        mcg->nlist->mnodes[mcg->nlist->nmnodes]->scc_index=lmcg->nlist->mnodes[i]->scc_index;//copy the scc index
        mcg->nlist->mnodes[mcg->nlist->nmnodes]->node_index=mcg->nlist->nmnodes;
        mcg->nlist->mnodes[mcg->nlist->nmnodes]->local_index=i;//copy the local index to the mcg
        mcg->nlist->mnodes[mcg->nlist->nmnodes]->used_tau=lmcg->nlist->mnodes[i]->used_tau;//copy the stored tau value
        lmcg->nlist->mnodes[i]->global_index=mcg->nlist->nmnodes;//copy the global index to the lmcg
        mcg->nlist->mnodes[mcg->nlist->nmnodes]->conley[0]=lmcg->nlist->mnodes[i]->conley[0];//copy the conley index
        mcg->nlist->mnodes[mcg->nlist->nmnodes]->conley[1]=lmcg->nlist->mnodes[i]->conley[1];//copy the conley index
        mcg->nlist->mnodes[mcg->nlist->nmnodes]->conley[2]=lmcg->nlist->mnodes[i]->conley[2];//copy the conley index

        mcg->nlist->mnodes[mcg->nlist->nmnodes]->priority=lmcg->nlist->mnodes[i]->priority;//Guoning for debugging 07/15/2010
        mcg->nlist->mnodes[mcg->nlist->nmnodes]->cancelled = false;

        /*
           Initialize the edge list
        */
        mcg->nlist->mnodes[mcg->nlist->nmnodes]->nedges = 0;
        mcg->nlist->mnodes[mcg->nlist->nmnodes]->graph_edges = nullptr;


        //mcg->nlist->mnodes[mcg->nlist->nmnodes]->boundaryE.clear();
        /*for(int j=0;j<lmcg->nlist->mnodes[i]->boundaryE.size();j++)
            mcg->nlist->mnodes[mcg->nlist->nmnodes]->boundaryE.push_back(lmcg->nlist->mnodes[i]->boundaryE[i]);
        lmcg->nlist->mnodes[i]->boundaryE.clear();*///empty

        //boundary
        //delete mcg->nlist->mnodes[mcg->nlist->nmnodes]->boundary;
        mcg->nlist->mnodes[mcg->nlist->nmnodes]->boundaryN=lmcg->nlist->mnodes[i]->boundaryN;//copy the number of the boundaries
        //copy the edges
        /*mcg->nlist->mnodes[mcg->nlist->nmnodes]->boundary= new int[mcg->nlist->mnodes[mcg->nlist->nmnodes]->boundaryN];
        for(int j=0;j<mcg->nlist->mnodes[mcg->nlist->nmnodes]->boundaryN;j++)
            mcg->nlist->mnodes[mcg->nlist->nmnodes]->boundary[j]=lmcg->nlist->mnodes[i]->boundary[j];*/
        mcg->nlist->mnodes[mcg->nlist->nmnodes]->boundary=lmcg->nlist->mnodes[i]->boundary;


        //(b)update the labelindex of each new mcg nodes according to the type of the nodes
        mcg->nlist->mnodes[mcg->nlist->nmnodes]->type=lmcg->nlist->mnodes[i]->type;

        switch(mcg->nlist->mnodes[mcg->nlist->nmnodes]->type)
        {
        case 0:mcg->nlist->mnodes[mcg->nlist->nmnodes]->labelindex=mcg->r_counter;
            mcg->r_counter++;
            break;
        case 1:mcg->nlist->mnodes[mcg->nlist->nmnodes]->labelindex=mcg->a_counter;
            mcg->a_counter++;
            break;
        case 2:mcg->nlist->mnodes[mcg->nlist->nmnodes]->labelindex=mcg->s_counter;
            mcg->s_counter++;
            break;
        }

        //(c)update the node_index of the corresponding SCC
        morse_decomp->scclist->scccomponents[mcg->nlist->mnodes[mcg->nlist->nmnodes]->scc_index]->node_index=mcg->nlist->nmnodes;

        //increase the mcg nodes #
        mcg->nlist->nmnodes++;
    }
    mcg->cur_mcgnode_index=mcg->nlist->nmnodes;
}

void RegionTauMap::MCGRefinedMorseCopy(int FillNode)
{
    //3. copy the new one to the old one's position, but keep the labelindex
    int FillSCC=-1;
    FillSCC=lmcg->nlist->mnodes[FillNode]->scc_index;
    mcg->nlist->mnodes[refined_morse]->scc_index=FillSCC;//copy the scc index
    mcg->nlist->mnodes[refined_morse]->local_index=FillNode;//copy the local index to the mcg
    mcg->nlist->mnodes[refined_morse]->type=lmcg->nlist->mnodes[FillNode]->type;//copy the type
    mcg->nlist->mnodes[refined_morse]->priority=lmcg->nlist->mnodes[FillNode]->priority;//copy the priority value
    //mcg->nlist->mnodes[refined_morse]->priority=1.;//Guoning for debugging 07/15/2010
    mcg->nlist->mnodes[refined_morse]->used_tau=lmcg->nlist->mnodes[FillNode]->used_tau;//copy the tau value
    mcg->nlist->mnodes[refined_morse]->conley[0]=lmcg->nlist->mnodes[FillNode]->conley[0];//copy the conley index
    mcg->nlist->mnodes[refined_morse]->conley[1]=lmcg->nlist->mnodes[FillNode]->conley[1];//copy the conley index
    mcg->nlist->mnodes[refined_morse]->conley[2]=lmcg->nlist->mnodes[FillNode]->conley[2];//copy the conley index
    //boundary
    delete mcg->nlist->mnodes[refined_morse]->boundary;
    mcg->nlist->mnodes[refined_morse]->boundaryN=lmcg->nlist->mnodes[FillNode]->boundaryN;//copy the number of the boundaries
    //copy the edges
    mcg->nlist->mnodes[refined_morse]->boundary= new int[mcg->nlist->mnodes[refined_morse]->boundaryN];
    for(int i=0;i<mcg->nlist->mnodes[refined_morse]->boundaryN;i++)
        mcg->nlist->mnodes[refined_morse]->boundary[i]=lmcg->nlist->mnodes[FillNode]->boundary[i];

    lmcg->nlist->mnodes[FillNode]->global_index=refined_morse;//copy the global index to the lmcg


    /*mcg->nlist->mnodes[refined_morse]->boundaryE.clear();
        for(int i=0;i<lmcg->nlist->mnodes[FillNode]->boundaryE.size();i++)
            mcg->nlist->mnodes[refined_morse]->boundaryE.push_back(lmcg->nlist->mnodes[FillNode]->boundaryE[i]);*/

    //lmcg->nlist->mnodes[FillNode]->boundaryE.clear();//empty



    //4. update the node_index of the corresponding SCC
    morse_decomp->scclist->scccomponents[FillSCC]->node_index=refined_morse;

}

void RegionTauMap::MCGNoSameType(void)
{
    //get the original type and current type
    unsigned char cur_type=lmcg->nlist->mnodes[0]->type;
    unsigned char ori_type=mcg->nlist->mnodes[refined_morse]->type;
    int ori_labelindex=mcg->nlist->mnodes[refined_morse]->labelindex;

    //copy the first MCG node to the refined morse set
    MCGRefinedMorseCopy(0);


    //change the type list of the MCG
    //delete the refined morse set from its type list
    //(a)assign the new refined morse the new type
    switch(cur_type)
    {
    case 0:mcg->nlist->mnodes[refined_morse]->labelindex=mcg->r_counter;
        mcg->r_counter++;
        break;
    case 1:mcg->nlist->mnodes[refined_morse]->labelindex=mcg->a_counter;
        mcg->a_counter++;
        break;
    case 2:mcg->nlist->mnodes[refined_morse]->labelindex=mcg->s_counter;
        mcg->s_counter++;
        break;
    }

    //(b)move following nodes to the previous position, because ori_type!= cur_type, so (a) won't affect the (b)
    for(int i=0;i<mcg->nlist->nmnodes;i++)
    {
        if(mcg->nlist->mnodes[i]->type==ori_type && mcg->nlist->mnodes[i]->labelindex > ori_labelindex)//the following nodes in the type list

            mcg->nlist->mnodes[i]->labelindex--;//shift to previous node
    }
    //(c)reduce the number of this type of nodes
    switch(ori_type)
    {
    case 0:
        mcg->r_counter--;
        break;
    case 1:
        mcg->a_counter--;
        break;
    case 2:
        mcg->s_counter--;
        break;
    }
}

void RegionTauMap::MCGNodeVanish(void)
{
    //cancel the node
    mcg->nlist->mnodes[refined_morse]->cancelled=true;

    //cancel the edges
    for(int i=0;i<mcg->nlist->mnodes[refined_morse]->nedges;i++)
    {
        int e=mcg->nlist->mnodes[refined_morse]->graph_edges[i];
        mcg->elist->edges[e]->cancel=true;
    }
    //build new edges
    vector<int>from;//record the nodes before the refined one
    vector<int>to;//record the nodes after the refined one
    //record from and to
    for(int i=0;i<mcg->nlist->mnodes[refined_morse]->nedges;i++)
    {
        int e=mcg->nlist->mnodes[refined_morse]->graph_edges[i];
        int n1=mcg->elist->edges[e]->node_index1;
        int n2=mcg->elist->edges[e]->node_index2;
        if(n1==refined_morse)
            to.push_back(n2);
        else
            from.push_back(n1);
    }

    //build the new edges ///// Removed by Guoning 07/06/2010
    //for(int i=0;i<from.size();i++)
    //{
    //	for(int j=0;j<to.size();j++)
    //	{
    //		//add edges
    //		mcg->add_to_mcg_edgelist(from[i],to[j]);
    //		mcg->add_edge_to_node(from[i], mcg->cur_mcgedge_index-1);
    //		mcg->add_edge_to_node(to[j],mcg->cur_mcgedge_index-1);
    //	}
    //}
}


//Try Edge

/*
    This MCG edge merging function is still buggy!
    Some time the corresponding SCC will point to a Morse set with negative index !
*/


void RegionTauMap::MergeMCGEdges(void)
{
    //initialize
    if(FromN.size())FromN.clear();
    if(ToN.size())ToN.clear();
    if(McgElist.size())McgElist.clear();

    //1.Add From and To nodes to the from list and to list, and delete the edges
    for(int i=0;i<mcg->nlist->mnodes[refined_morse]->nedges;i++)
    {
        int curE=mcg->nlist->mnodes[refined_morse]->graph_edges[i];
        if(mcg->elist->edges[curE]->node_index1==refined_morse)
        {
            //Add to To Nlist
            ToN.push_back(mcg->elist->edges[curE]->node_index2);
            //delete the edge from the other node elist
            RemoveMCGEFromNode(curE,mcg->elist->edges[curE]->node_index2);
        }
        if(mcg->elist->edges[curE]->node_index2==refined_morse)
        {
            //Add to From Nlist
            FromN.push_back(mcg->elist->edges[curE]->node_index1);
            //delete the edge from the other node elist
            RemoveMCGEFromNode(curE,mcg->elist->edges[curE]->node_index1);
        }
        //delete the edge
        mcg->elist->edges[curE]->cancel=true;
    }

    //2. delete edges from current node
    free(mcg->nlist->mnodes[refined_morse]->graph_edges);
    mcg->nlist->mnodes[refined_morse]->nedges=0;

    //print out all the nodes
    //FILE* grow=fopen("grow.txt","w");
    //fprintf(grow,"From:\t");
    //for(int i=0;i<FromN.size();i++)
    //	fprintf(grow,"%d\t",FromN[i]);
    //fprintf(grow,"\nlmcg:\t");
    //for(int i=0;i<lmcg->nlist->nmnodes;i++)
    //	fprintf(grow,"%d\t",lmcg->nlist->mnodes[i]->global_index);
    //fprintf(grow,"\nTo:\t");
    //for(int i=0;i<ToN.size();i++)
    //	fprintf(grow,"%d\t",ToN[i]);
    //
    //fclose(grow);


    //3.Grow Edges
    //3.1 From
    //for(int i=0;i<FromN.size();i++)
    //	for(int j=0;j<lmcg->nlist->nmnodes;j++)
    //	{
    //		if(ExistEdge(FromN[i], lmcg->nlist->mnodes[j]->global_index,true))
    //		{
    //			mcg->add_to_mcg_edgelist(FromN[i],lmcg->nlist->mnodes[j]->global_index);
    //			mcg->add_edge_to_node(FromN[i], mcg->cur_mcgedge_index-1);
    //			mcg->add_edge_to_node(lmcg->nlist->mnodes[j]->global_index,mcg->cur_mcgedge_index-1);
    //		}
    //	}
    ////3.2 To
    //for(int i=0;i<ToN.size();i++)
    //	for(int j=0;j<lmcg->nlist->nmnodes;j++)
    //	{
    //		if(ExistEdge( lmcg->nlist->mnodes[j]->global_index,ToN[i],false))
    //		{
    //			mcg->add_to_mcg_edgelist(lmcg->nlist->mnodes[j]->global_index,ToN[i]);
    //			mcg->add_edge_to_node(lmcg->nlist->mnodes[j]->global_index, mcg->cur_mcgedge_index-1);
    //			mcg->add_edge_to_node(ToN[i],mcg->cur_mcgedge_index-1);
    //		}
    //	}
    //
    for(int i=0;i<mcg->nlist->nmnodes;i++)
        for(int j=0;j<lmcg->nlist->nmnodes;j++)
        {
            if(IsInLMCGN(i))continue;
            if(ExistEdge(i, lmcg->nlist->mnodes[j]->global_index,true))
            {
                mcg->add_to_mcg_edgelist(i,lmcg->nlist->mnodes[j]->global_index);
                mcg->add_edge_to_node(i, mcg->cur_mcgedge_index-1);
                mcg->add_edge_to_node(lmcg->nlist->mnodes[j]->global_index,mcg->cur_mcgedge_index-1);
            }

            if(ExistEdge( lmcg->nlist->mnodes[j]->global_index,i,false))
            {
                mcg->add_to_mcg_edgelist(lmcg->nlist->mnodes[j]->global_index,i);
                mcg->add_edge_to_node(lmcg->nlist->mnodes[j]->global_index, mcg->cur_mcgedge_index-1);
                mcg->add_edge_to_node(i,mcg->cur_mcgedge_index-1);
            }
        }


    //4.Add all the new edges
    //4.1 Add local mcg edges to the McgElist
    CopyLocalEdge();
    //4.2 Add McgElist to mcg elist
    AddRest();

    mcg->remove_redundant_edges();


}




void RegionTauMap::CopyLocalEdge(void)
{
    for(int i=0;i<lmcg->nlist->nmnodes;i++)
        for(int j=0;j<lmcg->nlist->mnodes[i]->nedges;j++)
        {
            int n1=lmcg->nlist->mnodes[lmcg->elist->edges[lmcg->nlist->mnodes[i]->graph_edges[j]]->node_index1]->global_index;
            int n2=lmcg->nlist->mnodes[lmcg->elist->edges[lmcg->nlist->mnodes[i]->graph_edges[j]]->node_index2]->global_index;

            TempEdge temp;
            McgElist.push_back(temp);
            McgElist[McgElist.size()-1].node_index1=n1;
            McgElist[McgElist.size()-1].node_index2=n2;
        }
}



bool RegionTauMap::IsInADD(int node1, int node2)
{
    for(int i=0;i<McgElist.size();i++)
    {
        if(McgElist[i].node_index1==node1 && McgElist[i].node_index2==node2)return true;
    }
    return false;
}



void RegionTauMap::AddRest(void)
{
    for(int i=0;i<McgElist.size();i++)
    {
        //if(McgElist[i].visited==0)
        {
            mcg->add_to_mcg_edgelist(McgElist[i].node_index1,McgElist[i].node_index2);
            mcg->add_edge_to_node(McgElist[i].node_index1, mcg->cur_mcgedge_index-1);
            mcg->add_edge_to_node(McgElist[i].node_index2,mcg->cur_mcgedge_index-1);
        }
    }
}


bool RegionTauMap::InList(int item, vector<int> v)
{
    int i=0;
    for(i=0;i<v.size();i++)
    {
        if(v[i]==item)return true;
    }
    return false;
}

void RegionTauMap::RemoveMCGEFromNode(int e,int n)
{
    //1. find the correct index of the MCG edge in the MCG node elist
    int index;
    for(int i=0;i<mcg->nlist->mnodes[n]->nedges;i++)
    {
        if(mcg->nlist->mnodes[n]->graph_edges[i]==e)
        {
            index=i;
            break;
        }
    }

    //2. move the following edges to the former position
    for(int i=index+1;i<mcg->nlist->mnodes[n]->nedges;i++)
        mcg->nlist->mnodes[n]->graph_edges[i-1]=mcg->nlist->mnodes[n]->graph_edges[i];
    mcg->nlist->mnodes[n]->nedges--;
}


bool RegionTauMap::ExistEdge(int n1, int n2,bool from)
{

    //initialize
    int scc1=mcg->nlist->mnodes[n1]->scc_index;
    int scc2=mcg->nlist->mnodes[n2]->scc_index;

    DirGraph_Node **sccnodes = morse_decomp->dg->nlist->dirnodes;
    Graph_Edge **sccedges = morse_decomp->dg->elist->edges;

    /*reset all the flags of the visited nodes*/
    for(int i=0;i<region.size();i++)
        sccnodes[region[i]]->visited = 0;

    if(region.size())region.clear();

    //use all the triangles in the region
    for(int i=0;i<morse_decomp->scclist->scccomponents[scc1]->nnodes;i++)
    {
        region.push_back(morse_decomp->scclist->scccomponents[scc1]->nodes[i]);
    }

    //start to grow
    for(int i=0;i<region.size();i++)
    {
        /*expand current node along its outgoing edges*/
        int cur_node = region[i];

        //to check if visited
        if(sccnodes[cur_node]->visited == 2)
            continue;
        sccnodes[cur_node]->visited = 2;

        //for each edge connecting to the cur_node
        for(int j = 0; j < sccnodes[cur_node]->nedges; j++)
        {
            //if cur_node is the second node of this edge, so there is no chance to grow out
            if(sccedges[sccnodes[cur_node]->edges[j]]->node_index2 == cur_node)
                continue;
            int other_node = sccedges[sccnodes[cur_node]->edges[j]]->node_index2;

            if(sccnodes[other_node]->sscomp_index==scc2)
                return true;
            if(sccnodes[other_node]->sscomp_index != scc1
                && sccnodes[other_node]->visited == 0)
            {
                if(mcg->IsMCGN(sccnodes[other_node]->sscomp_index)
                    && morse_decomp->scclist->scccomponents[scc1]->node_index>=0
                    && morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index>=0)
                {
                    if(!mcg->is_connected_mcg(morse_decomp->scclist->scccomponents[scc1]->node_index,
                                               morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index)
                        && !mcg->nlist->mnodes[morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index]->cancelled
                        )
                    {
                        if(mcg->is_valid_link(morse_decomp->scclist->scccomponents[scc1]->node_index,
                                               morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->node_index))
                        {
                            if(sccnodes[other_node]->sscomp_index==scc2)
                                return true;
                        }


                        /*add all the nodes in this Morse sets into the array*/
                        for(int k = 0; k < morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->nnodes; k++)
                        {
                            if(sccnodes[morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->nodes[k]]->visited
                                == 1)
                                continue;
                            region.push_back(morse_decomp->scclist->scccomponents[sccnodes[other_node]->sscomp_index]->nodes[k]);

                            /*set the flag*/
                            sccnodes[region[region.size()-1]]->visited = 1;
                        }
                    }
                }
                else
                {
                    region.push_back(other_node);
                    sccnodes[region[region.size()-1]]->visited = 1;
                }
            }
        }
    }

    //if all the nodes are growing and there is no node hit the morse set, then there is no edge
    return false;
}


bool RegionTauMap::IsInLMCGN(int node)
{
    for(int i=0;i<lmcg->nlist->nmnodes;i++)
    {
        if(lmcg->nlist->mnodes[i]->global_index==node)return true;
    }
    return false;
}



////////////////////////////////
//refine connection regions
////////////////////////////////
//refine region of one mcg edge
void RegionTauMap::RefineConnMCGE(int node1, int node2,double tau)
{
    //1.find the MCG edge
    int e=FindMCGEdgeFromNodes(node1, node2);

    //2.refine, update the flow combinatorlialization
    //2.1 get the trianglelist
    //yeah, just use this morse_nlist, because a lot of work are done based on this array in various functions
    if(morse_nlist.size())
        morse_nlist.clear();
    for(int i=0;i<mcg->elist->edges[e]->ntris;i++)
        morse_nlist.push_back(mcg->elist->edges[e]->triangles[i]);

    //2.2 clean connection region edges
    //2.2.1 clean the edges in the region
    CleanConnectionRegion(node1,node2);
    //2.2.2 delete the tris in mcg edges
    free(mcg->elist->edges[e]->triangles);
    mcg->elist->edges[e]->ntris=0;

    //2.3 Flow Combinatorialization
    //2.3.1 Get tau value from the two morse sets
    double tau_;
    //tau_=mcg->nlist->mnodes[node1]->used_tau;
    //if(mcg->nlist->mnodes[node2]->used_tau < tau_)  // Modified by Guoning: always use the smallest tau value
    //	tau_=mcg->nlist->mnodes[node2]->used_tau;
    //else
    tau_=tau;
    int ntris1=morse_decomp->scclist->scccomponents[mcg->nlist->mnodes[node1]->scc_index]->nnodes;
    int ntris2=morse_decomp->scclist->scccomponents[mcg->nlist->mnodes[node2]->scc_index]->nnodes;


    //tau_=tau;
    //tau_=(min(ntris1,ntris2))*tau_/10;


    //tau_=ConnectionTau();
    //2.3.2 compute the new edges

    local_flow_comb(tau_);

    //2.4 Add the new edges to
    AddConnectionRegionEdges(node1, node2);

    //3. Compute the region for the mcg edges after refinement.
    mcg->cal_mcgedge_regions();


    /*  record the tau value used to refine the triangles inside this region  */
    for(int nodeID=0;nodeID<morse_nlist.size();nodeID++)
    {
        int scc = morse_decomp->dg->nlist->dirnodes[morse_nlist[nodeID]]->sscomp_index;

        if (morse_decomp->scclist->scccomponents[scc]->node_index<0)
            object->tlist.tris[morse_nlist[nodeID]]->used_tau = tau;
    }


}
//refine all the regions, will call the function RefineConnMCGE
void RegionTauMap::RefineConnRegion(double tau)
{
    //refine all the regions by tau
    for(int i=0;i<mcg->elist->nedges;i++)
    {

        if(mcg->elist->edges[i]->ntris!=0)
        {
            /*  obtain the maximum tau value used to refine these two Morse sets  */
            double tau1 = mcg->nlist->mnodes[mcg->elist->edges[i]->node_index1]->used_tau;
            double tau2 = mcg->nlist->mnodes[mcg->elist->edges[i]->node_index2]->used_tau;
            double tau = max(tau1, tau2);
            tau = tau/4.;
            if (tau == 0) tau = 0.1;
            RefineConnMCGE(mcg->elist->edges[i]->node_index1,mcg->elist->edges[i]->node_index2,tau);
        }
    }
    //recompute the region
    mcg->cal_mcgedge_regions();
}

int RegionTauMap::FindMCGEdgeFromNodes(int node1, int node2)
{
    //search for the edges of node1
    int e;
    for(int i=0;i<mcg->nlist->mnodes[node1]->nedges;i++)
    {
        e=mcg->nlist->mnodes[node1]->graph_edges[i];
        if(mcg->elist->edges[e]->node_index2==node2)
            return e;
    }
}

void RegionTauMap::CleanConnectionRegion(int node1, int node2)
{
    //clean the edges in the region
    DirGraph_Node **sccnodes = morse_decomp->dg->nlist->dirnodes;
    Graph_Edge **sccedges = morse_decomp->dg->elist->edges;
    int num_sccedges=morse_decomp->dg->elist->nedges;

    //vector to store the M1 and M2
    vector<int>M1M2;


    //second, M1
    int scc1=mcg->nlist->mnodes[node1]->scc_index;
    for(int i=0;i<morse_decomp->scclist->scccomponents[scc1]->nnodes;i++)
        M1M2.push_back(morse_decomp->scclist->scccomponents[scc1]->nodes[i]);
    //third, M2
    int scc2=mcg->nlist->mnodes[node2]->scc_index;
    for(int i=0;i<morse_decomp->scclist->scccomponents[scc2]->nnodes;i++)
        M1M2.push_back(morse_decomp->scclist->scccomponents[scc2]->nodes[i]);

    //for each triangle in the
    for(int nodeID=0;nodeID<morse_nlist.size();nodeID++)
    {
        //get the triangle
        int cur_node=morse_nlist[nodeID];
        //searching for the edges to delete, well, it seems that all the edges need to be deleted
        for(int edgeID=0;edgeID<sccnodes[cur_node]->nedges;edgeID++)
        {
            num_sccedges = morse_decomp->dg->elist->nedges;  //modified by Guoning 09/24/09

            int cur_edge=sccnodes[cur_node]->edges[edgeID];
            //get two ends of the edge, these two variable should not be changed!
            int node1ID=sccedges[cur_edge]->node_index1;
            int node2ID=sccedges[cur_edge]->node_index2;


            if(InList(node1ID,M1M2)||InList(node2ID,M1M2))continue;

            //delete edge from sccedges
            sccedges[cur_edge]->cancel=sccedges[num_sccedges-1]->cancel;
            sccedges[cur_edge]->node_index1=sccedges[num_sccedges-1]->node_index1;
            sccedges[cur_edge]->node_index2=sccedges[num_sccedges-1]->node_index2;
            sccedges[cur_edge]->visited=sccedges[num_sccedges-1]->visited;

            //change the node.edge
            int n1c=sccedges[cur_edge]->node_index1;
            int n2c=sccedges[cur_edge]->node_index2;

            for(int edgeN1=0;edgeN1<sccnodes[n1c]->nedges;edgeN1++)
            {
                if(sccnodes[n1c]->edges[edgeN1]==num_sccedges-1)
                {
                    sccnodes[n1c]->edges[edgeN1]=cur_edge;
                    break;
                }
            }

            for(int edgeN2=0;edgeN2<sccnodes[n2c]->nedges;edgeN2++) // corrected by Guoning 09/24/09
            {
                if(sccnodes[n2c]->edges[edgeN2]==num_sccedges-1)
                {
                    sccnodes[n2c]->edges[edgeN2]=cur_edge;
                    break;
                }
            }
            morse_decomp->dg->elist->nedges--;

            //delete cur_edge from its original two end nodes
            for(int j=0;j<sccnodes[node1ID]->nedges;j++)
            {

                if(sccnodes[node1ID]->edges[j]==cur_edge)
                {
                    if(j!=sccnodes[node1ID]->nedges-1)
                        sccnodes[node1ID]->edges[j]=sccnodes[node1ID]->edges[sccnodes[node1ID]->nedges-1];
                    sccnodes[node1ID]->nedges--;

                    break;
                }
            }

            for(int j=0;j<sccnodes[node2ID]->nedges;j++)
            {

                if(sccnodes[node2ID]->edges[j]==cur_edge)
                {
                    if(j!=sccnodes[node2ID]->nedges-1)
                        sccnodes[node2ID]->edges[j]=sccnodes[node2ID]->edges[sccnodes[node2ID]->nedges-1];
                    sccnodes[node2ID]->nedges--;

                    break;
                }
            }
            edgeID--;
        }
    }
}

void RegionTauMap::AddConnectionRegionEdges(int mnode1, int mnode2)
{
    int node1ID;
    int node2ID;
    vector<int>M1M2R;//for store the triangles in M1 union M2 union R

    //Before Add the edges, we need to get the region list
    //first, the region R
    for(int i=0;i<morse_nlist.size();i++)
        M1M2R.push_back(morse_nlist[i]);
    //second, M1
    int scc1=mcg->nlist->mnodes[mnode1]->scc_index;
    for(int i=0;i<morse_decomp->scclist->scccomponents[scc1]->nnodes;i++)
        M1M2R.push_back(morse_decomp->scclist->scccomponents[scc1]->nodes[i]);
    //third, M2
    int scc2=mcg->nlist->mnodes[mnode2]->scc_index;
    for(int i=0;i<morse_decomp->scclist->scccomponents[scc2]->nnodes;i++)
        M1M2R.push_back(morse_decomp->scclist->scccomponents[scc2]->nodes[i]);


    for(int i=0;i<num_tedges;i++)
    {
        //get two ends of the edge
        node1ID=tedges[i].node_index1;
        node2ID=tedges[i].node_index2;

        //check that if the temp edge is in the region M1 union M2 union R, if in, just fill it in
        if(InList(node1ID,M1M2R)&&InList(node2ID,M1M2R))
        {
            morse_decomp->dg->add_to_edgelist(node1ID, node2ID, morse_decomp->dg->elist->nedges);
            morse_decomp->dg->add_edge_to_node(node1ID, morse_decomp->dg->elist->nedges-1);
            morse_decomp->dg->add_edge_to_node(node2ID, morse_decomp->dg->elist->nedges-1);
        }
    }
}



double RegionTauMap::TriangleTau(int tri)
{
    double tau;
    Triangle *t=object->tlist.tris[tri];
    //find the longest edge of the triangle
    double length=0.0;
    for(int i=0;i<3;i++)
    {
        if(length<t->edges[i]->length)
            length=t->edges[i]->length;
    }
    //find the shortest vector
    double vec=1000000;
    for(int i=0;i<3;i++)
    {
        double vx=t->verts[i]->g_vec.entry[0];
        double vy=t->verts[i]->g_vec.entry[1];
        double vz=t->verts[i]->g_vec.entry[2];
        double v=sqrt(vx*vx+vy*vy+vz*vz);
        if(vec>v)
            vec=v;
    }
    tau=length/vec;
    return tau;
}

double RegionTauMap::ConnectionTau()
{
    double tau=0.0;
    for(int i=0;i<morse_nlist.size();i++)
    {
        if(tau<TriangleTau(morse_nlist[i]))
            tau=TriangleTau(morse_nlist[i]);
    }
    return tau;
}




////////////////////////////////////////////////////////////////////////////////////
//// Store the previous tracing result to further improve the performance

void
RegionTauMap::init_tau_pts_at_verts()  // initialize the sample points at vertices
{
    unsigned int i;

    Trajectory *temp = new Trajectory(-1,1);
    for (i=0; i<object->vlist.nverts; i++)
    {
        Vertex *v=object->vlist.verts[i];
        int tri_id=-1;
        temp->pass_vertex(v->index, tri_id, 0);

        if(tri_id < 0)
        {
            tri_id = v->corners[0]->t;
        }

        v->tauPt_f.x = v->x;
        v->tauPt_f.y = v->y;
        v->tauPt_f.z = v->z;
        v->tauPt_f.which_tri = tri_id;
        v->tauPt_f.tau = 0.;

        temp->pass_vertex(v->index, tri_id, 1);

        if(tri_id < 0)
        {
            tri_id = v->corners[0]->t;
        }

        v->tauPt_b.x = v->x;
        v->tauPt_b.y = v->y;
        v->tauPt_b.z = v->z;
        v->tauPt_b.which_tri = tri_id;
        v->tauPt_b.tau = 0.;
    }

    delete temp;

}

void
RegionTauMap::init_tau_pts_at_tris()   // initialize the sample points at the centers of triangles
{
    unsigned int i, j;

    for (i=0; i<object->tlist.ntris; i++)
    {
        Triangle *t=object->tlist.tris[i];
        float center[3]={0.};
        for (j=0; j<t->nverts; j++)
        {
            center[0] += t->verts[j]->x;
            center[1] += t->verts[j]->y;
            center[2] += t->verts[j]->z;
        }

        center[0] = center[0]/3.;
        center[1] = center[1]/3.;
        center[2] = center[2]/3.;

        // initialize the forward tracer
        t->tauPt_f.x=center[0];
        t->tauPt_f.y=center[1];
        t->tauPt_f.z=center[2];
        t->tauPt_f.which_tri=i;
        t->tauPt_f.tau=0.;

        // initialize the backward tracer

        t->tauPt_b.x=center[0];
        t->tauPt_b.y=center[1];
        t->tauPt_b.z=center[2];
        t->tauPt_b.which_tri=i;
        t->tauPt_b.tau=0.;
    }
}


void
RegionTauMap::build_local_connected_graph_2(double tau)
{

    ///***1. first trace backward*/
    trace_local_Verts_b(tau);//rewrite

    ///*1.2. build the edges according to the result*/
    build_edges_local_Vers_b();

    ///*hard coded sampling*/
    trace_all_centers_tris_build_edges_b(tau);
    trace_all_edges_build_di_edges_adp(tau, 1);

    ///***2. perform forward tracing*/
    trace_local_Verts_f(tau);

    ///*2.2. build the edges according to the result*/
    build_edges_local_Vers_f();

    ///*hard coded sampling*/
    trace_all_centers_tris_build_edges_f(tau);
    trace_all_edges_build_di_edges_adp(tau, 0);
}





void
RegionTauMap::trace_local_Verts_f(double tau)
{
    int i, k;
    int tri_id, end_tris;
    Vertex *v;
    icVector3 stP, newP;

    double lp[2];
    double alpha[3];
    Triangle *face;
    Triangle *t;

    Trajectory *temp= new Trajectory(-1, 1);

    /////////////////////////////////////////////
    //for(i = 0; i < object->vlist.nverts; i++)
    //	object->vlist.verts[i]->visited=false;

    for(int j=0;j<morse_nlist.size();j++)
    {
        t=object->tlist.tris[morse_nlist[j]];
        for(k=0;k<t->nverts;k++)
            t->verts[k]->visited = false;
    }
    /////////////////////////////////////////////


    /////////////////////////////////////////////
    for(int j=0;j<morse_nlist.size();j++)
    {
        /////////////////////////////////////////
        t=object->tlist.tris[morse_nlist[j]];
        for(int n=0;n<t->nverts;n++)
        {
            v=t->verts[n];
            if(v->visited==false)
            {
                /////////////////////////////////////////
                v->visited=true;


                stP.entry[0] = v->tauPt_f.x;
                stP.entry[1] = v->tauPt_f.y;
                stP.entry[2] = v->tauPt_f.z;

                tri_id = v->tauPt_f.which_tri;
                morse_decomp->gt_tau = tau-v->tauPt_f.tau;


                morse_decomp->trace_Ver_f(tri_id, stP.entry, newP.entry, end_tris, tau);

                v->img_tau[0] = newP.entry[0];
                v->img_tau[1] = newP.entry[1];
                v->img_tau[2] = newP.entry[2];

                v->tauPt_f.x = newP.entry[0];
                v->tauPt_f.y = newP.entry[1];
                v->tauPt_f.z = newP.entry[2];

                v->tauPt_f.which_tri = end_tris;
                v->tauPt_f.tau = morse_decomp->gt_tau;


                v->imgtri = end_tris;
            }
        }
        ///////////////
    }

    delete temp;
}



void
RegionTauMap::trace_local_Verts_b(double tau)
{
    int i, k;
    int tri_id, end_tris;
    Vertex *v;
    icVector3 stP, newP;

    double lp[2];
    double alpha[3];
    Triangle *face;
    Triangle *t;

    Trajectory *temp= new Trajectory(-1, 1);

    /////////////////////////////////////////////
    //for(i = 0; i < object->vlist.nverts; i++)
    //	object->vlist.verts[i]->visited=false;

    for(int j=0;j<morse_nlist.size();j++)
    {
        t=object->tlist.tris[morse_nlist[j]];
        for(k=0;k<t->nverts;k++)
            t->verts[k]->visited = false;
    }
    /////////////////////////////////////////////


    /////////////////////////////////////////////
    for(int j=0;j<morse_nlist.size();j++)
    {
        /////////////////////////////////////////
        t=object->tlist.tris[morse_nlist[j]];
        for(int n=0;n<t->nverts;n++)
        {
            v=t->verts[n];
            if(v->visited==false)
            {
                /////////////////////////////////////////
                v->visited=true;


                stP.entry[0] = v->tauPt_b.x;
                stP.entry[1] = v->tauPt_b.y;
                stP.entry[2] = v->tauPt_b.z;

                tri_id = v->tauPt_b.which_tri;
                morse_decomp->gt_tau = tau-v->tauPt_b.tau;

                morse_decomp->trace_Ver_b(tri_id, stP.entry, newP.entry, end_tris, tau);

                v->img_tau[0] = newP.entry[0];
                v->img_tau[1] = newP.entry[1];
                v->img_tau[2] = newP.entry[2];

                v->tauPt_b.x = newP.entry[0];
                v->tauPt_b.y = newP.entry[1];
                v->tauPt_b.z = newP.entry[2];

                v->tauPt_b.which_tri = end_tris;
                v->tauPt_b.tau = morse_decomp->gt_tau;

                v->imgtri = end_tris;
            }
        }
        ///////////////
    }

    delete temp;
}


void
RegionTauMap::trace_all_centers_tris_build_edges_f(double tau)
{


    for(int j=0;j<morse_nlist.size();j++)
    {
        trace_center_tris_build_edge_f(morse_nlist[j], tau);

    }
}

void
RegionTauMap::trace_all_centers_tris_build_edges_b(double tau)
{


    for(int j=0;j<morse_nlist.size();j++)
    {
        trace_center_tris_build_edge_b(morse_nlist[j], tau);

    }
}


void
RegionTauMap::trace_center_tris_build_edge_f(int tri, double tau)
{
    double center[3] = {0., 0., 0.};

    Triangle *t=object->tlist.tris[tri];

    int i, endtri;

    center[0]=t->tauPt_f.x;
    center[1]=t->tauPt_f.y;
    center[2]=t->tauPt_f.z;

    tri = t->tauPt_f.which_tri;
    morse_decomp->gt_tau = tau-t->tauPt_f.tau;


    /*start tracing here*/
    morse_decomp->trace_Ver_f(tri, center, center, endtri, tau);

    if(endtri < 0)
        return;

    t->tauPt_f.x = center[0];
    t->tauPt_f.y = center[1];
    t->tauPt_f.z = center[2];

    t->tauPt_f.which_tri = endtri;
    t->tauPt_f.tau = morse_decomp->gt_tau;


    build_one_edge(tri, endtri);
}

void
RegionTauMap::trace_center_tris_build_edge_b(int tri, double tau)
{
    double center[3] = {0., 0., 0.};

    Triangle *t=object->tlist.tris[tri];

    int i, endtri;

    center[0]=t->tauPt_b.x;
    center[1]=t->tauPt_b.y;
    center[2]=t->tauPt_b.z;

    tri = t->tauPt_b.which_tri;
    morse_decomp->gt_tau = tau-t->tauPt_b.tau;


    /*start tracing here*/
    morse_decomp->trace_Ver_b(tri, center, center, endtri, tau);

    if(endtri < 0)
        return;

    t->tauPt_b.x = center[0];
    t->tauPt_b.y = center[1];
    t->tauPt_b.z = center[2];

    t->tauPt_b.which_tri = endtri;
    t->tauPt_b.tau = morse_decomp->gt_tau;

    /*backward tracing*/
    build_one_edge(endtri, tri);
}
