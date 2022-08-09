#pragma once
#ifndef REGIONTAUMAP_H
#define REGIONTAUMAP_H

#include "Analysis/MorseDecomp.h"
#include <vector>
using namespace std;

class RegionTauMap
{
public:
    RegionTauMap(void);
    ~RegionTauMap(void);
    //ConleyIndex *conley;
    void RefineRegion(int id, double tau);
    void Auto_Refine(double tau_max,double ini_tau,double min_pri, int max_iter);
    void Auto_Refine(double ini_tau,double min_pri, int ntrials);  // DO NOT set the limit of maximal tau value and maximal iterations
    bool RefineRegion_forAuto(int id, double tau);

    void Auto_Refine_2(double tau_max);

    void init_local_graph_2(void);


private:
    void CleanRegion(int id);
    void local_flow_comb(double tau);
    void init_Region(void);
    void build_local_connected_graph(double tau);
    void add_temp_to_edgelist(void);

    typedef struct TempEdge
    {
        int edge_index;                   //Edge index (unique)
        int node_index1, node_index2;     //Indices of the two nodes this edge connects
        unsigned char visited;                      //For graph searching
        //int trajID;                       //the trajetory index of the corresponding separatrix 2/23/06
        //double flow_length;               //store the flow length of the corresponding separatrix 08/10/06
        ////double geo_length;				  //store the geometry length of the two nodes 08/10/06
        bool cancelled;
    } TempEdge;

    TempEdge *tedges;
    vector<TempEdge>CleanElist;
    vector<TempEdge>McgElist;
    /*vector<TempEdge>AddEdge;*/
    //vector<int>FromNE;
    //vector<int>ToNE;


    int curMaxNumTempEdges;
    int num_tedges;

    int cur_SCC_id;

    void trace_local_Verts(double tau, int backward);
    void build_edges_local_Vers(int backward);
    void trace_all_centers_tris_build_edges(double tau, int backward);
    void trace_all_edges_build_di_edges_adp(double tau, int backward);


    int cur_end_tedge;

    void build_edges_Ver(int vertid, int backward);
    void build_one_edge(int node_from, int node_to);
    void Region_AddToEdge(int node1, int node2, int& cur_index);
    void trace_center_tris_build_edge(int tri, double tau, int backward);
    void init_local_Edges(void);
    void trace_an_edge_build_di_edges_adp(double st1[3], double st2[3], int t1, int t2, int tri, int neighbor_tri, double tau, int backward);
    void trace_recursive_an_edge(double v1[3], double v2[3], int& t1, int tri, int neighbor_tri, double tau, int backward, int& level);
    /*   non-recursive edge sampling  */
    stack <point3> stack_pts;
    void adp_edge_sampling(double st1[3], double st2[3],
                           int t1, int t2, int tri,
                           int neighbor_tri, double tau, int backward);
    void adp_edge_sampling_2(double st1[3], double st2[3],
                             int t1, int t2, int tri,
                             int neighbor_tri, double tau, int backward);

    void BuildLocalRegion(void);
    void init_local_graph(void);
    void AddLocalEdges(void);
    void AddLocalNodes(void);
    void LocalToGlobalSCC(void);
    void CopySCC(int l, int g);
    void CopySCCTemp(int i);

    void BuildLocalRegion_2(void);  // Guoning: we try to re-use the previously allocated space and avoid repeated freeing and allocating memory
    //void AddLocalEdges_2(void);
    //void AddLocalNodes_2(void);
    void reset_local_region_2(void);



    void BuildLocalMCG(void);
    void MergeMCG(void);
    void build_mcg_local(void);
    void assign_mcgnodes_local(void);
    void scc_mcg_nodes(int scc_index);
    void MergeMCGNodes(void);
    void MergeMCGEdges(void);




    Graph_EdgeList *telist;
    vector<int> repell_region;
    vector<int>growNodes;
    vector<int>addNodes;
    vector<vector<int>> tmcgN;
    int* repell_region_lmcg;
    int ntris_repell_region;


    int cur_tedges;


    void CopyLocalEdge(void);

    bool IsInADD(int node1, int node2);

    void grow_repeller_region_graph_lmcg(int scc_index, int inverse);
    bool IsInMcgN(int scc);

    vector<int>FromN;
    vector<int>ToN;
    vector<int>ConN;//nodes connecting to the refined morse set
    vector<int>oldE;//edges that in refined morse set.
    vector<int>FromNE;
    vector<int>ToNE;


public:

    void AddRest(void);

    void build_morse_nlist(int scc);
    vector<int> morse_nlist;

    //new mcg node merge
    void MCGNtoN(/*int,*/int);
    void MCGRefinedMorseCopy(int FillNode);
    void MCGNoSameType(void);
    void MCGNodeVanish(void);
    bool InList(int item, vector<int> v);

    /*  Guoning 02/24/2010  */
    void set_exclude_MorseSet(int id);

    void set_MorseSet_id_tri(int id);

    /*  Qingqing 02/24/2010  */
    void RefineConnMCGE(int node1, int node2,double tau);
    void RefineConnRegion(double tau);
    int  FindMCGEdgeFromNodes(int node1, int node2);
    void CleanConnectionRegion(int node1, int node2);
    void AddConnectionRegionEdges(int mnode1, int mnode2);


    /**  Qingqing 03/07/2010  **/
    void RemoveMCGEFromNode(int e,int n);
    bool ExistEdge(int n1, int n2,bool from);
    vector<int>region;//nodes region to record the nodes
    bool IsInLMCGN(int node);


    /**  Guoning 07/06/2010  **/
    void init_tau_pts_at_verts();  // initialize the sample points at vertices
    void init_tau_pts_at_tris();   // initialize the sample points at the centers of triangles
    void build_local_connected_graph_2(double tau);

    void trace_all_centers_tris_build_edges_f(double tau);
    void trace_all_centers_tris_build_edges_b(double tau);
    void trace_center_tris_build_edge_f(int tri, double tau);
    void trace_center_tris_build_edge_b(int tri, double tau);
    void trace_local_Verts_f(double tau);
    void trace_local_Verts_b(double tau);


    void build_edges_Ver_f(int vertid);
    void build_edges_Ver_b(int vertid);
    void build_edges_local_Vers_f();
    void build_edges_local_Vers_b();


    /*
       Added by Guoning 07/14/2010
    */
    void trace_local_Verts_f2(double tau);
    void trace_local_Verts_b2(double tau);

    void trace_recursive_an_edge_f(double v1[3], double v2[3], int& t1, int tri, int neighbor_tri, double tau, int& level);
    void trace_recursive_an_edge_b(double v1[3], double v2[3], int& t1, int tri, int neighbor_tri, double tau, int& level);
    void trace_all_edges_build_di_edges_adp_f(double tau);
    void trace_all_edges_build_di_edges_adp_b(double tau);
    //void trace_recursive_an_edge_f(double v1[3], double v2[3], int& t1, int tri, int neighbor_tri, double tau, int& level);
    //void trace_recursive_an_edge_b(double v1[3], double v2[3], int& t1, int tri, int neighbor_tri, double tau, int& level);
    void trace_an_edge_build_di_edges_adp_f(double st1[3], double st2[3], int t1, int t2, int tri, int neighbor_tri, double tau);
    void trace_an_edge_build_di_edges_adp_b(double st1[3], double st2[3], int t1, int t2, int tri, int neighbor_tri, double tau);

    double TriangleTau(int tri);

    double ConnectionTau();
};
#endif // REGIONTAUMAP_H
