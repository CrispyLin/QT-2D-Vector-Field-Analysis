/*
This header file contains the class of Morse decomposition (finding strongly connected components
in the obtained directed graph).

Created and Modified by Guoning Chen
        copyright @2007
    */

#pragma once
#ifndef __MORSEDECOMP_H__
#define __MORSEDECOMP_H__

#include "VField.h"
#include <stack>
using namespace std;

class MorseDecomp{
public:
    /*we need a data structure to store the directed graph*/
    DirGraph *dg;

    SCCList *scclist;

    stack <point3> stack_pts;

    double tau_max, tau_min;  // for the auto-refinement 03/03/2010



//GraphEdge *sccedges;
//int num_sccedges = 0;
//int curMaxNumDirGraphEdges;
//GraphNode2 *sccnodes;
//int *cur_nodes_order;

    /*we need to build the directed graph*/
    MorseDecomp()
    {
        dg = NULL;
        scclist = NULL;

        tau_max = tau_min = 0;
    }

    ~MorseDecomp()
    {
        if(dg != NULL)
            delete dg;
        if(scclist != NULL)
            delete scclist;
    }

    void init_graph();
    void build_DirGraph_Tri_no_Tau(int tri, int backward);
    void build_direted_graph();

    /*using the idea of tau maps*/
    void build_multivalued_graph(double tau); /*use the idea of \tau map, will appear later*/
    void trace_all_Verts(double tau, int backward);
    void trace_Ver(int tri, double st[2], double en[2], int &en_tri, double, int);
    void build_edges_all_Vers(int backward);
    void build_edges_Ver(int, int);
    void build_one_edge(int, int);
    void trace_all_edges_build_edges(double tau, int backward);
    void trace_center_tri_build_edge(int tri, double tau, int backward);
    void trace_all_centers_tris_build_edges(double tau, int backward);
    void trace_all_edges_build_di_edges_adp(double tau, int backward);
    void init_all_edges();
    bool are_close_neighbors(int, int);
    void trace_an_edge_build_di_edges_adp(double st1[3], double st2[3], int t1, int t2, int tri,
                                      int neighbor_tri, double tau, int backward, int edge_id);
    void trace_recursive_an_edge(double v1[3], double v2[3], int &t1, int tri, int neighbor_tri,
                             double tau, int backward, int &level, int edge_id);
    void morse_decomp_tau(double);

    double gt_tau;

    /*we need to compute the strongly connected components*/
    void build_SCCElemList();
    bool is_valid_SCC(int scc_index);
    void mark_all_valid_SCCS();

    void morse_decomp();

    /*calculate the separation and attachment points for the valid SCCs*/
    void cal_sepsandatts_valid_SCCs();

    /*************************************************************/
    /*   non-recursive edge sampling  */
    void adp_edge_sampling(double st1[3], double st2[3],
                            int t1, int t2, int tri,
                            int neighbor_tri, double tau, int backward, int edge_id);
    void show_tri_mapping (int);

    bool obtain_connected_imgs(int);     // judge whether the image and preimage of a given triangle is connected  02/18/2010
    bool obtain_connected_imgs_reg(int *, int); // this will check whether the local updated graph will provide connected image for each triangle inside it



    /********************************************************
     02/11/2010
    */
    void build_direted_graph_2();
    void build_DirEdges_edge(int);

    bool mark_all_valid_SCCS(int *, int);  // for local region refinement
    bool is_valid_scc(int *, int, int);

    /********************************************************
    02/26/2010
    for local refinement tracing
    */
    void trace_Ver_local(int tri, double st[3], double en[3], int &en_tri, double tau, int backward, int scc_id);

    /********************************************************
    02/28/2010
     Save the computation results of Morse decompositions for better visualization
    */

    void save_dirGraph_FG(char *filename);
    void load_dirGraph_FG(char *filename);

    /******/
    // 07/21/2010
    void out_sub_FG(char *filename);
    void load_sub_FG(char *filename);


    /******************************************************/
    /*  03/03/2010 for evaluation
    */
    void compare_with_ECG ();  // we verify to see whether the triangles containing the detected fixed points
                               // and periodic orbits are included in the Morse sets or not
                               // we can simply check a flag to see whether this is true or not
    void mark_cur_MorseSets();

    //void compare_two_MCGs();   // we compare the difference of two MCGs in terms of the triangles included in these MCGs

    void get_t_max_min();

    /*
        Added by Guoning at 06/29/2010
    */
    void build_SCCElemList_local();


    /*
        Added by Guoning at 07/15/2010
    */
    void trace_Ver_f(int tri, double st[2], double en[2], int &en_tri, double);
    void trace_Ver_b(int tri, double st[2], double en[2], int &en_tri, double);
    void trace_all_Verts_f(double tau);
    void trace_all_Verts_b(double tau);
    void trace_all_edges_build_edges_f(double tau);
    void trace_all_edges_build_edges_b(double tau);
    void trace_center_tri_build_edge_f(int tri, double tau);
    void trace_center_tri_build_edge_b(int tri, double tau);

    /*
       Compute rotational sum using the framework of Morse decomposition
       02/07/2013
    */
    void init_rot_sum();
    void trace_all_Verts_rot_sum(double tau, int backward);
    void trace_Ver_f_rot_sum(int tri,
                         double st[3],
                         double en[3],
                         int &en_tri,
                         double tau,
                         float &rot_sum);

    void trace_Ver_b_rot_sum(int tri,
                         double st[3],
                         double en[3],
                         int &en_tri,
                         double tau,
                         float &rot_sum);

    void trace_Ver_f_smooth(int tri,
                         double st[3],
                         double en[3],
                         int &en_tri,
                         double tau,
                         float &rot_sum,
                         int step);

    void trace_Ver_b_smooth(int tri,
                         double st[3],
                         double en[3],
                         int &en_tri,
                         double tau,
                         float &rot_sum,
                         int step);

    void trace_all_Verts_smooth(double tau, int backward, int step);
    void smooth_A_field();
    void LIC_like_smoothing();
    void find_out_min_max_rot_sum();

};

bool is_connected_reg(int *tris, int ntris);

#endif /* __MORSEDECOMP_H__ */
