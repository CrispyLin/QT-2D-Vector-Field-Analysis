#pragma once
#ifndef CONLEYINDEX_H
#define CONLEYINDEX_H

#include <vector>
using namespace std ;

class ConleyIndex
{
public:
    ConleyIndex(void);
    ~ConleyIndex(void);

    //morse set node list
    vector<int> morse_nlist;
    void build_morse_nlist(int scc,int local);
    void compute_conley_index(int local);
    //
    //int type;
    int priority;

public:
    void GetConley(int scc,int local);
    int cur_scc;
    bool isLocal;

    int XM;
    int XL;
    int classification;
    int conleyB0;
    int conleyB1;
    int conleyB2;

    int variance;


    //Compute X(M)
    int computeXM(void);

    //Compute X(L)
    void computeXL(void);
    void BoundaryEdgeExtraction(void);
    vector<int>boundary_edgelist;
    vector<int>boundary_tlist;
    vector<int>exit_elist;
    vector<int>exit_vlist;

    void get_exit_elist(void);
    void get_exit_vlist(void);
    int GetDir_Geo(int edge);
    int GetDir_Tau(int edge);
    enum Edge_Character{EXIT,ENTRANCE,MIX};

    //Compute Betti Number
    void computeBetti(void);

    bool IsInScc(int node);
    bool IsInBoundary(int node);

    //morse set properties
    typedef struct SCCProperties{
        int XM;
        int XL;
        int classification;
        int conleyB0;
        int conleyB1;
        int conleyB2;
        int size;//# of tris
    } SCCProperties;
    vector<SCCProperties> pro;
public:
    bool IsInLocal(int node);

    //int select_morse(double min_priority);
    double compute_variance_vector(void);
    //double compute_min_vec(void);
    //double compute_ave_edge_length(void);
    //double get_tau(int ID);

    int GetMCGType(int scc);
    int* GetMCGBoundary(int&);
    double GetMCGPriority(void);
    void GetMCGConley(int  conley[]);
};


#endif // CONLEYINDEX_H
