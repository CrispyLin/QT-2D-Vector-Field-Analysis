#pragma once
#ifndef VFIELD_H
#define VFIELD_H

/*VField.h
Vector field related data structures
Created and modified by Guoning Chen
copyright @2007
*/
#include <QDebug>
#include <stack>

#include "BuildGeom/Geometry.h"
#include "ConleyIndex.h"
#include "Others/common_routines.h"


using namespace std;

class ECG_Node;
extern class Trajectory;
extern class TriangleList;
extern class EdgeList;


class PeriodicOrbit{
public:
    unsigned char nnodes;
    unsigned char type;
    int index;                                        //the index of the limit cycle
    //int   trajID;                                     //the index of trajectory storing the closed orbit
    int *cellcycle;									  //store the triangle cell cycle
    int  ntris;
    Trajectory *traj;                                 //the closed streamline
        //we may consider to put it to the separatrix list!

    /*-------for ECG-graph --------if I can find a way to build the graph during analysis, we can safely remove these variables */
    int   *connected_fixedpts;                         //the list of singularities directly connect with it
    int   num_connect_fixedpts;                          //the number of singularities in the list
    int   *connected_POs;                  //the list of other limit cycles directly connect //with it
    int   num_connect_POs;                //the number of limit cycles in the list
    int   node_index;                                 //the index of corresponding node in C-graph
    ECG_Node **nodes;                                 //the list of nodes connecting with the singularity in ECG

    icVector3  handle_center;                         //the global coordinates of the handle

    /*---- Member functions---*/
    PeriodicOrbit()
    {
        cellcycle = NULL;
        ntris = 0;
        traj = NULL;
        nodes = NULL;
    }

    ~PeriodicOrbit()
    {
        if(cellcycle!=NULL)
            delete [] cellcycle;
        if(traj!=NULL)
            delete traj;
        if(nodes!=NULL)
            delete [] nodes;
    }

    inline void reset()
    {
        if(cellcycle!=NULL)
            delete [] cellcycle;
        if(traj!=NULL)
            delete traj;
        if(nodes!=NULL)
            delete [] nodes;
        ntris = 0;
        nnodes = 0;
    }

    void update_list_to_fixedpt(int singID);
    void update_list_to_cycle(int cycle);

    friend class PeriodicOrbit_Detector;
}; // end of LimitCycle class


class PeriodicOrbitList{
public:
    PeriodicOrbit **polist;
    int nporbits;
    int curMaxNumPOrbits;


    // The similar list operations
    PeriodicOrbitList(int initsize = 100) //construction
    {
        polist = (PeriodicOrbit **)malloc(sizeof(PeriodicOrbit *)*initsize);
        curMaxNumPOrbits = initsize;
        nporbits = 0;

        if(polist == NULL)
        {
            char rout[256], var[256];
            sprintf(rout, "%s", "PeriodicOrbitList Constructor");
            sprintf(var, "%s", "polist");

            //write_mem_error(rout, var, 0);
            curMaxNumPOrbits = 0;
            exit(-1);
        }

        for(int i = 0; i < initsize; i++)
            polist[i] = NULL;
    }

    ~PeriodicOrbitList()
    {
        if(polist != NULL)
        {
            for(int i = 0; i < curMaxNumPOrbits; i++)
            {
                if(polist[i] != NULL)
                    delete polist[i];
            }
        }
    }

    inline void reset_all_porbits()
    {
        if(polist != NULL)
        {
            for(int i = 0; i < curMaxNumPOrbits; i++)
            {
                if(polist[i] != NULL)
                    polist[i]->reset();
            }
        }
    }

    //add a new vertex to the end of the list, if it succeeds, return true
    inline bool append(PeriodicOrbit *s)
    {
        if(isFull ())
            if(!extend(50))
                return false;             //if not enough memory available, return false
        polist[nporbits] = s;
        //copyElem(s, polist[nporbits]);
        nporbits++;
        return true;
    }

    inline bool del_End() //delete the vertex at the end of the list
    {
        if(isEmpty())  return false;
        nporbits --;
        return true;
    }

    inline void copy_Elem(PeriodicOrbit *s, PeriodicOrbit *d)
    {
    }

    //delete the corresponding  vertex, if it succeeds, return true
    inline bool del_Node(PeriodicOrbit *s)
    {
        if(isEmpty())  return false;

        //find the vertex, if find it, delete and move the following vertices forward
        //otherwise, return false;

        int i, pos = -1;

        for(i = 0; i < nporbits; i++)
        {
            if(polist[i] == s)
            {
                pos = i;
                break;
            }
        }

        if(pos == -1) return false;

        //delete it
        for(i = pos; i < nporbits-1; i++)
        {
            //we need a copy function
            copy_Elem(polist[i], polist[i+1]);
        }

        nporbits--;

        return true;
    }

    inline bool isEmpty()  //judge whether the list is empty
    {
        if(nporbits == 0)   return true;
        return false;
    }

    inline bool isFull()
    {
        if(nporbits >= curMaxNumPOrbits) return true;
        return false;
    }

    //extend the original list, if it succeeds, return true
    inline bool extend(int step = 100)
    {
        //FILE *fp = fopen("extend_polist.txt", "w");
        //fprintf(fp, "start extending the space for the list.\n");
        //fclose(fp);

        PeriodicOrbit **temp = polist;
        polist = (PeriodicOrbit **) malloc(sizeof(PeriodicOrbit *) * (curMaxNumPOrbits + step));
        if( polist == NULL)
        {
            //fail
            polist = temp;
            char rout[256], var[256];
            sprintf(rout, "%s", "PeriodicOrbitList::extend");
            sprintf(var, "%s", "polist");

            //write_mem_error(rout, var, 1);
            curMaxNumPOrbits = 0;
            exit(-1);
            return false;
        }


        /*copy the pointers*/
        int i;
        for(i = 0; i < curMaxNumPOrbits; i++)
            polist[i] = temp[i];

        for(i = curMaxNumPOrbits; i < curMaxNumPOrbits+step; i++)
            polist[i] = NULL;
        curMaxNumPOrbits += step;

        //fp = fopen("extend_polist.txt", "w");
        //fprintf(fp, "start allocating memeory for new periodic orbits\n");
        //fprintf(fp, "Current maximum number of POs in the list is %d.\n",
        //	nporbits);
        //fclose(fp);

        free(temp);
        return true;
    }

    inline void reset()
    {
        nporbits = 0;
    }

};


class Basin{
public:
    int  index;
    int *Repellers;           //the indices of the repelling nodes in the ECG graph
    int *Attractors;          //the indices of the attracting nodes in the ECG graph
    int *Saddles;             //the indices of the saddle nodes in the ECG graph
    int *trajs;                   //the separatrices involving in the basin
    unsigned char num_saddle;
    unsigned char num_attractors;
    unsigned char num_repellers;
    unsigned char num_trajs;
}; // end of Basin class




class BasinList{
public:
    Basin **basins;
    int nbasins;
    int curMaxNumBasins;

    // The similar list operations
    BasinList(int initsize = 1000) //construction
    {
        basins = (Basin **)malloc(sizeof(Basin *)*initsize);
        curMaxNumBasins = initsize;
        nbasins = 0;

        if(basins == nullptr)
        {
            char rout[256], var[256];
            sprintf(rout, "%s", "BasinList Constructor");
            sprintf(var, "%s", "basins");

            //write_mem_error(rout, var, 0);
            curMaxNumBasins = 0;
            exit(-1);
        }
    }

    //add a new vertex to the end of the list, if it succeeds, return true
    inline bool append(Basin *s)
    {
        if(isFull ())
            if(!extend(100))
                return false;             //if not enough memory available, return false
        basins[nbasins] = s;
        //copyElem(s, polist[nporbits]);

        nbasins++;
        return true;
    }

    inline bool del_End() //delete the vertex at the end of the list
    {
        if(isEmpty())  return false;
        nbasins --;
        return true;
    }

    inline void copy_Elem(Basin *s, Basin *d)
    {
    }

    //delete the corresponding  vertex, if it succeeds, return true
    inline bool del_Node(Basin *s)
    {
        if(isEmpty())  return false;

        //find the vertex, if find it, delete and move the following vertices forward
        //otherwise, return false;

        int i, pos = -1;

        for(i = 0; i < nbasins; i++)
        {
            if(basins[i] == s)
            {
                pos = i;
                break;
            }
        }

        if(pos == -1) return false;

        //delete it
        for(i = pos; i < nbasins-1; i++)
        {
            //we need a copy function
            copy_Elem(basins[i], basins[i+1]);
        }

        nbasins--;

        return true;
    }

    inline bool isEmpty()  //judge whether the list is empty
    {
        if(nbasins == 0)   return true;
        return false;
    }

    inline bool isFull()
    {
        if(nbasins == curMaxNumBasins) return true;
        return false;
    }

    //extend the original list, if it succeeds, return true
    inline bool extend(int step = 100)
    {
        Basin **temp = basins;
        basins = (Basin **) malloc(sizeof(Basin *) * (curMaxNumBasins + step));
        if( basins == nullptr)
        {
            //fail
            char rout[256], var[256];
            sprintf(rout, "%s", "BasinList::extend");
            sprintf(var, "%s", "basins");

            //write_mem_error(rout, var, 1);
            curMaxNumBasins = 0;
            basins = temp;
            exit(-1);
            return false;
        }

        int i;
        for(i = 0; i < curMaxNumBasins; i++)
            basins[i] = temp[i];

        for(i=curMaxNumBasins; i<curMaxNumBasins+step; i++)
            basins[i] = nullptr;

        curMaxNumBasins += step;

        free(temp);
        return true;
    }

    inline void reset()
    {
        nbasins = 0;
    }
}; //end of BasinList class



/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//Data structures for topology graph of a vector field

class ECG_Node{
public:
    bool visited;
    bool cancelled;
    unsigned char nedges;
    unsigned char type;                //0~repeller; 1~attractor; 2~saddle
    int node_index;                   //Node index (unique)
    int fixedPtID;                //The fixed point ID that the node refers to
    int periodicorbitID;               //The periodic orbit ID that the node refers to
    int *graph_edges;               //An array of the edges incident to the node
    int labelindex;                     //for display the numeric labels of the node
    icVector2 pos;                    //the positions of the node in the display window
}; //end of ECG_Node class


class ECG_NodeList{
public:
    ECG_Node **enodes;
    int nenodes;
    int curMaxNumENodes;

    // The similar list operations
    ECG_NodeList(int initsize = 1000) //construction
    {
        enodes = (ECG_Node **)malloc(sizeof(ECG_Node *)*initsize);
        curMaxNumENodes = initsize;
        nenodes = 0;

        if(enodes == NULL)
        {
            char rout[256], var[256];
            sprintf(rout, "%s", "ECG_NodeList Constructor");
            sprintf(var, "%s", "enodes");

            //write_mem_error(rout, var, 0);
            curMaxNumENodes = 0;
            exit(-1);
        }

        for(int i = 0; i < initsize; i++)
        {
            enodes[i] = (ECG_Node *)malloc(sizeof(ECG_Node));

            if(enodes[i] == NULL)
            {
                char rout[256], var[256];
                sprintf(rout, "%s", "ECG_NodeList Constructor");
                sprintf(var, "%s", "enodes[i]");

                //write_mem_error(rout, var, 0);
                curMaxNumENodes = 0;
                exit(-1);
            }

            enodes[i]->graph_edges = NULL;
            enodes[i]->nedges = 0;
            enodes[i]->cancelled = false;
            enodes[i]->visited = false;
        }
    }

    ~ECG_NodeList()
    {
        if(enodes!=NULL)
        {
            for(int i=0; i<curMaxNumENodes; i++)
            {
                if(enodes[i] != NULL)
                    free(enodes[i]);
            }
            free(enodes);
        }
        enodes = NULL;
    }

    //add a new vertex to the end of the list, if it succeeds, return true
    inline bool append(ECG_Node *s)
    {
        if(isFull ())
            if(!extend(100))
                return false;             //if not enough memory available, return false
        enodes[nenodes] = s;
        //copyElem(s, polist[nporbits]);
        nenodes++;
        return true;
    }

    inline bool del_End() //delete the vertex at the end of the list
    {
        if(isEmpty())  return false;
        nenodes --;
        return true;
    }

    inline void copy_Elem(ECG_Node *s, ECG_Node *d)
    {
    }

    //delete the corresponding  vertex, if it succeeds, return true
    inline bool del_Node(ECG_Node *s)
    {
        if(isEmpty())  return false;

        //find the vertex, if find it, delete and move the following vertices forward
        //otherwise, return false;

        int i, pos = -1;

        for(i = 0; i < nenodes; i++)
        {
            if(enodes[i] == s)
            {
                pos = i;
                break;
            }
        }

        if(pos == -1) return false;

        //delete it
        for(i = pos; i < nenodes-1; i++)
        {
            //we need a copy function
            copy_Elem(enodes[i], enodes[i+1]);
        }

        nenodes--;

        return true;
    }

    inline bool isEmpty()  //judge whether the list is empty
    {
        if(nenodes == 0)   return true;
        return false;
    }

    inline bool isFull()
    {
        if(nenodes == curMaxNumENodes) return true;
        return false;
    }

    //extend the original list, if it succeeds, return true
    inline bool extend(int step = 100)
    {
        ECG_Node **temp = enodes;
        enodes = (ECG_Node **) malloc(sizeof(ECG_Node *) * (curMaxNumENodes + step));
        if( enodes == NULL)
        {
            //fail
            char rout[256], var[256];
            sprintf(rout, "%s", "ECG_NodeList::extend");
            sprintf(var, "%s", "enodes");

            //write_mem_error(rout, var, 1);
            //exit(-1);

            enodes = temp;
            return false;
        }

        int i;
        for(i = 0; i < curMaxNumENodes; i++)
            enodes[i] = temp[i];
        for(i = curMaxNumENodes; i < curMaxNumENodes+step; i++)
            enodes[i] = NULL;
        curMaxNumENodes += step;

        free(temp);
        return true;
    }

    inline void reset()
    {
        nenodes = 0;
    }

}; //end of ECG_NodeList class



/* Graph edge */
class Graph_Edge{
public:
    bool cancel;
    bool visited;
    int edge_index;                   //Edge index (unique)
    int node_index1, node_index2;     //Indices of the two nodes this edge connects
    //unsigned int traj;                          //the trajectory index of the corresponding separatrix

    /*the variables to represent the region corresponds to this edge 08/29/2007*/
    int *triangles;
    int ntris;

    Graph_Edge()
    {
        triangles = nullptr;
        ntris = 0;
    }

    ~Graph_Edge()
    {
        if (triangles != nullptr)
            delete [] triangles;
    }
};




class Graph_EdgeList{
public:
    Graph_Edge **edges;
    int nedges;
    int curMaxNumGedges;

    // The similar list operations as VertexList
    // The similar list operations
    Graph_EdgeList(int initsize = 1000) //construction
    {
        //edges = (Graph_Edge **)malloc(sizeof(Graph_Edge *)*initsize);
        edges = new Graph_Edge *[initsize]; // modified on 07/09/2010
        curMaxNumGedges = initsize;
        nedges = 0;

        if(edges == NULL)
        {
            char rout[256], var[256];
            sprintf(rout, "%s", "Graph_EdgeList constructor");
            sprintf(var, "%s", "edges");

            //write_mem_error(rout, var, 0);
            curMaxNumGedges = 0;
            exit(-1);
        }

        for(int i = 0; i < initsize; i++)
            edges[i] = NULL;

    }

    ~Graph_EdgeList()
    {
        if(edges != NULL)
        {
            int i;
            for(i = 0; i < curMaxNumGedges/*nedges*/; i++)
            {
                if(edges[i]!=NULL)
                    //free(edges[i]);
                    delete edges[i];
            }
            //free(edges);
            delete [] edges;
        }
    }

    //add a new vertex to the end of the list, if it succeeds, return true
    inline bool append(Graph_Edge *s)
    {
        if(isFull ())
            if(!extend(100))
                return false;             //if not enough memory available, return false
        edges[nedges] = s;
        //copyElem(s, polist[nporbits]);
        nedges++;
        return true;
    }

    inline bool del_End() //delete the vertex at the end of the list
    {
        if(isEmpty())  return false;
        nedges --;
        return true;
    }

    inline void copy_Elem(Graph_Edge *s, Graph_Edge *d)
    {
    }

    //delete the corresponding  vertex, if it succeeds, return true
    inline bool del_Node(Graph_Edge *s)
    {
        if(isEmpty())  return false;

        //find the vertex, if find it, delete and move the following vertices forward
        //otherwise, return false;

        int i, pos = -1;

        for(i = 0; i < nedges; i++)
        {
            if(edges[i] == s)
            {
                pos = i;
                break;
            }
        }

        if(pos == -1) return false;

        //delete it
        for(i = pos; i < nedges-1; i++)
        {
            //we need a copy function
            copy_Elem(edges[i], edges[i+1]);
        }

        nedges--;

        return true;
    }

    inline bool isEmpty()  //judge whether the list is empty
    {
        if(nedges == 0)   return true;
        return false;
    }

    inline bool isFull()
    {
        if(nedges == curMaxNumGedges) return true;
        return false;
    }

    //extend the original list, if it succeeds, return true
    inline bool extend(int step = 100)
    {
        Graph_Edge **temp = edges;
        //edges = (Graph_Edge **) malloc(sizeof(Graph_Edge *) * (curMaxNumGedges + step));
        edges = new Graph_Edge *[curMaxNumGedges + step];
        if( edges == NULL)
        {
            //fail
            char rout[256], var[256];
            sprintf(rout, "%s", "Graph_Edge::extend");
            sprintf(var, "%s", "edges");

            //write_mem_error(rout, var, 1);
            //exit(-1);

            edges = temp;
            return false;
        }

        int i;
        for(i = 0; i < curMaxNumGedges; i++)
            edges[i] = temp[i];

        for(i = curMaxNumGedges; i < curMaxNumGedges+step; i++)
            edges[i] = NULL;

        curMaxNumGedges += step;

        //free(temp);
        delete [] temp;

        return true;
    }

    inline void reset()
    {
        nedges = 0;
    }

}; //end of Graph_EdgeList class







class DirGraph_Node{
public:
    int node_index;                   //Node index (unique)
    int *edges;                       //An array of the edges incident to the node
    int nedges;
    int levels;                       //the level during DFS
    int parent;                       //store the parent of current nodes for backward tracking to build the SCC
    int sscomp_index;
    unsigned char visited;
    int global_index;                 //Node index in the global 12/03/2009 by Qingqing
    bool rev_visited;

    DirGraph_Node()
    {
        edges = nullptr;
        nedges = 0;
    }
};

class DirGraph_NodeList{
public:
    DirGraph_Node **dirnodes;
    int ndirnodes;
    int curMaxNumENodes;

    // The similar list operations
    DirGraph_NodeList(int initsize = 1000) //construction
    {
        //dirnodes = (DirGraph_Node **)malloc(sizeof(DirGraph_Node *)*initsize);
        dirnodes = new DirGraph_Node *[initsize];
        curMaxNumENodes = initsize;
        ndirnodes = 0;

        if(dirnodes == nullptr)
        {
            char rout[256], var[256];
            sprintf(rout, "%s", "DirGraph_NodeList Constructor");
            sprintf(var, "%s", "dirnodes");

            //write_mem_error(rout, var, 0);
            curMaxNumENodes = 0;
            exit(-1);
        }

        for(int i = 0; i < initsize; i++)
            dirnodes[i] = nullptr;
    }

    ~DirGraph_NodeList()
    {
        if(dirnodes != nullptr)
        //free(dirnodes);
        {
            //  we may need to get inside to delete each element
            for (int i=0; i<ndirnodes; i++)
                delete dirnodes[i];
            delete [] dirnodes;
        }

    }

    //add a new vertex to the end of the list, if it succeeds, return true
    inline bool append(DirGraph_Node *s)
    {
        if(isFull ())
            if(!extend(100))
                return false;             //if not enough memory available, return false
        dirnodes[ndirnodes] = s;
        //copyElem(s, polist[nporbits]);
        ndirnodes++;
        return true;
    }

    inline bool del_End() //delete the vertex at the end of the list
    {
        if(isEmpty())  return false;
        ndirnodes --;
        return true;
    }

    inline void copy_Elem(DirGraph_Node *s, DirGraph_Node *d)
    {
        d->node_index = s->node_index;
        d->levels = s->levels;
        d->parent = s->parent;
        d->sscomp_index = s->sscomp_index;
        d->visited = s->visited;
        d->nedges = s->nedges;

        d->edges = new int[d->nedges];
        int i;
        for(i = 0; i < d->nedges; i++)
            d->edges[i] = s->edges[i];
    }

    //delete the corresponding  vertex, if it succeeds, return true
    inline bool del_Node(DirGraph_Node *s)
    {
        if(isEmpty())  return false;

        //find the vertex, if find it, delete and move the following vertices forward
        //otherwise, return false;

        int i, pos = -1;

        for(i = 0; i < ndirnodes; i++)
        {
            if(dirnodes[i] == s)
            {
                pos = i;
                break;
            }
        }

        if(pos == -1) return false;

        //delete it
        for(i = pos; i < ndirnodes-1; i++)
        {
            //we need a copy function
            copy_Elem(dirnodes[i], dirnodes[i+1]);
        }

        ndirnodes--;

        return true;
    }

    inline bool isEmpty()  //judge whether the list is empty
    {
        if(ndirnodes == 0)   return true;
        return false;
    }

    inline bool isFull()
    {
        if(ndirnodes == curMaxNumENodes) return true;
        return false;
    }

    //extend the original list, if it succeeds, return true
    inline bool extend(int step = 100)
    {
        DirGraph_Node **temp = dirnodes;
        //dirnodes = (DirGraph_Node **) malloc(sizeof(DirGraph_Node *) * (curMaxNumENodes + step));
        dirnodes = new DirGraph_Node *[curMaxNumENodes + step];  // modified on 07/05/2010
        if( dirnodes == nullptr)
        {
            //fail
            char rout[256], var[256];
            sprintf(rout, "%s", "DirGraph_NodeList::extend");
            sprintf(var, "%s", "dirnodes");

            //write_mem_error(rout, var, 0);
            exit(-1);

            dirnodes = temp;
            return false;
        }

        int i;

        for(i = 0; i < curMaxNumENodes; i++)
            dirnodes[i] = temp[i];

        for(i=curMaxNumENodes; i<curMaxNumENodes+step; i++)
            dirnodes[i]=nullptr;

        curMaxNumENodes += step;

        //free(temp);
        delete [] temp;

        return true;
    }

    inline void reset()
    {
        ndirnodes = 0;
    }

}; //end of DirGraph_NodeList class


/**Directed graph for encoding the dynamics of the flow **/
class DirGraph{
public:
    DirGraph_NodeList *nlist;
    Graph_EdgeList *elist;

    int *cur_nodes_order; /*use for finding strongly connected components*/
    int cur_sccnode_index;
    int num_sccomps;

    stack<int> dfs_stack;



    /**Member functions**/
    DirGraph(int nnodes = 0, int nedges = 0)
    {

        if(nnodes == 0)
            nlist = nullptr;
        else
        {
            nlist = new DirGraph_NodeList(nnodes);
            nlist->curMaxNumENodes = nnodes;
            nlist->ndirnodes = nnodes;
        }

        if(nedges == 0)
            elist = nullptr;
        else
        {
            elist = new Graph_EdgeList(nedges);
            elist->curMaxNumGedges = nedges;
            elist->nedges = 0;
        }

        //dfs_stack = nullptr;
        while(!dfs_stack.empty())
            dfs_stack.pop();
    }

    ~DirGraph()
    {
        if(nlist != nullptr)
            delete nlist;
        if(elist != nullptr)
            delete elist;
    }

    void init_graph()
    {
        int i;
        for(i = 0; i < nlist->curMaxNumENodes; i++)
        {
            if(nlist->dirnodes[i] == nullptr)
                nlist->dirnodes[i] = new DirGraph_Node();
        }

        for(i = 0; i < elist->curMaxNumGedges; i++)
        {
            if(elist->edges[i] == nullptr)
                elist->edges[i] = new Graph_Edge();
        }
    }

    /*For the directed graph construction*/
    bool is_repeated_edge(int from, int to);
    bool is_repeated_edge_2(int from, int to);
    void add_to_edgelist(int node1, int node2, int &cur_index);
    void add_edge_to_node(int node, int edgeindex);

    //Graph related algorithms
    //Dijkstra
    //BFS
    //DFS
    //MST

    /*Find the strongly connected components*/
    void find_SCCS(); /*finding all SCCs using double DFS visiting*/
    void DFS_visit_nonrecursive(int u, int scc_index, int inverse); /*non-recursive DFS*/
    void DFS_visit(int u, int &time, int inverse); /*recursive DFS*/
    void add_to_dfsstack(int elem);

    void DFS(int inverse); /*DFS traverse all the nodes of given directed graph*/
    void reverse_edges();
    //void build_SCCElemList();

};



/** ECG Graph **/
class ECG_Graph{
public:
    ECG_NodeList *nlist;
    Graph_EdgeList *elist;
    int r_counter, a_counter, s_counter;
    int cur_ecgnode_index, cur_ecgedge_index;

    ECG_Graph(int nnodes = 50, int nedges = 50)
    {
        nlist = new ECG_NodeList(nnodes);
        elist = new Graph_EdgeList(nedges);

        //boundary_nodes = nullptr;
        //nboundarytris = 0;
    }

    ~ECG_Graph()
    {
        delete nlist;
        delete elist;
    }

    void init_ECG()
    {
        int i;
        for(i = 0; i < elist->curMaxNumGedges; i++)
            elist->edges[i] = new Graph_Edge();
    }

    /*we now need the routines to find all the possible connections here*/
    void trace_from_saddle_to_PO(int saddle); /*link between saddles and periodic orbits*/
    void trace_from_PO_to_PO(double s[3], int triangle, int cycleID, int type);
    void trace_from_nonsaddle_to_PO(int singID);
    void trace_and_build_connection(double s[3], int triangle, int singID, int type);
    void trace_connection_cycle(int cycle);

    /*member functions for MCG construction and modification*/
    /*build the mcg according to the obtained Morse sets*/
    void build_ecg();
    void assign_ecgnodes();
    void layout_ecg();
    void build_ecg_edges();
    void add_to_ecg_nodelist(int type, int &index, int fixedorcycle);
    void add_to_ecg_edgelist(int node1, int node2);
    void add_edge_to_node(int node, int edgeindex);

    bool has_interval(int node1, int node2);
    bool is_connected(int node1, int node2);

    void build_connections();

    //Graph related algorithms
    //Dijkstra
    //BFS
    //DFS
    //MST
}; //end of ECG_Graph class


/* MCG_Node class */
class MCG_Node{
public:
    bool visited;                     //for graph operations
    bool cancelled;                   //for cancellation operation
    unsigned char type;               //0-source_like, 1-sink_like, 2-saddle_like
    int node_index;
    //int basin_index;
    int scc_index;                    // the corresponding SCC index
    int *graph_edges;                 //An array of the edges incident to the node
    int nedges;
    int labelindex;                   //for display the numeric labels of the node
    icVector2 pos;                    //the positions of the node in the display window
    int local_index;				  //index of local mcg, added 12/07/2009
    int global_index;				  //index of global mcg, added 12/07/2009



    //conley
    double priority;//priority value in the priority queue when selecting the morse set to refine, added 2/2/2010
    vector<int> boundaryE;
    //for the boundaries
    int boundaryN;//number of boundary edges
    int *boundary;//boundary edges index
    int conley[3];//conley index 0,1,2

    /*  we record the t value used to obtain this Morse set 02/24/2010 */
    double used_tau;

}; //end of MCG_Node class


class MCG_NodeList{
public:
    MCG_Node **mnodes;
    int nmnodes;
    int curMaxNumMNodes;

    // The similar list operations
    MCG_NodeList(int initsize = 1000) //construction
    {
        mnodes = (MCG_Node **)malloc(sizeof(MCG_Node *)*initsize);
        curMaxNumMNodes = initsize;
        nmnodes = 0;

        if(mnodes == nullptr)
        {
            char rout[256], var[256];
            sprintf(rout, "%s", "MCG_NodeList Constructor");
            sprintf(var, "%s", "mnodes");

            //write_mem_error(rout, var, 0);
            exit(-1);
        }

        for(int i = 0; i < initsize; i++)
        {
            mnodes[i] = (MCG_Node *)malloc(sizeof(MCG_Node));

            if(mnodes[i] == nullptr)
            {
                char rout[256], var[256];
                sprintf(rout, "%s", "MCG_NodeList Constructor");
                sprintf(var, "%s", "mnodes[i]");

                //write_mem_error(rout, var, 0);
                exit(-1);
            }

            mnodes[i]->graph_edges = nullptr;
            mnodes[i]->nedges = 0;
            mnodes[i]->cancelled = false;
            mnodes[i]->visited = false;
        }
    }

    ~MCG_NodeList()
    {
        if(mnodes != nullptr)
        {
            int i;
            for(i = 0; i < curMaxNumMNodes; i++)
            {
                if(mnodes[i] != nullptr)
                    free(mnodes[i]);
            }
            free(mnodes);
        }
    }

    //add a new vertex to the end of the list, if it succeeds, return true
    inline bool append(MCG_Node *s)
    {
        if(isFull ())
            if(!extend(100))
                return false;             //if not enough memory available, return false
        mnodes[nmnodes] = s;
        //copyElem(s, polist[nporbits]);
        nmnodes++;
        return true;
    }

    inline bool del_End() //delete the vertex at the end of the list
    {
        if(isEmpty())  return false;
        nmnodes --;
        return true;
    }

    inline void copy_Elem(MCG_Node *s, MCG_Node *d)
    {
    }

    //delete the corresponding  vertex, if it succeeds, return true
    inline bool del_Node(MCG_Node *s)
    {
        if(isEmpty())  return false;

        //find the vertex, if find it, delete and move the following vertices forward
        //otherwise, return false;

        int i, pos = -1;

        for(i = 0; i < nmnodes; i++)
        {
            if(mnodes[i] == s)
            {
                pos = i;
                break;
            }
        }

        if(pos == -1) return false;

        //delete it
        for(i = pos; i < nmnodes-1; i++)
        {
            //we need a copy function
            copy_Elem(mnodes[i], mnodes[i+1]);
        }

        nmnodes--;

        return true;
    }

    inline bool isEmpty()  //judge whether the list is empty
    {
        if(nmnodes == 0)   return true;
        return false;
    }

    inline bool isFull()
    {
        if(nmnodes == curMaxNumMNodes) return true;
        return false;
    }

    //extend the original list, if it succeeds, return true
    inline bool extend(int step = 100)
    {
        MCG_Node **temp = mnodes;
        mnodes = (MCG_Node **) malloc(sizeof(MCG_Node *) * (curMaxNumMNodes + step));
        if( temp == nullptr)
        {
            //fail
            char rout[256], var[256];
            sprintf(rout, "%s", "MCG_NodeList::extend");
            sprintf(var, "%s", "mnodes");

            //write_mem_error(rout, var, 1);
            exit(-1);
            mnodes = temp;
            return false;
        }

        int i;
        for(i = 0; i < curMaxNumMNodes; i++)
            mnodes[i] = temp[i];

        for(i = curMaxNumMNodes; i < curMaxNumMNodes+step; i++)
            mnodes[i] = nullptr;
        curMaxNumMNodes += step;

        free(temp);
        return true;
    }

    inline void reset()
    {
        nmnodes = 0;
    }

}; //end of MCG_NodeList class

class MorseDecomp;
class MCG_Graph{
public:
    MCG_NodeList *nlist;
    Graph_EdgeList *elist;
    int r_counter, a_counter, s_counter;
    int cur_mcgnode_index, cur_mcgedge_index;
    int *boundary_nodes;
    int nboundarytris;

    ConleyIndex *conley;
    vector<int>NeedRefined;
    void AddQueue(int morse);
    bool NeedRefinement(int morse);

    //Graph related algorithms
    //Dijkstra
    //BFS
    //DFS
    //MST

    MCG_Graph(int nnodes = 50, int nedges = 50)
    {
        nlist = new MCG_NodeList(nnodes);
        elist = new Graph_EdgeList(nedges);
        conley= new ConleyIndex();

        boundary_nodes = nullptr;
        nboundarytris = 0;
    }

    ~MCG_Graph()
    {
        delete nlist;
        delete elist;
        delete conley;
    }

    void init_MCG()
    {
        int i;
        for(i = 0; i < elist->curMaxNumGedges; i++)
            elist->edges[i] = new Graph_Edge();
    }

    void reset_sccnodes_flags();

    /*member functions for MCG construction and modification*/
    /*build the mcg according to the obtained Morse sets*/
    void build_mcg();
    void assign_mcgnodes();
    void layout_mcg();
    void add_to_mcg_nodelist(int, int &);
    void add_to_mcg_edgelist(int node1, int node2);
    void add_edge_to_node(int node, int edgeindex);

    int get_region_type(int scc_index);
    int make_correction_saddle(int scc_index);

    int get_boundary_nodes(int scc_index);
    int get_triangle_type(int tri, int scc_index);

    void grow_repeller_region_graph(int scc_index, int inverse);
    bool has_longer_path_in_mcg(int n1, int n2, int edgeindex);

    bool is_connected_mcg(int node1, int node2);
    bool is_valid_link(int from, int to);

    void grow_all_mcgnodes();
    void remove_redundant_edges();

    /*
        new algorithm to find out the MCG edges
    */
    void grow_from_MorseSet(int MorseID, bool back_forth, MorseDecomp *);
    DynList_Int *grow_reg_from_to (int scc1, int scc2, bool back_forth, MorseDecomp *morse_decomp, DynList_Int *array_);
    DynList_Int *get_intersect_reg (DynList_Int *reg1, DynList_Int *reg2, DirGraph *dg, DynList_Int *intersect_reg);

    int *get_intersect_reg_int (int *reg1, int nreg1,
                               int *reg2, int nreg2, DirGraph *dg, int *intersect_reg, int &n_intersect_tris);
    void add_to_mcg_edge_list (int from, int to, DynList_Int *con_reg = nullptr);
    void grow_saddle_mcgnodes();
    void grow_other_mcgnodes();

    void remove_indir_mcgedges(); // we remove the indirecte connection MCG edges rather than the redundant edges
    bool has_dir_path_between_MorseSets (int from, int to, MorseDecomp *morse_decomp);
    stack <int> dir_path_search;  // a stack for direct path searching

    ///////////////////////////////////////////////////////////////////////////////////


    void build_mcg_edges_graph();
    bool IsMCGN(int scc);

    void cal_mcgedge_regions(void);

    void get_saddle_two_regions(int scc_index, 	int *forward_region, int &nforward_tris, int &curMaxNumForward, int *backward_region, int &nbackward_tris, 	int &curMaxNumBackward);

    void cal_forward_region_MFIG(int scc_index, int inverse,  int *repell_region, int &ntris_repell_region, int &curMaxNumTrisRepellRegion);

    //void get_saddle_two_regions(int scc_index, 	int *forward_region, int &nforward_tris, int &curMaxNumForward, int *backward_region, int &nbackward_tris, 	int &curMaxNumBackward);
    void intersect_twoRegions(int *source1, int num1, int *source2, int num2,  int *dest, int &num_dest);
    void remove_redundant(int *source1, int &num1, int *source2, int num2);

    void ReverseEdges();
    bool IsRepeated(int *a, int b, int num);

    //auto refine, select morse
    int select_morse(double min_pri);

    /*
       Build local mcg given a list of the triangles inside the region
       Since only the local mcg will call this function, the Morse sets in current MCG are already inside the local region
    */
    void grow_all_mcgnodes_local(MorseDecomp *local_decomp);
    void grow_repeller_region_graph_local(MorseDecomp *local_decomp, int scc_index, int inverse);
    void build_mcg_edges_graph_local(/*int *reg_tris, int ntris*/MorseDecomp *local_decomp);

    /**********************************************************************
      02/28/2010
      save the obtained global MCG
    */
    void save_MCG(char *filename);
    void load_MCG(char *filename);


}; //end of ECG_Graph class




/** The SCC component class **/
class SCComponent{
public:
    int *nodes;          //the list of triangle lists belonging to this component
    int nnodes;
    int node_index;      //the index of Morse sets
    bool valid;          //true梫alid SCC, need to detect limit cycle; false梚nvalid
    int nseppts;
    int nattpts;
    int nboundaries;
    int *singular_tri;    //the list of triangles containing singularity
    int nfixedpoints;
    int *periodicorbits;
    int nperiodicorbits;

    //conley index
    //added 12/03/2009 by Qingqing
    int XM;//Euler Characteristics of Mesh
    int XL;//Euler Characteristics of Boundary
    int Conley0;//Conley index 0
    int Conley1;//Conley index 1
    int Conley2;//Conley index 2
    int classification;//0=trival;1=source;2=sink;3=saddle
    double priority;//for auto refine
    double variance_vector;//vector variance

    //corresponding
    int global_SCC;//the SCC# of global scclist
    int local_SCC;//the SCC# of local scclist

    //function for finding a SCC component

    /*calculate all the separation and attachment points of a strongly connected component*/
    void cal_sep_attp_pts();

    void reset_edge_intersections();  /*reset the intersection information for all the edges
                                      of the strongly connected component*/

    bool is_connected();

};  // end of SCComponent class



class SCCList{
public:
    SCComponent **scccomponents;
    int nsccs;
    int curMaxNumSCCS;

    // The similar list operations
    SCCList(int initsize = 1000) //construction
    {
        scccomponents = (SCComponent **)malloc(sizeof(SCComponent *)*initsize);
        curMaxNumSCCS = initsize;
        nsccs = 0;

        if(scccomponents == nullptr)
        {
            char rout[256], var[256];
            sprintf(rout, "%s", "SCCList Constructor");
            sprintf(var, "%s", "scccomponents");

            //write_mem_error(rout, var, 0);
            curMaxNumSCCS = 0;
            exit(-1);
        }

        for(int i = 0; i < initsize; i++)
            scccomponents[i] = nullptr;
    }

    ~SCCList()
    {
        if(scccomponents != nullptr)
        {
            int i;
            for(i = 0; i < curMaxNumSCCS; i++)
            {
                if(scccomponents[i] != nullptr)
                {
                    free(scccomponents[i]->nodes);
                    free(scccomponents[i]);
                }
            }

            free(scccomponents);
        }
    }

    //add a new vertex to the end of the list, if it succeeds, return true
    inline bool append(SCComponent *s)
    {
        if(isFull ())
            if(!extend(100))
                return false;             //if not enough memory available, return false
        scccomponents[nsccs] = s;
        //copyElem(s, polist[nporbits]);
        nsccs++;
        return true;
    }

    inline bool del_End() //delete the vertex at the end of the list
    {
        if(isEmpty())  return false;
        nsccs --;
        return true;
    }

    /*not fully implemented the copy yet*/
    inline void copy_Elem(SCComponent *s, SCComponent *d)
    {
    }

    //delete the corresponding  vertex, if it succeeds, return true
    inline bool del_Node(SCComponent *s)
    {
        if(isEmpty())  return false;

        //find the vertex, if find it, delete and move the following vertices forward
        //otherwise, return false;

        int i, pos = -1;

        for(i = 0; i < nsccs; i++)
        {
            if(scccomponents[i] == s)
            {
                pos = i;
                break;
            }
        }

        if(pos == -1) return false;

        //delete it
        for(i = pos; i < nsccs-1; i++)
        {
            //we need a copy function
            copy_Elem(scccomponents[i], scccomponents[i+1]);
        }

        nsccs--;

        return true;
    }

    inline bool isEmpty()  //judge whether the list is empty
    {
        if(nsccs == 0)   return true;
        return false;
    }

    inline bool isFull()
    {
        if(nsccs == curMaxNumSCCS) return true;
        return false;
    }

    //extend the original list, if it succeeds, return true
    inline bool extend(int step = 100)
    {
        SCComponent **temp = scccomponents;
        scccomponents = (SCComponent **) malloc(sizeof(SCComponent *) * (curMaxNumSCCS + step));
        if( temp == nullptr)
        {
            //fail
            char rout[256], var[256];
            sprintf(rout, "%s", "SCCList::extend");
            sprintf(var, "%s", "sccomponents");

            //write_mem_error(rout, var, 1);
            curMaxNumSCCS = 0;
            exit(-1);

            scccomponents = temp;
            return false;
        }

        int i;
        for(i = 0; i < curMaxNumSCCS; i++)
            scccomponents[i] = temp[i];
        for(i = curMaxNumSCCS; i < curMaxNumSCCS+step; i++)
            scccomponents[i] = nullptr;
        curMaxNumSCCS += step;

        free(temp);
        return true;
    }

    inline void reset()
    {
        nsccs = 0;
    }
}; //end of SCCList class




/* Line segment data structure */
typedef struct LineSeg{
    int Triangle_ID;            //which triangle this line segment locates in
    icVector2 start, end;    //local coordinates for start and end points
    icVector3 gstart, gend;  //global coordinates for start and end points
    double length;              //we may need to store the length of the line segment
} LineSeg;


/*For temporary trajetory*/
typedef struct CurvePoints{
    int triangleid;
    double lpx, lpy;
    double gpx, gpy, gpz;
    double length;
}CurvePoints;


/** Trajectory **/
extern class Edge;
class Trajectory{
public:
    int index;
    int  nlinesegs;
    int  curMaxNumLinesegs;
    LineSeg *linesegs;

    int saddleID;          /*which saddle this trajectory belongs to*/

    double eulerstep_scalar;

    /*Construct the trajectory*/
    Trajectory(int index, int curMaxNum = 200);

    ~Trajectory()
    {
        if(curMaxNumLinesegs > 0)
            free(linesegs);
    }

    /*extend the line segment list if there is not enough space left*/
    bool extend_line_segments(int add_size);

    //get the length of the trajectory
    double get_length();

    //get any line segment according to the input index
    LineSeg *get_line_seg(int index);

    /*------------ routines for calculating streamlines -------------------------------*/
    //all the tracing calculation should be fulfilled here
    bool cal_next_point_euler1(double first[2], double second[2], int &face_id, double alpha[3], int type);
    bool cal_nextpt_2ndeuler(double first[2], double second[2], int &face_id, double alpha[3], int type);
    bool cal_nextpt_RK4(double first[2], double second[2], int &face_id, double alpha[3], int type);

    bool cal_nextpt_RK2_nonNorm(double first[2], double second[2], int &face_id, double alpha[3], int type);
    bool get_next_pt(double first[2], double second[2], int &face_id, double alpha[3], int type, unsigned char opt);

    bool store_to_global_line_segs(CurvePoints *temp, int num);
    void local_To_global(int faceid, double locpos[2], icVector3 &glpos);
    int trace_in_triangle(int &face_id, double globalp[3], int type, int &flag);
    void get_next_triangle(int &face_id, double pre[2], double cur[2], double param_t[2], int type,
                           int &PassVertornot, double alpha[3]);
    void cross_a_vertex(int &face_id, double cur_p[2], double pre_p[2],int type, int &passornot);
    void pass_edge(int &face_id, int which_edge);
    void pass_vertex(int vert_id, int &theone, int type);
    void pass_vertex_2(int vert_id, int &theone, int type);
    //void cross_boundary(double pre[2], double cur[2], int face_id, int &which_edge, double t[2]);
    void cross_boundary(double pre[2], double cur[2],
                        int face_id, double alpha[3],
                        int &which_edge, double t[2]);
    int get_intersection(double PointA[2], double PointB[2], double PointC[2], double PointD[2], double t[2]);
    void cal_one_traj(int face_id, double x, double y, double z, int type);

    /**------------- routines for the tracing of periodic orbit detection -------------**/
    bool trace_to_find_recurrence(double g[3], int &triangle, int type,
                                  int scc_index, Edge *chosen_edge,  int &flag);
    //int trace_in_triangle_po_detect(double g[3], int &face_id, int type,
    //						 Edge *cur_e, icVector3 &intersect, double &seg_length, int &flag);

    /*---------------------------------------------------------------------------------*/
    /**----------- routines for the tracing of finding the connection ------------**/
    int trace_in_triangle_connection(double g[3], int &face_id, int type, int &flag);

    /*---------------------------------------------------------------------------------*/
    /**----------- routines for the tracing of computing the multi-valued map ------------**/
    int trace_in_triangle_tau(int &face_id, double globalp[3], int type,
                              double tau, double &gt_tau, int &flag);
    bool cal_next_point_euler1_tau(double first[2], double second[2], int &face_id, double alpha[3], int type);
    bool cal_next_point_euler1_tau_nonNorm(double first[2], double second[2], int &face_id, double alpha[3], int type);
    double cur_vec_mag;   //the magnitude of the integral vector (the speed of the flow)
    double move_dist;

    bool cal_next_point_2ndEuler_tau_nonNorm(double first[2], double second[2], int &face_id,
                                             double alpha[3], int type);
    bool cal_next_point_2ndEuler_tau_Norm(double first[2], double second[2], int &face_id,
                                          double alpha[3], int type);


    bool cal_next_pt_RK4_tau(double first[2], double second[2], int &face_id,
                             double alpha[3], int type);

    bool cal_next_pt_RK4_tau_f(double first[2], double second[2], int &face_id, double alpha[3]);
    bool cal_next_pt_RK4_tau_b(double first[2], double second[2], int &face_id, double alpha[3]);

    bool get_next_pt_tau(double first[2], double second[2], int &face_id, double alpha[3], int type, unsigned char opt);
    bool get_next_pt_tau_f(double first[2], double second[2], int &face_id, double alpha[3], unsigned char opt);
    bool get_next_pt_tau_b(double first[2], double second[2], int &face_id, double alpha[3], unsigned char opt);
    bool cal_next_point_2ndEuler_tau_Norm_f(double first[2], double second[2], int &face_id,
                                            double alpha[3]);
    bool cal_next_point_2ndEuler_tau_Norm_b(double first[2], double second[2], int &face_id,
                                            double alpha[3]);


    bool cal_next_point_euler1_tau_nonNorm_f(double first[2], double second[2], int &face_id, double alpha[3]);
    bool cal_next_point_euler1_tau_nonNorm_b(double first[2], double second[2], int &face_id, double alpha[3]);

    bool cal_next_point_euler1_tau_f(double first[2], double second[2], int &face_id,
                                     double alpha[3]);
    bool cal_next_point_euler1_tau_b(double first[2], double second[2], int &face_id,
                                     double alpha[3]);

    int trace_in_triangle_tau_f(int &face_id, double globalp[3],
                                double tau, double &gt_tau, int &flag);
    int trace_in_triangle_tau_b(int &face_id, double globalp[3],
                                double tau, double &gt_tau, int &flag);

    int trace_in_triangle_tau_b_rot_sum(int &face_id, double globalp[3],
                                        double tau, double &gt_tau, float &rot_sum, int &flag);
    int trace_in_triangle_tau_f_rot_sum(int &face_id, double globalp[3],
                                        double tau, double &gt_tau, float &rot_sum, int &flag);
    int trace_in_triangle_tau_f_smooth(int &face_id, double globalp[3],
                                       double tau, double &gt_tau, float &rot_sum, int &flag);
    int trace_in_triangle_tau_b_smooth(int &face_id, double globalp[3],
                                       double tau, double &gt_tau, float &rot_sum, int &flag);

    bool get_next_pt_rotsum_f(double first[2], double second[2], int &face_id, double alpha[3], unsigned char opt);
    bool get_next_pt_rotsum_b(double first[2], double second[2], int &face_id, double alpha[3], unsigned char opt);
    bool cal_next_point_2ndEuler_Norm_f(double first[2], double second[2], int &face_id,
                                        double alpha[3]);
    bool cal_next_point_2ndEuler_Norm_b(double first[2], double second[2], int &face_id, double alpha[3]);
    bool cal_next_point_euler1_f(double first[2], double second[2], int &face_id,
                                 double alpha[3]);
    bool cal_next_point_euler1_b(double first[2], double second[2], int &face_id,
                                 double alpha[3]);
}; //end of Trajectory class




class TrajectoryList{
public:
    Trajectory **trajs;         //the trajectory list
    int ntrajs;                  //current number of existing trajectories
    int curMaxNumTrajs;          //maximum number of trajectories can be stored
    double length;                //the flow length of the trajectory

    // The similar list operations
    TrajectoryList(int initsize = 1000) //construction
    {
        trajs = (Trajectory **)malloc(sizeof(Trajectory *)*initsize);
        curMaxNumTrajs = initsize;
        ntrajs = 0;

        if(trajs == nullptr)
        {
            char rout[256], var[256];
            sprintf(rout, "%s", "TrajectoryList Constructor");
            sprintf(var, "%s", "trajs");

            //write_mem_error(rout, var, 0);
            curMaxNumTrajs = 0;
            exit(-1);
        }

        for(int i = 0; i < initsize; i++)
            trajs[i] = nullptr;
        curMaxNumTrajs = initsize;
    }

    ~TrajectoryList()
    {
        int i, j;

        for(i = 0; i < curMaxNumTrajs; i++)
        {
            if(trajs[i] != nullptr)
            {
                free(trajs[i]->linesegs);
            }
        }

        free(trajs);
    }

    //add a new vertex to the end of the list, if it succeeds, return true
    inline bool append(Trajectory *s)
    {
        if(isFull ())
            if(!extend(100))
                return false;             //if not enough memory available, return false
        trajs[ntrajs] = s;
        //copyElem(s, polist[nporbits]);
        ntrajs++;
        return true;
    }

    inline bool del_End() //delete the vertex at the end of the list
    {
        if(isEmpty())  return false;
        ntrajs --;
        return true;
    }

    inline void copy_Elem(Trajectory *s, Trajectory *d)
    {
    }

    //delete the corresponding  vertex, if it succeeds, return true
    inline bool del_Node(Trajectory *s)
    {
        if(isEmpty())  return false;

        //find the vertex, if find it, delete and move the following vertices forward
        //otherwise, return false;

        int i, pos = -1;

        for(i = 0; i < ntrajs; i++)
        {
            if(trajs[i] == s)
            {
                pos = i;
                break;
            }
        }

        if(pos == -1) return false;

        //delete it
        for(i = pos; i < ntrajs-1; i++)
        {
            //we need a copy function
            copy_Elem(trajs[i], trajs[i+1]);
        }

        ntrajs--;

        return true;
    }

    inline bool isEmpty()  //judge whether the list is empty
    {
        if(ntrajs == 0)   return true;
        return false;
    }

    inline bool isFull()
    {
        if(ntrajs == curMaxNumTrajs) return true;
        return false;
    }

    //extend the original list, if it succeeds, return true
    inline bool extend(int step = 100)
    {
        Trajectory **temp = trajs;
        trajs = (Trajectory **) malloc(sizeof(Trajectory *) * (curMaxNumTrajs + step));
        if( trajs == nullptr)
        {
            //fail
            char rout[256], var[256];
            sprintf(rout, "%s", "TrajectoryList::extend");
            sprintf(var, "%s", "trajs");

            //write_mem_error(rout, var, 1);
            curMaxNumTrajs = 0;
            trajs = temp;
            exit(-1);

            return false;
        }

        int i;
        for(i = 0; i < curMaxNumTrajs; i++)
            trajs[i] = temp[i];
        for(i = curMaxNumTrajs; i < curMaxNumTrajs+step; i++)
            trajs[i] = nullptr;

        curMaxNumTrajs += step;

        free(temp);
        return true;
    }

    inline void reset()
    {
        ntrajs = 0;
    }

    /*we now put the separatrix calculation here */
    void cal_startpt_sep(int triangleID, icVector3 sep_vector, double saddle_ce[3], double newpos[3]);
    void cal_separatrices();

}; //end of TrajectoryList class



/** Seed point structure **/
typedef struct Seed{
    double pos[3];                  //the coordinates of the seed point
    int triangle;                      //the triangle contains the seed point
    unsigned char state;
}Seed; // end of Seed structure

class SeedList{
public:
    Seed **seeds;
    int nseeds;
    int curMaxNumSeeds;

    // The similar list operations
    SeedList(int initsize = 3000) //construction
    {
        seeds = (Seed **)malloc(sizeof(Seed *)*initsize);
        curMaxNumSeeds = initsize;
        nseeds = 0;

        if(seeds == nullptr)
        {
            char rout[256], var[256];
            sprintf(rout, "%s", "SeedList Constructor");
            sprintf(var, "%s", "seeds");

            //write_mem_error(rout, var, 0);
            curMaxNumSeeds = 0;
            exit(-1);

        }

        for(int i = 0; i < initsize; i++)
            seeds[i] = nullptr;
    }

    ~SeedList()
    {
        for(int i = 0; i < curMaxNumSeeds; i++)
            if(seeds[i] != nullptr)
                free(seeds[i]);
        free(seeds);
    }

    //add a new vertex to the end of the list, if it succeeds, return true
    inline bool append(Seed *s)
    {
        if(isFull ())
            if(!extend(100))
                return false;             //if not enough memory available, return false
        seeds[nseeds] = s;
        //copyElem(s, polist[nporbits]);
        nseeds++;
        return true;
    }

    inline bool del_End() //delete the vertex at the end of the list
    {
        if(isEmpty())  return false;
        nseeds --;
        return true;
    }

    inline void copy_Elem(Seed *s, Seed *d)
    {
    }

    //delete the corresponding  vertex, if it succeeds, return true
    inline bool del_Node(Seed *s)
    {
        if(isEmpty())  return false;

        //find the vertex, if find it, delete and move the following vertices forward
        //otherwise, return false;

        int i, pos = -1;

        for(i = 0; i < nseeds; i++)
        {
            if(seeds[i] == s)
            {
                pos = i;
                break;
            }
        }

        if(pos == -1) return false;

        //delete it
        for(i = pos; i < nseeds-1; i++)
        {
            //we need a copy function
            copy_Elem(seeds[i], seeds[i+1]);
        }

        nseeds--;

        return true;
    }

    inline bool isEmpty()  //judge whether the list is empty
    {
        if(nseeds == 0)   return true;
        return false;
    }

    inline bool isFull()
    {
        if(nseeds == curMaxNumSeeds) return true;
        return false;
    }

    //extend the original list, if it succeeds, return true
    inline bool extend(int step = 100)
    {
        Seed **temp = seeds;
        seeds = (Seed **) malloc(sizeof(Seed *) * (curMaxNumSeeds + step));
        if( temp == nullptr)
        {
            //fail
            char rout[256], var[256];
            sprintf(rout, "%s", "SeedList::extend");
            sprintf(var, "%s", "seeds");

            //write_mem_error(rout, var, 1);
            curMaxNumSeeds = 0;
            seeds = temp;
            //MessageBox(nullptr, "failed to extend memory for seeds!", "Error", MB_OK);
            qDebug() << "failed to extend memory for seeds!";
            exit(-1);
            return false;
        }

        int i;
        for(i = 0; i < curMaxNumSeeds; i++)
            seeds[i] = temp[i];
        for(i = curMaxNumSeeds; i < curMaxNumSeeds+step; i++)
            seeds[i] = nullptr;
        curMaxNumSeeds += step;

        free(temp);
        return true;
    }

    inline void reset()
    {
        nseeds = 0;
    }
}; //end of SeedList class





/** Sample point structure **/
typedef struct SamplePt{
    double gpt[3];
    int triangle;
    int traj;                                //which trajectory this sample falls on
}SamplePt; //end of SamplePt class





class SamplePtList{
public:
    SamplePt **samples;
    int nsamples;
    int curMaxNumSamplePts;

    SamplePtList(int initsize = 1000) //construction
    {
        samples = (SamplePt **)malloc(sizeof(SamplePt *)*initsize);
        curMaxNumSamplePts = initsize;
        nsamples = 0;

        if(samples == nullptr)
        {
            char rout[256], var[256];
            sprintf(rout, "%s", "SamplePtList Constructor");
            sprintf(var, "%s", "samples");

            //write_mem_error(rout, var, 0);
            curMaxNumSamplePts = 0;
            exit(-1);
        }

        for(int i = 0; i < initsize; i++)
            samples[i] = nullptr;
    }

    ~SamplePtList()
    {
        if(samples!= nullptr)
        {
            for(int i = 0; i < curMaxNumSamplePts; i++)
            {
                if(samples[i] != nullptr)
                    free(samples[i]);
            }
            free(samples);
        }
    }

    //add a new vertex to the end of the list, if it succeeds, return true
    inline bool append(SamplePt *s)
    {
        if(isFull ())
            if(!extend(100))
                return false;             //if not enough memory available, return false
        samples[nsamples] = s;
        //copyElem(s, polist[nporbits]);
        nsamples++;
        return true;
    }

    inline bool del_End() //delete the vertex at the end of the list
    {
        if(isEmpty())  return false;
        nsamples --;
        return true;
    }

    inline void copy_Elem(SamplePt *s, SamplePt *d)
    {
    }

    //delete the corresponding  vertex, if it succeeds, return true
    inline bool del_Node(SamplePt *s)
    {
        if(isEmpty())  return false;

        //find the vertex, if find it, delete and move the following vertices forward
        //otherwise, return false;

        int i, pos = -1;

        for(i = 0; i < nsamples; i++)
        {
            if(samples[i] == s)
            {
                pos = i;
                break;
            }
        }

        if(pos == -1) return false;

        //delete it
        for(i = pos; i < nsamples-1; i++)
        {
            //we need a copy function
            copy_Elem(samples[i], samples[i+1]);
        }

        nsamples--;

        return true;
    }

    inline bool isEmpty()  //judge whether the list is empty
    {
        if(nsamples == 0)   return true;
        return false;
    }

    inline bool isFull()
    {
        if(nsamples == curMaxNumSamplePts) return true;
        return false;
    }

    //extend the original list, if it succeeds, return true
    inline bool extend(int step = 100)
    {
        SamplePt **temp = samples;
        samples = (SamplePt **) malloc(sizeof(SamplePt *) * (curMaxNumSamplePts + step));
        if( temp == nullptr)
        {
            //fail
            char rout[256], var[256];
            sprintf(rout, "%s", "SamplePtList::extend");
            sprintf(var, "%s", "samples");

            //write_mem_error(rout, var, 1);
            curMaxNumSamplePts = 0;
            samples = temp;
            exit(-1);
            return false;
        }

        int i;
        for(i = 0; i < curMaxNumSamplePts; i++)
            samples[i] = temp[i];
        for(i = curMaxNumSamplePts; i < curMaxNumSamplePts+step; i++)
            samples[i] = nullptr;
        curMaxNumSamplePts += step;

        free(temp);
        return true;
    }

    inline void reset()
    {
        nsamples = 0;
    }
}; //end of SamplePtList class


/** Data structure for vector field design **/
//Use for editing singular element

typedef struct EditBox{
    icVector2 p1, p2, p3, p4, Up;  //stands for the 4 points of the edit box
}EditBox;



//** Singular element **
class SingularElem {
public:
    //member variables for editing
    int           index;                               //the index of the singular element in the list
    int           TriangleID;                       //the triangle containing the singular element
    double        rotang;                             //rotation angle around the singular element
    double        sx, sy;                              //asymmetric scales
    double        s;                                      //uniform scale
    icVector3     gpos;                                 //global coordinates of singular element
    icMatrix3x3   transform_matrix;            //transform matrix for both editbox and Jacobian
    EditBox       editbox, cur_editbox;      //the edit box
    unsigned char type;

    //we should have some functions to access(get or set) those values
    //set the transformation
    //set the editbox
    //apply transformation

}; //end of SingularElem class






//** Regular element **
class RegularElem{
    //member variables
public:
    int            index;                             //index of the element in the list
    int            base_triangle;                //the triangle containing the base point of the element
    icVector3      base;                              //the base position of the element
    icVector3      Dir;                                //the direction of the element
    icMatrix3x3    transform_matrix;         //the transformation matrix of the element
    double         rotang;                           //rotation angle around the base
    double         s;                                   //uniform scale only
    unsigned char  type;                              //the type of the element

    //set the basis
    //set the direction
    //set the tranformation matrix

}; //end of RegularElem class




/** Control point list **/

class ControlPtList
{
    Point3Dd *ctrlpoints;
    int num_pts;
    int MaxNumControlPts;

    // The similar list operations as VertexList
}; //end of ControlPtList class


//---------------------------------------------------------------------------------------

/* Data structure for vector field editing */

/*
The region consists of the indices of the triangles inside it
*/
class Region_Tri_Index{
public:
    int index;
    int *tris;
    int ntris;
    int curMaxNumTris;

    Region_Tri_Index(int index = 0, int initsize = 1000)
    {
        tris = new int[initsize];
        curMaxNumTris = initsize;
        ntris = 0;

        this->index = index;
    }

    Region_Tri_Index(int *i_tris, int i_ntris, int index = 0)
    {
        tris = new int[i_ntris];

        int i;
        for(i = 0; i < i_ntris; i++)
        {
            tris[i] = i_tris[i];
        }
        curMaxNumTris = ntris = i_ntris;
    }

    ~Region_Tri_Index()
    {
        delete [] tris;
    }

    int cal_euler_value();  /*calculate the Euler characteristics value*/

    bool contain_fixedpts(); /*judge whether this region contains fixed points or not*/

    int count_num_fixedpts(); /*count the number of fixed points in this region*/


    ////Routines for region validation
    bool in_edge_list(Edge **edges, int num, Edge *anedge);
    bool in_vert_list(int *verts, int num, int vert_index);
};


//** Region List **
class RegionList{
    TriangleList* regionlist = nullptr;
    int num_regions;
    int MaxNumRegions;
    // The similar list operations as VertexList
}; // end of RegionList class


//** Boundary **
class Boundary{
public:
   EdgeList* edges = nullptr;

    //CornerList corners;
    int index;
}; //end of Boundary class



class BoundaryList{
public:
    Boundary *boundarylist = nullptr;
    int num_boundarys;
    int MaxNumBoundaries;

    // The similar list operations (it is one pointer list)
}; //end of BoundaryList class


//** Repeller/Attractor List **
class NodeList{
    int *nodes = nullptr;                      //the indices of the nodes in the ECG/MCG graph
    int num_nodes;
    int MaxNumNodes;

    // The similar list operations (it is one pointer list)
}; //end of NodeList class



#endif // VFIELD_H
