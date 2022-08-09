#pragma once
#ifndef GRAPH2D_H
#define GRAPH2D_H

#include <vector>
using namespace std;

/*
    The data structure of the contour
*/

typedef struct ContourVertex
{
    int index;

    icVector3 pos;

    std::vector<int> edges;
    bool visited;
    bool stop;

    int at_edge;  // the corresponding mesh vertex id
}ContourVertex;

typedef struct ContourEdge
{
    int verts[2];
    int index;
    float length;
    bool visited;
    int at_triangle;
}ContourEdge;


class Contour_Graph{
public:

    std::vector<ContourVertex> vlist;
    std::vector<ContourEdge> elist;
    double iso_val;

    // Later, we may need to get the different connected components


    Contour_Graph(int init_vert = 0, int init_edge = 0)
    {
    }

    ~Contour_Graph()
    {
        vlist.clear();
        elist.clear();
    }

    void add_one_node(ContourVertex &node)
    {
        node.index = vlist.size();
        node.edges.clear();
        vlist.push_back(node);
    }

    void add_one_edge(ContourEdge &e)
    {
        e.index=elist.size();
        elist.push_back(e);

        // remember to udpate the edge lists of the two nodes connected by e
        vlist[e.verts[0]].edges.push_back(e.index);
        vlist[e.verts[1]].edges.push_back(e.index);
    }

    void extract_iso_contour_from_distance_field(double iso_val);

    void display();

};


#endif // GRAPH2D_H
