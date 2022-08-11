/*
For geometry data structures

Created and Modified by Guoning Chen
copyright @2007
*/
#pragma once
#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include <math.h>
#include <conio.h>

#include "ExternalDependencies/icVector.h"
#include "ExternalDependencies/icMatrix.h"
#include "Point.h"
#include "Predefined.h"
#include "FileLoader/PlyLoader.h"
#include "FileLoader/vfply_io.h"
#include "VField.h"

extern char *elem_names[];

extern PlyProperty vert_props[];

extern PlyProperty face_props[];
/* forward declarations */

class Edge;
class Corner;
class Triangle;
class Trajectory;
/* the vertex structure for 3D surface vector field */
class Vertex{
public:
    bool angle_deficit;         //1—has angle deficit; 0—no
    bool visited;    //use to calculate the geodesic distance (3D)
    bool OnBoundary;         //flag to mark whether the vertex is a boundary vertex

    unsigned char nedges;
    unsigned char ncorners;
    unsigned char ntris;
    int index;                  //the index of the vertex
    int max_tris;
    int imgtri;            //for tau-map guided Morse decomposition  08/30/2007

    Edge **edges;             //the edges incident to the vertex
    Triangle **tris;           //the triangles sharing the vertex
    Corner **corners;       //the corners sharing the vertex
    /* other information needs to be store */

    float mag_rgb[3];
    double x, y, z;             //the coordinates of the vertex
    double t_angle;           //the angle of the vector value in tangential local frame
    double distance;          //Geodesic distance (for 3D surface)
    icVector3 g_vec;
    icVector3 t_vec;          //the vector tangential to the tangent plane of the vertex (3D )
                                //for 2D plane, we use icVector2
    icVector3 T, B, normal;       //the local coordinate frame defined according to the tangent plane (for 3D only)
    icVector2 texture_coord;  //texture coordinates for flow-like texture mapping
    icVector2 back_Tex;

    double img_tau[3];          //for tau-map guided Morse decomposition  08/30/2007

    bool in_out_let;            // mark the vertex if it is in the part of the inlet and outlet
    icVector3 raw_vec;

    /*
       We still need to store its original vector magnitude for extracting the iso-contours 06/13/2011
    */
    float vf_mag;

    float total_rot_sum;  // to record the total rotation of an integral curve starting from this vertex
    float ave_sum;
    /*
       The sample associated with it  07/06/2010
    */
    tauPt3 tauPt_f, tauPt_b;

    void *other_props;

    Vertex(double, double, double, double, double, double);
    Vertex(double ix, double iy, double iz,
               double vx, double vy, double vz, int i);
    ~Vertex();
};   //end of class Vertex


class VertexList{
public:
    Vertex **verts;
    int nverts;
    int curMaxNumVerts;

    VertexList(int initsize = 1000) //construction
    {
        verts = (Vertex **)malloc(sizeof(Vertex *)*initsize);
        curMaxNumVerts = initsize;
        nverts = 0;
    }

    ~VertexList()
    {
        int i;

        for(i = 0; i < nverts; i++)
        {
            free(verts[i]);
        }
        free(verts);
    }

    //add a new vertex to the end of the list, if it succeeds, return true
    inline bool append(Vertex *v)
    {
        if(isFull ())
            if(!extend(100))
                    return false;             //if not enough memory available, return false
        verts[nverts] = v;
        nverts++;
        return true;
    }

    inline bool del_End() //delete the vertex at the end of the list
    {
        if(isEmpty())  return false;
        nverts --;
        return true;
    }

    inline void copy_Elem(Vertex *s, Vertex *d)
    {
        d->x = s->x;
        d->y = s->y;
        d->z = s->z;

        d->angle_deficit = s->angle_deficit;
        d->B = s->B;
        d->T = s->T;
        d->normal = s->normal;

        d->g_vec = s->g_vec;

        d->t_vec = s->t_vec;
        d->texture_coord = s->texture_coord;

        d->nedges = s->nedges;

        int i;

        for(i = 0; i < d->nedges; i++)
            d->edges[i] = s->edges[i];

        d->ncorners = s->ncorners;

        for(i = 0; i < d->ncorners; i++)
            d->corners[i] = s->corners[i];

        d->ntris = s->ntris;

        for(i = 0; i < d->ntris; i++)
            d->tris[i] = s->tris[i];

        d->distance = s->distance;
        d->OnBoundary = s->OnBoundary;
        d->index = s->index;
        d->visited = s->visited;

    }

    //delete the corresponding  vertex, if it succeeds, return true
    inline bool del_Node(Vertex *v)
    {
        if(isEmpty())  return false;

        //find the vertex, if find it, delete and move the following vertices forward
        //otherwise, return false;

        int i, pos = -1;

        for(i = 0; i < nverts; i++)
        {
            if(verts[i] == v)
            {
                pos = i;
                break;
            }

        }

        if(pos == -1) return false;

        //delete it
        for(i = pos; i < nverts-1; i++)
        {
            //we need a copy function
            copy_Elem(verts[i], verts[i+1]);
        }

        nverts--;

        return true;
    }

    inline bool isEmpty()  //judge whether the list is empty
        {
        if(nverts == 0)   return true;
        return false;
    }

    inline bool isFull()
    {
        if(nverts == curMaxNumVerts) return true;
        return false;
    }

    //extend the original list, if it succeeds, return true
    inline bool extend(int step = 100)
    {
        Vertex **temp = verts;
        verts = (Vertex **) malloc(sizeof(Vertex *) * (curMaxNumVerts + step));
        if( verts == nullptr)
        {
            //fail
            verts = temp;
            return false;
        }

        int i;
        for(i = 0; i < curMaxNumVerts; i++)
            verts[i] = temp[i];
        for(i = curMaxNumVerts; i < curMaxNumVerts+step; i++)
            verts[i] = (Vertex*)malloc(sizeof(Vertex));
        curMaxNumVerts += step;

        free(temp);
        return true;
    }

    inline void reset()
    {
        nverts = 0;
    }

    /* the other important member functions for vertex list */

    void update_Texcoord() ;   //update the texture coordinates for IBFV(s)

    void get_LocalFrame();        //get the T, N, B for each vertex (for 3D only)

    void get_Curvature();  //calculate the curvature for each vertex (optional)

    void sort_Corner();  //sort the corners share the vertex

    void allocate_AngsForCorners() ;  //allocate angle for the corners incident to the vertex

    void cal_Geodesic(Vertex *);        //calculate the geodesic distance from source vertex

}; //end of VertexList class


typedef struct SampleListInTriangle{
    int which_traj;
    int which_sample;
}SampleListInTriangle;


/** Triangle for 3D surface **/
class Triangle{
public:
    bool orient;
    bool visited;                       //for SCC searching
    //bool inDesignCell;
    unsigned char num_dir_vecs; //number of local vectors
    unsigned char nverts;       //number of vertices
    int index;                         //index of the triangle
    /* variables for analysis */
    int singularityID;            //the singularity id inside the triangle if exists
    int num_samplepts;

    Vertex *verts[3];             //the pointers to the three vertices
    Edge *edges[3];              //the pointers to the three edges
    icVector3 LX, LY;           //the local frame defined in the triangle
    icVector3 normal;           //the normal of the triangle
    icVector2 *dir_vec;         //the local vector values at the three vertices
    double x1, x2, y2;          //the local coordinates of the three vertices under Lx, Ly
    double area;
    icMatrix2x2 Jacobian;    //the Jacobian of the field inside the triangle
    int local_index;	//index for local refinement morse set 12/03/2009 by Qingqing
    bool cur_morse;// in currently refining morse set

    bool in_img;

    bool has_zero_vec;         // added by guoning 03/09/2010

    double used_tau;         // the tau value that is used to track the image of this triangle

    bool is_part_of_MCG;     // true -- if the triangle is somehow included in a Morse set or connection region, false otherwise
    bool diff_ECG_MCG;

    unsigned char counter;

    double euler_step;

    //
    unsigned char MS_type;   // set the label of the triangle based on the type of the Morse set that it belongs to
                             // 0 - non-MS; 1 - source-like; 2 - sink-like; 3 - saddle-like

    /*
       The sample associated with it 07/06/2010
    */
    tauPt3 tauPt_f, tauPt_b;

    SampleListInTriangle *samplepts; /*for even placement of streamlines*/

    /* variables for design*/

    void *other_props;

    /*Qingqing Add July 20,2009*/
    bool exclude;
    /*Qingqing Add July 20,2009*/


    /* The functions used in Triangle class */
    inline void Reset()
    {
        //inDesignCell = false;
        visited = false;
        Jacobian.set(0.);

        for(int i = 0; i < num_dir_vecs; i++)
            dir_vec[i].set(0, 0);
    }

    inline void get_center(double c[3])
    {
        c[0] = (verts[0]->x+verts[1]->x+verts[2]->x)/3.;
        c[1] = (verts[0]->y+verts[1]->y+verts[2]->y)/3.;
        c[2] = (verts[0]->z+verts[1]->z+verts[2]->z)/3.;
    }

    inline void cal_a_close_pt(double old[3], double out[3])
    {
        icVector3 VP;

        /* we first try the center of the triangle */
        out[0] = 0.333*this->verts[0]->x + 0.333*this->verts[1]->x + 0.334*this->verts[2]->x;
        out[1] = 0.333*this->verts[0]->y + 0.333*this->verts[1]->y + 0.334*this->verts[2]->y;
        out[2] = 0.333*this->verts[0]->z + 0.333*this->verts[1]->z + 0.334*this->verts[2]->z;

        VP.entry[0] = out[0] - old[0];
        VP.entry[1] = out[1] - old[1];
        VP.entry[2] = out[2] - old[2];

        ////if it is too close to the input point, we move further to close to one of the edge
        if(length(VP) < 1e-6)
        {
            out[0] = 0.01*verts[0]->x + 0.495*verts[1]->x + 0.495*verts[2]->x;
            out[1] = 0.01*verts[0]->y + 0.495*verts[1]->y + 0.495*verts[2]->y;
            out[2] = 0.01*verts[0]->z + 0.495*verts[1]->z + 0.495*verts[2]->z;
        }
    }

    inline void get_center(icVector3 &c)
    {
        c.entry[0] = (verts[0]->x+verts[1]->x+verts[2]->x)/3.;
        c.entry[1] = (verts[0]->y+verts[1]->y+verts[2]->y)/3.;
        c.entry[2] = (verts[0]->z+verts[1]->z+verts[2]->z)/3.;
    }

    void cal_Jacobian();     //calculate Jacobian for teh triangle

    /*Given the input triangle neighborhood with the coordinates and vector values of the three vertices
    calculate the Jacobian inside it*/
    void cal_Jacobian2x2_neighborhood(icVector2 v[3], icVector2 vec[3], icMatrix2x2 &mat);

    void cal_localvec_with_angledeficit(int v1, int v2, int deficitvert, icVector2 &local_vec);

    void cal_localvec_of_vertex_in_tri(int vertID1, int vertID2, icVector2 &local_vec);

    void compute_SingularityPos();  //get the position of the singularity

    int get_SingularityType();  //get the type of the singularity (can be moved to singularity class)

    void get_EigenVector();  //calculate the eigen vectors of the Jacobian if exists

    void cal_LocalFrame();   //calculate Lx, Ly;

    void cal_LocalXY();       //calculate locate coordinates for the vertices

    void get_DirectVecs();    //get the vector values defined by the local frame of the triangle

    void getVectorAtPoint(double a, double b);
    //get the vector values inside a triangle given the locate coordinates of the point
    //it may need other functions to get the correct vector at particular point

    void get_Barycentric(double a, double b);

    void get_Barycentric(double x, double y, double z);

    ///////////////////////////////////////////
    bool fall_in_the_tri(double pos[3]);

    bool contain_fixedpt();  /*judge whether this triangle contains a fixed point or not*/

    /*
    Extend the sample point list inside one triangle, this could be a bottleneck!!! 04/30/07
    */
    SampleListInTriangle *extend_sampleList()
    {
        SampleListInTriangle *temp = samplepts;
        samplepts = (SampleListInTriangle *)malloc(sizeof(SampleListInTriangle)*(num_samplepts + 1));
        if(num_samplepts > 0)
        {
            for(int i = 0; i < num_samplepts; i++)
            {
                samplepts[i].which_sample = temp[i].which_sample;
                samplepts[i].which_traj = temp[i].which_traj;
            }

            free (temp);
        }
        return samplepts;
    }

    void reset_sampleList()
    {
        if(samplepts != nullptr)
            free(samplepts);
        samplepts = nullptr;
        num_samplepts = 0;
    }

    friend class Trajectory;

}; // end of Triangle class


class TriangleList{
public:
    Triangle **tris;
    int ntris;
    int curMaxNumTris;

    //Triangle **singulartris;
    int nsingulartris;
    //int curMaxNumSingularTris;

    // The similar list operations as VertexList

    TriangleList(int initsize = 1000) //construction
    {
        tris = (Triangle **)malloc(sizeof(Triangle *)*initsize);
        curMaxNumTris = initsize;
        ntris = 0;

        //singulartris = nullptr;
        nsingulartris = 0;
        //curMaxNumSingularTris = 0;
    }

    ~TriangleList()
    {
        int i;

        for(i = 0; i < ntris; i++)
        {
            free(tris[i]);
        }
        free(tris);
    }

    //add a new vertex to the end of the list, if it succeeds, return true
    inline bool append(Triangle *t)
    {
        if(isFull ())
            if(!extend(100))
                    return false;             //if not enough memory available, return false
        tris[ntris] = t;
        ntris++;
        return true;
    }

    inline bool del_End() //delete the vertex at the end of the list
    {
        if(isEmpty())  return false;
        ntris --;
        return true;
    }

    inline void copy_Elem(Triangle *s, Triangle *d)
    {
    }

    //delete the corresponding  vertex, if it succeeds, return true
    inline bool del_Node(Triangle *t)
    {
        if(isEmpty())  return false;

        //find the vertex, if find it, delete and move the following vertices forward
        //otherwise, return false;

        int i, pos = -1;

        for(i = 0; i < ntris; i++)
        {
            if(tris[i] == t)
            {
                pos = i;
                break;
            }

        }

        if(pos == -1) return false;

        //delete it
        for(i = pos; i < ntris-1; i++)
        {
            //we need a copy function
            copy_Elem(tris[i], tris[i+1]);
        }

        ntris--;

        return true;
    }

    inline bool isEmpty()  //judge whether the list is empty
    {
        if(ntris == 0)   return true;
        return false;
    }

    inline bool isFull()
    {
        if(ntris == curMaxNumTris) return true;
        return false;
    }

    //extend the original list, if it succeeds, return true
    inline bool extend(int step = 100)
    {
        Triangle **temp = tris;
        tris = (Triangle **) malloc(sizeof(Triangle *) * (curMaxNumTris + step));
        if( tris == nullptr)
        {
            //fail
            tris = temp;
            return false;
        }

        int i;
        for(i = 0; i < curMaxNumTris; i++)
            tris[i] = temp[i];
        for(i = curMaxNumTris; i < curMaxNumTris+step; i++)
            tris[i] = (Triangle*)malloc(sizeof(Triangle));
        curMaxNumTris += step;
        free(temp);
        return true;
    }

    inline void reset()
    {
        ntris = 0;
    }


    /* the other important member functions for triangle list */

    void capture_Singularities(); //capture all the singularities in the mesh

    friend class Singularity;
    friend class SingularityList;


}; //end of TriangleList class


/** Edge **/
class Edge{
public:

    /*Member variables*/
    bool visited;
    bool OnBoundary;
    bool valid;
    bool find_attp, find_sep;
    unsigned char ntris;
    unsigned char att_visit, sep_visit;    //0 -- alive; 1--visited/dead;  2--too close/pending  07/23/06
    unsigned char num_intersections;	/**----- For recording the intersections -----**/
    int index;                       //the index of the edge
    int triangle;             //for which triangle whose Jacobian is calculated to approximate the Jacobian of the edge
    Vertex *verts[2];           //the two ending points of the edge
    Triangle **tris;           //the two triangles sharing the edge
    Point3Dd *attp, *sep;        //the coordinates of the possible attachment and separation //points, if attp or sep are not equal to nullptr, it means we //find the attachment point or separation point. These are only for 3D surface, on plane, we can use icVector2
    bool boundary;               //boundary of morse set 12/04/2009 added by Qingqing
    /*----- To store the special points on the edge (two at most) ----*/
    /***---- For finding the separation and attachment points ----***/


    double length;               //the length of the edge(optional)
    icVector2 normal_2d;
    icVector3 intersections[2];  //we store only the recent two intersections
    icVector2 evec[2];
    icMatrix2x2 Jacobian;     //for calculate the decomposition of the Jacobian


    /*
       Middle point ID for extracting the iso-contour of the vector field magnitude field
    */
    int MidPointID;

    /* Function to calculate the attachment or separation point */

    void cal_neighborhood_with_vectors(int faceid, double alpha[3],
                                double sincx, double sincy,
                                icVector2 Coord[3], icVector2 vec[3]);
    void cal_Jacobian_edge(int edge_tri);
    void cal_eigenvals_edge();

    bool cal_attachment_pt(int triangle, double sep[3]);
    bool cal_separation_pt(int triangle, double sep[3]);

    bool is_valid_sep(int tri, int sccindex);
    bool is_valid_attp(int tri, int sccindex);

    /*routines for periodic orbit extraction*/
    void reset_intersections();
    bool approx_tangent(double seg_length, icVector3 new_intersect);
    bool is_good_intersections(int triangle, int scc_index);
    int *get_disc_coarse(double p[3], int triangle,
             int *NearbyTriangles, int &num_triangles);

    /*
    Judge the type of the edge
    */
    bool is_mixed_edge(int triangle);
    bool is_exit_edge(int);
    void get_vecs_withangledeficit(Vertex *v1, Vertex *v2, Triangle *face,
                             int deficit_index, icVector2 vec[2]);
    void get_vecs_withangledeficit_2(Vertex *v1, Vertex *v2, Triangle *face,
                             int deficit_index, icVector2 vec[2]);
    void get_vecs_at_edge_endingpts(Edge *cur_edge, int triangle, icVector2 vec[2]);

    /*
    calculate the edge normal based on one of the triangles sharing the edge
    */
    void cal_normal_at_edge(Edge *cur_edge, Triangle *face, int which_face);
    void cal_normal_through_edge(Edge *cur_edge, Triangle *face, int which_face);
    int edge_vertex_case(Edge *cur_edge, Triangle *face);

    bool is_pure_exit_edge(int tri);
    bool is_pure_entrance_edge(int tri);

}; //end of Edge class


class EdgeList{
public:
    Edge **edges;
    int nedges;
    int curMaxNumEdges;

    // The similar list operations as VertexList
    EdgeList(int initsize = 1000) //construction
    {
        edges = (Edge **)malloc(sizeof(Edge *)*initsize);
        curMaxNumEdges = initsize;
        nedges = 0;
    }

    ~EdgeList()
    {
        if(edges == nullptr)
            return;

        int i;

        for(i = 0; i < nedges; i++)
        {
            free(edges[i]);
        }
        free(edges);
    }

    //add a new vertex to the end of the list, if it succeeds, return true
    inline bool append(Edge *e)
    {
        if(isFull ())
            if(!extend(100))
                return false;             //if not enough memory available, return false
        edges[nedges] = e;

        nedges ++;
        return true;
    }

    inline bool del_End() //delete the vertex at the end of the list
    {
        if(isEmpty())  return false;
        nedges --;
        return true;
    }

    inline void copy_Elem(Edge *s, Edge *d)
    {
    }

    //delete the corresponding  vertex, if it succeeds, return true
    inline bool del_Node(Edge *e)
    {
        if(isEmpty())  return false;

        //find the vertex, if find it, delete and move the following vertices forward
        //otherwise, return false;

        int i, pos = -1;

        for(i = 0; i < nedges; i++)
        {
            if(edges[i] == e)
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
        if(nedges == curMaxNumEdges) return true;
        return false;
    }

    //extend the original list, if it succeeds, return true
    inline bool extend(int step = 100)
    {
        Edge **temp = edges;
        edges = (Edge **) malloc(sizeof(Edge *) * (curMaxNumEdges + step));
        if( edges == nullptr)
        {
            //fail
            edges = temp;
            return false;
        }

        int i;
        for(i = 0; i < curMaxNumEdges; i++)
            edges[i] = temp[i];
        for(i = curMaxNumEdges; i < curMaxNumEdges+step; i++)
            edges[i] = nullptr;

        curMaxNumEdges += step;
        free(temp);
        return true;
    }

    inline void reset()
    {
        nedges = 0;
    }


}; //end of EdgeList class


/** Corner **/
class Corner{
public:
    unsigned char Edge_count;   //special variable for edges search 1/21/05
    int index;
    int v;            //the ID of the vertex of the corner
    int n;            //the next corner according to the orientation
    int p;            //the previous corner according to the orientation
    int t;            //the triangle the corner belongs to
    int ot;           //the index of its opposite triangle for traversal 2/9/05
    int o;            //the index of its opposite corner
    Edge *e;          //the opposite edge of the corner
    //double angle;       //the angle of the corner
    double BeginAng, EndAng;  //for correct angle allocation of the corner around vertex v
    double r;
    bool orient;

    /////////////////////////////////////////////////
    /* Optional variables */
    Edge *edge[2];    //two edges associated with this corner

    double get_Angle();
};


class CornerList{
public:
    Corner **corners;
    int ncorners;
    int curMaxNumCorners;

    // The similar list operations as VertexList
    CornerList(int initsize = 1000) //construction
    {
        corners = (Corner **)malloc(sizeof(Corner *)*initsize);
        curMaxNumCorners = initsize;
        ncorners = 0;
    }

    ~CornerList()
    {
        int i;

        for(i = 0; i < ncorners; i++)
        {
            free(corners[i]);
        }
        free(corners);
    }

    //add a new vertex to the end of the list, if it succeeds, return true
    inline bool append(Corner *c)
    {
        if(isFull ())
            if(!extend(100))
                return false;             //if not enough memory available, return false
        corners[ncorners] = c;

        ncorners++;
        return true;
    }

    inline bool del_End() //delete the vertex at the end of the list
    {
        if(isEmpty())  return false;
        ncorners --;
        return true;
    }

    inline void copy_Elem(Corner *s, Corner *d)
    {
    }

    //delete the corresponding  vertex, if it succeeds, return true
    inline bool del_Node(Corner *c)
    {
        if(isEmpty())  return false;

        //find the vertex, if find it, delete and move the following vertices forward
        //otherwise, return false;

        int i, pos = -1;

        for(i = 0; i < ncorners; i++)
        {
            if(corners[i] == c)
            {
                pos = i;
                break;
            }

        }

        if(pos == -1) return false;

        //delete it
        for(i = pos; i < ncorners-1; i++)
        {
            //we need a copy function
            copy_Elem(corners[i], corners[i+1]);
        }

        ncorners--;

        return true;
    }

    inline bool isEmpty()  //judge whether the list is empty
    {
        if(ncorners == 0)   return true;
        return false;
    }

    inline bool isFull()
    {
        if(ncorners == curMaxNumCorners) return true;
        return false;
    }

    //extend the original list, if it succeeds, return true
    inline bool extend(int step = 100)
    {
        Corner **temp = corners;
        corners = (Corner **) malloc(sizeof(Corner *) * (curMaxNumCorners + step));
        if( temp == nullptr)
        {
            //fail
            corners = temp;
            return false;
        }

        int i;

        for(i = 0; i < curMaxNumCorners; i++)
            corners[i] = temp[i];

        for(i = curMaxNumCorners; i < curMaxNumCorners+step; i++)
            corners[i] = nullptr;

        curMaxNumCorners += step;
        free(temp);
        return true;
    }

    inline void reset()
    {
        ncorners = 0;
    }


    //void build_Corners();    //build the corner table

}; // end of CornerList class


class ECG_Node;

class Singularity{
public:
    bool connected;
    unsigned char nnodes;
    unsigned char  num_connect_cycles;
    unsigned char type;                    //the type of the singularity
    int  index;                            //the index of the singularity in the list
    int  TriangleID;                       //the triangle containing the singularity
    int *connect_cycles;
    int separatrices;                      //the index of the separatrix group belongs to the saddle 10/13/05
    int  node_index;                      //the index of the node in ECG graph

    icVector3  gpos;                       //global coordinates of singularity
    icVector2  lpos;                       //local coordinates of singularity
    icVector2  lincoming, loutgoing;       //the local e-vectors of the singularity (esp. saddle)
    icMatrix2x2 Jacobian;                  //the jacobian associated with the fixed point

    /*-------for ECG-graph --------if I can find a way to build the graph during analysis, we can safely remove these variables */
    //Singularity **connect_singulars;         //singularities connecting with the singularity
    //unsigned char num_connect_singulars;
    //LimitCycle **connect_cycles;             //limit cycles connecting with the singularity
    ECG_Node **nodes;                     //the list of nodes connecting with the singularity in ECG

    /*member functions for updating the connection information involved this fixed point*/
    void update_list_to_PO(int limitcycle);

}; // end of Singularity class


class SingularityList{
public:
    Singularity **slist;
    int nsingularities;
    int curMaxNumSingularities;

    // The similar list operations
    SingularityList(int initsize = 200) //construction
    {
        slist = (Singularity **)malloc(sizeof(Singularity *)*initsize);
        curMaxNumSingularities = initsize;
        nsingularities = 0;
    }

    //add a new vertex to the end of the list, if it succeeds, return true
    inline bool append(Singularity *s)
    {
        if(isFull ())
            if(!extend(100))
                    return false;             //if not enough memory available, return false
        slist[nsingularities] = s;
        nsingularities++;
        return true;
    }

    inline bool del_End() //delete the vertex at the end of the list
    {
        if(isEmpty())  return false;
        nsingularities --;
        return true;
    }

    inline void copy_Elem(Singularity *s, Singularity *d)
    {
    }

    //delete the corresponding  vertex, if it succeeds, return true
    inline bool del_Node(Singularity *s)
    {
        if(isEmpty())  return false;

        //find the vertex, if find it, delete and move the following vertices forward
        //otherwise, return false;

        int i, pos = -1;

        for(i = 0; i < nsingularities; i++)
        {
            if(slist[i] == s)
            {
                pos = i;
                break;
            }

        }

        if(pos == -1) return false;

        //delete it
        for(i = pos; i < nsingularities-1; i++)
        {
            //we need a copy function
            copy_Elem(slist[i], slist[i+1]);
        }

        nsingularities--;

        return true;
    }

    inline bool isEmpty()  //judge whether the list is empty
        {
        if(nsingularities == 0)   return true;
        return false;
    }

    inline bool isFull()
    {
        if(nsingularities == curMaxNumSingularities) return true;
        return false;
    }

    //extend the original list, if it succeeds, return true
    inline bool extend(int step = 100)
    {
        Singularity **temp = slist;
        slist = (Singularity **) malloc(sizeof(Singularity *) * (curMaxNumSingularities + step));
        if( temp == nullptr)
        {
            //fail
            slist = temp;
            return false;
        }

        int i;
        for(i = 0; i < curMaxNumSingularities; i++)
            slist[i] = temp[i];

        for(i = curMaxNumSingularities; i < curMaxNumSingularities+step; i++)
            slist[i] = (Singularity*)malloc(sizeof(Singularity));

        curMaxNumSingularities += step;
        return true;
    }

    inline void reset()
    {
        nsingularities = 0;
    }

    friend class TriangleList;


}; //end of singularityList class


/*
/** Polyhedron **/
class Polyhedron{
public:
    VertexList vlist;          //the vertex list
    TriangleList tlist;        //the triangle list
    EdgeList elist;            //the edge list
    CornerList clist;          //the corner list
    SingularityList slist;     //the singularity list
    icVector3 center;
    int orient;                //0--counter clockwise; 1--clockwise
    bool with_vec;

    double radius;
    double area;
    double shortest_edgelength;

    PlyOtherProp *vert_other,*face_other;


    icVector3 rot_center;

    PlyLoader myplyloader;    //for ply file loading

    int *singular_tri;
    int nfixedpts;


    double min_vfmagnitude, max_vfmagnitude, mean_vfmagnitude, min_VecLenVisualized, max_VecLenVisualized;

public:


    Polyhedron();

    Polyhedron(FILE *, int);

    void build_Edges();         //initialize all the edges based on vertex and triangle information
    void create_Edge(Vertex *, Vertex *);
    Triangle *find_common_Edge(Triangle *, Vertex *, Vertex *);
    Triangle *other_Triangle(Edge *, Triangle *);
    void create_Pointers();
    void vertex_to_tri_ptrs();
    void order_Vertex_to_Tri_Ptrs(Vertex *v);
    void add_Edge_to_Vertex(Edge *, Vertex *);

    void preprocess_Vertex(void);
    void average_Normals();
    void calc_Bounding_Sphere();
    void calc_Face_Normals_and_Area();
    void calc_Face_Normals_and_Area_2();
    void calc_Edge_Length();
    void get_vertex_normal();

    void assign_color_VFmagnitude();
    void assign_color_VFmagnitude2();
    void normalize_Field();
    void project_to_TangentPlane();

    void build_Corners();  //build the corner table
    void add_Corner_to_Vertex(Corner *, Vertex *v);
    void find_Opposite();
    void sort_Corner_on_Vertex();
    int get_First_Corner(Vertex *p);
    void alloc_Corner_Angle();
    double get_Angle(Corner *c);
    int get_orientation(Vertex *p, double cur_ang, Edge *e1);
    void find_all_corners_around_vertex(int vertid);

    void init_Local_Frame();
    void init_Trace(bool);
    double get_Phi(int vertid, int faceid, double vk[3]);
    icVector2 get_Miu(int vertid, int faceid, double vk[3]);
    void recal_vec_in_triangle();

    void initialize();
    void finalize();

    void cal_TexCoord(void);


    void get_2D_Barycentric_Facters(int faceid, double a, double b, double alpha[3]);
    icVector2 get_Vector_At_Point(int faceid, double vk[3], double alpha[3], double a, double b);
    icVector2 get_A_2D_Point(int faceid, double alpha[3]);
    void local_To_Global(int faceid, icVector2 locpos, icVector3 &glpos);
    bool is_Valid_Alpha(double alpha[]);


    void capture_Singularities(); //capture all the singularities in the mesh
    void compute_FixedPts();
    void get_FixedPt_Coord(int TriangleID, int switch_flag,
                      double vx[3], double vy[3],
                      icVector2& VM, icVector2 &WQM,
                      icVector2 tempV1, icVector2 tempV2,
                      double Q[2], double M[2]);
    void get_M_onRS(double R[2], double S[2],
                double Q[2], double M[2],
                icVector2 &VM, icVector2 &WQM,
                icVector2 Vr, icVector2 Vs);
    void cal_SingularCoord_without_angledeficit(double px[3], double py[3],
                                          double vx[3], double vy[3],
                                          double locl_p[2]);
    int get_Type_of_FixedPt(int TriangleID, double alpha[3]);
    void get_Neighbor_Triangle_Vectors(int faceid, double alpha[3],
                                double sincx, double sincy,
                                icVector2 Coord[3], icVector2 vec[3]);
    void add_To_SingularityList(int Triangle_ID, icVector3 gPos, double local_x, double local_y);



    void cal_separatrices();


    void cal_all_tri_Jacobians();



    void detect_PeriodicOrbit();   //call limit cycle detector

    void simplify_VectorField();  //automatic simplification, need to call the function from vector field simplifier


    //void subdive_Mesh();
    //void remeshing();
    //void deforme_Mesh();
    void smooth_VectorField();
    void smooth_Surface();


    void mark_in_outlets();


    void XAxisVF(void);

    void write_ply(FILE *file);

    friend class Trajectory;
    friend class PeriodicOrbit_Detector;

    float min_rot_sum, max_rot_sum;

};  //end of Mesh class

void  HsvRgb( float hsv[3], float rgb[3] );

#endif /* __GEOMETRY_H__ */
