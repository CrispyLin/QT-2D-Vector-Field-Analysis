/*
This header file contains the class declaration of periodic orbit detector

Created and Modified by Guoning Chen
        copyright @2007
    */

#pragma once
#ifndef __PERIODICORBITEXTRACTION_H__
#define __PERIODICORBITEXTRACTION_H__

#include "VField.h"
#include "Others/common_routines.h"
#include "MorseDecomp.h"

    class PeriodicOrbit_Detector{
public:
    Edge *chosen_edge;
    //int *cellcycle;
    //int ntris_incell;
    DynList_Int *cellcycle;
    Trajectory *detect_traj;

    /* Member Functions */
    PeriodicOrbit_Detector(int init_cellcycle = 300)
    {
        chosen_edge = NULL;
        cellcycle = new DynList_Int(init_cellcycle);  /*initialize the cellcycle list*/
    }

    ~PeriodicOrbit_Detector()
    {
        delete [] cellcycle;
    }

    void detect_periodicorbit_morsedecomp();

    void start_regular(int scc_index, int type);
    void start_seporatt_pts(int scc_index, int type);
    bool trace_to_find_recurrence(double g[3], int &triangle, int type,
                                  int scc_index, Edge *chosen_edge,  int &flag);
    int trace_in_triangle_po_detect(double globalp[3], int &face_id, int type,
                                    Edge *cur_e, icVector3 &intersect, double &seg_length, int &flag);

    bool cal_closedstreamline(Edge *the_e, int triangle, int type, int scc_index);

    void  store_current_cellcycle(int *acellcycle, int num);
    void store_cur_streamline(LineSeg *streamline, int num);

    /*this routine should go with each SCC*/
    void mark_neighbor_pts(double p[3], int triangle, int scc_index, int type);
    int *get_disc_coarse(double p[3], int triangle,
                         int *NearbyTriangles, int &num_triangles);

    /*compare current obtained cell cycle with existing periodic orbits
    and see whether they are matched*/
    bool is_existing_porbit(int triangle, int &pre_limit_index, unsigned char type);
};


#endif /* __PERIODICORBITEXTRACTION_H__ */
