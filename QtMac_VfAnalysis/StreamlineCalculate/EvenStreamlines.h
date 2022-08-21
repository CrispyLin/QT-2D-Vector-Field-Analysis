/*
This file contains the declaration of class EvenStreamlinePlace

Created and modified by Guoning Chen
        copyright @2007
    */
#pragma once
#ifndef EVENSTREAMLINES_H
#define EVENSTREAMLINES_H

#include <QString>

#include "VField.h"
#include "Others/common_routines.h"
#include "Edit/RotateReflectField.h"
#include "GlView.h"
#include "CGraphView.h"
#include "Windows/MainWindow.h"

extern const int DebugOn;

//#include "common_routines.h"
class EvenStreamlinePlace
{
public:
    TrajectoryList *evenstreamlines;
    SeedList *seedpts;
    SamplePtList **samplepts;
    int cur_traj; /*the index of current streamline*/
    CurvePoints *tracing_points /*= (CurvePoints*) malloc(sizeof(CurvePoints) * 805)*/;
    int num_tracingpoints ;
    DynList_Int *trianglelist;
    EdgeList *geo_boundary;

    /*important parameters to control the even placement of the streamlines*/
    double streamlinelength;
    double dsep;
    double percentage_dsep;
    double discsize;
    double sample_interval;
    int every_nsample;
    double loopdsep;
    double dist2sing;
    double seeddist;
    double minstartdist;

    EvenStreamlinePlace()
    {
        evenstreamlines = new TrajectoryList(300);
        samplepts = (SamplePtList **)malloc(sizeof(SamplePtList*)*evenstreamlines->curMaxNumTrajs);

        //FILE *fp = fopen("streamline_err_log.txt", "a");
        //fprintf(fp, "the address for evenstreamlines is %d.\n", evenstreamlines);
        //fprintf(fp, "the number of current streamlines is %d\n", evenstreamlines->curMaxNumTrajs);
        //fprintf(fp, "the address of sample list is %d.\n", samplepts);
        //fclose(fp);

        if(samplepts == NULL)
        {
            char rout[256], var[256];
            sprintf(rout, "%s", "EvenStreamlinePlace Constructor");
            sprintf(var, "%s", "samplepts");

            //write_mem_error(rout, var, 0);
            exit(-1);
        }

        for(int i = 0; i < evenstreamlines->curMaxNumTrajs; i++)
        {
            samplepts[i] = new SamplePtList(100);

            if(samplepts[i] == NULL)
            {
                char rout[256], var[256];
                sprintf(rout, "%s", "EvenStreamlinePlace Constructor");
                sprintf(var, "%s", "samplepts[i]");

                //write_mem_error(rout, var, 0);
                exit(-1);
            }
        }

        trianglelist = NULL;
        geo_boundary = NULL;
    }

    ~EvenStreamlinePlace()
    {
        delete evenstreamlines;

        for(int i = 0; i < evenstreamlines->curMaxNumTrajs; i++)
        {
            if(samplepts[i] != NULL)
                free(samplepts[i]);
        }
        free(samplepts);
    }

    /*initialize the element inside the lists*/
    inline void init()
    {
        int i, j;
        for(i = 0; i < evenstreamlines->curMaxNumTrajs; i++)
        {
            evenstreamlines->trajs[i] = new Trajectory(i);

            if(evenstreamlines->trajs[i] == NULL)
            {
                char rout[256], var[256];
                sprintf(rout, "%s", "EvenStreamlinePlace::init");
                sprintf(var, "%s", "evenstreamlines->trajs[i]");

                //write_mem_error(rout, var, 0);
                exit(-1);
            }
        }

        for(i = 0; i < evenstreamlines->curMaxNumTrajs; i++)
        {
            for(j = 0; j < samplepts[i]->curMaxNumSamplePts; j++)
            {
                //samplepts[i]->samples[j] = new SamplePt[1];
                samplepts[i]->samples[j] = (SamplePt *)malloc(sizeof(SamplePt));
                if(samplepts[i]->samples[j] == NULL)
                {
                    char rout[256], var[256];
                    sprintf(rout, "%s", "EvenStreamlinePlacement::init");
                    sprintf(var, "%s", "samplepts[i]->samples[j]");

                    //write_mem_error(rout, var, 0);
                    exit(-1);
                }
            }
        }
    }

    void reset_placement();

    //void place_streamlines(int num_initial, double dsep, double percentage_dsep,
    //			  double discsize, double sample_interval, int every_nsample,
    //			  double loopdsep, double dist2sing, double streamlinelength,
    //			  double seeddist, double minstartdist, int flag);
    void place_streamlines(int flag);

    void cal_init_streamlines(int num, double streamlinelength, double dtest, double discsize,
                              double Sample_interval, double loopdsep, double dist2sing);
    /*
    compute the initial set of streamlines considering separatrices and periodic orbits
    */
    void cal_init_streamlines_enhanced(double streamlinelength, double dtest, double discsize,
                                       double sample_interval, double dist2sing);


    /*
    Trace one streamline and judgement the density during tracing
    */
    bool grow_a_streamline(double seed_p[3], int triangle, double dtest, double discsize,
                           double Sample_interval, double loopdsep, double dist2sing, double streamlinelength);

    void cal_seeds(int traj, double dsep, int every_nsample);
    bool get_a_seed(double sample_p[3], int begin_triangle, int cur_traj, int type,
                    double end_p[3], int &end_triangle, double dsep);
    /*
    Trace inside a triangle for the even streamline placement
    */
    int trace_in_triangle(int &face_id, double globalp[2], int type,
                          double dtest, double loopsep, double dist2sing,
                          double sample_interval, double discsize, int &flag);

    /*
    Get a set of sample points during tracing a new streamline
    This sample points are temporary and for preventing the closed loop of a streamline
    It will keep finding all the samples from 'cur_line' line segment till the end of
    current streamline
    Note that we use local coordinates to get one sample point each time,
    after that, we need to transform it back to 3D global coordinates
    */
    SamplePt **cal_samplepts_when_tracing(int traj, double interval, int &cur_line, int &movetonext, double &cur_length,
                                          SamplePt **samples, int &num_samples);

    /*
    Reverse the order of line segments of a specific streamline
    */
    void reverse_streamline(int streamlineid);

    /*
    Update the smaple point list of all the triangles associated with the input sample list
    */
    void update_samples_in_triangle(int traj, SamplePt **samples, int num_samples);

    /*
    Perform a short local tracing from each sampling point to find a seed point,
    we still need to store the temporaylist to a global line segment list
    */
    int trace_in_triangle_seed(int &face_id, double globalp[3], int type,
                               int &flag, double pre_p[3], double cur_p[3],
                               double dsep, double &cur_length, int cur_traj);

    /*
    During seed point finding, we use local tracing with some constant step size
    so, we need to get the seed point having the exact distance to the streamline
    */
    void cal_exact_seed(double cur_p[3], double pre_p[3], int triangle,
                        double extra_length, double exact_seed[3]);


    /*distance judgement routines*/
    /*
    Judge whether a given point is too close to existing fixed points except for the singid
    */
    bool close_to_fixedPt_except(double p[3], int triangle, int singid, double threshold, double discsize);
    /*
    Judge whether current point is too close to existing streamline except for the set traj
    */
    bool close_to_cur_streamline(double p[3], int triangle, int *Except_trajs, int num_trajs,
                                 double separate_dist, double discsize, int Init);
    /*Using the similar streamline placement for 2D planar vector field*/
    bool close_to_cur_samplePt(double p[3], int triangle, SamplePt **samples, int num_samples,
                               double separate_dist, double discsize, double sample_interval);

    /*
    Get one sample point on the regular streamline, not periodic orbit
    Different from planar case, we use local coordinates here
    */
    bool cal_a_sample_of_streamline(int traj, int &cur_lineindex, int &movetonext,
                                    double curpt[2], double interval, double &cur_length);

    /*
    Add one sample to the particular triangle
    */
    void add_sample_to_triangle(int triangle, int which_traj, int which_sample);

    /*
    save the placement results into a file
    */
    void save_to_a_file(const char* filename);
    void load_from_a_file(const char *filename);

    void save_to_a_file(const char* filename, int flag);
    void load_from_a_file_enhanced(const char *filename);


    /*The following routines can be classified into the Geodesic calculation class*/

    int *cal_geodesic_dist(int triangle, double globalp[3], double dsep, double discsize,
                           int *trianglelist, int &num_triangles);
    void cal_geodesic_dist_2(int triangle, double globalp[3], double dsep, double discsize,
                             DynList_Int *trianglelist);
    void reset_dist(int *NearbyTriangles, int num);
    void reset_dist(DynList_Int *trianglelist);
    void get_unfold_rotMat(int org_triangle, int edgeid, int other_f, double rotmat[16]);
    int get_thirdVer_of_triangle(int v1, int v2, int triangle);

    void update_one_vertex(int verta, int vertb, int vertc, int triangle);
    void set_boundary_edges(int *trianglelist, int num_triangles);
    void get_dis_for_obtuse_triangle(int v1, int v2, int v3, int origin_triangle, double lena, double lenb);
    /*
    During getting the shadow region, calculate the angles and the local frame associated with v3
    */
    void get_shadow_region(int &v1, int &v2, int v3, int origin_triangle, icVector3 &vec1, icVector3 &vec2,
                           double &theta1, double &theta2, icVector3 &LX, icVector3 &LY);
    /*
    Get the 4th vertex through unfolding a series of nearby triangles
    */
    bool unfold_to_find_ver(int v1, int v2, int v3, int origin_triangle, icVector3 vec1, icVector3 vec2,
                            icVector3 LX, icVector3 LY, double theta1, double theta2, int &theVer, double p[3]);

    int get_edgeIndex_of_triangle(int v1, int v2, int curtriangle);

    bool is_fall_in_shadow(double theta1, double theta2, icVector3 LX, icVector3 LY, double p[3], int origin_v3,
                           int &local_v1, int &local_v2, int local_v3);


};

#endif
