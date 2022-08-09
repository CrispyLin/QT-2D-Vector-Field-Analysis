/*
This file contains the routines for even streamline placement.

Created and modified by Guoning Chen
        copyright @2007
    */

#include "StreamlineCalculate/EvenStreamlines.h"


extern Polyhedron *object;
extern PeriodicOrbitList *periodic_orbits;
extern TrajectoryList *separatrices;


/*temporary solution for calculating the rotation matrix using OpenGL*/
extern CGlView *g_pclGlView;
extern CGraphView *g_pCgraphView;

EvenStreamlinePlace *evenplace = NULL;
extern int DebugOn;
extern clock_t g_start, g_finish;

double StreamlinePlace_slength;
double StreamlinePlace_dsep;


/*for testing the problem of the geodesic distance calculation*/
double g_dist = 0.;
double g_separate_dist = 0.;
int g_sample_tri = -1;

extern int ndisplay_trajs;


void RotateReflectField::rotate_wholefield_angle(double angle)
{
    int i;
    double rotm[4][4];
    icVector3 v_norm;

    double s[3] = {0.};

    for(i = 0; i < the_obj->vlist.nverts; i++)
    {
        v_norm = the_obj->vlist.verts[i]->normal;

        s[0] = the_obj->vlist.verts[i]->g_vec.entry[0];
        s[1] = the_obj->vlist.verts[i]->g_vec.entry[1];
        s[2] = the_obj->vlist.verts[i]->g_vec.entry[2];

        ////First, we need to get the rotation matrix for each vertex
        //glMatrixMode(GL_MODELVIEW_MATRIX);
        //glPushMatrix();
        //glLoadIdentity();
        //glRotatef(angle, v_norm.entry[0], v_norm.entry[1], v_norm.entry[2]);
        //glGetFloatv(GL_MODELVIEW_MATRIX, (float*)rotm);
        //glPopMatrix();

        //get_rotation(angle, v_norm, (double *)rotm);
        //g_pclGlView->get_rotation(angle, v_norm, (double *)rotm);
        //g_pCgraphView->get_rotation(angle, v_norm, (double *)rotm);

        get_rotation(angle*3.1415926/180, v_norm, (double*) rotm);

        the_obj->vlist.verts[i]->g_vec.entry[0] = rotm[0][0]*s[0] + rotm[0][1]*s[1] + rotm[0][2]*s[2];
        the_obj->vlist.verts[i]->g_vec.entry[1] = rotm[1][0]*s[0] + rotm[1][1]*s[1] + rotm[1][2]*s[2];
        the_obj->vlist.verts[i]->g_vec.entry[2] = rotm[2][0]*s[0] + rotm[2][1]*s[1] + rotm[2][2]*s[2];
    }

    the_obj->project_to_TangentPlane();
}


/*Implementation of member functions of class EvenStreamlinePlacement*/
void EvenStreamlinePlace::reset_placement()
{
    int i;
    Triangle *face;

    for(i = 0; i < object->tlist.ntris; i++)
    {
        face = object->tlist.tris[i];
        face->visited = false;
        //face->inDesignCell = false;

        //if(face->samplepts != NULL)
        //{
        //	free(face->samplepts);
        //	face->samplepts = NULL;

        //}
        //face->num_samplepts = 0;
        face->reset_sampleList();
    }

    for(i = 0; i < object->vlist.nverts; i++)
        object->vlist.verts[i]->distance = 1e49;
}


/*
Use sample points to control the termination of streamlines
*/
//void EvenStreamlinePlace::place_streamlines(int num_initial, double dsep, double percentage_dsep,
//				  double discsize, double sample_interval, int every_nsample,
//				  double loopdsep, double dist2sing, double streamlinelength,
//				  double seeddist, double minstartdist, int flag)
void EvenStreamlinePlace::place_streamlines(int flag)
{
    int i;
    int locate_preseed = 0;
    int begin_sample = 0;

    /*-----Set a set of default initial values for testing 5/2/06-----*/
    /*--- this setting should be provided by user ---*/

    //dsep = 0.024 * object->radius;    //using the radius of the object instead of the edge
    dsep = StreamlinePlace_dsep * object->radius;    //using the radius of the object instead of the edge
    percentage_dsep = 0.55;
    discsize = 2.;
    sample_interval = 0.4*dsep;
    every_nsample = 3;
    loopdsep = 0.8*dsep;
    dist2sing = 0.9*dsep;
    //streamlinelength = 0.03*object->radius;
    streamlinelength = StreamlinePlace_slength*object->radius;
    seeddist = 1.*dsep;
    minstartdist = 0.8*dsep;

    /*******************************************************************/
    /*allocate space for seed points*/
    seedpts = new SeedList();

    cur_traj = 0; //set the first streamline as current streamline

    ndisplay_trajs = 0;

    /*Record execution information into file*/
    //FILE *fp;

    //fp = fopen("streamline_err_log.txt", "a");
    //fprintf(fp, "start initializing the first streamline!\n");
    //fprintf(fp, "the flag is %d...\n", flag);
    //fclose(fp);

    //Grow a set of initial streamlines
    if(flag == 0)
    {

        evenstreamlines->ntrajs = 0;
        cal_init_streamlines(1, streamlinelength, percentage_dsep*dsep,
                             discsize, sample_interval, loopdsep, dist2sing);
    }

    else{
        cal_init_streamlines_enhanced(streamlinelength, percentage_dsep*dsep, discsize, sample_interval,
                                      dist2sing);
    }

    ////Initialize some flags
    /*---------------------------------*/

    seedpts->nseeds = 0; //reset the seedpts list

    //fp = fopen("streamline_err_log.txt", "a");
    //fprintf(fp, "start calculating seed points!\n");
    //fclose(fp);

    cal_seeds(cur_traj, seeddist, every_nsample);

    //fp = fopen("streamline_err_log.txt", "a");
    //fprintf(fp, "%d seed points have been found!\n", seedpts->nseeds);
    //fclose(fp);

    locate_preseed = 0;
    while(1 /*&& evenstreamlines->ntrajs < 2*/) //if there still are streamlines not being processed
    {
        ////at each step, we grow streamlines for all the seed points associated with current streamline
        for(i = locate_preseed; i < seedpts->nseeds; i++)
        {
            ////Judge whether this is a valid seed or not
            if(!close_to_cur_streamline(seedpts->seeds[i]->pos,
                                         seedpts->seeds[i]->triangle, &cur_traj, 1, minstartdist, discsize, 0))
            {
                ////if we find a valid seed point, break
                locate_preseed = i+1;
                break;
            }

            seedpts->seeds[i]->state = 2;  //reject before starting
        }

        if(i >= seedpts->nseeds ) // it means no more seeds available for current streamline
        {
            locate_preseed = seedpts->nseeds;
            cur_traj ++;  //let next streamline as current streamline

            if(cur_traj >= evenstreamlines->ntrajs) //no more streamlines available
            {
                /*release the seed point list*/
                delete seedpts;
                return;
            }

            //Get seeds for the current streamline
            cal_seeds(cur_traj, seeddist, every_nsample);
            continue;
        }

        ////else, we find a valid seed point, then grow a new streamline from this seed

        //fp = fopen("streamline_err_log.txt", "a");
        //fprintf(fp, "choose seed point %d: (%f, %f, %f) at triangle %d.\n", i,
        //	seedpts->seeds[i]->pos[0], seedpts->seeds[i]->pos[1], seedpts->seeds[i]->pos[2],
        //	seedpts->seeds[i]->triangle);
        //fclose(fp);

        if(grow_a_streamline(seedpts->seeds[i]->pos, seedpts->seeds[i]->triangle,
                              percentage_dsep*dsep, discsize, sample_interval, loopdsep, dist2sing, streamlinelength))
        {
            seedpts->seeds[i]->state = 1;

            if(evenstreamlines->isFull())
            {
                int oldnum = evenstreamlines->curMaxNumTrajs;
                if(!evenstreamlines->extend())
                {
                    /*record the error information*/
                    char rout[256], var[256];
                    sprintf(rout, "%s", "EvenStreamlinePlacement::place_streamlines");
                    sprintf(var, "%s", "evenstreamlines");

                    //write_mem_error(rout, var, 1);
                    delete seedpts;
                    //exit(-1);

                    /*release the seed point list*/
                    return;
                }

                /*allocate memeory for the elements*/
                for(i = oldnum; i < evenstreamlines->curMaxNumTrajs; i++)
                {
                    evenstreamlines->trajs[i] = new Trajectory(i);
                    if(evenstreamlines->trajs[i] == NULL)
                    {
                        char rout[256], var[256];
                        sprintf(rout, "%s", "EvenStreamlinePlacement::place_streamlines");
                        sprintf(var, "%s", "evenstreamlines->trajs[i]");

                        //write_mem_error(rout, var, 1);
                        return;
                    }
                }

                /*increase the sample point lists as well*/
                SamplePtList **temp = samplepts;
                samplepts = (SamplePtList **)malloc(sizeof(SamplePtList*)*evenstreamlines->curMaxNumTrajs);
                if(samplepts == NULL)
                {
                    /*report error*/
                    if(samplepts == NULL)
                    {
                        char rout[256], var[256];
                        sprintf(rout, "%s", "EvenStreamlinePlacement::place_streamlines");
                        sprintf(var, "%s", "samplepts");

                        //write_mem_error(rout, var, 1);
                        return;
                    }

                    /*release the seed point list*/
                    delete  seedpts;
                    return;
                }

                /*copy the information from the old list*/
                for(i = 0; i < oldnum; i++)
                    samplepts[i] = temp[i];
                delete [] temp;

                //FILE *fp;
                //	fp = fopen("traj_realloc.txt", "w");
                //	fprintf(fp, "assign the sample list\n");
                //	fclose(fp);

                /*allocate the pointers for the real samples*/
                for(i = oldnum; i < evenstreamlines->curMaxNumTrajs; i++)
                {
                    samplepts[i] = new SamplePtList(100);
                    if(samplepts[i] == NULL)
                    {
                        char rout[256], var[256];
                        sprintf(rout, "%s", "EvenStreamlinePlacement::place_streamlines");
                        sprintf(var, "%s", "samplepts[i]");

                        //write_mem_error(rout, var, 1);
                        return;
                    }

                    for(int j = 0; j < samplepts[i]->curMaxNumSamplePts; j++)
                    {
                        samplepts[i]->samples[j] = (SamplePt *)malloc(sizeof(SamplePt));

                        if(samplepts[i]->samples[j] == NULL)
                        {
                            char rout[256], var[256];
                            sprintf(rout, "%s", "EvenStreamlinePlacement::place_streamlines");
                            sprintf(var, "%s", "samplepts[i]->samples[j]");

                            //write_mem_error(rout, var, 1);
                            return;
                        }
                    }
                }
                //fp = fopen("traj_realloc.txt", "w");
                //fprintf(fp, "finish assigning the sample list\n");
                //fclose(fp);

            }
        }

        else
        {
            seedpts->seeds[i]->state = 3;   //the seed rejected due to the short streamline
        }

    }

    /*release the seed point list*/
    delete  seedpts;

    /*release the sample point list*/
    for(i = 0; i < evenstreamlines->curMaxNumTrajs; i++)
        delete samplepts[i];
    delete [] samplepts;
    reset_placement();

}

/*
calculate a set of initial streamlines without considering the obtained separatrices and periodic orbits
*/
void EvenStreamlinePlace::cal_init_streamlines(int num, double streamlinelength, double dtest, double discsize,
                                               double Sample_interval, double loopdsep, double dist2sing)
{
    double begin_p[3] = {0.};
    int pre_triangle = 0;
    int triangle;

    int cur_traj_index;

    triangle = 10;
    //triangle = 1;  /*need to avoid the triangle containing fixed point*/
    //triangle = (rand())%object->tlist.ntris;  //we should use farthest scheme to generate

    /////////////////////////

    while(1)
    {
        if(evenstreamlines->ntrajs > 0)
            triangle = (pre_triangle + rand())%object->tlist.ntris;  //we should use farthest scheme to generate

        object->tlist.tris[triangle]->get_center(begin_p);

        cur_traj_index = evenstreamlines->ntrajs;

        //FILE *fp;

        //fp = fopen("streamline_err_log.txt", "a");
        //fprintf(fp, "3 start initializing the first streamline!\n");
        //fclose(fp);

        if(evenstreamlines->ntrajs > 0 && close_to_cur_streamline(begin_p, triangle, &cur_traj_index, 1, 10*dtest, discsize, 0))
        {
            pre_triangle = triangle;
            continue;
        }

        FILE *fp;

        fp = fopen("streamline_err_log.txt", "w");
        fprintf(fp, "start initializing the first streamline!\n");
        fprintf(fp, "from (%f, %f, %f) at triangle %d.\n",
                begin_p[0], begin_p[1], begin_p[2], triangle);
        fclose(fp);


        if(grow_a_streamline(begin_p, triangle,
                              dtest, discsize, Sample_interval, loopdsep, dist2sing, streamlinelength))
        {
            pre_triangle = triangle;
        }

        fp = fopen("streamline_err_log.txt", "a");
        fprintf(fp, "finish initializing the first streamline!\n");
        fclose(fp);

        if(evenstreamlines->ntrajs >= num)
            return;
    }
}


void  EvenStreamlinePlace::cal_init_streamlines_enhanced(double streamlinelength, double dtest, double discsize,
                                                        double Sample_interval, double dist2sing)
{
    int i, j, k;
    int traj;
    Triangle *face;

    int pre_triangle = -1;
    int first_triangle = -1;
    int second_triangle = -1;
    int meetcount = 0;

    int cur_line = 0;
    double cur_length = 0;
    int movetonext = 0;

    evenstreamlines->ntrajs = 0;

    /*testing code*/
    FILE *fp;

    /*we need to store the precomputed periodic orbits first
    we require the periodic orbits being extracted first
    */
    for(i = 0; i < periodic_orbits->nporbits; i++)
    {

        ////Update according to the limit cycles at the same time

        /*-------------------------------------------------------------------*/
        first_triangle = periodic_orbits->polist[i]->traj->linesegs[0].Triangle_ID;

        ////find the second triangle to stop the repeating drawing of limit cycle 3/15/06
        for(j = 0; j < periodic_orbits->polist[i]->traj->nlinesegs; j++)
        {
            if(periodic_orbits->polist[i]->traj->linesegs[j].Triangle_ID != first_triangle)
            {
                second_triangle = periodic_orbits->polist[i]->traj->linesegs[j].Triangle_ID;
                break;
            }
        }
        meetcount = 0;
        /*-------------------------------------------------------------------*/

        for(j = 0; j < periodic_orbits->polist[i]->traj->nlinesegs; j++)
        {
            copy_array_Double(periodic_orbits->polist[i]->traj->linesegs[j].gstart.entry,
                              evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[j].gstart.entry, 3);

            copy_array_Double(periodic_orbits->polist[i]->traj->linesegs[j].gend.entry,
                              evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[j].gend.entry, 3);

            copy_array_Double(periodic_orbits->polist[i]->traj->linesegs[j].start.entry,
                              evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[j].start.entry, 2);

            copy_array_Double(periodic_orbits->polist[i]->traj->linesegs[j].end.entry,
                              evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[j].end.entry, 2);

            evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[j].Triangle_ID=
                periodic_orbits->polist[i]->traj->linesegs[j].Triangle_ID;

            //update the length of the line segment

            evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[j].length=
                periodic_orbits->polist[i]->traj->linesegs[j].length;

            evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs++;

            ////extend the space for each trajectory
            if(evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs >=
                evenstreamlines->trajs[evenstreamlines->ntrajs]->curMaxNumLinesegs)
            {

                if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->extend_line_segments(100))
                {
                    return ;
                }
            }


            if(pre_triangle == periodic_orbits->polist[i]->traj->linesegs[j].Triangle_ID)
                continue;

            ///*--------------------------3/15/06---------------------------*/
            if(second_triangle == periodic_orbits->polist[i]->traj->linesegs[j].Triangle_ID)
            {
                if(meetcount == 1)
                    break;
                else
                    meetcount ++;
            }
            ///*-------------------------------------------------------------------*/
            pre_triangle = periodic_orbits->polist[i]->traj->linesegs[j].Triangle_ID;
        }

        pre_triangle = periodic_orbits->polist[i]->traj->nlinesegs =
            evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs;

        /*Build the sample point list for the periodic orbit */

        //////Resample the streamline after the reversion
        cur_line = 0;
        cur_length = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].length;
        samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[0] = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart.entry[0];
        samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[1] = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart.entry[1];
        samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[2] = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart.entry[2];
        samplepts[evenstreamlines->ntrajs]->samples[0]->triangle = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].Triangle_ID;
        samplepts[evenstreamlines->ntrajs]->nsamples = 1;

        movetonext = 0;
        cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
                                   samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);

        /*Update the sample point lists inside the corresponding triangles*/
        update_samples_in_triangle(evenstreamlines->ntrajs, samplepts[evenstreamlines->ntrajs]->samples,
                                   samplepts[evenstreamlines->ntrajs]->nsamples);

        /*since we keep the memory of periodic orbits, we don't need to save them into
        the evenstreamline list*/
        evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs = 0;
        delete evenstreamlines->trajs[evenstreamlines->ntrajs];
        evenstreamlines->trajs[evenstreamlines->ntrajs] = NULL;
        evenstreamlines->ntrajs++;

        /*if we don't have enough memory for more streamlines, extend the list*/
        if(evenstreamlines->isFull())
        {
            int oldnum = evenstreamlines->curMaxNumTrajs;
            if(!evenstreamlines->extend())
            {
                /*record the error information*/

                return;
            }

            /*allocate memeory for the elements*/
            for(i = oldnum; i < evenstreamlines->curMaxNumTrajs; i++)
            {
                evenstreamlines->trajs[i] = new Trajectory(i);
            }

            /*increase the sample point lists as well*/
            SamplePtList **temp = samplepts;
            samplepts = (SamplePtList **)malloc(sizeof(SamplePtList*)*evenstreamlines->curMaxNumTrajs);
            if(samplepts == NULL)
            {
                /*report error*/
                char rout[256], var[256];
                sprintf(rout, "%s", "EvenStreamlinePlacement::cal_init_streamlines_enhanced");
                sprintf(var, "%s", "samplepts");

                //write_mem_error(rout, var, 1);
                return;

                return;
            }

            /*copy the information from the old list*/
            for(i = 0; i < oldnum; i++)
                samplepts[i] = temp[i];
            delete [] temp;

            /*allocate the pointers for the real samples*/
            for(i = oldnum; i < evenstreamlines->curMaxNumTrajs; i++)
            {
                samplepts[i] = new SamplePtList(100);
                if(samplepts[i] == NULL)
                {
                    char rout[256], var[256];
                    sprintf(rout, "%s", "EvenStreamlinePlacement::cal_init_streamlines_enhanced");
                    sprintf(var, "%s", "samplepts[i]");

                    //write_mem_error(rout, var, 1);
                    return;
                }

                for(int j = 0; j < samplepts[i]->curMaxNumSamplePts; j++)
                {
                    samplepts[i]->samples[j] = (SamplePt *)malloc(sizeof(SamplePt));
                    if(samplepts[i]->samples[j] == NULL)
                    {
                        char rout[256], var[256];
                        sprintf(rout, "%s", "EvenStreamlinePlacement::cal_init_streamline_enhanced");
                        sprintf(var, "%s", "samplepts[i]->samples[j]");

                        //write_mem_error(rout, var, 1);
                        return;
                    }
                }
            }

        }
    }


    int except_trajs[4] = {0};
    ////Update according to the separatrices
    for(i = 0; i < separatrices->ntrajs; i+=4)
    {
        except_trajs[0] = i+periodic_orbits->nporbits/*i*/;
        except_trajs[1] = i+1+periodic_orbits->nporbits/*i+1*/;
        except_trajs[2] = i+2+periodic_orbits->nporbits/*i+2*/;
        except_trajs[3] = i+3+periodic_orbits->nporbits/*i+3*/;

        for(j = 0; j < 4; j++)
        {
            traj = i+j;  /*get the index of the trajectory in the separatrix list*/

            ////Build the sample point list for the periodic orbit
            cur_line = 0;
            if(separatrices->trajs[traj]->nlinesegs == 0)
                continue;

            cur_length = separatrices->trajs[traj]->linesegs[0].length;
            samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[0] = separatrices->trajs[traj]->linesegs[0].gstart.entry[0];
            samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[1] = separatrices->trajs[traj]->linesegs[0].gstart.entry[1];
            samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[2] = separatrices->trajs[traj]->linesegs[0].gstart.entry[2];
            samplepts[evenstreamlines->ntrajs]->samples[0]->triangle = separatrices->trajs[traj]->linesegs[0].Triangle_ID;
            samplepts[evenstreamlines->ntrajs]->nsamples = 1;
            ////Search the last linesegment that may be close to other existing streamline
            ////or existing fixed point

            except_trajs[0] = traj;

            /*we still need to copy it to the evenstreamline list for calculating the sample points*/
            for(k = 0; k < separatrices->trajs[traj]->nlinesegs; k++)
            {
                copy_array_Double(separatrices->trajs[traj]->linesegs[k].gstart.entry,
                                  evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[k].gstart.entry, 3);

                copy_array_Double(separatrices->trajs[traj]->linesegs[k].gend.entry,
                                  evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[k].gend.entry, 3);

                copy_array_Double(separatrices->trajs[traj]->linesegs[k].start.entry,
                                  evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[k].start.entry, 2);

                copy_array_Double(separatrices->trajs[traj]->linesegs[k].end.entry,
                                  evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[k].end.entry, 2);

                evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[k].Triangle_ID=
                    separatrices->trajs[traj]->linesegs[k].Triangle_ID;

                evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[k].length=
                    separatrices->trajs[traj]->linesegs[k].length;

                evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs++;

                ////extend the space for each trajectory
                if(evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs >=
                    evenstreamlines->trajs[evenstreamlines->ntrajs]->curMaxNumLinesegs)
                {

                    if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->extend_line_segments(100))
                    {
                        return ;
                    }
                }
            }

            cur_line = 0;
            movetonext = 0;

            for(k = 30; k < separatrices->trajs[traj]->nlinesegs; k++)
            {
                ////Truncate the separatrix if it is too close to existing other streamline
                ////or too close to a fixed point
                ////for obeying the even placement rule  4/25/06

                /*Important change 06/18/07*/
                evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs = k;

                cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
                                           samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);


                if(close_to_cur_streamline(separatrices->trajs[traj]->linesegs[k].gstart.entry,
                                            separatrices->trajs[traj]->linesegs[k].Triangle_ID, except_trajs, 1,
                                            dtest, discsize, 1)
                    || close_to_fixedPt_except(separatrices->trajs[traj]->linesegs[k].gstart.entry,
                                               separatrices->trajs[traj]->linesegs[k].Triangle_ID, separatrices->trajs[traj]->saddleID,
                                               dist2sing, discsize)
                    || close_to_cur_samplePt(/*separatrices->trajs[traj]->linesegs[k].gend.entry,*/
                                             separatrices->trajs[traj]->linesegs[k].gstart.entry,
                                             separatrices->trajs[traj]->linesegs[k].Triangle_ID,
                                             samplepts[evenstreamlines->ntrajs]->samples,
                                             samplepts[evenstreamlines->ntrajs]->nsamples,
                                             dtest, discsize, sample_interval))

                {
                    separatrices->trajs[traj]->nlinesegs = k;
                    evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs = k;
                    break;
                }
            }


            /*Update the sample point lists inside the corresponding triangles*/
            update_samples_in_triangle(evenstreamlines->ntrajs, samplepts[evenstreamlines->ntrajs]->samples,
                                       samplepts[evenstreamlines->ntrajs]->nsamples);

            /*since we keep the memory for the separatrix list, we don't really need to copy it to the
            evenstreamline list*/
            evenstreamlines->trajs[evenstreamlines->ntrajs] = 0;
            delete evenstreamlines->trajs[evenstreamlines->ntrajs];
            evenstreamlines->trajs[evenstreamlines->ntrajs] = NULL;
            evenstreamlines->ntrajs++;

            /*if we don't have enough memory for more streamlines, extend the list*/
            if(evenstreamlines->isFull())
            {
                int oldnum = evenstreamlines->curMaxNumTrajs;
                if(!evenstreamlines->extend())
                {
                    /*record the error information*/
                    char rout[256], var[256];
                    sprintf(rout, "%s", "EvenStreamlinePlacement::cal_init_streamlines_enhanced");
                    sprintf(var, "%s", "evenstreamlines");

                    //write_mem_error(rout, var, 1);

                    return;
                }

                /*allocate memeory for the elements*/
                for(i = oldnum; i < evenstreamlines->curMaxNumTrajs; i++)
                {
                    evenstreamlines->trajs[i] = new Trajectory(i);
                    if(evenstreamlines->trajs[i] == NULL)
                    {
                        char rout[256], var[256];
                        sprintf(rout, "%s", "EvenStreamlinePlacement::cal_init_streamlines_enhanced");
                        sprintf(var, "%s", "evenstreamlines->trajs[i]");

                        //write_mem_error(rout, var, 1);
                        return;
                    }
                }

                /*increase the sample point lists as well*/
                SamplePtList **temp = samplepts;
                samplepts = (SamplePtList **)malloc(sizeof(SamplePtList*)*evenstreamlines->curMaxNumTrajs);
                if(samplepts == NULL)
                {
                    /*report error*/
                    char rout[256], var[256];
                    sprintf(rout, "%s", "EvenStreamlinePlacement::cal_init_streamlines_enhanced");
                    sprintf(var, "%s", "samplepts");

                    //write_mem_error(rout, var, 1);

                    return;
                }

                /*copy the information from the old list*/
                for(i = 0; i < oldnum; i++)
                    samplepts[i] = temp[i];
                delete [] temp;

                /*allocate the pointers for the real samples*/
                for(i = oldnum; i < evenstreamlines->curMaxNumTrajs; i++)
                {
                    samplepts[i] = new SamplePtList(100);
                    if(samplepts[i] == NULL)
                    {
                        char rout[256], var[256];
                        sprintf(rout, "%s", "EvenStreamlinePlacement::cal_init_streamlines_enhanced");
                        sprintf(var, "%s", "samplepts[i]");

                        //write_mem_error(rout, var, 1);
                        return;
                    }

                    for(int j = 0; j < samplepts[i]->curMaxNumSamplePts; j++)
                    {
                        samplepts[i]->samples[j] = (SamplePt *)malloc(sizeof(SamplePt));
                        if(samplepts[i]->samples[j] == NULL)
                        {
                            char rout[256], var[256];
                            sprintf(rout, "%s", "EvenStreamlinePlacement::cal_init_streamlines_enhanced");
                            sprintf(var, "%s", "samplepts[i]->samples[j]");

                            //write_mem_error(rout, var, 1);
                            return;
                        }
                    }
                }

            }


        }
    }

}

bool EvenStreamlinePlace::grow_a_streamline(double seed_p[3], int triangle, double dtest, double discsize,
                                            double Sample_interval, double loopdsep, double dist2sing, double streamlinelength)
{
    int i;
    int flag = -1;

    int pre_face, cur_face;
    double globalp[3] = {0.};
    int cur_line = 0;
    double cur_length = 0;
    int movetonext = 0;

    if(triangle < 0 || triangle > object->tlist.ntris)
        return false;

    pre_face = cur_face = triangle;
    globalp[0] = seed_p[0];
    globalp[1] = seed_p[1];
    globalp[2] = seed_p[2];

    FILE *fp;


    evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs = 0;


    samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[0] = seed_p[0];
    samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[1] = seed_p[1];
    samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[2] = seed_p[2];
    samplepts[evenstreamlines->ntrajs]->samples[0]->triangle = triangle;
    samplepts[evenstreamlines->ntrajs]->nsamples=1;


    cur_line = 0;
    cur_length = 0;
    //////////////////////////////////////////////////////////////////////////


    ////Backward tracing
    int NUMTRACETRIS = (int)sqrt((double)object->tlist.ntris);

    for(i = 0; i < 3*NUMTRACETRIS/*NUMTRACINGTRIANGLE*/; i++)
    {
        ////The triangle does not exist. Something is wrong!
        if(cur_face < 0 || cur_face > object->tlist.ntris)
        {
            break;
        }

        pre_face = cur_face;
        cur_face = trace_in_triangle(cur_face, globalp, 1, dtest, loopdsep, dist2sing,
                                     Sample_interval, discsize, flag);

        //if this is the first line segment, there may be a case that the start point is on the edge !!!
        if(pre_face == triangle && evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs>0)
        {
            cur_length = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].length;
        }


        ////We need to select the sampling points from current trajectory  4/22/06
        //samplepts[evenstreamlines->ntrajs]->samples =
        cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
                                   samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);


        if(flag == 1 || flag == 2 || flag == 3 || pre_face == cur_face )
        {
            //fp1 = fopen("streamline_err_log.txt", "a");
            //fprintf(fp1, "stoping reason: flag %d\n", flag);
            //fprintf(fp1, "# of triangles being travelled:  %d\n\n", i);
            //fclose(fp1);

            break;
        }

    }

    ////Reverse the order of the obtained line segments
    if(evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs>0)
    {
        reverse_streamline(evenstreamlines->ntrajs);

        //////Resample the streamline after the reversion
        cur_line = 0;
        cur_length = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].length;
        samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[0] = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart.entry[0];
        samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[1] = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart.entry[1];
        samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[2] = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart.entry[2];
        samplepts[evenstreamlines->ntrajs]->samples[0]->triangle = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].Triangle_ID;
        samplepts[evenstreamlines->ntrajs]->nsamples = 1;

        movetonext = 0;
        //samplepts[evenstreamlines->ntrajs]->samples =
        cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
                                   samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);
    }

    //////////Forward tracing
    pre_face = cur_face = triangle;
    globalp[0] = seed_p[0];
    globalp[1] = seed_p[1];
    globalp[2] = seed_p[2];
    flag = -1;

    for(i = 0; i < 3*NUMTRACETRIS/*NUMTRACINGTRIANGLE*/; i++)
    {

        if(cur_face < 0)
        {
            break;
        }

        pre_face = cur_face;
        cur_face = trace_in_triangle(cur_face, globalp, 0, dtest, loopdsep, dist2sing,
                                     Sample_interval, discsize, flag);


        ////We need to select the sampling points from current trajectory  3/9/06
        //samplepts[evenstreamlines->ntrajs]->samples =
        cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
                                   samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);


        if(flag == 1 || flag == 2 || flag == 3|| pre_face == cur_face )
        {
            //fp1 = fopen("streamline_err_log.txt", "a");
            //fprintf(fp1, "stoping reason: flag %d\n", flag);
            //fprintf(fp1, "# of triangles being travelled:  %d\n\n", i);
            //fclose(fp1);

            //if(i >= 2 && i < 4)
            //{
            //fp1 = fopen("streamline_err_log.txt", "a");
            //fprintf(fp1, "gdist: %f\n", g_dist);
            //fprintf(fp1, "gseparate_dist: %f\n", g_separate_dist);
            //fclose(fp1);
            //}

            break;
        }

        //if(i >= 2 && i < 4)
        //{
        //fp1 = fopen("streamline_err_log.txt", "a");
        //fprintf(fp1, "gdist: %f\n", g_dist);
        //fprintf(fp1, "gseparate_dist: %f\n", g_separate_dist);
        //fclose(fp1);
        //}

    }


    //fp1 = fopen("streamline_err_log.txt", "a");
    //fprintf(fp1, "finish calculating the %dth streamline.\n", evenstreamlines->ntrajs);
    //for(i=0; i<samplepts[evenstreamlines->ntrajs]->nsamples; i++)
    //{
    //	fprintf(fp1, "sample %d : (%f, %f, %f) at triangle %d.\n", i,
    //		samplepts[evenstreamlines->ntrajs]->samples[i]->gpt[0],
    //		samplepts[evenstreamlines->ntrajs]->samples[i]->gpt[1],
    //		samplepts[evenstreamlines->ntrajs]->samples[i]->gpt[2],
    //		samplepts[evenstreamlines->ntrajs]->samples[i]->triangle);
    //}
    //fclose(fp1);

    ///*write down the line segment information*/
    //fp1 = fopen("streamline_err_log.txt", "a");
    //fprintf(fp1, "\n");
    //for(i=0; i<evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs; i++)
    //{
    //	fprintf(fp1,"linesegment %d:\n", i);
    //	fprintf(fp1,"gstart (%f, %f, %f)\n",
    //		evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[i].gstart.entry[0],
    //		evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[i].gstart.entry[1],
    //		evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[i].gstart.entry[2]);
    //}
    //fclose(fp1);


    // calculate the whole length of the streamline, if it less than some threshold, discard it

    if(evenstreamlines->ntrajs == 0 ||
        evenstreamlines->trajs[evenstreamlines->ntrajs]->get_length() > streamlinelength)
    {
        update_samples_in_triangle(evenstreamlines->ntrajs, samplepts[evenstreamlines->ntrajs]->samples,
                                   samplepts[evenstreamlines->ntrajs]->nsamples);

        evenstreamlines->ntrajs ++;

        ndisplay_trajs = evenstreamlines->ntrajs-1;

        /*Testing codes*/
        //FILE *fp = fopen("num_trajs_chamber3.txt","w");

        FILE *fp;
        if(DebugOn == 1){
            //fp = fopen("num_trajs_chamber3_2.txt","w");
            //fp = fopen("num_trajs_chamber2.txt","w");
            fp = fopen("num_trajs_streamline_log.txt","w");
            //FILE *fp = fopen("num_trajs_cooling.txt","w");
            //FILE *fp = fopen("num_trajs_swirl.txt","w");
            fprintf(fp, "current number of streamlines is %d\n", evenstreamlines->ntrajs);
            fprintf(fp, "current being dealt streamline is %d\n", cur_traj);
            fprintf(fp, "current number of seed points is %d\n", seedpts->nseeds);
            for(i = 0; i < evenstreamlines->ntrajs; i++)
            {
                if(evenstreamlines->trajs[i] == NULL)
                {
                    fprintf(fp,"traj %d: #samples (%d)\n", i,
                            samplepts[i]->nsamples);
                }
                else
                    fprintf(fp,"traj %d: #lines (%d), #samples (%d)\n", i,
                            evenstreamlines->trajs[i]->nlinesegs, samplepts[i]->nsamples);
            }
            fclose(fp);
        }

        QString str;
        //int global_loc;
        int loc;

        if(DebugOn == 1)
        {
            /*
            g_finish = clock();
            str.Format(_T("\n----Streamline placement info-----\n"));
            loc=g_ctrlDebugList->GetCount();
            g_ctrlDebugList->InsertString(loc, str);
            str.Format(_T("Current number of streamlines is %d\n"), evenstreamlines->ntrajs);
            loc=g_ctrlDebugList->GetCount();
            g_ctrlDebugList->InsertString(loc, str);
            str.Format(_T("current being dealt streamline is %d\n"), cur_traj);
            loc=g_ctrlDebugList->GetCount();
            g_ctrlDebugList->InsertString(loc, str);
            str.Format(_T("current number of seed points is %d\n"), seedpts->nseeds);
            loc=g_ctrlDebugList->GetCount();
            g_ctrlDebugList->InsertString(loc, str);
            str.Format(_T("time for even placement till now is %f seconds\n"),(double)(g_finish - g_start)/CLOCKS_PER_SEC);
            loc=g_ctrlDebugList->GetCount();
            g_ctrlDebugList->InsertString(loc, str);
            str.Format(_T("see the log file for more details\n"));
            loc=g_ctrlDebugList->GetCount();
            g_ctrlDebugList->InsertString(loc, str);
            str.Format("--------------------------------------");
            loc=g_ctrlDebugList->GetCount();
            g_ctrlDebugList->InsertString(loc, str);
            */
        }

        return true;
    }

}

void EvenStreamlinePlace::cal_seeds(int traj, double dsep, int every_nsample)
{
    int i;
    double sample_p[3], newseed[3];
    int triangle, new_triangle;
    double ang = 90.;

    //Rotate current field with 90 degree
    RotateReflectField *rot = new RotateReflectField(object);
    //rot->save_old_field();
    rot->rotate_wholefield_angle(ang);

    object->recal_vec_in_triangle();

    //Get the seed points
    for(i = 0; i < samplepts[traj]->nsamples; i++)
    {
        if(i % every_nsample != 0)
        {
            if(i == samplepts[traj]->nsamples-1) //we consider the last sample point
                goto LL;
            else
                continue;
        }

    LL:		sample_p[0] = samplepts[traj]->samples[i]->gpt[0];
        sample_p[1] = samplepts[traj]->samples[i]->gpt[1];
        sample_p[2] = samplepts[traj]->samples[i]->gpt[2];

        triangle = samplepts[traj]->samples[i]->triangle;

        if(triangle < 0 || triangle >= object->tlist.ntris)
            continue;

        ////Get the forward seed, we scale the separate distance 4/18/06
        if(get_a_seed(sample_p, triangle, traj, 0, newseed, new_triangle, dsep))
        {
            if(new_triangle < 0 || new_triangle > object->tlist.ntris)
                continue;

            ////Add to the seed point list

            //seedpts->seeds[seedpts->nseeds]->pos[0] = newseed[0];
            //seedpts->seeds[seedpts->nseeds]->pos[1] = newseed[1];
            //seedpts->seeds[seedpts->nseeds]->pos[2] = newseed[2];
            //
            //seedpts->seeds[seedpts->nseeds]->triangle = new_triangle;

            //seedpts->seeds[seedpts->nseeds]->state = 0; //set it is active

            //seedpts->nseeds++;

            /*use new method to create and link the element of the seeds*/
            Seed *s = (Seed *) malloc(sizeof(Seed));
            if(s == NULL)
            {
                char rout[256], var[256];
                sprintf(rout, "%s", "EvenStreamlinePlacement::cal_seed");
                sprintf(var, "%s", "a seed");

                //write_mem_error(rout, var, 0);
                return;
            }

            s->pos[0] = newseed[0];
            s->pos[1] = newseed[1];
            s->pos[2] = newseed[2];

            s->triangle = new_triangle;
            s->state = 0; //set it is active

            seedpts->append(s);
        }

        ////Get the backward seed, we scale the separate distance
        if(get_a_seed(sample_p, triangle, traj, 1, newseed, new_triangle, dsep))
        {
            if(new_triangle < 0 || new_triangle > object->tlist.ntris)
                continue;

            ////Add to the seed point list

            //seedpts->seeds[seedpts->nseeds]->pos[0] = newseed[0];
            //seedpts->seeds[seedpts->nseeds]->pos[1] = newseed[1];
            //seedpts->seeds[seedpts->nseeds]->pos[2] = newseed[2];
            //
            //seedpts->seeds[seedpts->nseeds]->triangle = new_triangle;

            //seedpts->seeds[seedpts->nseeds]->state = 0; //set it is active

            //seedpts->nseeds++;

            /*use new method to create and link the element of the seeds*/
            Seed *s = (Seed *) malloc(sizeof(Seed));
            if(s == NULL)
            {
                char rout[256], var[256];
                sprintf(rout, "%s", "EvenStreamlinePlacement::cal_seeds");
                sprintf(var, "%s", "a seed");

                //write_mem_error(rout, var, 0);
                return;
            }

            s->pos[0] = newseed[0];
            s->pos[1] = newseed[1];
            s->pos[2] = newseed[2];

            s->triangle = new_triangle;
            s->state = 0; //set it is active

            seedpts->append(s);
        }

        ////Extend the space for seeds if needed

        if(seedpts->isFull())
        {
            int oldnum = seedpts->curMaxNumSeeds;
            if(!seedpts->extend())
            {
                /*save the error information*/
                rot->restore_cur_field(object);
                delete rot;
                return;
            }

            /*allocate the space for the elements, do we need that?*/
            //for(int i = oldnum; i < seedpts->curMaxNumSeeds; i++)
            //	seedpts->seeds[i] = new Seed[1];
        }
    }

    //Rotate the field with -90 degree
    //rot->rotate_wholefield_angle(-ang);
    rot->restore_cur_field(object);
    delete rot;

    object->recal_vec_in_triangle();
}


bool EvenStreamlinePlace::get_a_seed(double sample_p[3], int begin_triangle, int cur_traj, int type,
                                     double end_p[3], int &end_triangle, double dsep)
{
    int i;
    int flag = -1;
    double globalp[3];
    int pre_face, cur_face;
    double pre_p[3] = {0.};
    double cur_p[3] = {0.};
    double cur_length = 0;
    double smallest_dist = dsep + 1.;
    int closest_traj = -1;

    pre_face = cur_face = begin_triangle;

    globalp[0] = sample_p[0];
    globalp[1] = sample_p[1];
    globalp[2] = sample_p[2];

    tracing_points = (CurvePoints*) malloc(sizeof(CurvePoints) * 805);
    if(tracing_points == NULL)
    {
        char rout[256], var[256];
        sprintf(rout, "%s", "EvenStreamlinePlacement::get_a_seed");
        sprintf(var, "%s", "tracing_points");

        //write_mem_error(rout, var, 0);
        return false;
    }

    num_tracingpoints = 0;

    int NUMTRACETRIS = (int)sqrt((double)object->tlist.ntris);

    for(i = 0; i < NUMTRACETRIS/*NUMTRACINGTRIANGLE*/; i++)
    {
        if(cur_face < 0)
        {
            free(tracing_points);
            return false;
        }

        pre_face = cur_face;
        cur_face = trace_in_triangle_seed(cur_face, globalp, type, flag,
                                          pre_p, cur_p, dsep, cur_length, cur_traj); ////0 means always forward tracing here


        if(flag == 2 || flag == 3 || cur_face < 0 /*|| cur_face == pre_face*/) //flag = 2--reach singularity;   flag = 3--reach maximum linesegments
        {
            free(tracing_points);
            return false;
        }

        if(flag == 1) //reach the threshold, which means the length >= dsep
        {
            ////Now we need to get the exact seed point
            double extra_length = cur_length - dsep;
            cal_exact_seed(cur_p, pre_p, cur_face, extra_length, end_p);

            //end_triangle = cur_face;
            end_triangle = pre_face;  //modified at 5/8/06

            if(end_triangle < 0)
            {
                free(tracing_points);
                return false;
            }

            free(tracing_points);
            return true;
        }
    }
}


int EvenStreamlinePlace::trace_in_triangle(int &face_id, double globalp[2], int type,
                                           double dtest, double loopsep, double dist2sing,
                                           double sample_interval, double discsize, int &flag)
{
    int i;
    double alpha[3];
    double cur_point[2], pre_point[2];
    double vert0[3];
    icVector3 VP, globalv;

    Triangle *face = object->tlist.tris[face_id];

    Triangle *pre_f = face;

    ////Temporary curve point array

    CurvePoints *temp_point_list = (CurvePoints*) malloc(sizeof(CurvePoints) * 70);

    if(temp_point_list == NULL)
    {
        char rout[256], var[256];
        sprintf(rout, "%s", "EvenStreamlinePlacement::trace_in_triangle");
        sprintf(var, "%s", "temp_point_list");

        //write_mem_error(rout, var, 0);
        exit(-1);
    }

    int NumPoints = 0;

    ////initialize
    VP.entry[0] = globalp[0] - face->verts[0]->x;
    VP.entry[1] = globalp[1] - face->verts[0]->y;
    VP.entry[2] = globalp[2] - face->verts[0]->z;

    pre_point[0] = cur_point[0] = dot(VP, face->LX);
    pre_point[1] = cur_point[1] = dot(VP, face->LY);

    vert0[0] = face->verts[0]->x;   ////for update the global point
    vert0[1] = face->verts[0]->y;
    vert0[2] = face->verts[0]->z;

    //globalface = face_id;

    ////////////////////////////////////////////////////
    for(i = 0; i < 70; i++)
    {
        ////1. calculate the barycentric coordinates for current point
        object->get_2D_Barycentric_Facters(face_id, cur_point[0], cur_point[1], alpha);

        ////2. if current point is inside current triangle
        if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-10) && alpha[0] <= 1
            && (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-10) && alpha[1] <= 1
            && (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-10) && alpha[2] <= 1)
        {
            ////store the point into the temp curve points list
            temp_point_list[NumPoints].gpx = globalp[0];
            temp_point_list[NumPoints].gpy = globalp[1];
            temp_point_list[NumPoints].gpz = globalp[2];
            temp_point_list[NumPoints].lpx = cur_point[0];
            temp_point_list[NumPoints].lpy = cur_point[1];
            temp_point_list[NumPoints].triangleid = face->index;
            NumPoints++;

            pre_point[0] = cur_point[0];
            pre_point[1] = cur_point[1];

            //if(evenstreamlines->trajs[evenstreamlines->ntrajs]->cal_next_point_euler1
            //	(pre_point, cur_point, face_id, alpha, type))
            if(evenstreamlines->trajs[evenstreamlines->ntrajs]->cal_nextpt_2ndeuler
                (pre_point, cur_point, face_id, alpha, type))/*use 2nd order euler*/
            {
                globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

                globalp[0] = vert0[0] + globalv.entry[0];
                globalp[1] = vert0[1] + globalv.entry[1];
                globalp[2] = vert0[2] + globalv.entry[2];

                ////using distance to judge whether a point close to current fixed points
                if(close_to_fixedPt_except(globalp, face_id, -1, dist2sing, discsize))
                {
                    flag = 1;
                    if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs(temp_point_list, NumPoints))
                    {
                        ////Not enough memory
                        flag = 2;
                        free(temp_point_list);
                        return face_id;
                    }

                    free(temp_point_list);
                    return face_id;
                }
            }

            else{
                ////Judge whether it really reaches a fixed point
                flag = 1;

                ////Store the record into global line segment array
                if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs(temp_point_list, NumPoints))
                {
                    ////Not enough memory
                    flag = 2;
                    free(temp_point_list);
                    return face_id;
                }

                free(temp_point_list);
                return face_id;
            }


            ////Judge whether it is too close to other existing streamlines
            if(evenstreamlines->ntrajs > 0 &&
                close_to_cur_streamline(globalp, face_id, &evenstreamlines->ntrajs, 1, dtest, discsize, 0)) //scale the separate distance
            {
                if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs(temp_point_list, NumPoints))
                {
                    ////Not enough memory
                    flag = 2;
                    free(temp_point_list);
                    return face_id;
                }

                free(temp_point_list);
                flag = 3;

                return face_id;
            }

            ////We may also need to compare the current point with the sampling point on itself!
            if(close_to_cur_samplePt(globalp, face_id, samplepts[evenstreamlines->ntrajs]->samples,
                                      samplepts[evenstreamlines->ntrajs]->nsamples, loopsep, discsize, sample_interval)) //scale the separate distance
            {
                if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs(temp_point_list, NumPoints))
                {
                    ////Not enough memory
                    flag = 2;
                    free(temp_point_list);
                    return face_id;
                }

                flag = 3; //form a loop! 4/18/06

                //FILE *fp1 = fopen("streamline_err_log.txt", "a");
                //fprintf(fp1, "Form a closed loop!\n", flag);
                //fclose(fp1);


                free(temp_point_list);
                return face_id;
            }
        }

        ////3. if the point is out of current triangle
        else{
            double t[2];

            int PassVertornot = 0;
            int which_edge = -1;

            evenstreamlines->trajs[evenstreamlines->ntrajs]->get_next_triangle(face_id, pre_point,
                                                                               cur_point, t, type, PassVertornot, alpha);

            ////update the global point here
            if(PassVertornot > 0)
            {
                //we first need to know which vertex it is in the new triangle 01/30/07
                Vertex* vertid = pre_f->verts[PassVertornot-1];
                Triangle *cur_f = object->tlist.tris[face_id];
                int vert_new = 0;
                for(int k = 0; k < 3; k++)
                {
                    if(cur_f->verts[k] == vertid)
                    {
                        vert_new = k;
                        break;
                    }
                }

                alpha[vert_new]=1-0.00001;
                alpha[(vert_new+1)%3]=0.000005;
                alpha[(vert_new+2)%3]=0.000005;


                /* Get the new cur_point */
                cur_point[0] = alpha[1]*cur_f->x1+alpha[2]*cur_f->x2;
                cur_point[1] = alpha[2]*cur_f->y2;

                globalv = cur_point[0] * cur_f->LX + cur_point[1] * cur_f->LY;

                globalp[0] = cur_f->verts[0]->x + globalv.entry[0];
                globalp[1] = cur_f->verts[0]->y + globalv.entry[1];
                globalp[2] = cur_f->verts[0]->z + globalv.entry[2];
                face=cur_f;
            }

            else{
                ////transfer it to the global coordinates
                globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;

                globalp[0] = vert0[0] + globalv.entry[0];
                globalp[1] = vert0[1] + globalv.entry[1];
                globalp[2] = vert0[2] + globalv.entry[2];
            }

            ////Add the intersection point to the temporary points' list
            temp_point_list[NumPoints].gpx = globalp[0];
            temp_point_list[NumPoints].gpy = globalp[1];
            temp_point_list[NumPoints].gpz = globalp[2];
            temp_point_list[NumPoints].lpx = cur_point[0];
            temp_point_list[NumPoints].lpy = cur_point[1];
            temp_point_list[NumPoints].triangleid = face->index;  ////cause problem 05/25/05
            NumPoints++;

            if(NumPoints > 1){
                ////Store the record into global line segment array
                if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs(temp_point_list, NumPoints))
                {   ////Not enough memory
                    flag = 2;
                    free(temp_point_list);
                    return face_id;
                }
            }

            free(temp_point_list);

            return face_id;
        }
    }

    if(NumPoints > 0)
        evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs(temp_point_list, NumPoints);

    free(temp_point_list);

    return face_id;
}


SamplePt **EvenStreamlinePlace::cal_samplepts_when_tracing(int traj, double interval, int &cur_line, int &movetonext, double &cur_length,
                                                           SamplePt **samples, int &num_samples)
{
    double cur_p[2] = {0.};
    Triangle *face;
    icVector3 globalv;

    while(cur_line < evenstreamlines->trajs[traj]->nlinesegs)
    {
        if(cal_a_sample_of_streamline(traj, cur_line, movetonext,
                                       cur_p, interval, cur_length))
        {

            FILE *fp;

            /*Transform back to global coordinates*/
            face = object->tlist.tris[evenstreamlines->trajs[traj]->linesegs[cur_line].Triangle_ID];
            globalv = cur_p[0]*face->LX + cur_p[1]*face->LY;

            samplepts[traj]->samples[num_samples]->gpt[0] = face->verts[0]->x + globalv.entry[0];
            samplepts[traj]->samples[num_samples]->gpt[1] = face->verts[0]->y + globalv.entry[1];
            samplepts[traj]->samples[num_samples]->gpt[2] = face->verts[0]->z + globalv.entry[2];

            samplepts[traj]->samples[num_samples]->triangle = face->index;

            if(samplepts[traj]->samples[num_samples]->triangle < 0
                || samplepts[traj]->samples[num_samples]->triangle >= object->tlist.ntris)
            {
                continue;
            }

            num_samples++;

            if(num_samples >= samplepts[traj]->curMaxNumSamplePts)
            {

                int oldnum = samplepts[traj]->curMaxNumSamplePts;

                if(!samplepts[traj]->extend(200))
                {
                    /*record the error*/
                    char rout[256], var[256];
                    sprintf(rout, "%s", "EvenStreamlinePlacement::cal_samplepts_when_tracing");
                    sprintf(var, "%s", "samplepts[traj]");

                    //write_mem_error(rout, var, 1);
                    return NULL;
                }

                /*allocate memory for the elements*/
                for(int i = oldnum ; i < samplepts[traj]->curMaxNumSamplePts; i++)
                {
                    samplepts[traj]->samples[i] = (SamplePt *)malloc(sizeof(SamplePt));
                }

            }
        }
    }

    cur_line--;

    if(cur_line < 0)
        cur_line = 0;

    return samplepts[traj]->samples;
}


void EvenStreamlinePlace::reverse_streamline(int streamlineid)
{
    ////
    int i;
    int num_lines = evenstreamlines->trajs[streamlineid]->nlinesegs;
    Trajectory *traj = evenstreamlines->trajs[streamlineid];

    LineSeg *temp = (LineSeg *)malloc(sizeof(LineSeg)*num_lines+1);

    if(temp == NULL)
    {
        char rout[256], var[256];
        sprintf(rout, "%s", "EvenStreamlinePlacement::reverse_streamline");
        sprintf(var, "%s", "temporary linesegment list");

        //write_mem_error(rout, var, 0);
        return;
    }


    int newnum_lines = 0;

    ////store the line segment in reversed order

    for(i = num_lines-1; i >= 0; i--)
    {
        if(traj->linesegs[i].Triangle_ID < 0
            || traj->linesegs[i].Triangle_ID > object->tlist.ntris
            || traj->linesegs[i].length < 0)
        {
            continue;
        }

        temp[newnum_lines].gstart = traj->linesegs[i].gend;

        temp[newnum_lines].gend = traj->linesegs[i].gstart;

        temp[newnum_lines].start = traj->linesegs[i].end;

        temp[newnum_lines].end = traj->linesegs[i].start;

        temp[newnum_lines].length = traj->linesegs[i].length;
        temp[newnum_lines].Triangle_ID = traj->linesegs[i].Triangle_ID;


        newnum_lines++;
    }

    ////Copy it back to the origin array
    for(i = 0; i < newnum_lines; i++)
    {
        traj->linesegs[i].gstart = temp[i].gstart;

        traj->linesegs[i].gend = temp[i].gend;

        traj->linesegs[i].start = temp[i].start;

        traj->linesegs[i].end = temp[i].end;

        traj->linesegs[i].length = temp[i].length;
        traj->linesegs[i].Triangle_ID = temp[i].Triangle_ID;
    }

    traj->nlinesegs = newnum_lines;

    free(temp);
}


void EvenStreamlinePlace::update_samples_in_triangle(int traj, SamplePt **samples, int num_samples)
{
    int i;

    for(i = 0; i < num_samples; i++)
    {
        add_sample_to_triangle(samples[i]->triangle, traj, i);
    }
}

int EvenStreamlinePlace::trace_in_triangle_seed(int &face_id, double globalp[3], int type,
                                                int &flag, double pre_p[3], double cur_p[3],
                                                double dsep, double &cur_length, int cur_traj)
{
    if(face_id < 0)
        return -1;

    int i;
    double alpha[3];
    double cur_point[2], pre_point[2];
    double vert0[3];
    icVector3 VP, globalv, dist;

    Triangle *face = object->tlist.tris[face_id];
    Triangle *pre_f = face;

    Trajectory *temp_traj = new Trajectory(-1);;

    double smallest_dist = dsep+1.;
    int closest_traj = -1;


    ////initialize
    VP.entry[0] = globalp[0] - face->verts[0]->x;
    VP.entry[1] = globalp[1] - face->verts[0]->y;
    VP.entry[2] = globalp[2] - face->verts[0]->z;

    pre_point[0] = cur_point[0] = dot(VP, face->LX);
    pre_point[1] = cur_point[1] = dot(VP, face->LY);

    vert0[0] = face->verts[0]->x;   ////for update the global point
    vert0[1] = face->verts[0]->y;
    vert0[2] = face->verts[0]->z;

    //globalface = face_id;

    ////////////////////////////////////////////////////
    for(i = 0; i < 60; i++)
    {
        ////1. calculate the barycentric coordinates for current point
        object->get_2D_Barycentric_Facters(face_id, cur_point[0], cur_point[1], alpha);

        ////2. if current point is inside current triangle
        if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1
            && (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1
            && (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
        {
            ////store the point into the temp curve points list

            tracing_points[num_tracingpoints].gpx = globalp[0];
            tracing_points[num_tracingpoints].gpy = globalp[1];
            tracing_points[num_tracingpoints].gpz = globalp[2];
            tracing_points[num_tracingpoints].lpx = cur_point[0];
            tracing_points[num_tracingpoints].lpy = cur_point[1];
            tracing_points[num_tracingpoints].triangleid = face->index;

            ////Get the length for each line segment
            dist.entry[0] = cur_point[0] - pre_point[0];
            dist.entry[1] = cur_point[1] - pre_point[1];
            cur_length += tracing_points[num_tracingpoints].length = length(dist); //sum the lengths
            num_tracingpoints++;

            if(cur_length > dsep) //if the total length of the tracing curve is larger than the threshold
            {
                flag = 1;
                globalv = pre_point[0] * face->LX + pre_point[1] * face->LY;
                pre_p[0] = vert0[0] + globalv.entry[0];
                pre_p[1] = vert0[1] + globalv.entry[1];
                pre_p[2] = vert0[2] + globalv.entry[2];

                globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;
                cur_p[0] = vert0[0] + globalv.entry[0];
                cur_p[1] = vert0[1] + globalv.entry[1];
                cur_p[2] = vert0[2] + globalv.entry[2];

                delete temp_traj;
                return face_id;
            }

            ////Testing codes 4/17/06
            if(num_tracingpoints>=800) //can not get enough length
            {
                flag = 3;

                delete temp_traj;
                return face_id;
            }

            pre_point[0] = cur_point[0];
            pre_point[1] = cur_point[1];

            //if(evenstreamlines->trajs[cur_traj]->cal_next_point_euler1(pre_point, cur_point, face_id, alpha, type))
            if(temp_traj->cal_next_point_euler1(pre_point, cur_point, face_id, alpha, type))
            {
                ////update the global point
                globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

                globalp[0] = vert0[0] + globalv.entry[0];
                globalp[1] = vert0[1] + globalv.entry[1];
                globalp[2] = vert0[2] + globalv.entry[2];

            }

            else{  ////the curve reach a singularity
                flag = 2;

                delete temp_traj;
                return face_id;
            }

        }

        ////3. if the point is out of current triangle
        else{
            double t[2];

            int PassVertornot = 0;
            //int which_edge;

            //evenstreamlines->trajs[cur_traj]->get_next_triangle
            //	(face_id, pre_point, cur_point, t, type, PassVertornot, alpha);
            temp_traj->get_next_triangle
                (face_id, pre_point, cur_point, t, type, PassVertornot, alpha);

            ////update the global point here
            if(PassVertornot > 0)
            {
                //////we should not directly use the vertex as next point!!
                //////we may move a little bit along the VF direction, but make sure it is still inside
                //////the new triangle
                Vertex* vertid = pre_f->verts[PassVertornot-1];
                Triangle *cur_f = object->tlist.tris[face_id];
                int vert_new = 0;
                for(int k = 0; k < 3; k++)
                {
                    if(cur_f->verts[k] == vertid)
                    {
                        vert_new = k;
                        break;
                    }
                }

                alpha[vert_new]=1-0.0001;
                alpha[(vert_new+1)%3]=0.00005;
                alpha[(vert_new+2)%3]=0.00005;


                /* Get the new cur_point */
                cur_point[0] = alpha[1]*cur_f->x1+alpha[2]*cur_f->x2;
                cur_point[1] = alpha[2]*cur_f->y2;

                globalv = cur_point[0] * cur_f->LX + cur_point[1] * cur_f->LY;

                globalp[0] = cur_f->verts[0]->x + globalv.entry[0];
                globalp[1] = cur_f->verts[0]->y + globalv.entry[1];
                globalp[2] = cur_f->verts[0]->z + globalv.entry[2];
                face=cur_f;
            }

            else{
                globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;

                globalp[0] = vert0[0] + globalv.entry[0];
                globalp[1] = vert0[1] + globalv.entry[1];
                globalp[2] = vert0[2] + globalv.entry[2];
            }

            ////Add the intersection point to the temporary points' list
            tracing_points[num_tracingpoints].gpx = globalp[0];
            tracing_points[num_tracingpoints].gpy = globalp[1];
            tracing_points[num_tracingpoints].gpz = globalp[2];
            tracing_points[num_tracingpoints].lpx = cur_point[0];
            tracing_points[num_tracingpoints].lpy = cur_point[1];
            tracing_points[num_tracingpoints].triangleid = pre_f->index;

            ////Get the length for each line segment
            dist.entry[0] = cur_point[0] - pre_point[0];
            dist.entry[1] = cur_point[1] - pre_point[1];
            cur_length += tracing_points[num_tracingpoints].length = length(dist);
            num_tracingpoints++;

            ////if the length of the tracing curve is larger than the threshold
            if(cur_length > dsep)
            {
                flag = 1;
                globalv = pre_point[0] * pre_f->LX + pre_point[1] * pre_f->LY;
                pre_p[0] = vert0[0] + globalv.entry[0];
                pre_p[1] = vert0[1] + globalv.entry[1];
                pre_p[2] = vert0[2] + globalv.entry[2];

                globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;
                cur_p[0] = vert0[0] + globalv.entry[0];
                cur_p[1] = vert0[1] + globalv.entry[1];
                cur_p[2] = vert0[2] + globalv.entry[2];
            }

            ////Testing codes 4/17/06
            if(num_tracingpoints>=800) //can not get enough length
            {
                flag = 3;

                delete temp_traj;
                return face_id;
            }

            delete temp_traj;
            return face_id;
        }

    }

    flag = 3;   //can not go outside of the triangle

    delete temp_traj;
    return face_id;
}


void EvenStreamlinePlace::cal_exact_seed(double cur_p[3], double pre_p[3], int triangle,
                                         double extra_length, double exact_seed[3])
{
    icVector3 line;
    line.entry[0] = cur_p[0] - pre_p[0];
    line.entry[1] = cur_p[1] - pre_p[1];
    line.entry[2] = cur_p[2] - pre_p[2];

    double t = 1-(extra_length/length(line)); //it seems that this can get more even separate distance

    exact_seed[0] = (1-t) * pre_p[0] + t * cur_p[0];
    exact_seed[1] = (1-t) * pre_p[1] + t * cur_p[1];
    exact_seed[2] = (1-t) * pre_p[2] + t * cur_p[2];
}


bool EvenStreamlinePlace::close_to_fixedPt_except(double p[3], int triangle,
                                                  int singid, double threshold, double discsize)
{
    int i;
    int *NearbyTriangles = NULL;
    int num_triangles = 0;
    Triangle *face;
    Vertex *v;
    icVector3 VP, v0;
    double dist, a, b, alpha[3];

    ////For unfolding
    int edgeid;
    double rotmat[16] = {0.};
    icVector3 pt;
    set_Indentity16(rotmat);

    ////Get the triangles that fall into the circle region with p as the center and 3*separate_dist as radius
    //NearbyTriangles = cal_geodesic_dist(triangle, p, threshold, discsize, NearbyTriangles, num_triangles);

    //if(trianglelist == NULL)
    trianglelist = new DynList_Int();
    trianglelist->nelems = 0;
    cal_geodesic_dist_2(triangle, p, threshold, discsize, trianglelist);
    NearbyTriangles = trianglelist->elems;
    num_triangles = trianglelist->nelems;

    for(i = 0; i < num_triangles; i++)
    {
        face = object->tlist.tris[NearbyTriangles[i]];

        if(face->singularityID < 0) //this triangle does not contain fixed point
            continue;

        if(face->singularityID == singid)
            continue;

        ////For the center 4 triangles we unfold them and use the Euclidean distance directly
        if(i > 0 && i <= 3)
        {
            edgeid = i-1; //according to the order of how we add the first 4 triangles
            get_unfold_rotMat(NearbyTriangles[0], edgeid, NearbyTriangles[i], rotmat);
        }

        if(i <= 3)
        {
            pt= object->slist.slist[face->singularityID]->gpos;
            transform_point3D(pt.entry, rotmat);

            VP.entry[0] = pt.entry[0] - p[0];
            VP.entry[1] = pt.entry[1] - p[1];
            VP.entry[2] = pt.entry[2] - p[2];

            dist = length(VP);

        }

        else
        {
            v0.entry[0] = face->verts[0]->x;
            v0.entry[1] = face->verts[0]->y;
            v0.entry[2] = face->verts[0]->z;

            VP = object->slist.slist[face->singularityID]->gpos - v0;

            ////Calculate the barycentric coordinates for the fixed point
            a = dot(VP, face->LX);
            b = dot(VP, face->LY);
            object->get_2D_Barycentric_Facters(NearbyTriangles[i], a, b, alpha);

            dist = alpha[0]*face->verts[0]->distance
                   + alpha[1]*face->verts[1]->distance
                   + alpha[2]*face->verts[2]->distance;
        }

        if(dist <= threshold)
        {
            //reset_dist(NearbyTriangles,num_triangles);
            //free(NearbyTriangles);

            reset_dist(trianglelist);
            delete trianglelist;
            return true;
        }
    }

    //reset_dist(NearbyTriangles,num_triangles);
    //free(NearbyTriangles);

    reset_dist(trianglelist);
    delete trianglelist;
    return false;
}


bool EvenStreamlinePlace::close_to_cur_streamline(double p[3], int triangle, int *Except_trajs, int num_trajs,
                                                  double separate_dist, double discsize, int Init)
{
    int *NearbyTriangles = NULL;
    int num_triangles = 0;

    int i, j, k;
    int traj;
    Triangle *face;
    Vertex *v1, *v2, *v3;
    double dis, alpha[3];
    int sampleid;

    icVector3 VP;

    ////For unfolding 5/16/06
    int edgeid;
    double rotmat[16] = {0.};
    double pt[3] = {0.};
    set_Indentity16(rotmat);

    ////first, we get the nearby triangles and propogate the distance to the nearby vertices
    //NearbyTriangles = cal_geodesic_dist(triangle, p, separate_dist, discsize, NearbyTriangles, num_triangles);

    //if(trianglelist == NULL)
    trianglelist = new DynList_Int();
    trianglelist->nelems = 0;
    cal_geodesic_dist_2(triangle, p, separate_dist, discsize, trianglelist);
    NearbyTriangles = trianglelist->elems;
    num_triangles = trianglelist->nelems;

    ////then we test all the curve points in the nearby triangles
    for(i = 0; i < num_triangles; i++)
    {
        face = object->tlist.tris[NearbyTriangles[i]];

        v1 = face->verts[0];
        v2 = face->verts[1];
        v3 = face->verts[2];

        if(face->num_samplepts <= 0)
            continue;

        ////For the center 4 triangles we unfold them and use the Euclidean distance directly
        if(i > 0 && i <= 3)
        {
            edgeid = i-1; //according to the order of how we add the first 4 triangles
            get_unfold_rotMat(NearbyTriangles[0], edgeid, NearbyTriangles[i], rotmat);
        }
        /////////

        for(j = 0; j < face->num_samplepts; j++)
        {
            traj = face->samplepts[j].which_traj;

            if(is_repeated_elem(Except_trajs, traj, num_trajs))
                continue;

            //if(Init == 1 /*&& SeparatrixorNot[traj] == 1*/) /*need to judge whether it is a separatrix*/
            //	continue;

            sampleid = face->samplepts[j].which_sample;

            ////For the center 4 triangles we unfold them and use the Euclidean distance directly

            if(i <= 3)/*for the 4 center triangles, we use the Euclidean distance directly*/
            {
                pt[0] = samplepts[traj]->samples[sampleid]->gpt[0];
                pt[1] = samplepts[traj]->samples[sampleid]->gpt[1];
                pt[2] = samplepts[traj]->samples[sampleid]->gpt[2];
                transform_point3D(pt, rotmat);

                VP.entry[0] = pt[0] - p[0];
                VP.entry[1] = pt[1] - p[1];
                VP.entry[2] = pt[2] - p[2];

                dis = length(VP);
            }


            else
            {
                ////Get the bary centric coordinates of the curve point
                VP.entry[0] = samplepts[traj]->samples[sampleid]->gpt[0] - v1->x;
                VP.entry[1] = samplepts[traj]->samples[sampleid]->gpt[1] - v1->y;
                VP.entry[2] = samplepts[traj]->samples[sampleid]->gpt[2] - v1->z;

                double a = dot(VP, face->LX);
                double b = dot(VP, face->LY);

                object->get_2D_Barycentric_Facters(NearbyTriangles[i], a, b, alpha);

                ////Get the distance through the distances on the three vertices
                dis = alpha[0]*v1->distance + alpha[1]*v2->distance + alpha[2]*v3->distance;
            }

            if(dis <= separate_dist) //scale the distance a little bit
            {
                //reset_dist(NearbyTriangles,num_triangles);
                //free(NearbyTriangles);

                reset_dist(trianglelist);
                delete trianglelist;

                return true;
            }
        }
    }

    if(NearbyTriangles != NULL)
    {
        //reset_dist(NearbyTriangles,num_triangles);
        //free(NearbyTriangles);

        reset_dist(trianglelist);
        delete trianglelist;
    }
    return false;
}


bool EvenStreamlinePlace::close_to_cur_samplePt(double p[3], int triangle, SamplePt **samples, int num_samples,
                                                double separate_dist, double discsize, double sample_interval)
{
    int i;
    double dist, alpha[3], a, b;
    int stop_locate = floor(separate_dist/sample_interval)+3;
    Triangle *face;
    Vertex *v;
    icVector3 VP, v0;

    int *NearbyTriangles = NULL;
    int num_triangles = 0;


    ////For unfolding 5/16/06
    int edgeid;
    double rotmat[16] = {0.};
    double pt[3] = {0.};
    set_Indentity16(rotmat);
    int unfold_flag = 0;

    /*really temporary testing here!!*/
    //if(num_samples >=13 && triangle == 101)
    //{
    //	int test = 0;
    //}

    ////Get the triangles that fall into the circle region with p as the center and 3*separate_dist as radius
    //NearbyTriangles = cal_geodesic_dist(triangle, p, separate_dist, discsize, NearbyTriangles, num_triangles);


    //if(trianglelist == NULL)
    trianglelist = new DynList_Int();
    trianglelist->nelems = 0;
    cal_geodesic_dist_2(triangle, p, separate_dist, discsize, trianglelist);
    NearbyTriangles = trianglelist->elems;
    num_triangles = trianglelist->nelems;


    for(i = 0 ; i < num_samples-stop_locate; i++)
    {
        unfold_flag = 0;

        /*really temporary testing*/
        //double pr[3]={0.610012, 0.124964, 0.015919};
        //icVector3 dis_p;
        //dis_p.entry[0] = pr[0] - p[0];
        //dis_p.entry[1] = pr[1] - p[1];
        //dis_p.entry[2] = pr[2] - p[2];
        //double len = length(dis_p);

        //if(i == 5 && num_samples >=13 && triangle == 101 /*&& length(dis_p) < 1e-8*/)
        //{
        //	FILE *ff = fopen("neighboring_triangles.txt", "w");
        //	fprintf(ff, "the point is (%f, %f, %f)\n",
        //		p[0], p[1], p[2]);
        //	for(int k=0; k<num_triangles; k++)
        //	{
        //		fprintf(ff, "%d\n", NearbyTriangles[k]);
        //		//for(int l=0; l<3; l++)
        //			fprintf(ff, "distances on vertices: (%f, %f, %f)\n",
        //			object->tlist.tris[NearbyTriangles[k]]->verts[0]->distance,
        //			object->tlist.tris[NearbyTriangles[k]]->verts[1]->distance,
        //			object->tlist.tris[NearbyTriangles[k]]->verts[2]->distance);
        //	}
        //	fclose(ff);
        //}

        if(!is_repeated_elem(NearbyTriangles, samples[i]->triangle, num_triangles))
            continue;

        face = object->tlist.tris[samples[i]->triangle];

        ////For the center 4 triangles we unfold them and use the Euclidean distance directly
        /*probably need to consider one-ring neighboring triangles*/
        if(samples[i]->triangle == NearbyTriangles[1])
        {
            edgeid = 0; //according to the order of how we add the first 4 triangles
            get_unfold_rotMat(NearbyTriangles[0], edgeid, NearbyTriangles[1], rotmat);
            unfold_flag = 1;
        }
        else if(samples[i]->triangle == NearbyTriangles[2])
        {
            edgeid = 1; //according to the order of how we add the first 4 triangles
            get_unfold_rotMat(NearbyTriangles[0], edgeid, NearbyTriangles[2], rotmat);
            unfold_flag = 1;
        }
        else if(samples[i]->triangle == NearbyTriangles[3])
        {
            edgeid = 2; //according to the order of how we add the first 4 triangles
            get_unfold_rotMat(NearbyTriangles[0], edgeid, NearbyTriangles[3], rotmat);
            unfold_flag = 1;
        }
        else if(samples[i]->triangle == NearbyTriangles[0])
        {
            set_Indentity16(rotmat);
            unfold_flag = 1;
        }

        if(unfold_flag == 1)
        {
            pt[0] = samples[i]->gpt[0];
            pt[1] = samples[i]->gpt[1];
            pt[2] = samples[i]->gpt[2];
            transform_point3D(pt, rotmat);

            VP.entry[0] = pt[0] - p[0];
            VP.entry[1] = pt[1] - p[1];
            VP.entry[2] = pt[2] - p[2];

            dist = length(VP);
        }

        else
        {
            v0.entry[0] = face->verts[0]->x;
            v0.entry[1] = face->verts[0]->y;
            v0.entry[2] = face->verts[0]->z;

            VP.entry[0] = samples[i]->gpt[0] - v0.entry[0];
            VP.entry[1] = samples[i]->gpt[1] - v0.entry[1];
            VP.entry[2] = samples[i]->gpt[2] - v0.entry[2];

            //Calculate the barycentric coordinates for the sample point
            a = dot(VP, face->LX);
            b = dot(VP, face->LY);
            object->get_2D_Barycentric_Facters(samples[i]->triangle, a, b, alpha);

            dist = alpha[0]*face->verts[0]->distance
                   + alpha[1]*face->verts[1]->distance
                   + alpha[2]*face->verts[2]->distance;
        }

        if(dist < separate_dist)
        {
            //reset_dist(NearbyTriangles,num_triangles);
            //free(NearbyTriangles);

            reset_dist(trianglelist);
            delete trianglelist;

            //FILE *fp1 = fopen("streamline_err_log.txt", "a");
            //fprintf(fp1, "at triangle %d of sample %d\n", face->index, i);
            //fprintf(fp1, "dist: %f, separate_dist: %f\n", dist, separate_dist);
            //for(int k=0; k<16; k++)
            //	fprintf(fp1, "rotMat: %f\n", rotmat[k]);
            //fclose(fp1);

            //g_sample_tri = face->index;
            //g_dist = dist;
            //g_separate_dist = separate_dist;

            return true;
        }
        //g_dist = dist;
        //g_separate_dist = separate_dist;
    }

    //reset_dist(NearbyTriangles,num_triangles);
    //free(NearbyTriangles);

    reset_dist(trianglelist);
    delete trianglelist;
    return false;
}


void EvenStreamlinePlace::reset_dist(int *NearbyTriangles, int num)
{
    int i, j;
    Triangle *face;
    for(i = 0; i < num; i++)
    {
        face = object->tlist.tris[NearbyTriangles[i]];
        face->visited = false;
        for(j = 0; j < 3; j++)
            face->verts[j]->distance = 1e49;
    }
}

void EvenStreamlinePlace::reset_dist(DynList_Int *trianglelist)
{
    int i, j;
    Triangle *face;
    for(i = 0; i < trianglelist->nelems; i++)
    {
        face = object->tlist.tris[trianglelist->elems[i]];
        face->visited = false;
        for(j = 0; j < 3; j++)
            face->verts[j]->distance = 1e49;
    }
}


bool EvenStreamlinePlace::cal_a_sample_of_streamline(int traj, int &cur_lineindex, int &movetonext,
                                                     double curpt[2], double interval, double &cur_length)
{
    int i;
    icVector2 len_vec;
    double alpha;

    int num_lines = evenstreamlines->trajs[traj]->nlinesegs;

    curpt[0] = curpt[1] = -1;

    if(cur_length >= interval)
    {
        alpha = (cur_length-interval)/evenstreamlines->trajs[traj]->linesegs[cur_lineindex].length;
        curpt[0] = alpha*evenstreamlines->trajs[traj]->linesegs[cur_lineindex].start.entry[0]
                   + (1-alpha)*evenstreamlines->trajs[traj]->linesegs[cur_lineindex].end.entry[0];
        curpt[1] = alpha*evenstreamlines->trajs[traj]->linesegs[cur_lineindex].start.entry[1]
                   + (1-alpha)*evenstreamlines->trajs[traj]->linesegs[cur_lineindex].end.entry[1];

        cur_length -= interval;
        return true;
    }

    else{
        cur_lineindex++;
        if(cur_lineindex >= num_lines)
        {
            return false;
        }
        cur_length += evenstreamlines->trajs[traj]->linesegs[cur_lineindex].length;
        return false;
    }
}

void EvenStreamlinePlace::add_sample_to_triangle(int triangle, int which_traj, int which_sample)
{
    if(triangle < 0 || triangle > object->tlist.ntris)
        return;

    Triangle *face = object->tlist.tris[triangle];

    face->samplepts = face->extend_sampleList();

    face->samplepts[face->num_samplepts].which_traj = which_traj;
    face->samplepts[face->num_samplepts].which_sample = which_sample;

    face->num_samplepts++;
}


int *EvenStreamlinePlace::cal_geodesic_dist(int triangle, double globalp[3], double dsep, double discsize,
                                            int *trianglelist, int &num_triangles)
{
    ////we suppose the global point always falls in the triangle
    Triangle *face = object->tlist.tris[triangle];
    Vertex *v1, *v2, *v3/*, *cur_v*/;
    Edge *cur_edge;
    double cur_dis;
    icVector3 dis_vec;
    int other_f, third_vert;

    ////initial the triangle list  (this could be the bottleneck 04/30/07)
    trianglelist = extend_link(trianglelist, num_triangles);
    trianglelist[num_triangles]=triangle;
    num_triangles++;

    int i, j, k;

    ////Reset the "distance" of all vertices
    //for(i = 0; i < object->vlist.nverts; i++)
    //{
    //	object->vlist.verts[i]->distance = 1e49;
    //}

    ////Calculate the distance for the 3 vertices of the triangle
    v1 = face->verts[0];
    v2 = face->verts[1];
    v3 = face->verts[2];

    dis_vec.entry[0] = v1->x - globalp[0];
    dis_vec.entry[1] = v1->y - globalp[1];
    dis_vec.entry[2] = v1->z - globalp[2];
    cur_dis = v1->distance = length(dis_vec);

    /**--------------------------------------------------------*/
    ////Testing codes
    //maxdis = mindis = cur_dis;

    dis_vec.entry[0] = v2->x - globalp[0];
    dis_vec.entry[1] = v2->y - globalp[1];
    dis_vec.entry[2] = v2->z - globalp[2];
    v2->distance = length(dis_vec);

    /**--------------------------------------------------------*/
    ////Testing codes
    //if(v2->distance < mindis) mindis = v2->distance;
    //if(v2->distance > maxdis) maxdis = v2->distance;

    dis_vec.entry[0] = v3->x - globalp[0];
    dis_vec.entry[1] = v3->y - globalp[1];
    dis_vec.entry[2] = v3->z - globalp[2];
    v3->distance = length(dis_vec);

    /**--------------------------------------------------------*/
    ////Testing codes
    //if(v2->distance < mindis) mindis = v3->distance;
    //if(v2->distance > maxdis) maxdis = v3->distance;

    ////Propogate the distance using "Fast Marching method" until the distance larger than some threshold

    int previous_position = 0;
    int prev_num_triangles = num_triangles;

    ////Grow the three neighboring triangles first in whatever cases 5/16/06
    face = object->tlist.tris[triangle];
    for(i = 0; i < 3; i++)
    {
        cur_edge = face->edges[i];

        if(cur_edge->tris[0]!=NULL)
        {
            other_f = cur_edge->tris[0]->index;
            if(other_f == triangle && cur_edge->tris[1] == NULL)
                continue;
            else if(other_f == triangle && cur_edge->tris[1] !=NULL)
                other_f = cur_edge->tris[1]->index;
        }
        else if(cur_edge->tris[1] !=NULL && cur_edge->tris[1]->index != triangle)
            other_f=cur_edge->tris[1]->index;
        else
            continue;

        trianglelist = extend_link(trianglelist, num_triangles);
        trianglelist[num_triangles]=other_f;
        num_triangles++;

        ////Get the two ending points of the current edge
        v1 = cur_edge->verts[0];
        v2 = cur_edge->verts[1];

        third_vert = get_thirdVer_of_triangle(v1->index,
                                              v2->index, other_f);

        ////Calculate the distance
        if(v1->distance < v2->distance)
        {
            update_one_vertex(v1->index, v2->index, third_vert, other_f);
        }
        else
        {
            update_one_vertex(v2->index, v1->index, third_vert, other_f);
        }
    }
    prev_num_triangles = num_triangles;

    ////Get the disc from the center 4 triangles
    while(previous_position < prev_num_triangles && num_triangles < object->tlist.ntris)
    {
        set_boundary_edges(trianglelist, num_triangles);

        ////Check all the boundary edges,
        ////if there is one edge whose vertex has a distance smaller than the threshold
        ////then, grow one step to its opposite triangle and record the new triangle

        for(i = previous_position; i < prev_num_triangles; i++)
        {
            face = object->tlist.tris[trianglelist[i]];

            for(j = 0; j < 3; j++)
            {
                cur_edge = face->edges[j];

                if(!cur_edge->OnBoundary)
                    continue;

                ////get the two ending points of the current edge
                v1 = cur_edge->verts[0];
                v2 = cur_edge->verts[1];

                if(v1->distance > dsep*discsize && v2->distance > dsep*discsize)
                    continue;

                if(cur_edge->tris[0]!=NULL)
                {
                    other_f = cur_edge->tris[0]->index;
                    if(other_f == trianglelist[i] && cur_edge->tris[1] == NULL)
                        continue;
                    else if(other_f == trianglelist[i] && cur_edge->tris[1] !=NULL)
                        other_f = cur_edge->tris[1]->index;
                }
                else if(cur_edge->tris[1] !=NULL && cur_edge->tris[1]->index != trianglelist[i])
                    other_f=cur_edge->tris[1]->index;
                else
                    continue;

                //other_f = cur_edge->tris[0]->index;

                //if(other_f == trianglelist[i])
                //	other_f = cur_edge->tris[1]->index;

                ////If the triangle is not inside the region, add it to the region
                //if(!is_repeated_elem(trianglelist, other_f, num_triangles))
                if(!object->tlist.tris[other_f]->visited)
                {
                    trianglelist = extend_link(trianglelist, num_triangles);
                    trianglelist[num_triangles]=other_f;
                    num_triangles++;
                    object->tlist.tris[other_f]->visited = true;
                }

                ////Get the index of the third point
                for(k = 0; k < 3; k++)
                {
                    third_vert = object->tlist.tris[other_f]->verts[k]->index;

                    if(third_vert != v1->index && third_vert != v2->index)
                        break;
                }

                ////In the following codes, we suppose the distance at third vertex is alwayes
                ////larger than v1's and v2's
                if(v1->distance < v2->distance)
                {
                    update_one_vertex(v1->index, v2->index, third_vert, other_f);
                }
                else
                {
                    update_one_vertex(v2->index, v1->index, third_vert, other_f);
                }

                /**--------------------------------------------------------*/
                ////Testing codes for visualization
                //if(object->vlist.verts[third_vert]->distance < mindis)
                //	mindis = object->vlist.verts[third_vert]->distance;
                //if(object->vlist.verts[third_vert]->distance > maxdis)
                //	maxdis = object->vlist.verts[third_vert]->distance;
            }
        }

        previous_position = prev_num_triangles;
        prev_num_triangles = num_triangles;
    }

////Try to make up the gap
LL:	set_boundary_edges(trianglelist, num_triangles);
    for(i = 0; i < num_triangles; i++)
    {
        face = object->tlist.tris[trianglelist[i]];

        for(j = 0; j < 3; j++)
        {
            cur_edge = face->edges[j];

            if(!cur_edge->OnBoundary)
                continue;

            //other_f = cur_edge->tris[0]->index;

            //if(other_f == trianglelist[i])
            //	other_f = cur_edge->tris[1]->index;

            if(cur_edge->tris[0]!=NULL)
            {
                other_f = cur_edge->tris[0]->index;
                if(other_f == trianglelist[i] && cur_edge->tris[1] == NULL)
                    continue;
                else if(other_f == trianglelist[i] && cur_edge->tris[1] !=NULL)
                    other_f = cur_edge->tris[1]->index;
            }
            else if(cur_edge->tris[1] !=NULL && cur_edge->tris[1]->index != trianglelist[i])
                other_f=cur_edge->tris[1]->index;
            else
                continue;

            for(k = 0; k < 3; k++)
            {
                v1 = object->tlist.tris[other_f]->verts[k];

                if(v1->distance > 1e47)
                    break;
            }

            if(k >= 3)
            {
                //if(!is_repeated_elem(trianglelist, other_f, num_triangles))
                if(!object->tlist.tris[other_f]->visited)
                {
                    trianglelist = extend_link(trianglelist, num_triangles);
                    trianglelist[num_triangles]=other_f;
                    num_triangles++;
                    object->tlist.tris[other_f]->visited = true;
                }
            }
        }
    }


    return trianglelist;
}



/*
Use new data structure to save the geodesic disc
*/
void EvenStreamlinePlace::cal_geodesic_dist_2(int triangle, double globalp[3], double dsep, double discsize,
                                              DynList_Int *trianglelist)
{
    ////we suppose the global point always falls in the triangle
    Triangle *face = object->tlist.tris[triangle];
    Vertex *v1, *v2, *v3;
    Edge *cur_edge;
    double cur_dis;
    icVector3 dis_vec;
    int other_f, third_vert;

    /*try to use more efficient method here*/
    trianglelist->nelems = 0;
    trianglelist->add_New_2(triangle);
    face->visited = true;

    int i, j, k;

    ////Reset the "distance" of all vertices
    //for(i = 0; i < object->vlist.nverts; i++)
    //{
    //	object->vlist.verts[i]->distance = 1e49;
    //}

    ////Calculate the distance for the 3 vertices of the triangle
    v1 = face->verts[0];
    v2 = face->verts[1];
    v3 = face->verts[2];

    dis_vec.entry[0] = v1->x - globalp[0];
    dis_vec.entry[1] = v1->y - globalp[1];
    dis_vec.entry[2] = v1->z - globalp[2];
    cur_dis = v1->distance = length(dis_vec);

    /**--------------------------------------------------------*/

    dis_vec.entry[0] = v2->x - globalp[0];
    dis_vec.entry[1] = v2->y - globalp[1];
    dis_vec.entry[2] = v2->z - globalp[2];
    v2->distance = length(dis_vec);

    /**--------------------------------------------------------*/

    dis_vec.entry[0] = v3->x - globalp[0];
    dis_vec.entry[1] = v3->y - globalp[1];
    dis_vec.entry[2] = v3->z - globalp[2];
    v3->distance = length(dis_vec);

    /**--------------------------------------------------------*/

    ////Propogate the distance using "Fast Marching method" until the distance larger than some threshold

    int previous_position = 0;
    int prev_num_triangles = trianglelist->nelems;

    ////Grow the three neighboring triangles first in whatever cases 5/16/06
    face = object->tlist.tris[triangle];
    for(i = 0; i < 3; i++)
    {
        cur_edge = face->edges[i];

        if(cur_edge->tris[0]!=NULL)
        {
            other_f = cur_edge->tris[0]->index;
            if(other_f == triangle && cur_edge->tris[1] !=NULL)
                other_f = cur_edge->tris[1]->index;
            else
                continue;
        }
        else if(cur_edge->tris[1] !=NULL && cur_edge->tris[1]->index != triangle)
            other_f=cur_edge->tris[1]->index;
        else
            continue;

        //trianglelist = extend_link(trianglelist, num_triangles);
        //trianglelist[num_triangles]=other_f;
        //num_triangles++;
        trianglelist->add_New_2(other_f);
        object->tlist.tris[other_f]->visited=true;

        ////Get the two ending points of the current edge
        v1 = cur_edge->verts[0];
        v2 = cur_edge->verts[1];

        third_vert = get_thirdVer_of_triangle(cur_edge->verts[0]->index,
                                              cur_edge->verts[1]->index, other_f);

        ////Calculate the distance
        if(v1->distance < v2->distance)
        {
            update_one_vertex(v1->index, v2->index, third_vert, other_f);
        }
        else
        {
            update_one_vertex(v2->index, v1->index, third_vert, other_f);
        }
    }
    prev_num_triangles = trianglelist->nelems;

    ////Get the disc from the center 4 triangles
    while(previous_position < prev_num_triangles && trianglelist->nelems < object->tlist.ntris)
    {
        set_boundary_edges(trianglelist->elems, trianglelist->nelems);

        ////Check all the boundary edges,
        ////if there is one edge whose vertex has a distance smaller than the threshold
        ////then, grow one step to its opposite triangle and record the new triangle

        for(i = previous_position; i < prev_num_triangles; i++)
        {
            face = object->tlist.tris[trianglelist->elems[i]];

            for(j = 0; j < 3; j++)
            {
                cur_edge = face->edges[j];

                if(!cur_edge->OnBoundary)
                    continue;

                ////get the two ending points of the current edge
                v1 = cur_edge->verts[0];
                v2 = cur_edge->verts[1];

                if(v1->distance > dsep*discsize && v2->distance > dsep*discsize)
                    continue;

                //other_f = cur_edge->tris[0]->index;

                //if(other_f == trianglelist->elems[i])
                //	other_f = cur_edge->tris[1]->index;

                if(cur_edge->tris[0]!=NULL)
                {
                    other_f = cur_edge->tris[0]->index;
                    if(other_f == trianglelist->elems[i] && cur_edge->tris[1] == NULL)
                        continue;
                    else if(other_f == trianglelist->elems[i] && cur_edge->tris[1] !=NULL)
                        other_f = cur_edge->tris[1]->index;
                }
                else if(cur_edge->tris[1] !=NULL && cur_edge->tris[1]->index !=
                                                          trianglelist->elems[i])
                    other_f=cur_edge->tris[1]->index;
                else
                    continue;

                ////If the triangle is not inside the region, add it to the region
                //if(!is_repeated_elem(trianglelist->elems, other_f, trianglelist->nelems))
                //	trianglelist->add_New(other_f);

                if(!object->tlist.tris[other_f]->visited){
                    trianglelist->add_New_2(other_f);
                    object->tlist.tris[other_f]->visited=true;
                }

                ////Get the index of the third point
                for(k = 0; k < 3; k++)
                {
                    third_vert = object->tlist.tris[other_f]->verts[k]->index;

                    if(third_vert != v1->index && third_vert != v2->index)
                        break;
                }

                ////In the following codes, we suppose the distance at third vertex is alwayes
                ////larger than v1's and v2's
                if(v1->distance < v2->distance)
                {
                    update_one_vertex(v1->index, v2->index, third_vert, other_f);
                }
                else
                {
                    update_one_vertex(v2->index, v1->index, third_vert, other_f);
                }

            }
        }

        previous_position = prev_num_triangles;
        prev_num_triangles = trianglelist->nelems;
    }

////Try to make up the gap
LL:	set_boundary_edges(trianglelist->elems, trianglelist->nelems);
    for(i = 0; i < trianglelist->nelems; i++)
    {
        face = object->tlist.tris[trianglelist->elems[i]];

        for(j = 0; j < 3; j++)
        {
            cur_edge = face->edges[j];

            if(!cur_edge->OnBoundary/* == false*/)
                continue;

            //other_f = cur_edge->tris[0]->index;

            //if(other_f == trianglelist->elems[i])
            //	other_f = cur_edge->tris[1]->index;

            if(cur_edge->tris[0]!=NULL)
            {
                other_f = cur_edge->tris[0]->index;
                if(other_f == trianglelist->elems[i] && cur_edge->tris[1] == NULL)
                    continue;
                else if(other_f == trianglelist->elems[i] && cur_edge->tris[1] !=NULL)
                    other_f = cur_edge->tris[1]->index;
            }
            else if(cur_edge->tris[1] !=NULL && cur_edge->tris[1]->index != trianglelist->elems[i])
                other_f=cur_edge->tris[1]->index;
            else
                continue;

            for(k = 0; k < 3; k++)
            {
                v1 = object->tlist.tris[other_f]->verts[k];

                if(v1->distance > 1e47)
                    break;
            }

            if(k >= 3)
            {
                //if(!is_repeated_elem(trianglelist->elems, other_f, trianglelist->nelems))
                //	trianglelist->add_New(other_f);

                if(!object->tlist.tris[other_f]->visited)
                {
                    trianglelist->add_New_2(other_f);
                    object->tlist.tris[other_f]->visited=true;
                }
            }
        }
    }

    //return trianglelist;

    /*Testing codes*/
    //FILE *fp = fopen("cooling_test.txt", "a");
    //fprintf(fp, "%d triangles have been found.\n", trianglelist->nelems);
    //fclose(fp);

    for(i=0; i<trianglelist->nelems; i++)
    {
        object->tlist.tris[trianglelist->elems[i]]->visited = false;
    }
}

void EvenStreamlinePlace::get_unfold_rotMat(int org_triangle, int edgeid, int other_f, double rotmat[16])
{
    Edge *cur_e;
    icVector3 rot_axis,  normal1, normal2, local_LX, local_LY;
    double base[3] = {0.};
    double rot_ang;

    cur_e = object->tlist.tris[org_triangle]->edges[edgeid];

    /*we need to make sure they share the same edge!*/
    int count = 1;
    //LL:	if(!(cur_e->tris[0]->index == org_triangle&&cur_e->tris[1]->index == other_f)
    //		&& !(cur_e->tris[0]->index==other_f && cur_e->tris[1]->index == org_triangle))
    //	{
    //		cur_e = object->tlist.tris[org_triangle]->edges[(edgeid+count)%3];
    //		count++;
    //		if(count > 2) return; //error!!!
    //		goto LL;
    //	}

    //Get rotation axis based on the edge 'cur_e'
    rot_axis.entry[0] = cur_e->verts[1]->x - cur_e->verts[0]->x;
    rot_axis.entry[1] = cur_e->verts[1]->y - cur_e->verts[0]->y;
    rot_axis.entry[2] = cur_e->verts[1]->z - cur_e->verts[0]->z;
    normalize(rot_axis);

    base[0] = cur_e->verts[0]->x;
    base[1] = cur_e->verts[0]->y;
    base[2] = cur_e->verts[0]->z;

    //Get the two normals
    normal1 = object->tlist.tris[org_triangle]->normal;
    normal2 = object->tlist.tris[other_f]->normal;

    //Get the angle between the normals of cur_triangle and other_f
    local_LX = normal1;
    local_LY = cross(rot_axis, normal1);

    double a = dot(normal2, local_LX);
    double b = dot(normal2, local_LY);

    rot_ang = atan2(b, a);

    //rot_ang = rot_ang * 180/3.1415926;

    //Get rotation matrix
    //get_rotation(base, rot_axis, -rot_ang, rotmat);
    //g_pclGlView->get_rotation(base, rot_axis, -rot_ang, rotmat);
    //g_pCgraphView->get_rotation(base, rot_axis, -rot_ang, rotmat);

    get_rotation(base, rot_axis, -rot_ang, rotmat);
}



int EvenStreamlinePlace::get_thirdVer_of_triangle(int v1, int v2, int triangle)
{
    Triangle *face = object->tlist.tris[triangle];

    for(int i = 0; i < face->nverts; i++)
    {
        if(face->verts[i]->index != v1 && face->verts[i]->index != v2)
            return face->verts[i]->index;
    }

    return -1; //something is wrong here
}



void EvenStreamlinePlace::update_one_vertex(int verta, int vertb, int vertc, int triangle)
{
    Vertex *va, *vb, *vc;
    va = object->vlist.verts[verta];
    vb = object->vlist.verts[vertb];
    vc = object->vlist.verts[vertc];
    icVector3 veca, vecb;

    double lena, lenb;
    veca.entry[0] = va->x - vc->x;
    veca.entry[1] = va->y - vc->y;
    veca.entry[2] = va->z - vc->z;
    lenb = length(veca);

    vecb.entry[0] = vb->x - vc->x;
    vecb.entry[1] = vb->y - vc->y;
    vecb.entry[2] = vb->z - vc->z;
    lena = length(vecb);

    double theta, sin_theta, cos_theta;
    normalize(veca);
    normalize(vecb);
    theta = acos(dot(veca, vecb));
    sin_theta = sin(theta);
    cos_theta = dot(veca, vecb);

    if(dot(veca, vecb) < 0)
    {
        get_dis_for_obtuse_triangle(verta, vertb, vertc, triangle, lenb, lena);
        return;
    }

    double u = vb->distance - va->distance; ////here we may suppose u>=0

    ////
    double A, B, C, delta, t1, t2;
    A = lena*lena + lenb*lenb - 2*lena*lenb*cos_theta;
    B = 2*lenb*u*(lena*cos_theta-lenb);
    C = lenb*lenb*(u*u - lena*lena*sin_theta*sin_theta);

    delta =  B*B - 4*A*C;

    if(delta < 0)
    {
        double possible_dis = min(lenb+va->distance, lena+vb->distance);
        vc->distance = min(vc->distance, possible_dis);
    }

    else{
        t1 = (-B + sqrt(delta))/(2*A);
        t2 = (-B - sqrt(delta))/(2*A);

        double pending = lenb*(t1 - u)/t1;

        ////we may need to worry about the case when theta = PI/2;
        if(t1 > 0 && (t1 > u) && pending > lena*cos_theta && pending < lena / cos_theta)
        {
            vc->distance = min(vc->distance, t1+va->distance);
        }
        else
        {
            double possible_dis = min(lenb+va->distance, lena+vb->distance);
            vc->distance = min(vc->distance, possible_dis);
        }
    }
}

/*
we need an extra structure to save the previous boundary edges
*/

void EvenStreamlinePlace::set_boundary_edges(int *trianglelist, int num_triangles)
{
    int i, j;
    Triangle *face;
    Edge *cur_e;
    int other_f;

    if(geo_boundary == NULL)
        geo_boundary = new EdgeList();  /*initially, we have 1000 edges*/

    ////Reset all the flags of the edges in the region
    for(i=0; i<geo_boundary->nedges; i++)
    {
        geo_boundary->edges[i]->visited = false;
        geo_boundary->edges[i]->OnBoundary = false;
    }

    //for(i = 0; i < num_triangles; i++)
    //{
    //	face = object->tlist.tris[trianglelist[i]];

    //	for(j = 0; j < 3; j++)
    //	{
    //		cur_e = face->edges[j];
    //		cur_e->OnBoundary = false;
    //		cur_e->visited = false;
    //	}
    //}

    geo_boundary->nedges = 0;

    ////Mark only boundary edges
    for(i = 0; i < num_triangles; i++)
    {
        face = object->tlist.tris[trianglelist[i]];

        for(j = 0; j < 3; j++)
        {
            cur_e = face->edges[j];
            if(cur_e->visited)
                continue;

            //other_f = cur_e->tris[0]->index;

            //if(other_f == trianglelist[i])
            //	other_f = cur_e->tris[1]->index;

            if(cur_e->tris[0]!=NULL)
            {
                other_f = cur_e->tris[0]->index;
                if(other_f == trianglelist[i] && cur_e->tris[1] == NULL)
                    continue;
                else if(other_f == trianglelist[i] && cur_e->tris[1] !=NULL)
                    other_f = cur_e->tris[1]->index;
            }
            else if(cur_e->tris[1] !=NULL && cur_e->tris[1]->index != trianglelist[i])
                other_f=cur_e->tris[1]->index;
            else
                continue;

            //if(!is_repeated_elem(trianglelist, other_f, num_triangles))
            if(!object->tlist.tris[other_f]->visited)
            {
                cur_e->OnBoundary = true;

                geo_boundary->append(cur_e);
            }

            cur_e->visited = true;  /*new added 05/01/07*/
        }
    }
}



void EvenStreamlinePlace::get_dis_for_obtuse_triangle(int v1, int v2, int v3, int origin_triangle, double lena, double lenb)
{
    int theVer;
    double p[3] = {0.};
    icVector3 vec1, vec2, dis_vec;
    double theta1, theta2;
    icVector3 LX, LY;

    //Get the shadow region which represented by v3 and two rays vec1, and vec2
    get_shadow_region(v1, v2, v3, origin_triangle, vec1, vec2, theta1, theta2, LX, LY);

    ////Calculate the distance
    if(unfold_to_find_ver(v1, v2, v3, origin_triangle, vec1, vec2, LX, LY,
                           theta1, theta2, theVer, p))
    {
        //calculat the Euclidean distance between v3 and p
        dis_vec.entry[0] = p[0] - object->vlist.verts[v3]->x;
        dis_vec.entry[1] = p[1] - object->vlist.verts[v3]->y;
        dis_vec.entry[2] = p[2] - object->vlist.verts[v3]->z;

        object->vlist.verts[v3]->distance =
            min(object->vlist.verts[v3]->distance, object->vlist.verts[theVer]->distance + length(dis_vec));

        ////Testing codes
        icVector3 VP;
        VP.entry[0] = p[0] - object->vlist.verts[v3]->x;
        VP.entry[1] = p[1] - object->vlist.verts[v3]->y;
        VP.entry[2] = p[2] - object->vlist.verts[v3]->z;

        double c = dot(VP, object->tlist.tris[origin_triangle]->normal);
    }

    else // use Dijkstra distance
    {
        double possible = min(object->vlist.verts[v1]->distance + lena,
                              object->vlist.verts[v2]->distance + lenb);
        object->vlist.verts[v3]->distance = min(object->vlist.verts[v3]->distance, possible);
    }
}


void EvenStreamlinePlace::get_shadow_region(int &v1, int &v2, int v3, int origin_triangle, icVector3 &vec1, icVector3 &vec2,
                                            double &theta1, double &theta2, icVector3 &LX, icVector3 &LY)
{
    int vert;
    icVector3 vec31, vec32, tmp1, tmp2, normal;
    double rot_mat1[16], rot_mat2[16], p[3];

    p[0] = object->vlist.verts[v3]->x;
    p[1] = object->vlist.verts[v3]->y;
    p[2] = object->vlist.verts[v3]->z;

    normal = object->tlist.tris[origin_triangle]->normal;

    ////Find the shadow region
    vec31.entry[0] = object->vlist.verts[v1]->x - object->vlist.verts[v3]->x;
    vec31.entry[1] = object->vlist.verts[v1]->y - object->vlist.verts[v3]->y;
    vec31.entry[2] = object->vlist.verts[v1]->z - object->vlist.verts[v3]->z;

    vec32.entry[0] = object->vlist.verts[v2]->x - object->vlist.verts[v3]->x;
    vec32.entry[1] = object->vlist.verts[v2]->y - object->vlist.verts[v3]->y;
    vec32.entry[2] = object->vlist.verts[v2]->z - object->vlist.verts[v3]->z;

    normalize(vec31);
    normalize(vec32);

    ////Get the rotation matrix
    //get_rotation(p, normal, 90., rot_mat1);
    //get_rotation(p, normal, -90., rot_mat2);
    //g_pclGlView->get_rotation(p, normal, 90., rot_mat1);
    //g_pclGlView->get_rotation(p, normal, -90., rot_mat2);

    //g_pCgraphView->get_rotation(p, normal, 90., rot_mat1);
    //g_pCgraphView->get_rotation(p, normal, -90., rot_mat2);


    get_rotation(p, normal, 90./180*3.1415926, rot_mat1);
    get_rotation(p, normal, -90./180*3.1415926, rot_mat2);

    ////Rotate the two vectors
    tmp1 = vec31; tmp2 = vec32;

    transform_vec3D(tmp1, rot_mat1);
    transform_vec3D(tmp2, rot_mat2);

    /* Use other transformation */
    if(dot(tmp1, vec32) < 0 || dot(tmp2, vec31) < 0)
    {
        tmp1 = vec31; tmp2 = vec32;

        transform_vec3D(tmp1, rot_mat2);
        transform_vec3D(tmp2, rot_mat1);

        if(dot(tmp1, tmp2) < 0)  //something is wrong here
        {
            int test = 0;
        }
    }
    normalize(tmp1);
    normalize(tmp2);

    vec1 = tmp1; vec2 = tmp2;

    //Get the two angles for vec1 and vec2 using atan2
    LX = vec31;
    LY = cross(normal, LX);

    double a, b, theta;
    a = dot(vec32, LX);
    b = dot(vec32, LY);

    theta = atan2(b, a);

    if(theta < 0)
    {
        LX = vec32;
        LY = cross(normal, LX);

        //swap the vec1 and vec2, v1 and v2
        icVector3 tmp = vec1;
        vec1 = vec2;
        vec2 = tmp;

        int tmpv = v1;
        v1 = v2;
        v2 = tmpv;
    }

    a = dot(vec1, LX);
    b = dot(vec1, LY);
    theta1 = atan2(b, a);

    a = dot(vec2, LX);
    b = dot(vec2, LY);
    theta2 = atan2(b, a);
}


/*
Get the 4th vertex through unfolding a series of nearby triangles
This is the most important subroutine for the distance calculation of obtuse case
in fast marching algorithm
*/
bool EvenStreamlinePlace::unfold_to_find_ver(int v1, int v2, int v3, int origin_triangle, icVector3 vec1, icVector3 vec2,
                                             icVector3 LX, icVector3 LY, double theta1, double theta2, int &theVer, double p[3])
{
    int counter = 0;
    int local_v1, local_v2, local_v3;
    Vertex *l_v1, *l_v2;
    int cur_triangle = origin_triangle;
    int other_f;
    int edgei;
    Edge *cur_e;
    icVector3 rot_axis,  normal1, normal2, local_LX, local_LY;
    double local_p1[3], local_p2[3];
    double old_rotmat[16], rotmat[16];
    double base[3] = {0.};
    double rot_ang;

    ////initial
    local_v1 = v1;
    local_v2 = v2;
    set_Indentity16(old_rotmat);
    set_Indentity16(rotmat);

    while ( counter <= 10 )
    {
        //Get the two vertices
        local_p1[0] = object->vlist.verts[local_v1]->x;
        local_p1[1] = object->vlist.verts[local_v1]->y;
        local_p1[2] = object->vlist.verts[local_v1]->z;

        local_p2[0] = object->vlist.verts[local_v2]->x;
        local_p2[1] = object->vlist.verts[local_v2]->y;
        local_p2[2] = object->vlist.verts[local_v2]->z;

        //Get the edge that consists of the two vertices
        edgei = get_edgeIndex_of_triangle(local_v1, local_v2, cur_triangle);

        cur_e = object->tlist.tris[cur_triangle]->edges[edgei];

        //Get rotation axis based on the edge 'cur_e'
        rot_axis.entry[0] = cur_e->verts[1]->x - cur_e->verts[0]->x;
        rot_axis.entry[1] = cur_e->verts[1]->y - cur_e->verts[0]->y;
        rot_axis.entry[2] = cur_e->verts[1]->z - cur_e->verts[0]->z;
        normalize(rot_axis);

        base[0] = cur_e->verts[0]->x;
        base[1] = cur_e->verts[0]->y;
        base[2] = cur_e->verts[0]->z;

        //Get the neighboring triangle according to 'cur_e'
        //other_f = cur_e->tris[0]->index;
        //if(other_f == cur_triangle)
        //	other_f = cur_e->tris[1]->index;

        if(cur_e->tris[0]!=NULL)
        {
            other_f = cur_e->tris[0]->index;
            if(other_f == cur_triangle && cur_e->tris[1] == NULL)
                continue;
            else if(other_f == cur_triangle && cur_e->tris[1] !=NULL)
                other_f = cur_e->tris[1]->index;
        }
        else if(cur_e->tris[1] !=NULL && cur_e->tris[1]->index != cur_triangle)
            other_f=cur_e->tris[1]->index;
        else
            continue;

        //Get the third vertex of other_f rather than local_v1, local_v2
        local_v3 = get_thirdVer_of_triangle(local_v1, local_v2, other_f);

        //Get the two normals
        normal1 = object->tlist.tris[cur_triangle]->normal;
        normal2 = object->tlist.tris[other_f]->normal;

        //Get the angle between the normals of cur_triangle and other_f
        local_LX = normal1;
        local_LY = cross(rot_axis, normal1);

        double a = dot(normal2, local_LX);
        double b = dot(normal2, local_LY);

        rot_ang = atan2(b, a);

        //rot_ang = rot_ang * 180/3.14159265;

        //Get rotation matrix
        //get_rotation(base, rot_axis, -rot_ang, rotmat);

        //g_pCgraphView->get_rotation(base, rot_axis, -rot_ang, rotmat);

        //FILE *fp;
        //fp = fopen("rot_cg.txt", "w");
        //for(int i=0; i<16; i++)
        //	fprintf(fp, "%f\n", rotmat[i]);
        //fclose(fp);

        /*test whether we get the same rotation matrix or not*/
        //double t_rotmat[16] = {0.};
        //double t_oldrot[16] = {0.};
        //set_Indentity16(t_rotmat);
        //for(int i=0; i<16; i++)
        //	t_oldrot[i] = old_rotmat[i];

        get_rotation(base, rot_axis, -rot_ang, rotmat);
        //fp = fopen("rot_cg_t.txt", "w");
        //for(int i=0; i<16; i++)
        //	fprintf(fp, "%f\n", t_rotmat[i]);
        //fclose(fp);

        //Update the old_rotmat by right multiplying rotmat
        //g_pCgraphView->rightMultiply16(old_rotmat, rotmat);

        rightMultiply16_2(old_rotmat, rotmat);



        //Transform the third vertex local_v3
        p[0] = object->vlist.verts[local_v3]->x;
        p[1] = object->vlist.verts[local_v3]->y;
        p[2] = object->vlist.verts[local_v3]->z;

        transform_point3D(p, old_rotmat);

        //Judge whether p falls in the shadow region or not
        if(is_fall_in_shadow(theta1, theta2, LX, LY, p, v3,
                              local_v1, local_v2, local_v3))
        {
            theVer = local_v3;
            if(object->vlist.verts[theVer]->distance > 1.e47) return false;
            return true;
        }

        // Update current triangle and the counter
        cur_triangle = other_f;
        counter++;
    }
    return false;
}


int EvenStreamlinePlace::get_edgeIndex_of_triangle(int v1, int v2, int curtriangle)
{
    Triangle *face = object->tlist.tris[curtriangle];
    Edge *cur_e;
    int i;

    for(i = 0; i < 3; i++)
    {
        cur_e = face->edges[i];

        if((cur_e->verts[0]->index == v1 && cur_e->verts[1]->index == v2)
            || cur_e->verts[0]->index== v2 && cur_e->verts[1]->index == v1)
            return i;
    }
}


bool EvenStreamlinePlace::is_fall_in_shadow(double theta1, double theta2, icVector3 LX, icVector3 LY, double p[3], int origin_v3,
                                            int &local_v1, int &local_v2, int local_v3)
{
    Vertex *v3 = object->vlist.verts[origin_v3];
    icVector3 v3p;
    v3p.entry[0] = p[0] - v3->x;
    v3p.entry[1] = p[1] - v3->y;
    v3p.entry[2] = p[2] - v3->z;
    normalize(v3p);

    double a = dot(v3p, LX);
    double b = dot(v3p, LY);

    double theta3 = atan2(b, a);


    if(theta3 >= theta2 && theta3 <= theta1)
        return true;

    else if(theta3 > theta1)
    {
        local_v2 = local_v3;
        return false;
    }

    else if(theta3 < theta2)
    {
        local_v1 = local_v3;
        return false;
    }
}


/*
save the placement results into a file
Format:
#number of trajectories
#number of points
traj:1
#number of line segments
x,y,z,triangle
......
traj:2
#number of line segments
x,y,z,triangle
*/
void EvenStreamlinePlace::save_to_a_file(const char* filename)
{
    FILE *fp = fopen(filename, "w");

    fprintf(fp,"#%d\n", evenstreamlines->ntrajs);
    fprintf(fp,"#%d\n", 0);

    int i, j;
    Trajectory *traj;
    for(i = 0; i < evenstreamlines->ntrajs; i++)
    {
        fprintf(fp,"traj:%d\n", i);
        traj = evenstreamlines->trajs[i];
        fprintf(fp,"#%d\n",traj->nlinesegs);
        for(j = 0; j < evenstreamlines->trajs[i]->nlinesegs; j++)
        {
            fprintf(fp,"%f,%f,%f,%d\n", traj->linesegs[j].gstart.entry[0],
                    traj->linesegs[j].gstart.entry[1],
                    traj->linesegs[j].gstart.entry[2],
                    traj->linesegs[j].Triangle_ID);

            //fprintf(fp,"%f,%f,%f,%d\n", traj->linesegs[j].gend.entry[0],
            //	traj->linesegs[j].gend.entry[1],
            //	traj->linesegs[j].gend.entry[2],
            //	traj->linesegs[j].Triangle_ID);
        }
    }
    fclose(fp);
}


/*
save the placement results into a file
Format:
#number of trajectories
#number of points
traj:1
#number of line segments
x,y,z,triangle
......
traj:2
#number of line segments
x,y,z,triangle
*/
void EvenStreamlinePlace::save_to_a_file(const char* filename, int flag)
{
    FILE *fp = fopen(filename, "w");

    fprintf(fp, "#%d\n", evenstreamlines->ntrajs);
    fprintf(fp, "periodic_orbits: %d\n", periodic_orbits->nporbits);
    fprintf(fp, "separatrices: %d\n", separatrices->ntrajs);
    //fprintf(fp,"#%d\n", 0);

    int i, j;
    Trajectory *traj;

    /*save the periodic orbits*/
    for(i = 0; i < periodic_orbits->nporbits; i++)
    {
        fprintf(fp, "PO: %d, type: %d, #%d\n", i, periodic_orbits->polist[i]->type,
                periodic_orbits->polist[i]->traj->nlinesegs);
        traj = periodic_orbits->polist[i]->traj;
        for(j = 0; j < traj->nlinesegs; j++)
        {
            fprintf(fp,"%f,%f,%f,%d\n", traj->linesegs[j].gstart.entry[0],
                    traj->linesegs[j].gstart.entry[1],
                    traj->linesegs[j].gstart.entry[2],
                    traj->linesegs[j].Triangle_ID);
        }
    }

    /*save the separatrices*/
    for(i = 0; i < separatrices->ntrajs; i++)
    {
        fprintf(fp, "Sep: %d, type: %d, #%d\n", i, i%2,
                separatrices->trajs[i]->nlinesegs); // 0--outgoing separatrix, 1--incoming one
        traj = separatrices->trajs[i];;
        for(j = 0; j < traj->nlinesegs; j++)
        {
            fprintf(fp,"%f,%f,%f,%d\n", traj->linesegs[j].gstart.entry[0],
                    traj->linesegs[j].gstart.entry[1],
                    traj->linesegs[j].gstart.entry[2],
                    traj->linesegs[j].Triangle_ID);
        }
    }

    /*save the other streamlines*/
    for(i = periodic_orbits->nporbits + separatrices->ntrajs; i < evenstreamlines->ntrajs; i++)
    {
        fprintf(fp,"traj:%d\n", i);
        traj = evenstreamlines->trajs[i];
        fprintf(fp,"#%d\n",traj->nlinesegs);
        for(j = 0; j < evenstreamlines->trajs[i]->nlinesegs; j++)
        {
            fprintf(fp,"%f,%f,%f,%d\n", traj->linesegs[j].gstart.entry[0],
                    traj->linesegs[j].gstart.entry[1],
                    traj->linesegs[j].gstart.entry[2],
                    traj->linesegs[j].Triangle_ID);

            //fprintf(fp,"%f,%f,%f,%d\n", traj->linesegs[j].gend.entry[0],
            //	traj->linesegs[j].gend.entry[1],
            //	traj->linesegs[j].gend.entry[2],
            //	traj->linesegs[j].Triangle_ID);
        }
    }
    fclose(fp);
}


void EvenStreamlinePlace::load_from_a_file(const char *filename)
{
    int ntrajs = 0;
    int totallines = 0;
    int cur_traj = 0;

    float x, y, z;

    int i, j;

    FILE *fp = fopen(filename, "r");

    if(fp == NULL)
    {
        fprintf(stderr, "can't open file %s\n", filename);
        return;
    }

    /*need to check the file format here, if it is not correct, return error information*/

    fscanf(fp, "#%d\n", &ntrajs);
    fscanf(fp, "#%d\n", &totallines);

    /*allocate memory for the streamlines*/

    /*first, release the old lists*/

    /*release the sample point list*/
    if(samplepts != NULL)
    {
        for(i = 0; i < evenstreamlines->curMaxNumTrajs; i++)
        {
            if(samplepts[i] != NULL)
                delete samplepts[i];
        }
        delete [] samplepts;
    }

    if(evenstreamlines != NULL)
        delete evenstreamlines;

    /*second, allocate the memory according to the number of streamlines read from the file*/
    evenstreamlines = new TrajectoryList(ntrajs);

    for(i = 0; i < ntrajs; i++)
    {
        fscanf(fp, "traj:%d\n", &cur_traj);
        fscanf(fp, "#%d\n", &totallines);

        /*allocate memory for this streamline*/
        evenstreamlines->trajs[cur_traj] = new Trajectory(cur_traj, totallines);
        for(j = 0; j < totallines; j++)
        {
            fscanf(fp, "%f,%f,%f,%d\n", &x, &y, &z, &evenstreamlines->trajs[cur_traj]->linesegs[j].Triangle_ID);
            evenstreamlines->trajs[cur_traj]->linesegs[j].gstart.entry[0] = x;
            evenstreamlines->trajs[cur_traj]->linesegs[j].gstart.entry[1] = y;
            evenstreamlines->trajs[cur_traj]->linesegs[j].gstart.entry[2] = z;

        }
        evenstreamlines->trajs[cur_traj]->nlinesegs = totallines;

        /*assign the gend*/
        for(j = 0; j < totallines-1; j++)
        {
            evenstreamlines->trajs[cur_traj]->linesegs[j].gend =
                evenstreamlines->trajs[cur_traj]->linesegs[j+1].gstart;
        }
        evenstreamlines->trajs[cur_traj]->linesegs[j].gend =
            evenstreamlines->trajs[cur_traj]->linesegs[j].gstart;
    }

    evenstreamlines->ntrajs = ntrajs;
    fclose(fp);
}


void EvenStreamlinePlace::load_from_a_file_enhanced(const char *filename)
{
    int ntrajs = 0;
    int npos = 0;
    int nseparatrices = 0;
    int totallines = 0;
    int cur_traj = 0;
    int po_type = 0;
    int po_id = 0;

    float x, y, z;

    int i, j;

    FILE *fp = fopen(filename, "r");

    if(fp == NULL)
    {
        fprintf(stderr, "can't open file %s\n", filename);
        return;
    }

    /*need to check the file format here, if it is not correct, return error information*/

    fscanf(fp, "#%d\n", &ntrajs);
    //fscanf(fp, "#%d\n", &totallines);

    fscanf(fp, "periodic_orbits: %d\n", &npos);
    fscanf(fp, "separatrices: %d\n", &nseparatrices);

    /*load the periodic orbits*/

    /*allocate space for periodic orbits*/
    if(periodic_orbits != NULL)
        delete periodic_orbits;
    periodic_orbits = new PeriodicOrbitList(npos);

    for(i=0; i<npos; i++)
    {
        periodic_orbits->polist[i]=new PeriodicOrbit();
        fscanf(fp,"PO: %d, type: %d, #%d\n", &po_id, &po_type, &totallines);
        periodic_orbits->polist[i]->type = po_type;
        periodic_orbits->polist[i]->traj = new Trajectory(po_id, totallines);

        /*read the gstart*/
        for(j=0; j<totallines; j++)
        {
            fscanf(fp, "%f,%f,%f,%d\n", &x, &y, &z, &periodic_orbits->polist[i]->traj->linesegs[j].Triangle_ID);
            periodic_orbits->polist[i]->traj->linesegs[j].gstart.entry[0] = x;
            periodic_orbits->polist[i]->traj->linesegs[j].gstart.entry[1] = y;
            periodic_orbits->polist[i]->traj->linesegs[j].gstart.entry[2] = z;
        }

        /*set the gend*/
        for(j=0; j<totallines-1; j++)
        {
            periodic_orbits->polist[i]->traj->linesegs[j].gend =
                periodic_orbits->polist[i]->traj->linesegs[j+1].gstart;
        }
        periodic_orbits->polist[i]->traj->linesegs[j].gend =
            periodic_orbits->polist[i]->traj->linesegs[j].gstart;

        periodic_orbits->polist[i]->traj->nlinesegs = totallines;
    }
    periodic_orbits->nporbits = npos;

    /*load the separatrices*/

    if(separatrices != NULL)
        delete separatrices;
    separatrices = new TrajectoryList(nseparatrices);

    for(i=0; i<nseparatrices; i++)
    {
        fscanf(fp, "Sep: %d, type: %d, #%d\n", &po_id, &po_type, &totallines);
        separatrices->trajs[i]=new Trajectory(i, totallines);

        /*read the gstart*/
        for(j=0; j<totallines; j++)
        {
            fscanf(fp, "%f,%f,%f,%d\n", &x, &y, &z, &separatrices->trajs[i]->linesegs[j].Triangle_ID);
            separatrices->trajs[i]->linesegs[j].gstart.entry[0] = x;
            separatrices->trajs[i]->linesegs[j].gstart.entry[1] = y;
            separatrices->trajs[i]->linesegs[j].gstart.entry[2] = z;
        }

        /*set the gend*/
        for(j=0; j<totallines-1; j++)
        {
            separatrices->trajs[i]->linesegs[j].gend =
                separatrices->trajs[i]->linesegs[j+1].gstart;
        }
        separatrices->trajs[i]->linesegs[j].gend =
            separatrices->trajs[i]->linesegs[j].gstart;

        separatrices->trajs[i]->nlinesegs=totallines;
    }
    separatrices->ntrajs = nseparatrices;

    /*allocate memory for the streamlines*/

    /*first, release the old lists*/

    /*release the sample point list*/
    if(samplepts != NULL)
    {
        for(i = 0; i < evenstreamlines->curMaxNumTrajs; i++)
        {
            if(samplepts[i] != NULL)
                delete samplepts[i];
        }
        delete [] samplepts;
    }

    if(evenstreamlines != NULL)
        delete evenstreamlines;

    /*second, allocate the memory according to the number of streamlines read from the file*/
    evenstreamlines = new TrajectoryList(ntrajs);

    for(i=0; i<npos+nseparatrices; i++)
        evenstreamlines->trajs[i] = NULL; /*Null for periodic orbits and separatrices*/

    for(i = npos+nseparatrices; i < ntrajs; i++)
    {
        fscanf(fp, "traj:%d\n", &cur_traj);
        fscanf(fp, "#%d\n", &totallines);

        /*allocate memory for this streamline*/
        evenstreamlines->trajs[cur_traj] = new Trajectory(cur_traj, totallines);
        for(j = 0; j < totallines; j++)
        {
            fscanf(fp, "%f,%f,%f,%d\n", &x, &y, &z, &evenstreamlines->trajs[cur_traj]->linesegs[j].Triangle_ID);
            evenstreamlines->trajs[cur_traj]->linesegs[j].gstart.entry[0] = x;
            evenstreamlines->trajs[cur_traj]->linesegs[j].gstart.entry[1] = y;
            evenstreamlines->trajs[cur_traj]->linesegs[j].gstart.entry[2] = z;

        }
        evenstreamlines->trajs[cur_traj]->nlinesegs = totallines;

        /*assign the gend*/
        for(j = 0; j < totallines-1; j++)
        {
            evenstreamlines->trajs[cur_traj]->linesegs[j].gend =
                evenstreamlines->trajs[cur_traj]->linesegs[j+1].gstart;
        }
        evenstreamlines->trajs[cur_traj]->linesegs[j].gend =
            evenstreamlines->trajs[cur_traj]->linesegs[j].gstart;
    }

    evenstreamlines->ntrajs = ntrajs;

    fclose(fp);
}
