/*
This file contains routines of class Trajectory and TrajectoryList that we use to compute
the streamlines for steady flow.
(trajectory for time-varying vector field in the future)
Note that Streamlines are not the inherent features of the surfaces, therefore, we separate
them from the surface geometry.

Created and Modified by Guoning Chen
        copyright @2007
    */


//#include "VField.h"

#include "Analysis/MorseDecomp.h"

#include "Analysis/PeriodicOrbitExtraction.h"
#include "ExternalDependencies/nr.h"

#include "Others/common_routines.h"

extern Polyhedron *object;
extern MorseDecomp *morse_decomp;
extern PeriodicOrbitList *periodic_orbits;
extern PeriodicOrbit_Detector *periodicorbit_detector;
extern const int DebugOn;
extern int Integrator_opt;

/*
Do we want a separate trajectory list from the separatrix list
or combine them together?
Now we separate them into two lists. Therefore, we need a global index for all trajectories
*/

TrajectoryList *separatrices = NULL;

int g_traj_index = 0;


Trajectory::Trajectory(int index, int curMaxNum)
{
    this->index = index;

    linesegs = (LineSeg *)malloc(sizeof(LineSeg)*curMaxNum);

    if(linesegs == NULL)
    {
        char rout[256], var[256];
        sprintf(rout, "%s", "Trajectory::Trajectory");
        sprintf(var, "%s", "linesegs");

        //write_mem_error(rout, var, 0);
        curMaxNumLinesegs = 0;
        exit(-1);
    }

    curMaxNumLinesegs = curMaxNum;
    nlinesegs = 0;

    //eulerstep_scalar = 	0.004 * object->radius; /*should consider the shortest edge of the mesh*/

    /*for coolingjacket and swirl motion data sets*/
    //eulerstep_scalar = 	min(object->shortest_edgelength, 0.004 * object->radius);
    eulerstep_scalar = 	0.0035 * object->radius;  //0.004 before modification

}


bool Trajectory::extend_line_segments(int add_size)
{
    //LineSeg *temp = linesegs;

    /*Can we trust realloc function here ? */
    //linesegs = (LineSeg *)realloc(linesegs, sizeof(LineSeg)*(curMaxNumLinesegs+add_size));

    //if(linesegs == NULL)
    //{
    //	char rout[256], var[256];
    //	sprintf(rout, "%s", "Trajectory::extend_line_segments");
    //	sprintf(var, "%s", "linesegs");

    //	write_mem_error(rout, var, 1);
    //	return false;
    //}

    //free(temp);
    //curMaxNumLinesegs += add_size;
    //return true;

    LineSeg *extendlist=(LineSeg*)malloc(sizeof(LineSeg)*(curMaxNumLinesegs+add_size));

    if(extendlist == NULL)
    {
        char rout[256], var[256];
        sprintf(rout, "%s", "Trajectory::extend_line_segments");
        sprintf(var, "%s", "linesegs");

        //write_mem_error(rout, var, 1);
        return false;
    }

    int i;
    for(i = 0; i < curMaxNumLinesegs; i++)
    {
        extendlist[i].end = linesegs[i].end;

        extendlist[i].start = linesegs[i].start;

        extendlist[i].gend = linesegs[i].gend;

        extendlist[i].gstart = linesegs[i].gstart;

        extendlist[i].length = linesegs[i].length;
        extendlist[i].Triangle_ID = linesegs[i].Triangle_ID;
    }

    free(linesegs);
    linesegs = extendlist;
    curMaxNumLinesegs += add_size;
    return true;
}



void Trajectory::local_To_global(int faceid, double locpos[2], icVector3 &glpos)
{
    Triangle *t = object->tlist.tris[faceid];

    icVector3 dir3D = locpos[0] * t->LX + locpos[1] * t->LY;

    glpos.entry[0] = t->verts[0]->x + dir3D.entry[0];
    glpos.entry[1] = t->verts[0]->y + dir3D.entry[1];
    glpos.entry[2] = t->verts[0]->z + dir3D.entry[2];
}



bool Trajectory::cal_next_point_euler1(double first[2], double second[2], int &face_id, double alpha[3], int type)
{

    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    icVector2 VecAtPoint = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if(length(VecAtPoint) < 1e-15 ) return false;

    //scale the vector??? 01/23/07 !!!!!!!!!!!!!Important modification
    normalize(VecAtPoint);

    Triangle *t = object->tlist.tris[face_id];
    VecAtPoint = 0.4*length(t->verts[0]->t_vec)*VecAtPoint;  // for the normalized vector field 11/19/2010
    //VecAtPoint = 0.004*length(t->verts[0]->t_vec)*VecAtPoint;  // for the non-normalized vector field 11/19/2010

    ////Using first order Euler method to test 12/11/05
    if(type == 0)
    {
        second[0] = first[0] + eulerstep_scalar * VecAtPoint.entry[0];
        second[1] = first[1] + eulerstep_scalar * VecAtPoint.entry[1];
    }
    else
    {
        second[0] = first[0] - eulerstep_scalar * VecAtPoint.entry[0];
        second[1] = first[1] - eulerstep_scalar * VecAtPoint.entry[1];
    }

    return true;
}


/*
second order Euler method
*/

//bool Trajectory::cal_nextpt_2ndeuler(double first[2], double second[2], int &face_id, double alpha[3], int type)
//{
//	////Using first order Euler method to get the next point
//	icVector3 t_gp;
//	local_To_global(face_id, first, t_gp);
//
//	icVector2 VecAtPoint = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);
//
//	if(length(VecAtPoint) < 1e-15) return false;
//
//	//scale the vector??? 01/23/07 !!!!!!!!!!!!!Important modification
//	normalize(VecAtPoint);
//
//	Triangle *t = object->tlist.tris[face_id];
//	VecAtPoint = 0.08*length(t->verts[0]->t_vec)*VecAtPoint; /*evaluate in very small step now*/
//
//	double temp[2] = {0.};
//
//	if(type == 0)
//	{
//		temp[0] = first[0] + /*2**/eulerstep_scalar*VecAtPoint.entry[0];
//		temp[1] = first[1] + /*2**/eulerstep_scalar*VecAtPoint.entry[1];
//	}
//	else
//	{
//		temp[0] = first[0] - /*2**/eulerstep_scalar*VecAtPoint.entry[0];
//		temp[1] = first[1] - /*2**/eulerstep_scalar*VecAtPoint.entry[1];
//	}
//
//	////get the vector at next point
//	double alpha1[3] = {0.};
//	local_To_global(face_id, temp, t_gp);
//	object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha1);
//
//	if ( (alpha1[0] >= 0 ||fabs(alpha1[0]) <= 1e-8) && alpha1[0] <= 1
//			&& (alpha1[1] >= 0 || fabs(alpha1[1]) <= 1e-8) && alpha1[1] <= 1
//			&& (alpha1[2] >= 0 || fabs(alpha1[2]) <= 1e-8) && alpha1[2] <= 1)
//
//	{
//		icVector2 VecAtPoint2 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha1, temp[0], temp[1]);
//		normalize(VecAtPoint);
//		normalize(VecAtPoint2);
//
//		icVector2 total_v;
//		total_v.entry[0] = (VecAtPoint.entry[0] +  VecAtPoint2.entry[0])/2;
//		total_v.entry[1] = (VecAtPoint.entry[1] +  VecAtPoint2.entry[1])/2;
//		total_v = length(t->verts[0]->t_vec)*VecAtPoint;
//
//		if(type == 0)
//		{
//			second[0] = first[0] + 2*eulerstep_scalar*total_v.entry[0];
//			second[1] = first[1] + 2*eulerstep_scalar*total_v.entry[1];
//		}
//		else
//		{
//			second[0] = first[0] - 2*eulerstep_scalar*total_v.entry[0];
//			second[1] = first[1] - 2*eulerstep_scalar*total_v.entry[1];
//		}
//	}
//
//	else
//	{
//		if(type == 0)
//		{
//			second[0] = first[0] + 2*eulerstep_scalar*VecAtPoint.entry[0];
//			second[1] = first[1] + 2*eulerstep_scalar*VecAtPoint.entry[1];
//		}
//		else
//		{
//			second[0] = first[0] - 2*eulerstep_scalar*VecAtPoint.entry[0];
//			second[1] = first[1] - 2*eulerstep_scalar*VecAtPoint.entry[1];
//		}
//	}
//
//	return true;
//
//}


bool Trajectory::cal_nextpt_2ndeuler(double first[2], double second[2], int &face_id, double alpha[3], int type)
{
    ////Using first order Euler method to get the next point
    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    icVector2 VecAtPoint = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if(length(VecAtPoint) < 1e-15) return false;

    //scale the vector??? 01/23/07 !!!!!!!!!!!!!Important modification
    normalize(VecAtPoint);

    Triangle *t = object->tlist.tris[face_id];
    VecAtPoint = 0.4*length(t->verts[0]->t_vec)*VecAtPoint; /*evaluate in very small step now*/
    //VecAtPoint = 0.004*length(t->verts[0]->t_vec)*VecAtPoint; /*evaluate in very small step now*/

    double temp[2] = {0.};
    double euler2nd_step = eulerstep_scalar;

    if(type == 0)
    {
        temp[0] = first[0] + /*2**/euler2nd_step/2.*VecAtPoint.entry[0];
        temp[1] = first[1] + /*2**/euler2nd_step/2*VecAtPoint.entry[1];
    }
    else
    {
        temp[0] = first[0] - /*2**/euler2nd_step/2*VecAtPoint.entry[0];
        temp[1] = first[1] - /*2**/euler2nd_step/2*VecAtPoint.entry[1];
    }

    ////get the vector at next point
    double alpha1[3] = {0.};
    local_To_global(face_id, temp, t_gp);
    object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha1);

    if ( (alpha1[0] >= 0 ||fabs(alpha1[0]) <= 1e-8) && alpha1[0] <= 1
        && (alpha1[1] >= 0 || fabs(alpha1[1]) <= 1e-8) && alpha1[1] <= 1
        && (alpha1[2] >= 0 || fabs(alpha1[2]) <= 1e-8) && alpha1[2] <= 1)

    {
        icVector2 VecAtPoint2 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha1, temp[0], temp[1]);
        normalize(VecAtPoint);
        normalize(VecAtPoint2);

        icVector2 total_v;
        total_v.entry[0] = (VecAtPoint.entry[0] +  VecAtPoint2.entry[0])/2;
        total_v.entry[1] = (VecAtPoint.entry[1] +  VecAtPoint2.entry[1])/2;
        //total_v = length(t->verts[0]->t_vec)*VecAtPoint;                // for normalized vector field
        normalize(total_v);
        //total_v = 0.004*length(t->verts[0]->t_vec)*total_v;          // for non-normalized vector field
        total_v = 0.4*length(t->verts[0]->t_vec)*total_v;          // for normalized vector field

        if(type == 0)
        {
            second[0] = first[0] + euler2nd_step/*2*eulerstep_scalar*/*total_v.entry[0];
            second[1] = first[1] + euler2nd_step/*2*eulerstep_scalar*/*total_v.entry[1];
        }
        else
        {
            second[0] = first[0] - euler2nd_step/*2*eulerstep_scalar*/*total_v.entry[0];
            second[1] = first[1] - euler2nd_step/*2*eulerstep_scalar*/*total_v.entry[1];
        }
    }

    else
    {
        if(type == 0)
        {
            second[0] = first[0] + euler2nd_step/*2*eulerstep_scalar*/*VecAtPoint.entry[0];
            second[1] = first[1] + euler2nd_step/*2*eulerstep_scalar*/*VecAtPoint.entry[1];
        }
        else
        {
            second[0] = first[0] - euler2nd_step/*2*eulerstep_scalar*/*VecAtPoint.entry[0];
            second[1] = first[1] - euler2nd_step/*2*eulerstep_scalar*/*VecAtPoint.entry[1];
        }
    }

    return true;
}


bool
Trajectory::cal_nextpt_RK2_nonNorm(double first[2], double second[2], int &face_id, double alpha[3], int type)
{
    ////Using first order Euler method to get the next point
    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    icVector2 VecAtPoint = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if(length(VecAtPoint) < 1e-10) return false;

    //scale the vector??? 01/23/07 !!!!!!!!!!!!!Important modification
    //normalize(VecAtPoint);

    Triangle *t = object->tlist.tris[face_id];
    //VecAtPoint = 0.08*length(t->verts[0]->t_vec)*VecAtPoint; /*evaluate in very small step now*/

    double temp[2] = {0.};

    if(type == 0)
    {
        temp[0] = first[0] + eulerstep_scalar/2*VecAtPoint.entry[0];
        temp[1] = first[1] + eulerstep_scalar/2*VecAtPoint.entry[1];
    }
    else
    {
        temp[0] = first[0] - eulerstep_scalar/2*VecAtPoint.entry[0];
        temp[1] = first[1] - eulerstep_scalar/2*VecAtPoint.entry[1];
    }

    ////get the vector at next point
    double alpha1[3] = {0.};
    local_To_global(face_id, temp, t_gp);
    object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha1);

    if ( (alpha1[0] >= 0 ||fabs(alpha1[0]) <= 1e-8) && alpha1[0] <= 1
        && (alpha1[1] >= 0 || fabs(alpha1[1]) <= 1e-8) && alpha1[1] <= 1
        && (alpha1[2] >= 0 || fabs(alpha1[2]) <= 1e-8) && alpha1[2] <= 1)

    {
        icVector2 VecAtPoint2 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha1, temp[0], temp[1]);
        //normalize(VecAtPoint);
        //normalize(VecAtPoint2);

        icVector2 total_v;
        total_v.entry[0] = (VecAtPoint.entry[0] +  VecAtPoint2.entry[0])/2;
        total_v.entry[1] = (VecAtPoint.entry[1] +  VecAtPoint2.entry[1])/2;
        //total_v = length(t->verts[0]->t_vec)*VecAtPoint;

        if(type == 0)
        {
            second[0] = first[0] + 2*eulerstep_scalar*total_v.entry[0];
            second[1] = first[1] + 2*eulerstep_scalar*total_v.entry[1];
        }
        else
        {
            second[0] = first[0] - 2*eulerstep_scalar*total_v.entry[0];
            second[1] = first[1] - 2*eulerstep_scalar*total_v.entry[1];
        }
    }

    else
    {
        if(type == 0)
        {
            second[0] = first[0] + 2*eulerstep_scalar*VecAtPoint.entry[0];
            second[1] = first[1] + 2*eulerstep_scalar*VecAtPoint.entry[1];
        }
        else
        {
            second[0] = first[0] - 2*eulerstep_scalar*VecAtPoint.entry[0];
            second[1] = first[1] - 2*eulerstep_scalar*VecAtPoint.entry[1];
        }
    }

    return true;

}


double RK4_step;
bool
Trajectory::cal_nextpt_RK4(double first[2], double second[2], int &face_id, double alpha[3], int type)
{
    ////Using 4th Runge-Kutta method to get the next point
    icVector2 t_v;
    double temp[2] = {0.};

    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    Triangle *t = object->tlist.tris[face_id];

    /*compute K1*/
    icVector2 V1 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if(length(V1) < 1e-14) return false;


    RK4_step = eulerstep_scalar/2.;  // for normalized
    //RK4_step = eulerstep_scalar/16.; // for diesel slice

    t_v = V1;
    normalize(t_v);
    //t_v = 0.08*length(t->verts[0]->t_vec)*t_v; /*evaluate in very small step now*/

    if(type == 0)
    {
        temp[0] = first[0] + RK4_step/2*t_v.entry[0];
        temp[1] = first[1] + RK4_step/2*t_v.entry[1];
    }
    else
    {
        temp[0] = first[0] - RK4_step/2*t_v.entry[0];
        temp[1] = first[1] - RK4_step/2*t_v.entry[1];
    }

    /*compute K2*/
    //Get2DBarycentricFacters(face_id, temp[0], temp[1], alpha);
    //icVector2 V2 = GetVectorAtPoints(face_id, alpha);
    local_To_global(face_id, temp, t_gp);
    object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha);
    icVector2 V2 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, temp[0], temp[1]);
    t_v = V2;
    normalize(t_v);
    //t_v = 0.08*length(t->verts[0]->t_vec)*t_v; /*evaluate in very small step now*/

    if(type == 0)
    {
        temp[0] = first[0] + RK4_step/2*t_v.entry[0];
        temp[1] = first[1] + RK4_step/2*t_v.entry[1];
    }
    else
    {
        temp[0] = first[0] - RK4_step/2*t_v.entry[0];
        temp[1] = first[1] - RK4_step/2*t_v.entry[1];
    }

    /*compute K3*/
    local_To_global(face_id, temp, t_gp);
    object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha);
    icVector2 V3 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, temp[0], temp[1]);
    t_v = V3;
    normalize(t_v);
    //t_v = 0.08*length(t->verts[0]->t_vec)*t_v; /*evaluate in very small step now*/

    if(type == 0)
    {
        temp[0] = first[0] + RK4_step*t_v.entry[0];
        temp[1] = first[1] + RK4_step*t_v.entry[1];
    }
    else
    {
        temp[0] = first[0] - RK4_step*t_v.entry[0];
        temp[1] = first[1] - RK4_step*t_v.entry[1];
    }

    /*compute K4*/
    local_To_global(face_id, temp, t_gp);
    object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha);
    icVector2 V4 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, temp[0], temp[1]);

    icVector2 total_v = 1./6.*(V1+2*V2+2*V3+V4);
    t_v = total_v;
    normalize(t_v);
    //t_v = 0.08*length(t->verts[0]->t_vec)*t_v; /*evaluate in very small step now*/

    if(type == 0)
    {
        second[0] = first[0] + RK4_step*t_v.entry[0];
        second[1] = first[1] + RK4_step*t_v.entry[1];
    }

    else
    {
        second[0] = first[0] - RK4_step*t_v.entry[0];
        second[1] = first[1] - RK4_step*t_v.entry[1];
    }

    return true;
}


bool
Trajectory::cal_next_pt_RK4_tau(double first[2], double second[2], int &face_id, double alpha[3], int type)
{
    ////Using 4th Runge-Kutta method to get the next point
    icVector2 t_v;
    double temp[2] = {0.};

    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    /*compute K1*/
    icVector2 V1 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if(length(V1) < 1e-10) return false;

    //double shortest = get_shortestedge_tri(face_id);
    //double scalar = 0.06;

    //if(Mag_Scheme == 0)
    //{
    //	RK4_step = scalar*shortest*0.5;
    //}
    //else if(Mag_Scheme == 1)
    //{
    //	RK4_step = scalar*shortest;
    //}
    //else
    //{
    //	RK4_step = scalar*shortest*2;
    //}

    RK4_step = eulerstep_scalar/4.;
    //RK4_step = eulerstep_scalar/16.; // for diesel slice

    t_v = V1;
    normalize(t_v);

    if(type == 0)
    {
        temp[0] = first[0] + RK4_step/2*t_v.entry[0];
        temp[1] = first[1] + RK4_step/2*t_v.entry[1];
    }
    else
    {
        temp[0] = first[0] - RK4_step/2*t_v.entry[0];
        temp[1] = first[1] - RK4_step/2*t_v.entry[1];
    }

    /*compute K2*/
    //Get2DBarycentricFacters(face_id, temp[0], temp[1], alpha);
    //icVector2 V2 = GetVectorAtPoints(face_id, alpha);
    local_To_global(face_id, temp, t_gp);
    object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha);
    icVector2 V2 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, temp[0], temp[1]);
    t_v = V2;
    normalize(t_v);

    if(type == 0)
    {
        temp[0] = first[0] + RK4_step/2*t_v.entry[0];
        temp[1] = first[1] + RK4_step/2*t_v.entry[1];
    }
    else
    {
        temp[0] = first[0] - RK4_step/2*t_v.entry[0];
        temp[1] = first[1] - RK4_step/2*t_v.entry[1];
    }

    /*compute K3*/
    local_To_global(face_id, temp, t_gp);
    object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha);
    icVector2 V3 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, temp[0], temp[1]);
    t_v = V3;
    normalize(t_v);

    if(type == 0)
    {
        temp[0] = first[0] + RK4_step*t_v.entry[0];
        temp[1] = first[1] + RK4_step*t_v.entry[1];
    }
    else
    {
        temp[0] = first[0] - RK4_step*t_v.entry[0];
        temp[1] = first[1] - RK4_step*t_v.entry[1];
    }

    /*compute K4*/
    local_To_global(face_id, temp, t_gp);
    object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha);
    icVector2 V4 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, temp[0], temp[1]);

    icVector2 total_v = 1./6.*(V1+2*V2+2*V3+V4);
    t_v = total_v;
    normalize(t_v);

    if(type == 0)
    {
        second[0] = first[0] + RK4_step*t_v.entry[0];
        second[1] = first[1] + RK4_step*t_v.entry[1];
    }

    else
    {
        second[0] = first[0] - RK4_step*t_v.entry[0];
        second[1] = first[1] - RK4_step*t_v.entry[1];
    }

    /*ready for the temporary tau  08/30/07*/
    icVector2 line_v;
    line_v.entry[0] = second[0]-first[0];
    line_v.entry[1] = second[1]-first[1];
    normalize(line_v);
    double proj_len = dot(line_v, t_v); /*here we use the dot product in the local frame*/
    cur_vec_mag = fabs(proj_len);

    return true;
}



bool
Trajectory::cal_next_pt_RK4_tau_f(double first[2], double second[2], int &face_id, double alpha[3])
{
    ////Using 4th Runge-Kutta method to get the next point
    icVector2 t_v;
    double temp[2] = {0.};

    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    /*compute K1*/
    icVector2 V1 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if(length(V1) < 1e-10) return false;


    RK4_step = eulerstep_scalar/4.;
    //RK4_step = eulerstep_scalar/16.; // for diesel slice

    t_v = V1;
    normalize(t_v);

    temp[0] = first[0] + RK4_step/2*t_v.entry[0];
    temp[1] = first[1] + RK4_step/2*t_v.entry[1];

    /*compute K2*/
    //Get2DBarycentricFacters(face_id, temp[0], temp[1], alpha);
    //icVector2 V2 = GetVectorAtPoints(face_id, alpha);
    local_To_global(face_id, temp, t_gp);
    object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha);
    icVector2 V2 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, temp[0], temp[1]);
    t_v = V2;
    normalize(t_v);

    temp[0] = first[0] + RK4_step/2*t_v.entry[0];
    temp[1] = first[1] + RK4_step/2*t_v.entry[1];

    /*compute K3*/
    local_To_global(face_id, temp, t_gp);
    object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha);
    icVector2 V3 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, temp[0], temp[1]);
    t_v = V3;
    normalize(t_v);

    temp[0] = first[0] + RK4_step*t_v.entry[0];
    temp[1] = first[1] + RK4_step*t_v.entry[1];


    /*compute K4*/
    local_To_global(face_id, temp, t_gp);
    object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha);
    icVector2 V4 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, temp[0], temp[1]);

    icVector2 total_v = 1./6.*(V1+2*V2+2*V3+V4);
    t_v = total_v;
    normalize(t_v);

    second[0] = first[0] + RK4_step*t_v.entry[0];
    second[1] = first[1] + RK4_step*t_v.entry[1];

    /*ready for the temporary tau  08/30/07*/
    icVector2 line_v;
    line_v.entry[0] = second[0]-first[0];
    line_v.entry[1] = second[1]-first[1];
    move_dist = length(line_v);
    //normalize(line_v);

    if (move_dist>0)
        line_v = 1./move_dist * line_v;
    double proj_len = dot(line_v, t_v); /*here we use the dot product in the local frame*/
    cur_vec_mag = fabs(proj_len);

    return true;
}


bool
Trajectory::cal_next_pt_RK4_tau_b(double first[2], double second[2], int &face_id, double alpha[3])
{
    ////Using 4th Runge-Kutta method to get the next point
    icVector2 t_v;
    double temp[2] = {0.};

    icVector3 t_gp;
    local_To_global(face_id, first, t_gp);

    /*compute K1*/
    icVector2 V1 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, first[0], first[1]);

    if(length(V1) < 1e-10) return false;


    RK4_step = eulerstep_scalar/4.;
    //RK4_step = eulerstep_scalar/16.; // for diesel slice

    t_v = V1;
    normalize(t_v);

    temp[0] = first[0] - RK4_step/2*t_v.entry[0];
    temp[1] = first[1] - RK4_step/2*t_v.entry[1];

    /*compute K2*/
    //Get2DBarycentricFacters(face_id, temp[0], temp[1], alpha);
    //icVector2 V2 = GetVectorAtPoints(face_id, alpha);
    local_To_global(face_id, temp, t_gp);
    object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha);
    icVector2 V2 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, temp[0], temp[1]);
    t_v = V2;
    normalize(t_v);

    temp[0] = first[0] - RK4_step/2*t_v.entry[0];
    temp[1] = first[1] - RK4_step/2*t_v.entry[1];

    /*compute K3*/
    local_To_global(face_id, temp, t_gp);
    object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha);
    icVector2 V3 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, temp[0], temp[1]);
    t_v = V3;
    normalize(t_v);

    temp[0] = first[0] - RK4_step*t_v.entry[0];
    temp[1] = first[1] - RK4_step*t_v.entry[1];

    /*compute K4*/
    local_To_global(face_id, temp, t_gp);
    object->get_2D_Barycentric_Facters(face_id, temp[0], temp[1], alpha);
    icVector2 V4 = object->get_Vector_At_Point(face_id, t_gp.entry, alpha, temp[0], temp[1]);

    icVector2 total_v = 1./6.*(V1+2*V2+2*V3+V4);
    t_v = total_v;
    normalize(t_v);

    second[0] = first[0] - RK4_step*t_v.entry[0];
    second[1] = first[1] - RK4_step*t_v.entry[1];

    /*ready for the temporary tau  08/30/07*/
    icVector2 line_v;
    line_v.entry[0] = second[0]-first[0];
    line_v.entry[1] = second[1]-first[1];
    move_dist = length(line_v);
    //normalize(line_v);

    if (move_dist>0)
        line_v = 1./move_dist * line_v;
    double proj_len = dot(line_v, t_v); /*here we use the dot product in the local frame*/
    cur_vec_mag = fabs(proj_len);

    return true;
}


bool
Trajectory::get_next_pt(double first[2],
                        double second[2],
                        int &face_id,
                        double alpha[3],
                        int type,
                        unsigned char opt)
{
    switch (opt)
    {
    case 0:
        return cal_next_point_euler1(first, second, face_id, alpha, type);
    case 1:
        return cal_nextpt_2ndeuler(first, second, face_id, alpha, type);
        //return cal_nextpt_RK2_nonNorm(first, second, face_id, alpha, type);
    case 2:
        return cal_nextpt_RK4(first, second, face_id, alpha, type);

    default:
        return false;
    }
}





////calculate the trajectory in a single triangle
int Trajectory::trace_in_triangle(int &face_id, double globalp[3], int type, int &flag)
{
    int i;
    double alpha[3];
    double cur_point[2], pre_point[2];
    double vert0[3];
    icVector3 VP, globalv;

    Triangle *face = object->tlist.tris[face_id];

    Triangle *pre_f = face;

    ////Temporary curve point array

    CurvePoints *temp_point_list = (CurvePoints*) malloc(sizeof(CurvePoints) * 100);
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
    for(i = 0; i < 100; i++)
    {
        ////1. calculate the barycentric coordinates for current point

        object->get_2D_Barycentric_Facters(face_id, cur_point[0], cur_point[1], alpha);

        ////2. if current point is inside current triangle
        if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1
            && (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1
            && (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
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

            //if(cal_next_point_euler1(pre_point, cur_point, face_id, alpha, type))
            //if(cal_nextpt_2ndeuler(pre_point, cur_point, face_id, alpha, type))
            //if(cal_nextpt_RK4(pre_point, cur_point, face_id, alpha, type))

            if (get_next_pt(pre_point, cur_point, face_id, alpha, type, (unsigned char)(Integrator_opt)))
            {
                ////update the global point

                globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

                //Get global coordinates of the point
                globalp[0] = vert0[0] + globalv.entry[0];
                globalp[1] = vert0[1] + globalv.entry[1];
                globalp[2] = vert0[2] + globalv.entry[2];

            }

            else{  ////the curve reach a singularity
                flag = 1;

                ////Store the record into global line segment array

                if(!store_to_global_line_segs(temp_point_list, NumPoints))
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

        ////3. if the point is out of current triangle
        else{
            double t[2] = {0.};

            int PassVertornot = 0;

            get_next_triangle(face_id, pre_point, cur_point, t, type, PassVertornot, alpha);

            ////update the globalpoint here (Modified on 01/30/07)
            if(PassVertornot > 0)
            {
                if (face_id < 0)
                {
                    free(temp_point_list);
                    return -1;
                }

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

                alpha[vert_new]=1-0.001;
                alpha[(vert_new+1)%3]=0.0005;
                alpha[(vert_new+2)%3]=0.0005;


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
                //// transfer it to the global coordinates
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
                if(!store_to_global_line_segs(temp_point_list, NumPoints))
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

    if(NumPoints > 1)
        store_to_global_line_segs(temp_point_list, NumPoints);

    free(temp_point_list);
    return face_id;
}



/*****************************************************************
store temp curve point array to global line segment array
*****************************************************************/

bool Trajectory::store_to_global_line_segs(CurvePoints *temp, int num)
{
    int i;
    int tempid = nlinesegs;
    icVector3 dis_vec;

    ////if the number of the line segements over the maximum number of the line segments each trajectory can store
    ////extend the space for each trajectory
    if(tempid + num - 1 >= curMaxNumLinesegs)
    {

        if(!extend_line_segments(100))
        {
            return false;
        }
    }

    /*save to the global list*/

    for( i = 0; i < num-1; i++)
    {
        ////Build the line segment
        linesegs[tempid+i].gstart.entry[0] = temp[i].gpx;
        linesegs[tempid+i].gstart.entry[1] = temp[i].gpy;
        linesegs[tempid+i].gstart.entry[2] = temp[i].gpz;
        linesegs[tempid+i].start.entry[0] = temp[i].lpx;
        linesegs[tempid+i].start.entry[1] = temp[i].lpy;

        linesegs[tempid+i].gend.entry[0] = temp[i+1].gpx;
        linesegs[tempid+i].gend.entry[1] = temp[i+1].gpy;
        linesegs[tempid+i].gend.entry[2] = temp[i+1].gpz;
        linesegs[tempid+i].end.entry[0] = temp[i+1].lpx;
        linesegs[tempid+i].end.entry[1] = temp[i+1].lpy;

        ////Use local coordinates to calculate the length
        dis_vec.entry[0] = temp[i+1].lpx - temp[i].lpx;
        dis_vec.entry[1] = temp[i+1].lpy - temp[i].lpy;
        dis_vec.entry[2] = 0;

        linesegs[tempid+i].length = length(dis_vec);

        linesegs[tempid+i].Triangle_ID = temp[i].triangleid;
    }

    nlinesegs = tempid + num - 1;
    return true;
}



/*****************************************************************
store temp curve point array to global array
Entry:
Output:
*****************************************************************/
void Trajectory::get_next_triangle(int &face_id, double pre[2], double cur[2], double param_t[2], int type,
                                   int &PassVertornot, double alpha[3])
{
    int which_edge = -1;

    int prev_face_id = face_id;

    Triangle *prev_face = object->tlist.tris[face_id];

    Vertex *vert = NULL;

    PassVertornot = 0;

    ////We should put pass vertex testing here before testing crossing edge
    int passornot = 0;
    cross_a_vertex(face_id, cur, pre, type, passornot);
    if(passornot > 0)
    {
        PassVertornot = passornot;
        return ;
    }

    face_id = prev_face_id;  //////added on 06/08/05

    //CrossBoundary2(pre, cur, face_id, which_edge, param_t); //06/29/06
    cross_boundary(pre, cur, face_id, alpha, which_edge, param_t);

    if(param_t[0] == -1 && param_t[1] == -1)
        face_id = prev_face_id;   ////something wrong here


    ////if not passing a vertex, judge which triangle it will enter later
    pass_edge(face_id, which_edge);

}



////New routine for crossing vertex testing
void Trajectory::cross_a_vertex(int &face_id, double cur_p[2], double pre_p[2],int type, int &passornot)
{
    //	double tx, ty;
    double vert[2];
    int alpha_index = 0;
    double max_alpha ;
    int newtriangleid = 0;
    Vertex* crossVert;

    Triangle *face = object->tlist.tris[face_id];


    ////New way to get the possible crossing vertex
    icVector2 test_dir;
    test_dir.entry[0] = cur_p[0] ;
    test_dir.entry[1] = cur_p[1] ;
    max_alpha = length(test_dir);
    alpha_index = 0;

    test_dir.entry[0] = cur_p[0] - face->x1;
    test_dir.entry[1] = cur_p[1] ;
    if(length(test_dir) < max_alpha)
    {
        max_alpha = length(test_dir);
        alpha_index = 1;
    }

    test_dir.entry[0] = cur_p[0] - face->x2;
    test_dir.entry[1] = cur_p[1] - face->y2;
    if(length(test_dir) < max_alpha)
    {
        max_alpha = length(test_dir);
        alpha_index = 2;
    }

    crossVert = face->verts[alpha_index];

    if(alpha_index == 0)
    {
        vert[0] = vert[1] = 0;
    }

    else if(alpha_index == 1)
    {
        vert[0] = face->x1;
        vert[1] = 0;
    }

    else
    {
        vert[0] = face->x2;
        vert[1] = face->y2;
    }


    double A, B, C;
    A = pre_p[1] - cur_p[1];
    B = cur_p[0] - pre_p[0];
    C = (pre_p[0]*cur_p[1] - cur_p[0]*pre_p[1]);

    double pending = A*vert[0] + B*vert[1] + C;

    if(fabs(pending) <= 1e-8) ////passing the vertex
    {
        pass_vertex(crossVert->index, newtriangleid, type);

        if (newtriangleid < 0)  // perform another test to locate the next triangle to continue...
            pass_vertex(crossVert->index, newtriangleid, type);

        face_id = newtriangleid;
        passornot = alpha_index+1;
        return;
    }

    passornot = 0;
}



void Trajectory::pass_vertex(int vert_id, int &theone, int type)
{
    Vertex *vert = object->vlist.verts[vert_id];

    int i;

    double vang;
    //Face *face;
    Corner *c;
    int NewTID = -1;

    int orient;

    ////for 3D surfaces, we should not use the following angle directly!
    vang = vert->t_angle;

    if(type == 1)
        vang += M_PI;

    if(vang < 0)
        vang += 2*M_PI;

    if(vang > 2*M_PI)
        vang -= 2*M_PI;

    ////Get the orientation of the angle allocation
    //orient = Object.clist[vert->Corners[0]]->orient;
    //orient = object->orient;
    //orient = object->get_orientation(vert, vang, vert->corners[0]->edge[1]);

    //orient = 0; // temporary solution at 07/01/2009! This assumes that the mesh has been oriented!

    if(vert->corners[0]->orient)
        orient=1;
    else
        orient=0;

    for( i = 0; i < vert->ncorners; i++)
    {

        c = vert->corners[i];

        ////first, we check which angle area the vector on the vertex will fall in
        if(orient > 0)
        {
            if(c->BeginAng > c->EndAng)
            {
                if((vang >= c->BeginAng && vang < 2*M_PI)|| (vang <= c->EndAng && vang >= 0))
                {
                    NewTID = i;
                    break;
                }
            }
            else{
                if(vang >= c->BeginAng && vang <= c->EndAng)
                {
                    NewTID = i;
                    break;
                }
            }
        }
        else{
            if(c->BeginAng < c->EndAng)
            {
                if((vang <= c->BeginAng && vang >= 0)|| (vang >= c->EndAng && vang < 2*M_PI))
                {
                    NewTID = i;
                    break;
                }
            }
            else{
                if(vang <= c->BeginAng && vang >= c->EndAng)
                {
                    NewTID = i;
                    break;
                }
            }
        }
    }

    if(NewTID < 0) //avoid crash!
    {
        theone = -1;
        return;
    }


    theone = vert->corners[NewTID]->t;

}


void Trajectory::pass_vertex_2(int vert_id, int &theone, int type)
{
    Vertex *vert = object->vlist.verts[vert_id];

    icVector3 v_vec = vert->raw_vec;
    if (type == 1)
        v_vec = -v_vec;

    int i;

    for (i=0; i<vert->ntris; i++)
    {
        Triangle *t = vert->tris[i];

        if (dot (vert->raw_vec, t->normal)==0) continue;

        // project the vector into the local frame
        icVector3 proj_vec = v_vec - dot(t->normal, v_vec) * t->normal;
        double a = dot(proj_vec, t->LX);
        double b = dot(proj_vec, t->LY);

        icVector2 tmp_vec(a,b);
        normalize(tmp_vec);

        icVector3 p_v(vert->x, vert->y, vert->z);

        double p_x = dot(p_v, t->LX);
        double p_y = dot(p_v, t->LY);

        double smallest_edge_len = 1.e8;
        for(int j=0; j<3; j++)
        {
            if (t->edges[j]->length < smallest_edge_len)
                smallest_edge_len = t->edges[j]->length;
        }

        // move along the direction of the project vector
        double tp_x = p_x + smallest_edge_len/4. * tmp_vec.entry[0];
        double tp_y = p_y + smallest_edge_len/4. * tmp_vec.entry[1];

        // determine whether the new point falls in the same triangle or not
        double alpha[3] = {0.};

        object->get_2D_Barycentric_Facters(t->index, tp_x, tp_y, alpha);
        if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1
            && (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1
            && (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
        {
            theone = t->index;
            return;
        }
    }
    theone = -1;
}

/********************************************************************
If we have already judge which edge the curve will cross
We can use the edge information to get next triangle
********************************************************************/
void Trajectory::pass_edge(int &face_id, int which_edge)
{
    ////Using edge information to get next triangle

    Triangle *face = object->tlist.tris[face_id];

    Edge *theone ;
    Vertex * vertindex;

    for(int i = 0; i < 3; i++)
    {
        vertindex = face->verts[which_edge];

        theone = face->edges[i];

        if(theone->verts[0] != vertindex && theone->verts[1] != vertindex)
            break;
    }

    if(theone->tris[0]!=NULL && theone->tris[0]->index != face_id)
        face_id = theone->tris[0]->index;
    else if(theone->tris[1]!=NULL && theone->tris[1]->index != face_id)
        face_id = theone->tris[1]->index;
    else
        face_id = -1;

}


//void Trajectory::cross_boundary(double pre[2], double cur[2], int face_id, int &which_edge, double t[2])
//{
//	Triangle *face = object->tlist.tris[face_id];
//	double theta1, theta2, theta3, theta;
//	icVector2 ang_vec;
//
//	which_edge = -1;
//
//	//calculate theta1
//	ang_vec.entry[0] = - pre[0];
//	ang_vec.entry[1] = - pre[1];
//	theta1 = atan2(ang_vec.entry[1], ang_vec.entry[0]);
//	if(theta1 < 0)
//		theta1 += 2*M_PI;
//
//	//calculate theta2
//	ang_vec.entry[0] = face->x1 - pre[0];
//	ang_vec.entry[1] = - pre[1];
//	theta2 = atan2(ang_vec.entry[1], ang_vec.entry[0]);
//	if(theta2 < 0)
//		theta2 += 2*M_PI;
//
//	//calculate theta3
//	ang_vec.entry[0] = face->x2 - pre[0];
//	ang_vec.entry[1] = face->y2 - pre[1];
//	theta3 = atan2(ang_vec.entry[1], ang_vec.entry[0]);
//	if(theta3 < 0)
//		theta3 += 2*M_PI;
//
//	//calculate theta
//	ang_vec.entry[0] = cur[0] - pre[0];
//	ang_vec.entry[1] = cur[1] - pre[1];
//	theta = atan2(ang_vec.entry[1], ang_vec.entry[0]);
//	if(theta < 0)
//		theta += 2*M_PI;
//
//	//let the angle begin from theta1
//	theta2 -= theta1;
//	if(theta2 < 0)
//		theta2 += 2*M_PI;
//
//	theta3 -= theta1;
//	if(theta3 < 0)
//		theta3 += 2*M_PI;
//
//	theta  -= theta1;
//	if(theta < 0)
//		theta += 2*M_PI;
//
//	theta1 = 0;
//
//	if(face->xy[2][1] < 0) //clockwise orientation
//	{
//		if(theta >= theta1 && theta < theta3)
//			which_edge = 1;  //v0v2 edge
//
//		else if(theta >= theta3 && theta < theta2)
//			which_edge = 0;  //v1v2 edge
//
//		else
//			which_edge = 2;  //v0v1 edge
//	}
//
//	else //counter clockwise orientation
//	{
//		if(theta >= theta1 && theta < theta2)
//			which_edge = 2; //v0v1 edge
//
//		else if(theta >= theta2 && theta < theta3)
//			which_edge = 0; //v1v2 edge
//
//		else
//			which_edge = 1; //v0v2 edge
//	}
//
//	if(which_edge == 0)
//	{
//		get_intersection(pre, cur, face->xy[1], face->xy[2], t);
//		//GetIntersection(pre, cur, face->xy[1], face->xy[2], t);
//		cur[0] = face->xy[1][0] + t[1] * (face->xy[2][0] - face->xy[1][0]);
//		cur[1] = face->xy[1][1] + t[1] * (face->xy[2][1] - face->xy[1][1]);
//	}
//
//	else if(which_edge == 1)
//	{
//		get_intersection(pre, cur, face->xy[2], face->xy[0], t);
//		//GetIntersection(pre, cur, face->xy[2], face->xy[0], t);
//		cur[0] = face->xy[2][0] + t[1] * (face->xy[0][0] - face->xy[2][0]);
//		cur[1] = face->xy[2][1] + t[1] * (face->xy[0][1] - face->xy[2][1]);
//	}
//
//	else if(which_edge == 2)
//	{
//		get_intersection(pre, cur, face->xy[0], face->xy[1], t);
//		//GetIntersection(pre, cur, face->xy[0], face->xy[1], t);
//		cur[0] = face->xy[0][0] + t[1] * (face->xy[1][0] - face->xy[0][0]);
//		cur[1] = face->xy[0][1] + t[1] * (face->xy[1][1] - face->xy[0][1]);
//	}
//}


void  Trajectory::cross_boundary(double pre[2], double cur[2],
                                int face_id, double alpha[3],
                                int &which_edge, double t[2])
{
    Triangle *face = object->tlist.tris[face_id];

    if(alpha[0] < 0 && alpha[1] < 0)
    {
        double p0[2] = {0, 0};
        double p1[2] = {face->x2, face->y2};
        if(get_intersection(pre, cur, p0, p1, t)==1)
        {
            which_edge = 1;
            cur[0] = p0[0] + t[1] * (p1[0] - p0[0]);
            cur[1] = p0[1] + t[1] * (p1[1] - p0[1]);
            return;
        }

        //Calculate the intersection with edge v1v2

        p0[0] = face->x1; p0[1] = 0;

        if(get_intersection(pre, cur, p0, p1, t)==1)
        {
            which_edge = 0;
            cur[0] = p0[0] + t[1] * (p1[0] - p0[0]);
            cur[1] = p0[1] + t[1] * (p1[1] - p0[1]);
            return;
        }
    }

    else if(alpha[0] < 0 && alpha[2] < 0)
    {
        //Calculate the intersection with edge v0v1
        double p0[2] = {0, 0};
        double p1[2] = {face->x1, 0};

        if(get_intersection(pre, cur, p0, p1, t)==1)
        {
            which_edge = 2;
            cur[0] = p0[0] + t[1] * (p1[0] - p0[0]);
            cur[1] = p0[1] + t[1] * (p1[1] - p0[1]);
            return;
        }

        //Calculate the intersection with edge v1v2
        p0[0] = p1[0]; p0[1] = p1[1];
        p1[0] = face->x2; p1[1] = face->y2;

        if(get_intersection(pre, cur, p0, p1, t)==1)
        {
            which_edge = 0;
            cur[0] = p0[0] + t[1] * (p1[0] - p0[0]);
            cur[1] = p0[1] + t[1] * (p1[1] - p0[1]);
            return;
        }
    }

    else if(alpha[1] < 0 && alpha[2] < 0)
    {
        //Calculate the intersection with edge v0v1
        double p0[2] = {0, 0};
        double p1[2] = {face->x1, 0};

        if(get_intersection(pre, cur, p0, p1, t)==1)
        {
            which_edge = 2;
            cur[0] = p0[0] + t[1] * (p1[0] - p0[0]);
            cur[1] = p0[1] + t[1] * (p1[1] - p0[1]);
            return;
        }

        //Calculate the intersection with edge v0v2
        p1[0] = face->x2; p1[1] = face->y2;

        if(get_intersection(pre, cur, p0, p1, t)==1)
        {
            which_edge = 1;
            cur[0] = p0[0] + t[1] * (p1[0] - p0[0]);
            cur[1] = p0[1] + t[1] * (p1[1] - p0[1]);
            return;
        }
    }

    else if(alpha[0] < 0)
    {
        double p0[2] = {face->x1, 0};
        double p1[2] = {face->x2, face->y2};

        which_edge = 0;
        get_intersection(pre, cur, p0, p1, t);
        cur[0] = p0[0] + t[1] * (p1[0] - p0[0]);
        cur[1] = p0[1] + t[1] * (p1[1] - p0[1]);
        return;
    }

    else if(alpha[1] < 0)
    {
        double p0[2] = {face->x2, face->y2};
        double p1[2] = {0, 0};
        which_edge = 1;
        get_intersection(pre, cur, p0, p1, t);
        cur[0] = p0[0] + t[1] * (p1[0] - p0[0]);
        cur[1] = p0[1] + t[1] * (p1[1] - p0[1]);
        return;
    }

    else if(alpha[2] < 0)
    {
        double p0[2] = {0, 0};
        double p1[2] = {face->x1, 0};
        which_edge = 2;
        get_intersection(pre, cur, p0, p1, t);
        cur[0] = p0[0] + t[1] * (p1[0] - p0[0]);
        cur[1] = p0[1] + t[1] * (p1[1] - p0[1]);
        return;
    }
}


/***************************************************************
New method to calculate the intersection of two line segments
****************************************************************/

/* meaning of return value
 0----Intersection dosn't exists
 1----Intersection exists.
 2----two line segments are parallel.
 3----two line segments are collinear, but not overlap.
 4----two line segments are collinear, and share one same end point.
 5----two line segments are collinear, and overlap.
*/

int Trajectory::get_intersection(double PointA[2], double PointB[2],
                                 double PointC[2], double PointD[2], double t[2])
{

    double delta;
    double t1,t2;
    double a,b,c,d;
    double xba,yba,xdc,ydc,xca,yca;

    xba=PointB[0]-PointA[0];    yba=PointB[1]-PointA[1];
    xdc=PointD[0]-PointC[0];    ydc=PointD[1]-PointC[1];
    xca=PointC[0]-PointA[0];    yca=PointC[1]-PointA[1];

    delta=xba*ydc-yba*xdc;
    t1=xca*ydc-yca*xdc;
    t2=xca*yba-yca*xba;

    if(delta!=0)
    {
        t[0]=t1/delta;   t[1]=t2/delta;
        /*two segments intersect (including intersect at end points)*/
        //if ( t[0]<=1 && t[0]>=0 && t[1]<=1 && t[1]>=0 ) return 1;
        if ( t[0]<=1 && (t[0]>=0 || fabs(t[0]) < 1e-8)
            && t[1]<=1 && (t[1]>=0 || fabs(t[1]) < 1e-8 ))
            return 1;
        else return 0;
    }

    else
    {
        /* AB & CD are parallel. */
        if ( (t1!=0) && (t2!=0) ) return 2;

        /* when AB & CD are collinear */

        /*if AB isn't a vertical line segment, project to x-axis */
        if(PointA[0]!=PointB[0])
        {
            a=MIN(PointA[0],PointB[0]); b=MAX(PointA[0],PointB[0]);
            c=MIN(PointC[0],PointD[0]); d=MAX(PointC[0],PointD[0]);

            if ( (d<a) || (c>b) ) return  3;
            else if( (d==a) || (c==b) ) return 4;
            else return 5;
        }

        else         /* if AB is a vertical line segment, project to y-axis */
        {

            a=MIN(PointA[1],PointB[1]); b=MAX(PointA[1],PointB[1]);
            c=MIN(PointC[1],PointD[1]); d=MAX(PointC[1],PointD[1]);

            if( (d<a) || (c>b) ) return  3;
            else if( (d==a) || (c==b) ) return 4;
            else return 5;
        }
    }
}

/*get the flow length of the streamline*/
double Trajectory::get_length()
{
    int i;
    double len = 0;
    for(i = 0 ; i < nlinesegs; i++)
        len += linesegs[i].length;
    return len;
}


/*Calculate a single trajectory
*/
/*------------------------------------------------------------
Calculate the trajectory using local frame
Need to extend to 3D surfaces
-------------------------------------------------------------*/
void Trajectory::cal_one_traj(int face_id, double x, double y, double z, int type)
{
    int i;
    int flag = 0;
    double globalp[3];
    int pre_face, cur_face;

    pre_face = cur_face = face_id;

    globalp[0] = x;   globalp[1] = y;    globalp[2] = z;
    int NUMTRACETRIS = (int)sqrt((double)object->tlist.ntris);

    for(i = 0; i < 10*NUMTRACETRIS/*NUMTRACINGTRIANGLE*/; i++)
    {

        if(cur_face == -1)
        {
            return;
        }

        pre_face = cur_face;
        cur_face = trace_in_triangle(cur_face, globalp, type, flag); ////0 means always forward tracing here

        if(flag > 0 || pre_face == cur_face )
        {
            return;
        }

    }
}


/****************************************************************
Fix the beginning position of the calculation of a separatrix
1/1/06
****************************************************************/

void TrajectoryList::cal_startpt_sep(int triangleID, icVector3 sep_vector, double saddle_ce[3], double newpos[3])
{
    int count = 2;
    Triangle *t = object->tlist.tris[triangleID];

    ////Get the first position of the beginning point
    newpos[0] = saddle_ce[0] + SEPARATRIXSTEP * sep_vector.entry[0];
    newpos[1] = saddle_ce[1] + SEPARATRIXSTEP * sep_vector.entry[1];
    newpos[2] = saddle_ce[2] + SEPARATRIXSTEP * sep_vector.entry[2];

    while(count <= 500)
    {
        if(t->fall_in_the_tri(newpos))
            return;

        newpos[0] = saddle_ce[0] + SEPARATRIXSTEP/(/*2**/count) * sep_vector.entry[0];
        newpos[1] = saddle_ce[1] + SEPARATRIXSTEP/(/*2**/count) * sep_vector.entry[1];
        newpos[2] = saddle_ce[2] + SEPARATRIXSTEP/(/*2**/count) * sep_vector.entry[2];

        count+=1;
    }


}

/*can we put the separatrices calculation here*/

void Polyhedron::cal_separatrices()
{
    int i;

    icVector3 g_out, g_in, LX, LY;
    double newpos[3] = {0.};
    int flag = -1;

    /*
    we may first estimate the number of trjactories based on the number of saddles
    */

    int predict_num = 0;

    for(i = 0; i < slist.nsingularities; i++)
    {
        if(slist.slist[i]->type == SADDLE)
            predict_num++;   /*count the number of saddles*/
    }

    /*
    Allocate the initial number of trajectories
    */
    if(separatrices != NULL)
        delete separatrices;

    //FILE *fp;
    //if(DebugOn == 1)
    //{
    //	fp = fopen("detect_process_log.txt", "a");
    //	fprintf(fp, "successfully release the old separatrices!\n");
    //	fclose(fp);
    //}

    separatrices = new TrajectoryList(predict_num*4); //each saddle will have four separatrices
    separatrices->curMaxNumTrajs = predict_num*4;
    separatrices->ntrajs = 0;

    //if(DebugOn == 1)
    //{
    //	fp = fopen("detect_process_log.txt", "a");
    //	fprintf(fp, "successfully reallocte the separatrices!\n");
    //	fprintf(fp, "current number of separatrices is %d!\n", predict_num);
    //	fclose(fp);
    //}

    /**/
    for(i = 0; i < slist.nsingularities; i++)
    {
        if(slist.slist[i]->type == SADDLE)
        {

            ////Get the local frame of the triangle
            LX = tlist.tris[slist.slist[i]->TriangleID]->LX;
            LY = tlist.tris[slist.slist[i]->TriangleID]->LY;

            ////transfer the eigenvectors to global space
            g_out = slist.slist[i]->loutgoing.entry[0] * LX
                    + slist.slist[i]->loutgoing.entry[1] * LY;
            g_in  = slist.slist[i]->lincoming.entry[0] * LX
                   + slist.slist[i]->lincoming.entry[1] * LY;

            /*we record the first separatrix of the saddle to the saddle data structure
            05/11/07
            */
            slist.slist[i]->separatrices = separatrices->ntrajs;

            /*1. Positive outgoing separatrix */
            separatrices->cal_startpt_sep(slist.slist[i]->TriangleID, g_out,
                                          slist.slist[i]->gpos.entry, newpos);

            if(separatrices->ntrajs >= separatrices->curMaxNumTrajs)
            {
                if(!separatrices->extend(4))
                {
                    char rout[256], var[256];
                    sprintf(rout, "%s", "TrajectoryList::cal_separatrices");
                    sprintf(var, "%s", "separatrices");

                    //write_mem_error(rout, var, 0);
                    return;
                }
            }

            separatrices->trajs[separatrices->ntrajs] = new Trajectory(separatrices->ntrajs, 100);
            //separatrices->trajs[separatrices->ntrajs] = new Trajectory(g_traj_index, 300);

            separatrices->trajs[separatrices->ntrajs]->cal_one_traj(slist.slist[i]->TriangleID,
                                                                    newpos[0],newpos[1], newpos[2], 0);

            separatrices->trajs[separatrices->ntrajs]->saddleID = i;
            separatrices->ntrajs++;



            /* 2. Positive incoming separatrix */
            separatrices->cal_startpt_sep(slist.slist[i]->TriangleID, g_in,
                                          object->slist.slist[i]->gpos.entry, newpos);

            if(separatrices->ntrajs >= separatrices->curMaxNumTrajs)
            {
                if(!separatrices->extend(4))
                {
                    char rout[256], var[256];
                    sprintf(rout, "%s", "TrajectoryList::cal_separatrices");
                    sprintf(var, "%s", "separatrices");

                    //write_mem_error(rout, var, 0);
                    return;
                }
            }

            separatrices->trajs[separatrices->ntrajs] = new Trajectory(separatrices->ntrajs, 100);
            //separatrices->trajs[separatrices->ntrajs] = new Trajectory(g_traj_index, 300);

            separatrices->trajs[separatrices->ntrajs]->cal_one_traj(slist.slist[i]->TriangleID,
                                                                    newpos[0],newpos[1], newpos[2], 1);

            separatrices->trajs[separatrices->ntrajs]->saddleID = i;
            separatrices->ntrajs++;


            /*3. Negative outgoing separatrix*/

            separatrices->cal_startpt_sep(slist.slist[i]->TriangleID, -g_out,
                                          slist.slist[i]->gpos.entry, newpos);

            if(separatrices->ntrajs >= separatrices->curMaxNumTrajs)
            {
                if(!separatrices->extend(4))
                {
                    char rout[256], var[256];
                    sprintf(rout, "%s", "TrajectoryList::cal_separatrices");
                    sprintf(var, "%s", "separatrices");

                    //write_mem_error(rout, var, 0);
                    return;
                }
            }

            separatrices->trajs[separatrices->ntrajs] = new Trajectory(separatrices->ntrajs, 100);
            //separatrices->trajs[separatrices->ntrajs] = new Trajectory(g_traj_index, 300);

            separatrices->trajs[separatrices->ntrajs]->cal_one_traj(slist.slist[i]->TriangleID,
                                                                    newpos[0],newpos[1], newpos[2], 0);

            separatrices->trajs[separatrices->ntrajs]->saddleID = i;
            separatrices->ntrajs++;


            /*4. Negative incoming separatrix*/
            separatrices->cal_startpt_sep(slist.slist[i]->TriangleID, -g_in,
                                          slist.slist[i]->gpos.entry, newpos);

            if(separatrices->ntrajs >= separatrices->curMaxNumTrajs)
            {
                if(!separatrices->extend(4))
                {
                    char rout[256], var[256];
                    sprintf(rout, "%s", "TrajectoryList::cal_separatrices");
                    sprintf(var, "%s", "separatrices");

                    //write_mem_error(rout, var, 0);
                    return;
                }
            }

            separatrices->trajs[separatrices->ntrajs] = new Trajectory(separatrices->ntrajs, 100);
            //separatrices->trajs[separatrices->ntrajs] = new Trajectory(g_traj_index, 300);

            separatrices->trajs[separatrices->ntrajs]->cal_one_traj(slist.slist[i]->TriangleID,
                                                                    newpos[0],newpos[1], newpos[2], 1);

            separatrices->trajs[separatrices->ntrajs]->saddleID = i;
            separatrices->ntrajs++;
        }
    }

}




/*------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------*/
/*Here are the routines of tracing for periodic orbit detection*/

//bool Trajectory::trace_to_find_recurrence(double g[3], int &triangle,
//													  int type, int scc_index, Edge *chosen_edge,  int &flag)
//{
//	int i;
//	flag = 0;
//	double globalp[3];
//	int pre_face, cur_face;
//
//	pre_face = cur_face = triangle;
//
//	globalp[0] = g[0];   globalp[1] = g[1]; globalp[2] = g[2];
//
//
//	Edge *the_e = NULL;
//	icVector3 intersect;
//	chosen_edge = NULL;
//	Edge *pre_e = NULL;
//
//	double seg_length = 0;
//
//	num_celltriangle = 0;
//	add_To_Cellcycle(triangle);
//
//	int NUMTRACETRIS = (int)sqrt(object->tlist.ntris);
//
//	for(i = 0; i < 10*NUMTRACETRIS && i < object->tlist.ntris; i++)
//	{
//
//		if(cur_face == -1)
//		{
//			flag = 1;       //fail, reaches boundary of the mesh, no cycle exists on this path
//			return false;
//		}
//
//		pre_face = cur_face;
//
//		////Here we need a sub routine that can return the intersection and the corresponding edge
//		////when it enters a new triangle
//
//		pre_e = the_e;
//
//		cur_face = trace_in_triangle_po_detect(globalp, cur_face, type, the_e, intersect, seg_length, flag);
//
//
//		if(flag == 1 || pre_face == cur_face)
//		{
//			flag = 1;      //fail, reach singularity or other error, no cycle has been found on this path
//			return false;
//		}
//
//		if(!is_repeated_elem(morse_decomp->scclist->scccomponents[scc_index]->nodes, cur_face,
//			morse_decomp->scclist->scccomponents[scc_index]->nnodes))
//		{
//			flag = 1;     //fail, goes out of current SCC
//			return false;
//		}
//
//		/* we need to judge whether it forms a loop! 02/08/07 */
//		if(!periodicorbit_detector->cellcycle->add_New(cur_face) && periodicorbit_detector->cellcycle->nelems)
//		{
//			////Move the repeated cycle judgement here 07/25/06
//			int pre_limit_index;
//
//			if(periodicorbit_detector->is_existing_porbit(cur_face, pre_limit_index, 1-type))
//			{
//				flag = 1;
//				return false;
//			}
//
//			//the_e = g_theedge;
//
//			if(the_e == NULL)  //cross a vertex
//				continue;
//
//			if(the_e->num_intersections > 1)
//			{
//				//do not consider the tangent point now
//				if(pre_e == the_e)
//					if(the_e->approx_tangent(seg_length, intersect))
//						continue;
//
//				//validate the intersections on the edge
//				the_e->intersections[1] = intersect;        //save the current intersection to the_e->intersections[1]
//				if(the_e->is_good_intersections(pre_face, scc_index))  //this judgement may not be very good for highly curled field
//				{
//					//g_chosenedge = chosen_edge = the_e;
//					chosen_edge = the_e;
//					chosen_edge->intersections[1] = intersect;
//					triangle = cur_face;
//					return true;                             //we find a cycle here 07/23/06
//				}
//			}
//
//			the_e->intersections[0] = intersect;  //save the intersection with the edge
//			the_e->num_intersections ++;
//		}
//	}
//
//	return false;  //can not converge!
//}


///*
//Trace inside a triangle for periodic orbit detection
//*/
//int Trajectory::trace_in_triangle_po_detect(double g[3], int &face_id, int type,
//							 Edge *cur_e, icVector3 &intersect, double &seg_length, int &flag)
//{
//
//	//similar to the regular tracing, except that you need to mark those special points too close to the cycle
//	int i;
//	double alpha[3];
//	double cur_point[2], pre_point[2];
//	double vert0[3];
//	icVector3 VP, globalv;
//
//	if(face_id < 0)
//		return -1;
//
//	Triangle *face = object->tlist.tris[face_id];
//	Triangle *pre_f = face;
//
//
//	icVector2 dis;
//	seg_length = 0;
//
//	////initialize
//	VP.entry[0] = globalp[0] - face->verts[0]->x;
//	VP.entry[1] = globalp[1] - face->verts[0]->y;
//	VP.entry[2] = globalp[2] - face->verts[0]->z;
//
//	pre_point[0] = cur_point[0] = dot(VP, face->LX);
//	pre_point[1] = cur_point[1] = dot(VP, face->LY);
//
//	vert0[0] = face->verts[0]->x;   ////for update the global point
//	vert0[1] = face->verts[0]->y;
//	vert0[2] = face->verts[0]->z;
//
//	//globalface = face_id;
//
//	////////////////////////////////////////////////////
//    for(i = 0; i < 60; i++)
//	{
//		////1. calculate the barycentric coordinates for current point
//
//		object->get_2D_Barycentric_Facters(face_id, cur_point[0], cur_point[1], alpha);
//
//		////2. if current point is inside current triangle
//		if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1
//			&& (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1
//			&& (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
//		{
//
//			pre_point[0] = cur_point[0];
//			pre_point[1] = cur_point[1];
//
//			if(cal_next_point_euler1(pre_point, cur_point, face_id, alpha, type))
//			{
//				////update the global point
//
//				globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;
//
//				//Get global coordinates of the point
//				globalp[0] = vert0[0] + globalv.entry[0];
//				globalp[1] = vert0[1] + globalv.entry[1];
//				globalp[2] = vert0[2] + globalv.entry[2];
//
//			}
//
//			else{  ////the curve reach a singularity
//				flag = 1;
//				return face_id;
//			}
//		}
//
//		////3. if the point is out of current triangle
//		else{
//			double t[2] = {0.};
//
//			int PassVertornot = 0;
//
//			int presave_face = face_id;
//
//			cross_a_vertex(face_id, cur_point, pre_point, type, PassVertornot);
//
//			//get_next_triangle(face_id, pre_point, cur_point, t, type, PassVertornot, alpha);
//
//			////update the globalpoint here (Modified on 01/30/07)
//			if(PassVertornot > 0)
//			{
//				//we first need to know which vertex it is in the new triangle 01/30/07
//				Vertex* vertid = pre_f->verts[PassVertornot-1];
//				Triangle *cur_f = object->tlist.tris[face_id];
//				int vert_new = 0;
//				for(int k = 0; k < 3; k++)
//				{
//					if(cur_f->verts[k] == vertid)
//					{
//						vert_new = k;
//						break;
//					}
//				}
//
//				alpha[vert_new]=1-0.0001;
//				alpha[(vert_new+1)%3]=0.00005;
//				alpha[(vert_new+2)%3]=0.00005;
//
//				/* Get the new cur_point */
//				cur_point[0] = alpha[1]*cur_f->x1+alpha[2]*cur_f->x2;
//				cur_point[1] = alpha[2]*cur_f->y2;
//
//				globalv = cur_point[0] * cur_f->LX + cur_point[1] * cur_f->LY;
//
//				globalp[0] = cur_f->verts[0]->x + globalv.entry[0];
//				globalp[1] = cur_f->verts[0]->y + globalv.entry[1];
//				globalp[2] = cur_f->verts[0]->z + globalv.entry[2];
//				face=cur_f;
//
//				cur_e = NULL;  //do not consider passing vertex cases now!!! 07/23/06
//			}
//
//			else{
//				face_id = presave_face;
//
//				int which_edge = -1;
//
//				cross_boundary(pre_point, cur_point, face_id, alpha, which_edge, t);
//
//				pass_edge(face_id, which_edge);
//
//				//// transfer it to the global coordinates
//				globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;
//
//				intersect.entry[0] = globalp[0] = vert0[0] + globalv.entry[0];
//				intersect.entry[1] = globalp[1] = vert0[1] + globalv.entry[1];
//				intersect.entry[2] = globalp[2] = vert0[2] + globalv.entry[2];
//
//				////Get the corresponding edge
//				for(int k = 0; k < 3; k++)
//				{
//					int vertindex = face->verts[which_edge];
//
//					cur_e = face->edges[k];
//					//g_theedge = cur_e;
//					periodicorbit_detector->chosen_edge = cur_e;
//
//					if(cur_e->verts[0]->index != vertindex && cur_e->verts[1]->index != vertindex)
//						break;
//				}
//			}
//
//			return face_id;
//		}
//
//	}
//
//	return face_id;
//}
