#pragma once
#ifndef COMMON_ROUTINES_H
#define COMMON_ROUTINES_H

#include <cmath>

#include "GL_LIB/glew.h"
#include "ExternalDependencies/icMatrix.h"
#include "ExternalDependencies/icVector.h"

int *extend_link(int *edge_link, int Num_edges);
double cal_determinant(icMatrix2x2 mat);
double cal_determinant2x2(double a[][2]);
bool cal_inverse2x2(double a[][2], double inverse_a[][2]);

void cal_eigenvectors(icMatrix2x2 mat, double evalues[2], icVector2 ev[2]);
void cal_eigenval_for_Matrix2x2(icMatrix2x2 mat, icVector2 evec[2]);
bool is_repeated_elem(int *a, int b, int num);
bool get_Cycle(int *triangles, int &num, int cur_tri);
void set_Indentity16(double a[]);
void rightMultiply16(double old_rotmat[16], double rot_mat[16]);
void rightMultiply16_2(double old_rotmat[16], double rot_mat[16]);
void rightMultiply4x4(double old_rotmat[][4], double rot_mat[][4]);
void transform_vec3D(icVector3 &vec, double rot_mat[16]);
void transform_point3D(double p[3], double rot_mat[16]);
void get_rotation(double p[3], icVector3 rot_axis, double rotang, double rotmatrix[16]);
//void get_rotation(double rotang, icVector3 rot_axis, double rotmatrix[16]);
void get_rotation(double rotang, icVector3 rot_axis, double rotmatrix[][4]);
void get_rotation(double rotang, icVector3 rot_axis, double *rotmatrix);
bool del_one_edge_from(int *edges, int &n, int del_e);
void copy_array_Int(int *source, int *dest, int num);
void copy_array_Double(double *source, double *dest, int num);

class DynList_Int{
public:
    //member variables
    int *elems;
    int  nelems;
    int  curMaxNum;
    int  extendstep;

    //member functions
    //constructions and destruction
    DynList_Int(int MaxNum = 500)
    {
        elems = new int [MaxNum];
        curMaxNum = MaxNum;
        nelems = 0;
    }
    ~DynList_Int()
    {delete [] elems; }

    inline bool extend(int step = 100)
    {
        int *temp = elems;
        elems = new int[curMaxNum + step];
        if(elems == NULL)
        {
            elems = temp;
            return false;
        }

        for(int i = 0; i < curMaxNum; i++)
            elems[i] = temp[i];

        curMaxNum += step;

        delete [] temp;
        return true;
    }

    inline bool add_New(int t)
    {
        if(is_Full())
        {
            if(!extend())
            {
                return false;
            }
        }

        if(!is_repeated_elem(elems, t, nelems))
        {
            elems[nelems] = t;
            nelems++;
            return true;
        }

        return false;
    }

    inline bool add_New_2(int t)
    {
        if(nelems>=curMaxNum)
        {
            if(!extend())
            {
                return false;
            }
        }

        //if(!is_repeated_elem(elems, t, nelems))
        //{
        //	elems[nelems] = t;
        //	nelems++;
        //	return true;
        //}
        elems[nelems] = t;
        nelems++;
        return true;

        return false;
    }

    inline bool is_Full()
    {
        return nelems>=curMaxNum;
    }

    inline bool  del_Elem(int t)
    {
        int i, pos = 0;
        for(i = 0; i < nelems; i++)
        {
            if(elems[i] == t)
            {
                pos = i;
                break;
            }
        }
        if(pos >= nelems)
            return false;

        for(i = pos; i < nelems-1; i++)
            elems[i] = elems[i+1];

        nelems--;
        return true;
    }


    inline bool  del_Last()
    {
        if(nelems <= 0) return false;
        nelems--;
        return true;
    }
};

#endif /* __COMMON_ROUTINES_H__ */
