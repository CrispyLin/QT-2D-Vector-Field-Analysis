/*
This file contains the declaration of class RotateReflectField used for the
rotation and reflection operations on the whole vector field.

Created and Modified by Guoning Chen
        */
#pragma once

#include "BuildGeom/Geometry.h"
#include "ExternalDependencies/icVector.h"

extern class Polyhedron;

class RotateReflectField
{
public:
    icVector3 *before;
    Polyhedron *the_obj;

    //RotateReflectField(int nverts)
    //{
    //	before = new icVector3[nverts];
    //}

    RotateReflectField(Polyhedron *object)
    {
        before = (icVector3*)malloc(sizeof(icVector3)*object->vlist.nverts);
        the_obj = object;

        save_old_field(object);
    }


    ~RotateReflectField()
    {
        free(before);
    }

    inline void save_old_field(Polyhedron *object)
    {
        int i;
        for(i = 0; i < object->vlist.nverts; i++)
        {
            before[i] = object->vlist.verts[i]->g_vec;
        }
    }

    inline void restore_cur_field(Polyhedron *object)
    {
        int i;
        for(i = 0; i < object->vlist.nverts; i++)
        {
            object->vlist.verts[i]->g_vec = before[i];
        }
        object->project_to_TangentPlane();
    }

    /*rotate the whole with angle */
    void rotate_wholefield_angle(double angle);
};
