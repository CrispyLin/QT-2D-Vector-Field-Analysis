/*

Data structures for I/O of polygonal models containing vector field information

Created and Modified by Guoning Chen

*/
#pragma once
#ifndef VFPLY_IO_H
#define VFPLY_IO_H

#include <stddef.h>
#include "PlyLoader.h"

struct Vertex_io {
    float x,y,z;
    int angle_deficit;
    float vx, vy, vz;
    void *other_props;       /* other properties */
};


struct Face_io {
    unsigned char nverts;    /* number of vertex indices in list */
    int *verts;              /* vertex index list */
    void *other_props;       /* other properties */
};

#endif
