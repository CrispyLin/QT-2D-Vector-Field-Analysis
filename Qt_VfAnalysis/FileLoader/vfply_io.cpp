#include "FileLoader/vfply_io.h"


char *elem_names[] = {
    (char*)"vertex", (char*)"face"
};


PlyProperty vert_props[] = { /* list of property information for a vertex */
                            {(char*)"x", Float32, Float32, offsetof(Vertex_io, x), 0, 0, 0, 0},
                            {(char*)"y", Float32, Float32, offsetof(Vertex_io, y), 0, 0, 0, 0},
                            {(char*)"z", Float32, Float32, offsetof(Vertex_io, z), 0, 0, 0, 0},
                            {(char*)"angle_deficit", Int32, Int32, offsetof(Vertex_io,angle_deficit), 0, 0, 0, 0},
                            {(char*)"nx", Float32, Float32, offsetof(Vertex_io, vx), 0, 0, 0, 0},
                            {(char*)"ny", Float32, Float32, offsetof(Vertex_io, vy), 0, 0, 0, 0},
                            {(char*)"nz", Float32, Float32, offsetof(Vertex_io, vz), 0, 0, 0, 0},
                            };

PlyProperty face_props[] = { /* list of property information for a face */
                            {(char*)"vertex_indices", Int32, Int32, offsetof(Face_io, verts), 1, Uint8, Uint8, offsetof(Face_io, nverts)},
                            };
