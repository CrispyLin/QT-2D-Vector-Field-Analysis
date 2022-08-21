#include "Windows/VectorFieldWindow.h"

extern Polyhedron* object;

//void output_tex_to_file(GLubyte tex[NPIX][NPIX][3], string tex_name){
//    using namespace std;
//    ofstream outfile;
//    string filename = "QT_" + tex_name + "_" + to_string(version) + ".txt";
//    remove(filename.c_str());
//    outfile.open(filename);
//    outfile << "QT " + tex_name + "\n";
//    for (int i = 0; i < NPIX; i++) {
//        for (int j = 0; j < NPIX; j++) {
//            outfile << i << ", " << j << ": ";
//            for (int k = 0; k < 3; k++) {
//                outfile << (int)tex[i][j][k] << " ";
//            }
//            outfile << "\n";
//        }
//        outfile << "\n\n\n\n\n\n\n\n\n";
//    }
//    outfile.close();

//}


//void    output_pattern_to_file(GLubyte pattern[NPN][NPN][4], string pattern_name){
//    using namespace std;
//    ofstream outfile;
//    string filename = "QT_" + pattern_name + "_" + to_string(version) + ".txt";
//    remove(filename.c_str());
//    outfile.open(filename);
//    outfile << "QT " + pattern_name + "\n";
//    for (int i = 0; i < NPN; i++) {
//        for (int j = 0; j < NPN; j++) {
//            outfile << i << ", " << j << ": ";
//            for (int k = 0; k < 4; k++) {
//                outfile << (int)pattern[i][j][k] << " ";
//            }
//            outfile << "\n";
//        }
//        outfile << "\n\n\n\n\n\n\n\n\n";
//    }
//    outfile.close();
//}


void    display_sel_tri(int tri)
{
    Triangle *t = object->tlist.tris[tri];
    int i;
    glDepthFunc(GL_LEQUAL);
    glBegin(GL_LINE_LOOP);
    for (i=0; i<t->nverts; i++)
    {
        glVertex3f(t->verts[i]->x, t->verts[i]->y, t->verts[i]->z);
    }
    glEnd();
    glDepthFunc(GL_LESS);
}


void    matrix_ident( float m[4][4])
{
    int i;

    for (i = 0; i <= 3; i++) {
        m[i][0] = 0.0;
        m[i][1] = 0.0;
        m[i][2] = 0.0;
        m[i][3] = 0.0;
        m[i][i] = 1.0;
    }
}
