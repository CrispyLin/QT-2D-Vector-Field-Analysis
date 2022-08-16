/*
This file contains some routines used by several classes. It is hard to put them into any class.

Created and modified by Guoning Chen
        copyright @2007
    */

#include "Others/common_routines.h"

/***************************************************************
Extend the number of the edge incident to a specific vertex
The function extend one space each time
***************************************************************/

/***************************************************************
Extend the number of the edge incident to a specific vertex
The function extend one space each time
***************************************************************/

int *extend_link(int *edge_link, int Num_edges)
{
    int *temp = edge_link;
    int *temp_edge_link = new int[Num_edges + 1];
    if( Num_edges > 0)
    {
        for(int i = 0; i < Num_edges; i++)
            temp_edge_link[i] = temp[i];
        delete [] temp;
    }

    return temp_edge_link;
}




double cal_determinant(icMatrix2x2 mat)
{
    return (mat.entry[0][0]*mat.entry[1][1]-mat.entry[0][1]*mat.entry[1][0]);
}

double cal_determinant2x2(double a[][2])
{
    return (a[0][0]*a[1][1]-a[0][1]*a[1][0]);
}

/*
compute the inverse matrix of a 2X2 matrix
*/
bool cal_inverse2x2(double a[][2], double inverse_a[][2])
{
    double det;
    if((det=cal_determinant2x2(a)) == 0)  //degenerate case
        return false;
    inverse_a[0][0] = a[1][1]/det;
    inverse_a[1][1] = a[0][0]/det;
    inverse_a[0][1] = -a[0][1]/det;
    inverse_a[1][0] = -a[1][0]/det;
    return true;
}


void cal_eigenvectors(icMatrix2x2 mat, double evalues[2], icVector2 ev[2])
{
    //first calculate the dominant of the matrix, if it is zero, return
    if(abs(cal_determinant(mat)) < 1e-30)
        return;

    ev[0].entry[0] = mat.entry[1][1] - evalues[0];
    ev[0].entry[1] = - mat.entry[1][0];

    ev[1].entry[0] = mat.entry[1][1] - evalues[1];
    ev[1].entry[1] = - mat.entry[1][0];

}


/* The routine of calculate eigen vectors for a given 2x2 matrix */

void cal_eigenval_for_Matrix2x2(icMatrix2x2 mat, icVector2 evec[2])
{
    double la, lb, lc, ld, A, B, C, delta;

    la = mat.entry[0][0];
    lb = mat.entry[0][1];
    lc = mat.entry[1][0];
    ld = mat.entry[1][1];

    double evalues[2] = {0.};

    A = 1;
    B = -(la + ld);
    C = (la * ld - lb * lc);

    delta =  B*B - 4*A*C;

    if(delta >= 0)
    {
        evalues[0] = (-B + sqrt(delta))/2;
        evalues[1] = (-B - sqrt(delta))/2;

        ////find the largest (absolute) eigen values, and store it to the first element
        //if(abs(evalues[0]) < abs(evalues[1]))
        //{
        double temp = evalues[0];
        evalues[0] = evalues[1];
        evalues[1] = temp;
        //}

        //for real eigen values, we calculate the eigen vectors of it
        cal_eigenvectors(mat, evalues, evec);

        normalize(evec[0]);
        normalize(evec[1]);
    }

    else
    {
        evec[0].set(0, 0);
        evec[1].set(0, 0);
    }
}


/*
Judge whether an element has already been stored in an integer array
This can be extend to a template!
*/
bool is_repeated_elem(int *a, int b, int num)
{
    int i;
    for(i = 0; i < num; i++)
    {
        if(a[i] == b)
            return true;
    }
    return false;
}

/*
The routine is used to get the cycle from the triangle strip
We may also need to get the streamline starting from the real cycle
New added 02/08/07
*/

bool get_Cycle(int *triangles, int &num, int cur_tri)
{
    if(num == 1)
        return false;

    /* find out cur_tri in the triangles */
    int i;
    for(i = 0; i < num; i++)
    {
        if(cur_tri == triangles[i])
            break;
    }

    if(i >= num)
        return false;

    int n_repeated = num - i;
    int *temp = (int*)malloc(sizeof(int)*(n_repeated+1));

    for(int j = i; j < num; j++)
    {
        temp[j-i] = triangles[j];
    }

    for(i = 0; i < n_repeated; i++)
    {
        triangles[i] = temp[i];
    }
    num = n_repeated;

    free(temp);
    return true;
}


void set_Indentity16(double mat[])
{
    for(int i = 0; i < 16; i++)
        mat[i] = 0.;

    //set diagonal entries
    mat[0] = mat[5] = mat[10] = mat[15] = 1.;
}

/*
Multiply two 4x4 matrices.
The result will be saved in old_rotmat
*/
void rightMultiply16(double old_rotmat[16], double rot_mat[16])
{
    //int i, j, k;
    //double result[16] = {0.};

    //for(i = 0; i < 4; i++)
    //{
    //	for(j = 0; j < 4; j++)
    //	{
    //		for(k = 0; k < 4; k++)
    //			result[i*4+j] += old_rotmat[i*4+k] * rot_mat[k*4+j];
    //	}
    //}

    ////Restore the result to old_rotmat
    //for(i = 0; i < 16; i++)
    //{
    //	old_rotmat[i] = result[i];
    //}

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMultMatrixd ((double*)old_rotmat);
    glMultMatrixd ((double*)rot_mat);
    glGetDoublev(GL_MODELVIEW_MATRIX, (double *)old_rotmat);
    glPopMatrix();
}


void rightMultiply16_2(double old_rotmat[16], double rot_mat[16])
{
    int i, j, k;
    double result[16] = {0.};

    for(i = 0; i < 4; i++)
    {
        for(j = 0; j < 4; j++)
        {
            for(k = 0; k < 4; k++)
                result[i*4+j] += old_rotmat[i*4+k] * rot_mat[k*4+j];
        }
    }

    //Restore the result to old_rotmat
    for(i = 0; i < 16; i++)
    {
        old_rotmat[i] = result[i];
    }

}


/*
Multiply two 4x4 matrices.
The result will be saved in old_rotmat
*/
void rightMultiply4x4(double old_rotmat[][4], double rot_mat[][4])
{
    int i, j, k;
    double result[4][4] = {0.};

    for(i = 0; i < 4; i++)
    {
        for(j = 0; j < 4; j++)
        {
            for(k = 0; k < 4; k++)
                result[i][j] += old_rotmat[i][k] * rot_mat[k][j];
        }
    }

    //save to the old_rotmat
    for(i = 0; i < 4; i++)
        for(j = 0; j < 4; j++)
            old_rotmat[i][j] = result[i][j];
}


/*
Transform a vector according to the given transformation matrix
*/
void transform_vec3D(icVector3 &vec, double rot_mat[16])
{
    icVector3 tmp;

    tmp.entry[0] = rot_mat[0] * vec.entry[0] + rot_mat[4] * vec.entry[1] + rot_mat[8] * vec.entry[2];
    tmp.entry[1] = rot_mat[1] * vec.entry[0] + rot_mat[5] * vec.entry[1] + rot_mat[9] * vec.entry[2];
    tmp.entry[2] = rot_mat[2] * vec.entry[0] + rot_mat[6] * vec.entry[1] + rot_mat[10] * vec.entry[2];

    vec = tmp;
}


/*
Transform a 3D point according to the given transformation matrix
*/
void transform_point3D(double p[3], double rot_mat[16])
{
    double tmp[3] = {0.};

    tmp[0] = rot_mat[0] * p[0] + rot_mat[4] * p[1] + rot_mat[8]  * p[2] + rot_mat[12];
    tmp[1] = rot_mat[1] * p[0] + rot_mat[5] * p[1] + rot_mat[9]  * p[2] + rot_mat[13];
    tmp[2] = rot_mat[2] * p[0] + rot_mat[6] * p[1] + rot_mat[10] * p[2] + rot_mat[14];

    p[0] = tmp[0];
    p[1] = tmp[1];
    p[2] = tmp[2];
}

/*
Use OpenGL matrix stack to get the rotation matrix
*/
void get_rotation(double p[3], icVector3 rot_axis, double rotang, double rotmatrix[16])
{

    /*set up the translation matrix*/
    double translate[4][4] = {0.};
    translate[0][0] = translate[1][1] = translate[2][2] = translate[3][3] = 1.;

    translate[3][0] = p[0]; translate[3][1] = p[1]; translate[3][2] = p[2];

    double t_rotmatrix[4][4] = {0.};
    double reverse_translate[4][4] = {0.};
    reverse_translate[0][0] = reverse_translate[1][1] = reverse_translate[2][2] = reverse_translate[3][3] = 1.;
    /*set up the reversed translate*/
    reverse_translate[3][0] = -p[0]; reverse_translate[3][1] = -p[1]; reverse_translate[3][2] = -p[2];

    get_rotation(rotang, rot_axis, t_rotmatrix);

    rightMultiply4x4(t_rotmatrix, translate);
    rightMultiply4x4(reverse_translate, t_rotmatrix);
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            rotmatrix[i*4+j] = reverse_translate[i][j];
}


void get_rotation(double rotang, icVector3 rot_axis, double rotmatrix[][4])
{
    //rotang = rotang * 3.14159265/180.;
    double c = cos(rotang);
    double s = sin(rotang);
    double t = 1-c;

    normalize(rot_axis);

    double x = rot_axis.entry[0];
    double y = rot_axis.entry[1];
    double z = rot_axis.entry[2];

    for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
            rotmatrix[i][j] = 0;

    rotmatrix[0][0] = t*x*x+c;
    rotmatrix[0][1] = t*x*y+s*z;
    rotmatrix[0][2] = t*x*z-s*y;
    rotmatrix[1][0] = t*x*y-s*z;
    rotmatrix[1][1] = t*y*y+c;
    rotmatrix[1][2] = t*y*z+s*x;
    rotmatrix[2][0] = t*x*z+s*y;
    rotmatrix[2][1] = t*y*z-s*x;
    rotmatrix[2][2] = t*z*z+c;
    rotmatrix[3][3] = 1;

}


void get_rotation(double rotang, icVector3 rot_axis, double *rotmatrix)
{
    double c = cos(rotang);
    double s = sin(rotang);
    double t = 1-c;

    normalize(rot_axis);

    double x = rot_axis.entry[0];
    double y = rot_axis.entry[1];
    double z = rot_axis.entry[2];

    for(int i = 0; i < 16; i++)
        rotmatrix[i] = 0;

    rotmatrix[0] = t*x*x+c;
    rotmatrix[1] = t*x*y+s*z;
    rotmatrix[2] = t*x*z-s*y;
    rotmatrix[4] = t*x*y-s*z;
    rotmatrix[5] = t*y*y+c;
    rotmatrix[6] = t*y*z+s*x;
    rotmatrix[8] = t*x*z+s*y;
    rotmatrix[9] = t*y*z-s*x;
    rotmatrix[10] = t*z*z+c;
    rotmatrix[15] = 1;

}



/*
delete one element from an integer list
can be used to delete one edge from the edge list
*/

bool del_one_edge_from(int *edges, int &n, int del_e)
{
    int i, pos;

    for(i = 0; i < n; i++)
    {
        if(edges[i] == del_e)
        {
            pos = i;
            break;
        }
    }

    if(i >= n)
        return false;

    for(i = pos; i < n-1; i++)
    {
        edges[i] = edges[i+1];
    }

    n--;
    return true;
}



/*Some region operations
Copy one region to a new array
dest <- source
*/
void copy_array_Int(int *source, int *dest, int num)
{
    int i;

    for(i = 0; i < num; i++)
    {
        dest[i] = source[i];
    }
}



/*Some region operations
Copy one region to a new array
dest <- source
*/
void copy_array_Double(double *source, double *dest, int num)
{
    int i;

    for(i = 0; i < num; i++)
    {
        dest[i] = source[i];
    }
}


template <class T>
T * extend_link(T *old, int num)
{
    T *temp = old;
    old = new T[num+1];
    if(num > 0)
    {
        for(int i = 0; i < num; i++)
            old[i] = temp[i];   /*need to overwrite "=" operator for complex objects*/
        delete [] temp;
    }
    return old;
}


template <class T>
bool is_repeated_elem(T *a, T b, int num)
{
    for(int i = 0; i < num; i++)
    {
        if(a[i] == b) /*need to overwrite "==" operator for complex objects*/
            return true;
    }
    return false;
}

template <class T>
T * copy_array(T *source, T *dest, int num)
{
    int i;

    for(i = 0; i < num; i++)
    {
        dest[i] = source[i];
    }
}




