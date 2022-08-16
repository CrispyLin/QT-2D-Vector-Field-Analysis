/*
This file contains routines of class Polyhedron that we use to extract fixed points

Created and Modified by Guoning Chen
        copyright @2007
    */


#include "VField.h"
#include "Numerical.h"

#define FILEDEBUG
extern const int DebugOn;

extern Polyhedron *object;
double x_cp, y_cp;
icMatrix2x2 FieldMatrix;
double a, b, c, d, e, f, r1, r2, i1, i2;


bool OppositeVector(icVector2 vec1, icVector2 vec2)
{
    normalize(vec1);
    normalize(vec2);

    if(fabs(dot(vec1, vec2)+1) < 1e-10)
        return true;

    return false;
}


void TiltVectorsForATriangle(int triangle)
{
    int i, j;

    icVector2 vec1, vec2;

    Triangle *face = object->tlist.tris[triangle];

    for(i = 0; i < face->num_dir_vecs; i++)
    {
        vec1 = face->dir_vec[i];

        for(j = i+1; j < face->num_dir_vecs; j++)
        {
            vec2 = face->dir_vec[j];

            if(OppositeVector(vec1, vec2))
            {
                ////tilt vec2 only
                icVector2 temp;
                temp.entry[0] = cos(0.05)*vec2.entry[0] - sin(0.05)*vec2.entry[1];
                temp.entry[1] = sin(0.05)*vec2.entry[0] + cos(0.05)*vec2.entry[1];

                //vec2 = temp;
                face->dir_vec[j] = temp;
            }
        }
    }
}

/*------------------------------------------------------------------------------------------*/
/**
The following implement the routines for the fixed point related functionality in the
the Triangle class
**/

/*
This routine judge whether a triangle contains a fixed point or not
*/
bool Triangle::contain_fixedpt()
{
    int j;
    Triangle *face = this;
    icVector2 VertsVector;
    double *Polar_ang, *angle, ang_sum = 0;
    icVector2 tempV1, tempV2;  // to store the two vectors for the vertex that has angle deficit
    int direct_id = 0;         // index to read the directional vectors
    double vx[3], vy[3];       //

    /*  remove the singularities at the boundary triangles  */
    bool is_bound_tri = false;
    for (j=0; j<3; j++)
    {
        if (face->edges[j]->tris[0] == nullptr || face->edges[j]->tris[1] == nullptr)
            is_bound_tri = true;
    }
    if (is_bound_tri)
        return false;

    Polar_ang = new double[face->num_dir_vecs];
    angle = new double[face->num_dir_vecs];

    ang_sum = 0;
    direct_id = 0;

    ////we need to judge whether the singularity is on the vertex or not  04/18/05
    for(j = 0; j < 3; j++)
    {
        if(face->verts[j]->angle_deficit)
        {
            ////get the two vectors on the two edges
            tempV1 = face->dir_vec[direct_id];
            direct_id ++;
            tempV2 = face->dir_vec[direct_id];
            direct_id ++;


            /* Do not consider the fixed points at vertices now  */
            if((length(tempV1)< 1e-10) || (length(tempV2) < 1e-10))////set some threshold instead of 0
            {
                delete [] Polar_ang;
                delete [] angle;
                return false;
            }

        }

        else{
            /////This is for the process of Vertex without angle deficit
            tempV1 = VertsVector = face->dir_vec[direct_id];
            direct_id ++;

            /*  Do not consider the fixed points at vertices now  */
            if(length(tempV1) < 1e-10)
            {
                delete [] Polar_ang;
                delete [] angle;
                return false;
            }

            vx[j] = VertsVector.entry[0];
            vy[j] = VertsVector.entry[1];

        }
    }

    ///////////////////////////////////////////////////////////////////////////////


    ////calculate the angle between two vectors
    for(j = 0; j < face->num_dir_vecs; j++)
    {
        tempV1 = face->dir_vec[j];
        normalize(tempV1);
        Polar_ang[j] = atan2(tempV1.entry[1], tempV1.entry[0]);
    }

    ang_sum = 0;
    for(j = 0; j < face->num_dir_vecs; j++)
    {

        if(face->y2 < 0) //clock wise orientation
            angle[j] = Polar_ang[(j-1+face->num_dir_vecs)%face->num_dir_vecs] - Polar_ang[j];
        else
            angle[j] = Polar_ang[j] - Polar_ang[(j-1+face->num_dir_vecs)%face->num_dir_vecs];


        if( angle[j] < -M_PI)
            angle[j] += 2 * M_PI;

        if( angle[j] > M_PI)
            angle[j] -= 2 * M_PI;

        ang_sum += angle[j];
    }


    if(fabs(ang_sum) >= (2 * M_PI - 0.0001))
    {
        delete [] Polar_ang;
        delete [] angle;
        return true;
    }

    delete [] Polar_ang;
    delete [] angle;
    return false;
}


/*------------------------------------------------------------------------------------------*/
/**
The following implement the routines for the fixed point detection functionality in the
the Polyhedron class
**/


void Polyhedron::capture_Singularities()
{
    int i, j;
    Triangle *face;
    icVector2 VertsVector;
    double *Polar_ang, *angle, ang_sum = 0;
    icVector2 tempV1, tempV2;  // to store the two vectors for the vertex that has angle deficit
    int direct_id = 0;         // index to read the directional vectors
    double vx[3], vy[3];       //

    double dotresult = 0;     ////testing variables

    ////Test total Poincare Index
    int poincare_index = 0;
    nfixedpts = 0;

    singular_tri = new int[slist.curMaxNumSingularities];


    slist.nsingularities = 0;

    //////Test codes here, write to file//////
    FILE *fp;
    if(DebugOn == 1){
        fp = fopen("sing_inf.txt","w");
    }
    int positive = 0, negative = 0;

    ////These flags are used for pair cancellation
    ////Initialize all the flags!!!
    for(i = 0; i < tlist.ntris; i++)
    {
        tlist.tris[i]->singularityID = -1;
    }

    /////here we should estimate each triangle under its local frame
    for (i=0; i<tlist.ntris; i++) {

        face = tlist.tris[i];

        /*  remove the singularities at the boundary triangles  */
        bool is_bound_tri = false;
        for (j=0; j<3; j++)
        {
            if (face->edges[j]->tris[0] == nullptr || face->edges[j]->tris[1] == nullptr)
                is_bound_tri = true;
        }
        if (is_bound_tri)
            continue;

        bool has_zero_vert = false;
        for (j=0; j<face->nverts; j++){
            if (length(face->verts[j]->t_vec) < 1.e-8)
            {
                has_zero_vert = true; break;
            }
        }
        if (has_zero_vert) continue;


        Polar_ang = new double[face->num_dir_vecs];
        angle = new double[face->num_dir_vecs];

        ang_sum = 0;
        direct_id = 0;

        //TiltVectorsForATriangle(i);

        ////we need to judge whether the singularity is on the vertex or not  04/18/05
        for(j = 0; j < 3; j++)
        {
            if(face->verts[j]->angle_deficit)
            {
                ////get the two vectors on the two edges
                tempV1 = face->dir_vec[direct_id];
                direct_id ++;
                tempV2 = face->dir_vec[direct_id];
                direct_id ++;


                if((length(tempV1)<= 1e-20) || (length(tempV2) <= 1e-20))////set some threshold instead of 0
                {
                    continue;  // we don't consider the fixed points on vertices right now

                    //if(slist.isFull())
                    if(nfixedpts >= slist.curMaxNumSingularities)
                    {
                        if(!slist.extend())
                        {
                            //MessageBox(nullptr, "fail to allocate memory for singularity list", "", MB_OK);
                            exit(-1);
                        }
                        singular_tri = (int*)realloc(singular_tri, sizeof(int)*(slist.curMaxNumSingularities));
                        if(singular_tri == nullptr)
                        {
                            //MessageBox(nullptr, "fail to allocate memory for singular triangle list", "", MB_OK);
                            exit(-1);
                        }
                    }

                    singular_tri[nfixedpts]=face->index;

                    //face->singularityID = slist.nsingularities;         ////Mark current triangle as one containing singularity
                    //slist.nsingularities ++;
                    face->singularityID = nfixedpts;        ////Mark current triangle as one containing singularity
                    nfixedpts ++;

                    positive ++;
                    goto LL;
                }

            }

            else{
                /////This is for the process of Vertex without angle deficit
                tempV1 = VertsVector = face->dir_vec[direct_id];
                direct_id ++;

                vx[j] = VertsVector.entry[0];
                vy[j] = VertsVector.entry[1];


                if(fabs(vx[j]) <= 1e-20 && fabs(vy[j]) <= 1e-20)  ////set some threshold instead of 0
                {
                    //if(slist.isFull())
                    if(nfixedpts >= slist.curMaxNumSingularities)
                    {
                        if(!slist.extend())
                        {
                            //MessageBox(nullptr, "fail to allocate memory for singularity list", "", MB_OK);
                            exit(-1);
                        }
                        singular_tri = (int*)realloc(singular_tri, sizeof(int)*(slist.curMaxNumSingularities));
                        if(singular_tri == nullptr)
                        {
                            //MessageBox(nullptr, "fail to allocate memory for singular triangle list", "", MB_OK);
                            exit(-1);
                        }
                    }

                    singular_tri[nfixedpts]=face->index;

                    face->singularityID = nfixedpts;        ////Mark current triangle as one containing singularity
                    nfixedpts ++;

                    positive ++;
                    goto LL;
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////////////


        ////calculate the angle between two vectors
        for(j = 0; j < face->num_dir_vecs; j++)
        {
            tempV1 = face->dir_vec[j];

            normalize(tempV1);

            Polar_ang[j] = atan2(tempV1.entry[1], tempV1.entry[0]);

        }

        ang_sum = 0;
        for(j = 0; j < face->num_dir_vecs; j++)
        {

            if(face->y2 < 0) //clock wise orientation
                angle[j] = Polar_ang[(j-1+face->num_dir_vecs)%face->num_dir_vecs] - Polar_ang[j];
            else
                angle[j] = Polar_ang[j] - Polar_ang[(j-1+face->num_dir_vecs)%face->num_dir_vecs];



            if( angle[j] < -M_PI)
                angle[j] += 2 * M_PI;

            if( angle[j] > M_PI)
                angle[j] -= 2 * M_PI;

            ang_sum += angle[j];
        }


    LL:			if(fabs(ang_sum) >= (2 * M_PI - 0.0001))
        {

            //if(slist.isFull())
            if(nfixedpts >= slist.curMaxNumSingularities)
            {
                if(!slist.extend())
                {
                    //MessageBox(nullptr, "fail to allocate memory for singularity list", "", MB_OK);
                    exit(-1);
                }
                singular_tri = (int*)realloc(singular_tri, sizeof(int)*(slist.curMaxNumSingularities));
                if(singular_tri == nullptr)
                {
                    //MessageBox(nullptr, "fail to allocate memory for singular triangle list", "", MB_OK);
                    exit(-1);
                }
            }

            singular_tri[nfixedpts]=face->index;


            face->singularityID = nfixedpts;        ////Mark current triangle as one containing singularity
            nfixedpts ++;


            ////Add the index (at present, we just consider first order critical points
            if(ang_sum > 0){
                poincare_index ++;
                positive ++;
            }
            else{
                poincare_index --;
                negative ++;
            }

            if(DebugOn == 1){
                fprintf(fp, "Triangle : %d,   Index: %f\n", i, ang_sum/(2*M_PI));
                fprintf(fp, "       The first vectors on the vertices: ");
                for( int k = 0; k < face->num_dir_vecs; k ++)
                    fprintf(fp, " %f, ", length(face->dir_vec[k]));

                fprintf(fp, "  \n\n ");
            }

        }

        delete [] Polar_ang;
        delete [] angle;

    }



    //#ifdef FILEDEBUG

    if(DebugOn == 1){

        fprintf(fp, "Positive: %d,   Negative: %d,  Total: %d\n", positive, negative, poincare_index);

        fprintf(fp, "the Euler Characteristic: %d\n", tlist.ntris+vlist.nverts-elist.nedges);

        fclose(fp);
    }
    //#endif


    ////Calculate the coordinates of being found singularities
    if(nfixedpts > 0)
        compute_FixedPts();

    free(singular_tri);
}


/***************************************************************
Compute the accurate coordinates after capture the ID of
the triangles that contain singularities using new method 1/16
***************************************************************/

void Polyhedron::compute_FixedPts(void)
{
    int i, j;
    Triangle *face;
    Vertex *verts;
    double vx[3], vy[3];//
    icVector2 tempV1, tempV2;  // to store the two vectors for the vertex that has angle deficit
    int direct_id = 0;        // index to read the directional vectors
    int switch_flag  = 0;     //switch flag to decide whether we need to switch vx[i] vy[i] or not
    double alpha[3];
    icVector3 gPos;
    icVector2 VertsVector;

    /////The following variables are for calling of Get_M_onRS() routine to locate the singularities
    double Q[2] = {0.};
    double M[2] = {0.};
    icVector2 VM, WQM;

    double epsilon;

    ////Initialize the index of singularity list
    slist.nsingularities = 0;

    /*--------------------------------------------------------------------*/
    ////Testing codes here 11/30/05
    FILE *fp;

    if(DebugOn == 1){
        fp = fopen("singpos.txt", "w");
    }

    /*--------------------------------------------------------------------*/

    ///////This coordinates should be calculated under 2D local frame
    //////Then, transform back to 3D global frame
    for(i = 0; i < nfixedpts; i++)
    {
        //For each being captured triangle, compute the coordinates of singularities inside it
        face = tlist.tris[singular_tri[i]];

        direct_id = 0;
        switch_flag = -1;   ////it may be possible that all three vertices have no angle deficit, so set it as -1

        //For each triangle, calculate the vector for all vertices
        ////For the angle deficit vertex, we need to store the two vectors
        for (j=0; j< face->nverts; j++) {

            verts = face->verts[j];
            ////we should use the interpolation to get the vector on those vertices that have angle deficits
            ////Note that in current model, we always treat the each triangle have one and only one corner that
            ////has angle deficit properties
            ////Because we have already stored the directional angles, so we can just reuse the angles

            if(verts->angle_deficit == 1)
            {
                ////get the two vectors on the two edges
                tempV1 = face->dir_vec[direct_id];
                direct_id ++;
                tempV2 = face->dir_vec[direct_id];
                direct_id ++;

                /////If the angle deficit not happen on vertex 0, we may need to do switch
                switch_flag = j;

                if((length(tempV1)<= 1e-20) || (length(tempV2) <= 1e-20))////set some threshold instead of 0
                {
                    ////This case means that the vertex itself is a singularity
                    ////Add this vertex into the singularity list
                    ////We need to try to avoid this case

                    //get the center coordinates
                    gPos.entry[0] = verts->x;
                    gPos.entry[1] = verts->y;
                    gPos.entry[2] = verts->z;

                    /*-------*/
                    double lx, ly;

                    if(j == 0)
                    {
                        lx = ly = 0;
                    }
                    else if(j == 1)
                    {
                        lx = face->x1; ly = 0;
                    }
                    else
                    {
                        lx = face->x2; ly = face->y2;
                    }

                    add_To_SingularityList(singular_tri[i], gPos, lx, ly);

                    switch_flag = -2;  ////In this case, the vertex that has angle deficit is the singularity
                        ////At present, the programe can not find other singularity!!!!
                    ////Except that we treat it as an vertex without angle deficit, is that correct?


                    if(DebugOn == 1){
                        fprintf(fp, "singularity %d at triangle %d : %f\n",
                                i, singular_tri[i], length(tempV1));
                    }
                }
            }

            else{
                /////This is for the process of Vertex without angle deficit
                VertsVector = face->dir_vec[direct_id];
                direct_id ++;

                vx[j] = VertsVector.entry[0];
                vy[j] = VertsVector.entry[1];

                if(fabs(vx[j]) <= 1e-20 && fabs(vy[j]) <= 1e-20)  ////set some threshold instead of 0
                {
                    ////This case means that the vertex itself is a singularity
                    ////Add this vertex into the singularity list
                    ////We need to try to avoid this case

                    //get the center coordinates
                    gPos.entry[0] = verts->x;
                    gPos.entry[1] = verts->y;
                    gPos.entry[2] = verts->z;


                    /*-------*/
                    double lx, ly;

                    if(j == 0)
                    {
                        lx = ly = 0;
                    }
                    else if(j == 1)
                    {
                        lx = face->x1; ly = 0;
                    }
                    else
                    {
                        lx = face->x2; ly = face->y2;
                    }

                    add_To_SingularityList(singular_tri[i], gPos, lx, ly);

                    if(DebugOn == 1){
                        fprintf(fp, "singularity %d at triangle %d : %f\n",
                                i, singular_tri[i], vx[j]);
                    }
                }
            }
        }

        get_FixedPt_Coord(singular_tri[i], switch_flag, vx, vy, VM, WQM, tempV1, tempV2, Q, M);

        if(switch_flag != -2)
        {
            /////Using the VM and WQM to calculate the singularity point
            double center_x;
            double center_y;

            if(switch_flag == -1)  ////There is no angle deficit
            {
                ///////Get the center coordinates under 3D global frame
                center_x = M[0];
                center_y = M[1];

            }

            else
            {
                /////Using the VM and WQM to calculate the singularity point
                epsilon = length(WQM)/length(VM);

                double m = epsilon / (epsilon + 1);
                center_x = (1-m) * Q[0] + m*M[0];
                center_y = (1-m) * Q[1] + m*M[1];

            }


            ///////Get the center coordinates under 3D global frame
            icVector3 t_gPos = center_x * face->LX + center_y * face->LY;

            gPos.entry[0] = face->verts[0]->x + t_gPos.entry[0];
            gPos.entry[1] = face->verts[0]->y + t_gPos.entry[1];
            gPos.entry[2] = face->verts[0]->z + t_gPos.entry[2];

            add_To_SingularityList(singular_tri[i], gPos, center_x, center_y);

            /* Update the result here 07/04/06 */
            center_x = x_cp;
            center_y = y_cp;

            t_gPos = center_x * face->LX + center_y * face->LY;

            slist.slist[slist.nsingularities-1]->gpos.entry[0] = face->verts[0]->x + t_gPos.entry[0];
            slist.slist[slist.nsingularities-1]->gpos.entry[1] = face->verts[0]->y + t_gPos.entry[1];
            slist.slist[slist.nsingularities-1]->gpos.entry[2] = face->verts[0]->z + t_gPos.entry[2];

            /*--------------------------------------------------------------------*/
            //Test the vectors on these captured singularities 11/30/05
            get_2D_Barycentric_Facters(singular_tri[i], center_x, center_y, alpha);
            icVector2 t_v = get_Vector_At_Point(singular_tri[i], gPos.entry, alpha, center_x, center_y);

            if(DebugOn == 1){
                fprintf(fp, "singularity %d at triangle %d (type %d): %f\n",
                        i, singular_tri[i], slist.slist[slist.nsingularities-1]->type, length(t_v));
            }
            /*--------------------------------------------------------------------*/

        }
    }

    /*--------------------------------------------------------------------*/

    if(DebugOn == 1){
        fclose(fp);
    }
}


void Polyhedron::get_FixedPt_Coord(int TriangleID, int switch_flag,
                                   double vx[3], double vy[3],
                                   icVector2& VM, icVector2 &WQM,
                                   icVector2 tempV1, icVector2 tempV2,
                                   double Q[2], double M[2])
{
    // to store the two vectors for the vertex that has angle deficit
    icVector2 Vr, Vs;
    double R[2], S[2]/*, Q[2], M[2]*/;

    Triangle *face = tlist.tris[TriangleID];

    //////The following codes use binary search to locate the singularities
    if(switch_flag == 0) ////verts[0] has angle deficit
    {
        Q[0] = 0;
        Q[1] = 0;

        //////Note that the order may be wrong!!!!!!!!!!!!!!!!????
        VM = tempV1;    //Get VQR
        WQM = tempV2;   //Get VQS


        R[0] = face->x2;
        R[1] = face->y2;
        Vr = face->dir_vec[3];

        S[0] = face->x1;
        S[1] = 0;
        Vs = face->dir_vec[2];

        get_M_onRS(R, S, Q, M, VM, WQM, Vr, Vs);

    }

    else if(switch_flag == 1) ////verts[1] has angle deficit
    {
        Q[0] = face->x1;
        Q[1] = 0;

        //////Note that the order may be wrong!!!!!!!!!!!!!!!!
        VM = tempV1;
        WQM = tempV2;

        R[0] = 0;
        R[1] = 0;
        Vr = face->dir_vec[0];

        S[0] = face->x2;
        S[1] = face->y2;
        Vs = face->dir_vec[3];

        get_M_onRS(R, S, Q, M, VM, WQM, Vr, Vs);
    }

    else if(switch_flag == 2){ ////verts[2] has angle deficit
        Q[0] = face->x2;
        Q[1] = face->y2;

        //////Note that the order may be wrong!!!!!!!!!!!!!!!!
        VM = tempV1;
        WQM = tempV2;

        R[0] = face->x1;
        R[1] = 0;
        Vr = face->dir_vec[1];

        S[0] = 0;
        S[1] = 0;
        Vs = face->dir_vec[0];

        get_M_onRS(R, S, Q, M, VM, WQM, Vr, Vs);
    }

    else if(switch_flag == -1){

        //all the three vertices have no angle deficit!!!!
        ///so we can use linear interpolation to calculate the coordinates of the singularities
        ///But for consistent, we just use the same way above

        ////Do not use this routine, refer to the calculation in 2D plane
        double px[3], py[3];
        px[0] = 0;
        py[0] = 0;
        px[1] = face->x1;
        py[1] = 0;
        px[2] = face->x2;
        py[2] = face->y2;

        cal_SingularCoord_without_angledeficit(px, py, vx, vy, M);  ////here M is the coordinates of the singularity
    }

    else{
        ////do nothing at this moment
    }
}



////Get the position of the singularity inside the triangle without angle deficit
////the return coordinates are under the local frame of the triangle
void Polyhedron::cal_SingularCoord_without_angledeficit(double px[3], double py[3],
                                                        double vx[3], double vy[3],
                                                        double locl_p[2])
{
    double temp_ratio = 0.;
    double a, b, c, a1, b1;
    double vdx, vdy;

    /////The following codes using bilinear method to get the center of new critical points
    if(vx[2]!=0){
        temp_ratio = vy[2]/vx[2];
    }
    else{
        ///Switch the index of vectors to make sure it is not divided by zero
    }

    ////Calculate the a1 first according to the following equations
    // vdx = a1*vx[0] + b1* vx[1];
    // vdy = a1*vy[0] + b1* vy[1];
    // a1 + b1 = 1;
    // vdy/vdx = temp_ratio;

    a1 = (vy[1] - temp_ratio * vx[1])/
         (vy[1] - vy[0] + temp_ratio*vx[0] - temp_ratio*vx[1]);
    b1 = 1 - a1;

    ///Getting vdx, vdy using a1 and b1
    vdx = a1 * vx[0] + b1 * vx[1];
    vdy = a1 * vy[0] + b1 * vy[1];

    ///Calculate c using following equation
    //(1 - c)*vdx + c*vx[2] = 0
    c = vdx / (vdx - vx[2]); //need to judge whether vdx == vx[2] or not!!!!

    ///Get the a, b
    a = (1 - c)*a1;
    b = (1 - c)*b1;

    ///Get the coordinates of this singularities
    if((a>=0)&&(a<=1)&&(b>=0)&&(b<=1)&&(c>=0)&&(c<=1)){
        locl_p[0] = a * px[0] + b * px[1] + c * px[2];  ////using the barycentric coefficients to calculate
        locl_p[1] = a * py[0] + b * py[1] + c * py[2];  ////the center of the singularity
    }

}


/********************************************************************
The following routine uses binary search to find the singularity
R, S, and Q store the 2D local coordinates
M  stores the 2D local coordinates of the point M we want on the edge RS
VM input is the vector on the edge QR for the Vertex Q with angle deficit
   output is the Vector on the point M on the edge RS
WQM input is the vector on the edge QS for the Vertex Q with angle deficit
   output is the Vector on the opposite direction to VM on the QM
Vr is the Vector on the Vertex R
Vs is the Vector on the Vertex S
********************************************************************/
void Polyhedron::get_M_onRS(double R[2], double S[2],
                            double Q[2], double M[2],
                            icVector2 &VM, icVector2 &WQM,
                            icVector2 Vr, icVector2 Vs)
{
    double beta;

    double OriginR[2], OriginS[2];
    OriginR[0] = R[0];
    OriginR[1] = R[1];
    OriginS[0] = S[0];
    OriginS[1] = S[1];

    icVector2 VQR = VM;
    icVector2 VQS = WQM; /////the vectors for the vertex that has angle deficit (Modified on 5/2/05)

    double ang_VQR = atan2(VQR.entry[1], VQR.entry[0]);
    double ang_VQS = atan2(VQS.entry[1], VQS.entry[0]);

    icVector2 temp_VM, temp_WQM;

    //// angles for edge vectors QR and QS for angled based interpolation
    double Ang_MQS, Ang_RQM, x;

    double theta, theta2, theta3;

    icVector2 QR, QM, QS;
    icVector2 RM, RS;

    RS.entry[0] = S[0] - R[0];
    RS.entry[1] = S[1] - R[1];

    QS.entry[0] = S[0] - Q[0];
    QS.entry[1] = S[1] - Q[1];
    QR.entry[0] = R[0] - Q[0];
    QR.entry[1] = R[1] - Q[1];

    normalize(QS);
    normalize(QR);

    double ang_RQS = acos(dot(QS, QR)); ////store the
    ////for angle based interpolation 4/11/06
    double cos_ang_RQS = dot(QS, QR);
    double cos_ang_MQS, cos_ang_RQM;


    while(1)
    {
        ////Get the middle point of R and S;
        M[0] = (R[0] + S[0])/2.;
        M[1] = (R[1] + S[1])/2.;

        VM = 0.5 * Vr + 0.5 * Vs;


        //////Calculate the angles between vectors QR, QM  and  QM, QS
        QM.entry[0] = M[0] - Q[0];
        QM.entry[1] = M[1] - Q[1];
        normalize(QM);

        Ang_RQM = acos(dot(QM, QR));
        Ang_MQS = acos(dot(QM, QS));
        beta = Ang_RQM / (Ang_RQM + Ang_MQS);

        //WQM = (1-beta) * VQR + beta * VQS;

        /*use angle based interpolation 04/22/07*/
        double t_s = sin(Ang_RQM+Ang_MQS);
        double t_s1 = sin((1-beta)*(Ang_RQM+Ang_MQS));
        double t_s2 = sin(beta*(Ang_RQM+Ang_MQS));
        WQM = t_s1/t_s*VQR+t_s2/t_s*VQS;

        /////remember, we can not normalize the original VM and WQM
        temp_VM = VM;
        temp_WQM = WQM;

        normalize(temp_VM);
        normalize(temp_WQM);

        //////Calculate the angle between the two vectors)
        theta = acos(dot(temp_VM, temp_WQM));

        if(fabs(M_PI - theta) < 1e-10)  ////set the threshold to stop the searching, the two vectors are opposite approximately
        {   ///the point M is what we want
            break;
        }

        ////if it is too small of line segment RS
        else  if( (fabs(M[0] - S[0]) < 1e-10 && fabs(M[1] - S[1]) < 1e-10)
                 || (fabs(M[0] - R[0]) < 1e-10 && fabs(M[1] - R[1]) < 1e-10))
        {
            break;
        }

        else{
            double temp[2];
            temp[0] = M[0];
            temp[1] = M[1];

            icVector2 tempV = VM;

            M[0] = (R[0] + temp[0])/2.;
            M[1] = (R[1] + temp[1])/2.;

            VM = 0.5 * Vr + 0.5 * tempV;


            //////Calculate the angles between vectors QR, QM  and  QM, QS
            QM.entry[0] = M[0] - Q[0];
            QM.entry[1] = M[1] - Q[1];
            normalize(QM);

            Ang_RQM = acos(dot(QM, QR));
            Ang_MQS = acos(dot(QM, QS));
            beta = Ang_RQM / (Ang_RQM + Ang_MQS);

            /////Get the interpolation Vector on the vertex with angle deficit

            //WQM = (1-beta) * VQR + beta * VQS;

            /*use angle based interpolation 04/22/07*/
            double t_s = sin(Ang_RQM+Ang_MQS);
            double t_s1 = sin((1-beta)*(Ang_RQM+Ang_MQS));
            double t_s2 = sin(beta*(Ang_RQM+Ang_MQS));
            WQM = t_s1/t_s*VQR+t_s2/t_s*VQS;

            /////remember, we can not normalize the original VM and WQM
            temp_VM = VM;
            temp_WQM = WQM;

            normalize(temp_VM);
            normalize(temp_WQM);

            //////Calculate the angle between the two vectors)
            theta2 = acos(dot(temp_VM, temp_WQM));

            if(fabs(M_PI - theta2) < 1e-10)  ////set the threshold to stop the searching, the two vectors are opposite approximately
            {   ///the point M is what we want
                break;
            }

            ////if it is too small of line segment RS
            else  if( (fabs(M[0] - S[0]) < 1e-10 && fabs(M[1] - S[1]) < 1e-10)
                     || (fabs(M[0] - R[0]) < 1e-10 && fabs(M[1] - R[1]) < 1e-10))
            {
                break;
            }


            ////calculate the second half part
            M[0] = (S[0] + temp[0])/2.;
            M[1] = (S[1] + temp[1])/2.;

            VM = 0.5 * Vs + 0.5 * tempV;

            QM.entry[0] = M[0] - Q[0];
            QM.entry[1] = M[1] - Q[1];
            normalize(QM);

            Ang_RQM = acos(dot(QM, QR));
            Ang_MQS = acos(dot(QM, QS));
            beta = Ang_RQM / (Ang_RQM + Ang_MQS);

            //New method to calculat the beta 1/13/06

            /////Get the interpolation Vector on the vertex with angle deficit

            //WQM = (1-beta) * VQR + beta * VQS;

            /*use angle based interpolation 04/22/07*/
            t_s = sin(Ang_RQM+Ang_MQS);
            t_s1 = sin((1-beta)*(Ang_RQM+Ang_MQS));
            t_s2 = sin(beta*(Ang_RQM+Ang_MQS));
            WQM = t_s1/t_s*VQR+t_s2/t_s*VQS;

            /////remember, we can not normalize the original VM and WQM
            temp_VM = VM;
            temp_WQM = WQM;

            normalize(temp_VM);
            normalize(temp_WQM);

            //////Calculate the angle between the two vectors)
            theta3 = acos(dot(temp_VM, temp_WQM));

            if(fabs(M_PI - theta3) < 1e-10)  ////set the threshold to stop the searching, the two vectors are opposite approximately
            {   ///the point M is what we want
                break;
            }

            ////if it is too small of line segment RS
            else  if( (fabs(M[0] - S[0]) < 1e-10 && fabs(M[1] - S[1]) < 1e-10)
                     || (fabs(M[0] - R[0]) < 1e-10 && fabs(M[1] - R[1]) < 1e-10))
            {
                break;
            }

            if( fabs(M_PI - theta2) < fabs(M_PI - theta3))
            {
                ////Keep to binary search the first half part
                S[0] = temp[0];
                S[1] = temp[1];
                Vs = tempV;
            }

            else{
                ////keep to binary search the last half part
                R[0] = temp[0];
                R[1] = temp[1];
                Vr = tempV;
            }
        }
    }
}


/****************************************************************
Add the found singularity into the singularity list
and get the type of the singularity
*****************************************************************/
void Polyhedron::add_To_SingularityList(int Triangle_ID, icVector3 gPos, double local_x, double local_y)
{
    Triangle *face = tlist.tris[Triangle_ID];
    double alpha[3];

    //Allocate the memory for current fixed point
    slist.slist[slist.nsingularities] = (Singularity*)malloc(sizeof(Singularity));

    //get the center coordinates
    slist.slist[slist.nsingularities]->gpos = gPos;

    slist.slist[slist.nsingularities]->TriangleID = Triangle_ID;

    /////we may use another way to get the type of the singularities
    get_2D_Barycentric_Facters(Triangle_ID, local_x, local_y, alpha);

    //slist.slist[slist.nsingularities]->alpha[0] = alpha[0];
    //slist.slist[slist.nsingularities]->alpha[1] = alpha[1];
    //slist.slist[slist.nsingularities]->alpha[2] = alpha[2];

    slist.slist[slist.nsingularities]->type
        = get_Type_of_FixedPt(Triangle_ID, alpha); //Judge the type of this singularity

    slist.slist[slist.nsingularities]->Jacobian.set(FieldMatrix); //need to save its Jacobian for cancellation

    if(slist.slist[slist.nsingularities]->type == SADDLE)
    {
        ////Calculate the eigenvectors of the saddle

        slist.slist[slist.nsingularities]->loutgoing.entry[0]\
            = slist.slist[slist.nsingularities]->Jacobian.entry[0][1];
        slist.slist[slist.nsingularities]->loutgoing.entry[1]\
            = r1 - slist.slist[slist.nsingularities]->Jacobian.entry[0][0];//outgoing direction

        slist.slist[slist.nsingularities]->lincoming.entry[0]\
            = r2 - slist.slist[slist.nsingularities]->Jacobian.entry[1][1];
        slist.slist[slist.nsingularities]->lincoming.entry[1]\
            = slist.slist[slist.nsingularities]->Jacobian.entry[0][1]; //incoming direction

        normalize(slist.slist[slist.nsingularities]->loutgoing);
        normalize(slist.slist[slist.nsingularities]->lincoming);

    }

    /* -- Initial the connection list for C-graph */
    slist.slist[slist.nsingularities]->connect_cycles = nullptr;
    slist.slist[slist.nsingularities]->num_connect_cycles = 0;

    slist.slist[slist.nsingularities]->connected = false;

    ////
    face->singularityID = slist.nsingularities;
    slist.nsingularities++;
}

/***********************************************************************************************
Get the types of the singularities by new method.
We now reuse the method from 2D planar vector field,
except that the current triangle is a very tiny triangle that quite close to the singularity
************************************************************************************************/

int Polyhedron::get_Type_of_FixedPt(int TriangleID, double alpha[3])
{
    Triangle *face;

    icVector2 VertsVector[3];

    icVector2 ApproxVec[3];
    icVector2 ApproxCoord[3];
    icVector2 SingCoord;

    double x[3], y[3];

    int i;

    ///These information has been obtained before
    ///Can we reuse them and enhance the performance
    face = tlist.tris[TriangleID];

    //////Using the approximate tiny triangle to judge the type of singularity

    SingCoord = get_A_2D_Point(TriangleID, alpha);  ////Get the 2D coordinates of the singularity

    ////Get the coordinates and vectors on the vertices of the approximate triangle

    ////11/22/05
    get_Neighbor_Triangle_Vectors(TriangleID, alpha,
                                  SingCoord.entry[0], SingCoord.entry[1], ApproxCoord, ApproxVec);

    /////normalize the vector to get more accurate judgement!!!

    ////initialize for calculate the Jacobbian of the piecewise linear vector field
    ////inside the approximate triangle
    ////It seems that we should not treat the singularity as one vertex of this small region
    x[0] = ApproxCoord[0].entry[0];
    y[0] = ApproxCoord[0].entry[1];
    x[1] = ApproxCoord[1].entry[0];
    y[1] = ApproxCoord[1].entry[1];
    x[2] = ApproxCoord[2].entry[0];
    y[2] = ApproxCoord[2].entry[1];

    /////According to the following equations and the new coordinates to get the
    ///components of Jacobian Matrix

    if( SingCoord.entry[0] == 0 && SingCoord.entry[1] == 0){
        c = 0;  //singularity, the vector is equal to 0
        f = 0;

        a = (ApproxVec[0].entry[0] - c)/x[1];
        d = (ApproxVec[0].entry[1] - f)/x[1];

        b = (ApproxVec[1].entry[0] - a*x[2])/y[2];
        e = (ApproxVec[1].entry[1] - d*x[2])/y[2];
    }

    else{

        ////How about to solve a linear system with 6 unknown variables
        ////The following codes use LU decomposition to solve the system  11/22/05
        double coord[3][3],  *inver_coord ;  //inver_coord[3][3];

        for(i = 0; i < 3; i++)
        {
            coord[0][i] = x[i];
            coord[1][i] = y[i];
            coord[2][i] = 1.;
        }


        inver_coord = MatrixOpp((double*)coord, 3, 3);

        icMatrix3x3 result, rightM;

        result.set(ApproxVec[0].entry[0], ApproxVec[1].entry[0], ApproxVec[2].entry[0],
                   ApproxVec[0].entry[1],ApproxVec[1].entry[1],ApproxVec[2].entry[1],
                   1,1,1);
        rightM.set(inver_coord[0], inver_coord[1], inver_coord[2],
                   inver_coord[3], inver_coord[4], inver_coord[5],
                   inver_coord[6], inver_coord[7], inver_coord[8]);


        result.rightMultiply(rightM);

        //FieldMatrix.set(result);

        double M[][2] = {result.entry[0][0], result.entry[0][1], result.entry[1][0], result.entry[1][1]};
        FieldMatrix.set(M);

        a = result.entry[0][0];
        b = result.entry[0][1];
        c = result.entry[0][2];
        d = result.entry[1][0];
        e = result.entry[1][1];
        f = result.entry[1][2];
    }

    //use to calculate the coordinates of the singularity 07/04/06
    x_cp = (f*b - c*e)/(a*e - d*b);
    y_cp = (c*d - a*f)/(a*e - b*d);


    //////////////////////////////////////////////////////////////////////////////
    /////Second way to calculate the a b c d e f values

    //Getting the eigenvalue and eigen vectors here

    double A, B, C, delta;

    A = 1;
    B = -(a + e);
    C = (a * e - b * d);

    delta =  B*B - 4*A*C;

    if( delta >= 0)
    {
        i1 = i2 = 0;
        r1 = (-B + sqrt(delta))/2;
        r2 = (-B - sqrt(delta))/2;
    }
    else
    {
        r1 = r2 = - B/2.;
        i1 = sqrt(-delta)/2;
        i2 = -sqrt(-delta)/2;
    }

    ////////////////////////////////////////////////////////////////
    if( r1 > 0 && r2 > 0 && i1 == 0 && i2 == 0)
        return 1; //it is a source

    else if( r1 < 0 && r2 < 0 && i1 == 0 && i2 == 0)
        return 2; //it is a sink

    else if( r1 * r2 < 0 && i1 == 0 && i2 == 0)
        return 3; //it is a saddle

    else if( r1 == 0 && r2 == 0 && i1 != 0 && i2 != 0)
        return 4; //it is a center

    else if( r1 < 0 && r2 < 0 && i1 != 0 && i2 != 0)
        return 6; //it is an attracting focus

    else if( r1 > 0 && r2 > 0 && i1 != 0 && i2 != 0)
        return 7; //it is a repelling focus

    else
        return 0; //Unknow, there must be some error here!

}


/* Local coordinates to global coordinates transformation */
//void Polyhedron::local_To_Global(int faceid, double locpos[2], icVector3 &glpos)
//{
//	Triangle *face = tlist.tris[faceid];
//
//	icVector3 dir3D = locpos[0] * face->LX + locpos[1] * face->LY;
//
//	glpos.entry[0] = face->verts[0]->x + dir3D.entry[0];
//	glpos.entry[1] = face->verts[0]->y + dir3D.entry[1];
//	glpos.entry[2] = face->verts[0]->z + dir3D.entry[2];
//
//}


/*************************************************************
This routine is used to calculate the vectors on the vertices
of the approximate tiny triangle for singularity type judgement
**************************************************************/

void Polyhedron::get_Neighbor_Triangle_Vectors(int faceid, double alpha[3],
                                               double sincx, double sincy,
                                               icVector2 Coord[3], icVector2 vec[3])
{

    icVector3 glPos;

    icVector2 vec2D;

    double t_alpha[3];
    Triangle *face = tlist.tris[faceid];

    /*---------------------------------------*/
    //Testing routine

    if(!is_Valid_Alpha(alpha))
    {
        int k = 0;
        return;
    }


    vec[0].set(0,0);

    ////the singularity is on one of the vertices
    if(fabs(alpha[0] - 1) < 1e-6 || fabs(alpha[1] - 1) < 1e-6 || fabs(alpha[2] - 1) < 1e-6)
    {
        if(fabs(alpha[0]-1) < 1e-6) ////on 1st vertex
        {
            Coord[0].entry[0] = 0;
            Coord[0].entry[1] = 0;

            ///calculate the other two points along the two edges incident to v0 respectively
            vec2D.entry[0] = face->x1;
            vec2D.entry[1] = 0;

            Coord[1] = Coord[0] + 0.5 * vec2D;  //get the point on the edge v0v1

            ////Get the vector associated with the point
            get_2D_Barycentric_Facters(faceid, Coord[1].entry[0], Coord[1].entry[1], t_alpha);
            local_To_Global(faceid, Coord[1], glPos);
            vec[1] = get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);

            ////
            vec2D.entry[0] = face->x2;
            vec2D.entry[1] = face->y2;

            Coord[2] = Coord[0] + 0.5 * vec2D;  //get the point on the edge v0v2

            ////Get the vector associated with the point
            get_2D_Barycentric_Facters(faceid, Coord[2].entry[0], Coord[2].entry[1], t_alpha);
            local_To_Global(faceid, Coord[2], glPos);
            vec[2] = get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);
        }

        else if(fabs(alpha[1]-1) < 1e-6) //// on 2nd vertex
        {
            Coord[0].entry[0] = face->x1;
            Coord[0].entry[1] = 0;

            ///calculate the other two points along the two edges incident to v0 respectively
            vec2D.entry[0] = face->x2 - face->x1;
            vec2D.entry[1] = face->y2;

            Coord[1] = Coord[0] + 0.5 * vec2D;  //get the point on the edge v2v1

            ////Get the vector associated with the point
            get_2D_Barycentric_Facters(faceid, Coord[1].entry[0], Coord[1].entry[1], t_alpha);
            local_To_Global(faceid, Coord[1], glPos);
            vec[1] = get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);

            ////
            vec2D.entry[0] = - face->x1;
            vec2D.entry[1] = 0;

            Coord[2] = Coord[0] + 0.5 * vec2D;  //get the point on the edge v1v0

            ////Get the vector associated with the point
            get_2D_Barycentric_Facters(faceid, Coord[2].entry[0], Coord[2].entry[1], t_alpha);
            local_To_Global(faceid, Coord[2], glPos);
            vec[2] = get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);
        }

        else //// on 3rd vertex
        {
            Coord[0].entry[0] = face->x2;
            Coord[0].entry[1] = face->y2;

            ///calculate the other two points along the two edges incident to v0 respectively
            vec2D.entry[0] = - face->x2;
            vec2D.entry[1] = - face->y2;

            Coord[1] = Coord[0] + 0.5 * vec2D;  //get the point on the edge v0v2

            ////Get the vector associated with the point
            get_2D_Barycentric_Facters(faceid, Coord[1].entry[0], Coord[1].entry[1], t_alpha);
            local_To_Global(faceid, Coord[1], glPos);
            vec[1] = get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);

            ////
            vec2D.entry[0] = face->x1 - face->x2;
            vec2D.entry[1] = - face->y2;

            Coord[2] = Coord[0] + 0.5 * vec2D;  //get the point on the edge v1v2

            ////Get the vector associated with the point
            get_2D_Barycentric_Facters(faceid, Coord[2].entry[0], Coord[2].entry[1], t_alpha);
            local_To_Global(faceid, Coord[2], glPos);
            vec[2] = get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);
        }
    }

    ////the singularity is on one of the edges
    else if(fabs(alpha[0] - 0) < 1e-6 || fabs(alpha[1] - 0) < 1e-6 || fabs(alpha[2] - 0) < 1e-6)
    {
        Coord[0].entry[0] = sincx;
        Coord[0].entry[1] = sincy;

        ////On the edge v1v2
        if(fabs(alpha[0]-0) < 1e-6)
        {
            ///calculate the other two points along the two edges incident to v0 respectively
            vec2D.entry[0] = face->x2 - sincx;
            vec2D.entry[1] = face->y2 - sincy;

            Coord[1] = Coord[0] + 0.5 * vec2D;

            ////Get the vector associated with the point
            get_2D_Barycentric_Facters(faceid, Coord[1].entry[0], Coord[1].entry[1], t_alpha);
            local_To_Global(faceid, Coord[1], glPos);
            vec[1] = get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);

            ////
            vec2D.entry[0] = - face->x1;
            vec2D.entry[1] = 0;

            Coord[2] = Coord[0] + 0.2 * vec2D;

            ////Get the vector associated with the point
            get_2D_Barycentric_Facters(faceid, Coord[2].entry[0], Coord[2].entry[1], t_alpha);
            local_To_Global(faceid, Coord[2], glPos);
            vec[2] = get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);
        }

        ////On the edge v2v0
        else if(fabs(alpha[1]-0) < 1e-6)
        {
            ///calculate the other two points along the two edges incident to v0 respectively
            vec2D.entry[0] = - sincx;
            vec2D.entry[1] = - sincy;

            Coord[1] = Coord[0] + 0.5 * vec2D;

            ////Get the vector associated with the point
            get_2D_Barycentric_Facters(faceid, Coord[1].entry[0], Coord[1].entry[1], t_alpha);
            local_To_Global(faceid, Coord[1], glPos);
            vec[1] = get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);

            ////
            vec2D.entry[0] = face->x1;
            vec2D.entry[1] = 0;

            Coord[2] = Coord[0] + 0.2 * vec2D;

            ////Get the vector associated with the point
            get_2D_Barycentric_Facters(faceid, Coord[2].entry[0], Coord[2].entry[1], t_alpha);
            local_To_Global(faceid, Coord[2], glPos);
            vec[2] = get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);
        }

        ////On the edge v0v1
        else
        {
            ///calculate the other two points along the two edges incident to v0 respectively
            vec2D.entry[0] = face->x1 - sincx;
            vec2D.entry[1] = - sincy;

            Coord[1] = Coord[0] + 0.5 * vec2D;

            ////Get the vector associated with the point
            get_2D_Barycentric_Facters(faceid, Coord[1].entry[0], Coord[1].entry[1], t_alpha);
            local_To_Global(faceid, Coord[1], glPos);
            vec[1] = get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);

            ////
            vec2D.entry[0] = face->x2;
            vec2D.entry[1] = face->y2;

            Coord[2] = Coord[0] + 0.2 * vec2D;

            ////Get the vector associated with the point
            get_2D_Barycentric_Facters(faceid, Coord[2].entry[0], Coord[2].entry[1], t_alpha);
            local_To_Global(faceid, Coord[2], glPos);
            vec[2] = get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);
        }
    }

    ////The general cases
    else
    {
        //Coord[0].entry[0] = sincx;
        //Coord[0].entry[1] = sincy;
        icVector2 sinpos;
        sinpos.entry[0] = sincx;
        sinpos.entry[1] = sincy;

        vec2D.entry[0] = face->x2 - sincx;
        vec2D.entry[1] = face->y2 - sincy;

        Coord[0] = sinpos + 0.5 * vec2D;

        ////Get the vector associated with the point
        get_2D_Barycentric_Facters(faceid, Coord[0].entry[0], Coord[0].entry[1], t_alpha);
        local_To_Global(faceid, Coord[0], glPos);
        vec[0] = get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);


        vec2D.entry[0] = - sincx;
        vec2D.entry[1] = - sincy;

        Coord[1] = sinpos + 0.5 * vec2D;

        ////Get the vector associated with the point
        get_2D_Barycentric_Facters(faceid, Coord[1].entry[0], Coord[1].entry[1], t_alpha);
        local_To_Global(faceid, Coord[1], glPos);
        vec[1] = get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);

        ////
        vec2D.entry[0] = face->x1 - sincx;
        vec2D.entry[1] = - sincy;

        Coord[2] = sinpos + 0.5 * vec2D;

        ////Get the vector associated with the point
        get_2D_Barycentric_Facters(faceid, Coord[2].entry[0], Coord[2].entry[1], t_alpha);
        local_To_Global(faceid, Coord[2], glPos);
        vec[2] = get_Vector_At_Point(faceid, glPos.entry, t_alpha, 0, 0);
    }
}



