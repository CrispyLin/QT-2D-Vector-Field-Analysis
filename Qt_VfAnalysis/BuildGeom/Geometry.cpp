/*
This file will handle the basic ply file loading and modeling the object

Created and Modified by Guoning Chen
copyright @2007
*/

#include "BuildGeom/Geometry.h"
#include "Graph2D.h"

static PlyFile *in_ply;

#define FILEDEBUG

#define EQFIELD

extern const int DebugOn;

Polyhedron *object = NULL;
unsigned char orientation;  // 0=ccw, 1=cw


#include "Graph2D.h"
std::vector<Contour_Graph> isolines;

Corner **extend_Link(Corner **corner, int ncorners)
{
    if(ncorners == 0)
    {
        corner = new Corner *[1];
        return corner;
    }

    int i;
    Corner **temp = corner;

    corner = new Corner *[ncorners+1];
    if(corner == NULL)
    {
        //MessageBox(NULL, "fail to allocate memory for corner list of a vertex, program will be terminated", "", MB_OK);
        exit(-1);
    }

    for(i = 0; i < ncorners; i++)
        corner[i] = temp[i];

    corner = temp;
    free(temp);
    return corner;
}


/* Implementation of the memeber function for Vertex class */

Vertex::Vertex(double ix=0, double iy=0, double iz=0,
               double vx=0, double vy=0, double vz=0)
{
    x = ix; y = iy;  z = iz;
    g_vec.entry[0] = vx;
    g_vec.entry[1] = vy;
    g_vec.entry[2] = vz;

    raw_vec.entry[0] = vx;
    raw_vec.entry[1] = vy;
    raw_vec.entry[2] = vz;
}


Vertex::Vertex(double ix, double iy, double iz,
               double vx, double vy, double vz, int i)
{
    x = ix; y = iy;  z = iz;
    t_vec.entry[0] = vx;
    t_vec.entry[1] = vy;
    t_vec.entry[2] = vz;

    raw_vec.entry[0] = vx;
    raw_vec.entry[1] = vy;
    raw_vec.entry[2] = vz;
}

Vertex::~Vertex()
{
}

/* Implementation of the memeber function for VertexList class */
//VertexList::VertexList()
//{
//}

//VertexList::~VertexList()
//{
//}


/* Implementation of the memeber function for Triangle class */
bool Triangle::fall_in_the_tri(double pos[3])
{
    double alpha[3] = {0.};

    icVector3 VP;

    VP.entry[0] = pos[0] - verts[0]->x;
    VP.entry[1] = pos[1] - verts[0]->y;
    VP.entry[2] = pos[2] - verts[0]->z;

    ///1. transfer to local coordinates
    double a = dot(VP, this->LX);
    double b = dot(VP, this->LY);

    ///2. calculate the barycentric coordinates of the point under the local frame of the triangle
    object->get_2D_Barycentric_Facters(this->index, a, b, alpha);

    ///3. judge whether the point falls into the triangle
    if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1
        && (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1
        && (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
        return true;
    else
        return false;
}



/* Implementation of the memeber function for TriangleList class */


void TriangleList::capture_Singularities()
{
    int i, j;
    Triangle *face;
    int *verts;
    icVector2 VertsVector;
    double *Polar_ang, *angle, ang_sum = 0;
    icVector2 tempV1, tempV2;  // to store the two vectors for the vertex that has angle deficit
    int direct_id = 0;         // index to read the directional vectors
    double vx[3], vy[3];       //

    double dotresult = 0;     ////testing variables

    ////Test total Poincare Index
    int poincare_index = 0;


    nsingulartris = 0;
    //cur_singularity_index = 0;

////Initialize the unknown singularities link
/////Building the triangle ID link
//if(singulartris == NULL)
//{
//	curMaxNumSingularTris = (int)(curMaxNumTris/2);
//	singulartris = new Triangle *[curMaxNumSingularTris];

//	if(singulartris == NULL)
//	{
//		MessageBox(NULL, "failed to allocate memory to detect singularities\n", "", MB_OK);
//		exit(-1);
//	}
//}


//////Test codes here, write to file//////
#ifdef FILEDEBUG
    FILE *fp = fopen("tt.txt","w");
#endif

    int positive = 0, negative = 0;

    ////These flags are used for pair cancellation
    ////Initialize all the flags!!!
    for(i = 0; i < ntris; i++)
    {
        tris[i]->singularityID = -1;
    }

    /////here we should estimate each triangle under its local frame
    for (i=0; i<ntris; i++) {
        face = tris[i];
        //verts = face->verts;

        Polar_ang = new double[face->num_dir_vecs];
        angle = new double[face->num_dir_vecs];

        ang_sum = 0;
        direct_id = 0;

        //int vertexsing = 0;
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

                if((length(tempV1)<= 1e-18) || (length(tempV2) <= 1e-18))////set some threshold instead of 0
                {
                    //if(nsingulartris >= curMaxNumSingularTris)
                    //{
                    //	TriangleList **temp = singulartris;
                    //	singulartris = new Triangle *[curMaxNumSingularTris+50];

                    //	if(singulartris == NULL)
                    //	{
                    //		MessageBox(NULL, "failed to allocate memory to detect singularities\n", "", MB_OK);
                    //		exit(-1);
                    //	}

                    //	for(int k = 0; k < curMaxNumSingularTris; k++)
                    //	{
                    //		singulartris[k] = temp[k];
                    //	}
                    //	curMaxNumSingularTris += 50;
                    //}

                    //singulartris[FindTriangleIndex] = i;

                    face->singularityID = nsingulartris;         ////Mark current triangle as one containing singularity
                    nsingulartris ++;

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

                if(fabs(vx[j]) <= 1e-18 && fabs(vy[j]) <= 1e-18)  ////set some threshold instead of 0
                {
                    //if(FindTriangleIndex >= MaxNumSingularities)
                    //{
                    //	MaxNumSingularities += 50;
                    //	MarkTriangleID = (int*) realloc(MarkTriangleID, sizeof(int) * MaxNumSingularities);
                    //	singularities = (Singularities*) realloc(singularities, sizeof(Singularities) * MaxNumSingularities);
                    //}

                    //MarkTriangleID[FindTriangleIndex] = i;
                    face->singularityID = nsingulartris;         ////Mark current triangle as one containing singularity
                    nsingulartris ++;

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


    LL:			if(fabs(ang_sum) >= (2 * M_PI - 0.01))
        {
            face->singularityID = nsingulartris;         ////Mark current triangle as one containing singularity
            nsingulartris ++;

            ////Add the index (at present, we just consider first order critical points
            if(ang_sum > 0){
                poincare_index ++;
                positive ++;
            }
            else{
                poincare_index --;
                negative ++;
            }
            //#ifdef FILEDEBUG
            if(DebugOn == 1){
                fprintf(fp, "Triangle : %d,   Index: %f\n", i, ang_sum/(2*M_PI));
                fprintf(fp, "       The first vectors on the vertices: ");
                for( int k = 0; k < face->num_dir_vecs; k ++)
                    fprintf(fp, " %f, ", length(face->dir_vec[k]));

                fprintf(fp, "  \n\n ");
            }
            //#endif
        }

    }

    //#ifdef FILEDEBUG
    if(DebugOn == 1){
        fprintf(fp, "Positive: %d,   Negative: %d,  Total: %d\n", positive, negative, poincare_index);
        fclose(fp);
    }
    //#endif
}



/* Implementation of the memeber function for Edge class */


/* Build edge list */

/* Implementation of the memeber function for EdgeList class */


/* Implementation of the memeber function for Corner class */


/* Implementation of the memeber function for CornerList class */


/* Implementation of the memeber function for Polyhedron class */

/******************************************************************************
Read in a polyhedron from a file.
******************************************************************************/

Polyhedron::Polyhedron(FILE *file, int withvec)
{
    int i,j;
    int elem_count;
    char *elem_name;

    /*** Read in the original PLY object ***/
    in_ply = myplyloader.read_ply (file);

    for (i = 0; i < in_ply->num_elem_types; i++)
    {
        /* prepare to read the i'th list of elements */
        elem_name = myplyloader.setup_element_read_ply (in_ply, i, &elem_count);

        if (myplyloader.equal_strings ((char*)"vertex", elem_name)) {

            /* create a vertex list to hold all the vertices */
            vlist.nverts = vlist.curMaxNumVerts = elem_count;
            vlist.verts = new Vertex *[vlist.nverts];

            /* set up for getting vertex elements */

            myplyloader.setup_property_ply (in_ply, &vert_props[0]);
            myplyloader.setup_property_ply (in_ply, &vert_props[1]);
            myplyloader.setup_property_ply (in_ply, &vert_props[2]);
            myplyloader.setup_property_ply (in_ply, &vert_props[3]);

            if(withvec == 1)
            {
                myplyloader.setup_property_ply (in_ply, &vert_props[4]);
                myplyloader.setup_property_ply (in_ply, &vert_props[5]);
                myplyloader.setup_property_ply (in_ply, &vert_props[6]);
            }

            vert_other = myplyloader.get_other_properties_ply (in_ply,
                                                              offsetof(Vertex_io,other_props));

            /* grab all the vertex elements */
            for (j = 0; j < vlist.nverts; j++) {
                Vertex_io vert;
                vert.angle_deficit = 0;  /*for the file without angle deficit information*/
                myplyloader.get_element_ply (in_ply, (void *) &vert);

                /* copy info from the "vert" structure */
                if(withvec == 0)
                    vlist.verts[j] = new Vertex (vert.x, vert.y, vert.z);
                else
                    vlist.verts[j] = new Vertex (vert.x, vert.y, vert.z, vert.vx, vert.vy, vert.vz);
                //vlist.verts[j] = new Vertex (vert.x, vert.y, vert.z, vert.vx, vert.vy, vert.vz, 0);

                /* set angle deficit */
                if(vert.angle_deficit == 1)
                    vlist.verts[j]->angle_deficit = true;
                else
                    vlist.verts[j]->angle_deficit = false;

                /* set the index */
                vlist.verts[j]->index = j;

                /* initial other member variables */
                vlist.verts[i]->edges = NULL;
                vlist.verts[i]->corners = NULL;
                vlist.verts[j]->ncorners = vlist.verts[j]->nedges = 0;
                vlist.verts[j]->ntris = 0;

                /* set other properties of the vertex */
                vlist.verts[j]->other_props = vert.other_props;

                vlist.verts[j]->distance = 1e49;

            }
        }
        else if (myplyloader.equal_strings ((char*)"face", elem_name))
        {
            /* create a list to hold all the face elements */
            tlist.ntris =tlist.curMaxNumTris = elem_count;
            tlist.tris = new Triangle *[tlist.ntris];

            /* set up for getting face elements */
            myplyloader.setup_property_ply (in_ply, &face_props[0]);
            face_other = myplyloader.get_other_properties_ply (in_ply, offsetof(Face_io,other_props));

            /* grab all the face elements */
            for (j = 0; j < elem_count; j++)
            {
                Face_io face;
                myplyloader.get_element_ply (in_ply, (void *) &face);

                if (face.nverts != 3) {
                    fprintf (stderr, "Face has %d vertices (should be three).\n",
                            face.nverts);
                    exit (-1);
                }

                /* copy info from the "face" structure */
                tlist.tris[j] = new Triangle();
                tlist.tris[j]->nverts = 3;
                tlist.tris[j]->verts[0] = (Vertex *) face.verts[0];
                tlist.tris[j]->verts[1] = (Vertex *) face.verts[1];
                tlist.tris[j]->verts[2] = (Vertex *) face.verts[2];
                tlist.tris[j]->other_props = face.other_props;

                /*important resetting*/
                tlist.tris[j]->samplepts = NULL;  /*intialize the sample list 05/01/07*/
                tlist.tris[j]->num_samplepts = 0;

                //tlist.tris[j]->inDesignCell = false;

                tlist.tris[j]->counter = 0;

                //cprintf("Current faces: %d\n", j);
                if (j==2269)
                {
                    int stop = 1;
                }
            }
        }
        else
            myplyloader.get_other_element_ply (in_ply);

    }

    /* close the file */
    myplyloader.close_ply (in_ply);

    /* fix up vertex pointers in triangles */
    for (i = 0; i < tlist.ntris; i++) {
        tlist.tris[i]->verts[0] = vlist.verts[(int) tlist.tris[i]->verts[0]];
        tlist.tris[i]->verts[1] = vlist.verts[(int) tlist.tris[i]->verts[1]];
        tlist.tris[i]->verts[2] = vlist.verts[(int) tlist.tris[i]->verts[2]];
    }

    /* get rid of triangles that use the same vertex more than once */

    for (i = tlist.ntris-1; i >= 0; i--) {

        Triangle *tri = tlist.tris[i];
        Vertex *v0 = tri->verts[0];
        Vertex *v1 = tri->verts[1];
        Vertex *v2 = tri->verts[2];

        if (v0 == v1 || v1 == v2 || v2 == v0) {
            free (tlist.tris[i]);
            tlist.ntris--;
            tlist.tris[i] = tlist.tris[tlist.ntris];
        }
    }

    if(withvec == 1)
    {
        with_vec = true;
        ///* normalize the vector values */
        //normalize_Field();
        ///* get the tangential vectors for all vertices if vector values have been read from file*/
        //project_to_TangentPlane();
    }

    else
        with_vec = false;

    /* initialize the vector field variables here */
    //slist.curMaxNumSingularities = 300;
}


extern void  HsvRgb( float hsv[3], float rgb[3] );

void Polyhedron::assign_color_VFmagnitude()
{
    int i;
    min_vfmagnitude = 1e20;
    max_vfmagnitude = -1;
    double vecLenSum = 0;

    for(i=0; i<vlist.nverts; i++)
    {
        double vecLen = length(vlist.verts[i]->g_vec);

        if(vecLen<min_vfmagnitude) min_vfmagnitude = vecLen;
        if(vecLen>max_vfmagnitude) max_vfmagnitude = vecLen;

        vecLenSum += vecLen;
    }

    mean_vfmagnitude = vecLenSum/vlist.nverts;

    //compute the standard deviation
    double sd = 0;
    for(i=0; i<vlist.nverts; i++)
    {
        double vecLen = length(vlist.verts[i]->g_vec);
        sd += (vecLen - mean_vfmagnitude)*(vecLen - mean_vfmagnitude);
    }

    sd = sd/vlist.nverts;
    sd = sqrt(sd);

    min_VecLenVisualized = mean_vfmagnitude - 0.5*sd;
    max_VecLenVisualized = mean_vfmagnitude + 0.45*sd;

    // assign colors according to the vector field magnitude
    float hsv[3]={0, 1, 1};
    float rgb[3];
    for(i=0; i<vlist.nverts; i++)
    {
        double vecLen = length(vlist.verts[i]->g_vec);

        if(vecLen < min_VecLenVisualized) hsv[0] = 256.;
        else if(vecLen > max_VecLenVisualized) hsv[0] = 0.;
        else{
            hsv[0] = 256.-256.*(vecLen-min_VecLenVisualized)/(max_VecLenVisualized-min_VecLenVisualized);
        }
        HsvRgb(hsv, rgb);

        vlist.verts[i]->mag_rgb[0] = rgb[0];
        vlist.verts[i]->mag_rgb[1] = rgb[1];
        vlist.verts[i]->mag_rgb[2] = rgb[2];
    }
}


void
Polyhedron::assign_color_VFmagnitude2()
{
    int i;
    min_vfmagnitude = 1e20;
    max_vfmagnitude = -1;
    double vecLenSum = 0;

    /*
        For the earthquake data, we need to project to the unit look direction (0.3310 0.0710 0.9410)
        06/13/2011
    */
    icVector3 look_dir(0.3310, 0.0710, 0.9410);


    for(i=0; i<vlist.nverts; i++)
    {
        double vecLen = length(vlist.verts[i]->g_vec);

#ifdef EQFIELD
        vecLen = dot (vlist.verts[i]->g_vec, look_dir);
#endif

        vlist.verts[i]->vf_mag = vecLen;

        if(vecLen<min_vfmagnitude) min_vfmagnitude = vecLen;
        if(vecLen>max_vfmagnitude) max_vfmagnitude = vecLen;

        vecLenSum += vecLen;
    }


    //////////////////////////////////////////////////////

    //mean_vfmagnitude = vecLenSum/vlist.nverts;  // for cooling jacket
    mean_vfmagnitude = (min_vfmagnitude+max_vfmagnitude)/2;

    //compute the standard deviation
    double sd = 0;
    for(i=0; i<vlist.nverts; i++)
    {
        double vecLen = length(vlist.verts[i]->g_vec);
        sd += (vecLen - mean_vfmagnitude)*(vecLen - mean_vfmagnitude);
    }

    sd = sd/vlist.nverts;
    sd = sqrt(sd);


    min_VecLenVisualized = mean_vfmagnitude - 0.5*sd;  // for projected magnitude
    max_VecLenVisualized = mean_vfmagnitude + 0.85*sd;

    // assign colors according to the vector field magnitude
    float hsv[3]={0, 1, 1};
    float rgb[3];
    for(i=0; i<vlist.nverts; i++)
    {
        //double vecLen = length(vlist.verts[i]->g_vec);
        double vecLen = vlist.verts[i]->vf_mag;

#ifdef	EQFIELD
        //hsv[0] = 240.-240.*(vecLen-min_vfmagnitude)*1/(max_vfmagnitude-min_vfmagnitude);
        if(vecLen < min_VecLenVisualized) hsv[0] = 240.;
        else if(vecLen > max_VecLenVisualized) hsv[0] = 0.;
        else
            hsv[0] = 240.-240.*(vecLen-min_VecLenVisualized)/(max_VecLenVisualized-min_VecLenVisualized);
#else
        if(vecLen < min_VecLenVisualized) hsv[0] = 240.;
        else if(vecLen > max_VecLenVisualized) hsv[0] = 0.;
        else{
            hsv[0] = 240.-240.*(vecLen-min_VecLenVisualized)/(max_VecLenVisualized-min_VecLenVisualized);
            //hsv[0] = 240.-240.*(vecLen-min_vfmagnitude)/(max_vfmagnitude-min_vfmagnitude);
        }
#endif

        HsvRgb(hsv, rgb);

        vlist.verts[i]->mag_rgb[0] = rgb[0];
        vlist.verts[i]->mag_rgb[1] = rgb[1];
        vlist.verts[i]->mag_rgb[2] = rgb[2];
    }

    /***************************************************************
        Let us compute some iso-contour here 06/13/2011
    */
    int num_contours = 20;

    isolines.clear();
    isolines.resize(num_contours);

#ifdef EQFIELD
    double interval = (max_vfmagnitude-min_vfmagnitude)/(num_contours-1);
    for (i=0; i<num_contours; i++)
    {
        double cur_val = /*0.001+*/min_vfmagnitude+(i-.5)*interval;

        isolines[i].extract_iso_contour_from_distance_field(cur_val);
    }
#endif
}


extern double dmax;  //parameter related to the IBFV method


void Polyhedron::normalize_Field()
{
    int i;
    double r;
    Vertex *cur_v;

    for(i = 0; i < vlist.nverts; i++)
    {
        cur_v = vlist.verts[i];
        r = length(cur_v->g_vec);
        r *= r;

        if (r < DistanceThreshold)
        {
            r = DistanceThreshold;
            cur_v->g_vec *= dmax/r;
        }

        r = length(cur_v->g_vec);
        r *= r;

        if (r > dmax*dmax) {
            r  = sqrt(r);
            cur_v->g_vec *= dmax/r;
        }
    }
}


//void Polyhedron::normalize_TangentField()
//{
//	int i;
//    double r;
//	Vertex *cur_v;
//
//	for(i = 0; i < vlist.nverts; i++)
//	{
//		cur_v = vlist.verts[i];
//		r = length(cur_v->t_vec);
//		r *= r;
//
//		if (r < DistanceThreshold)
//		{
//			r = DistanceThreshold;
//			cur_v->t_vec *= dmax/r;
//		}
//
//		r = length(cur_v->t_vec);
//		r *= r;
//
//		if (r > dmax*dmax) {
//			r  = sqrt(r);
//			cur_v->t_vec *= dmax/r;
//		}
//	}
//}

/******************************************************************************
Find out if there is another face that shares an edge with a given face.

Entry:
        f1    - face that we're looking to share with
        v1,v2 - two vertices of f1 that define edge

                Exit:
                       return the matching face, or NULL if there is no such face
                ******************************************************************************/

            Triangle *Polyhedron::find_common_Edge(Triangle *f1, Vertex *v1, Vertex *v2)
{
    int i,j;
    Triangle *f2;
    Triangle *adjacent = NULL;

    /* look through all faces of the first vertex */

    for (i = 0; i < v1->ntris; i++) {
        f2 = v1->tris[i];
        if (f2 == f1)
            continue;
        /* examine the vertices of the face for a match with the second vertex */
        for (j = 0; j < f2->nverts; j++)
        {
            /* look for a match */
            if (f2->verts[j] == v2) {
                /* if we've got a match, return this face */
                return (f2);
            }
        }
    }

    return (adjacent);
}


/******************************************************************************
Create various face and vertex pointers.
******************************************************************************/

void Polyhedron::create_Pointers()
{
    int i;

    /* index the vertices and triangles */

    for (i = 0; i < vlist.nverts; i++)
        vlist.verts[i]->index = i;

    for (i = 0; i < tlist.ntris; i++)
        tlist.tris[i]->index = i;

    /* create pointers from vertices to triangles */
    vertex_to_tri_ptrs();

    /* make edges */
    build_Edges();

    //print out the edges that are adjacent to vertex 3729
    //for (i=0; i<elist.nedges; i++)
    //{
    //	if (elist.edges[i]->verts[0]->index == 3729
    //		|| elist.edges[i]->verts[1]->index == 3729)
    //	{
    //		Edge *e = elist.edges[i];
    //		_cprintf("Edge (%d, %d): tri1 = %d, tri2 = %d.\n",
    //			e->verts[0]->index, e->verts[1]->index,
    //			e->tris[0]->index, e->tris[1]!=NULL? e->tris[1]->index: -1);
    //	}
    //}


    /* order the pointers from vertices to faces */
    for (i = 0; i < vlist.nverts; i++){
        order_Vertex_to_Tri_Ptrs(vlist.verts[i]);

    }

    /* index the edges */

    for (i = 0; i < elist.nedges; i++){
        elist.edges[i]->index = i;
    }

}


/******************************************************************************
Create pointers from vertices to faces.
******************************************************************************/

void Polyhedron::vertex_to_tri_ptrs()
{
    int i,j;
    Triangle *f;
    Vertex *v;

    /* zero the count of number of pointers to faces */

    for (i = 0; i < vlist.nverts; i++)
        vlist.verts[i]->max_tris = 0;

    /* first just count all the face pointers needed for each vertex */

    for (i = 0; i < tlist.ntris; i++) {
        f = tlist.tris[i];
        for (j = 0; j < f->nverts; j++)
            f->verts[j]->max_tris++;
    }

    /* allocate memory for face pointers of vertices */

    for (i = 0; i < vlist.nverts; i++) {
        vlist.verts[i]->tris = (Triangle **)
            malloc (sizeof (Triangle *) * vlist.verts[i]->max_tris);
        vlist.verts[i]->ntris = 0;
    }

    /* now actually create the face pointers */

    for (i = 0; i < tlist.ntris; i++) {
        f = tlist.tris[i];
        for (j = 0; j < f->nverts; j++) {
            v = f->verts[j];
            v->tris[v->ntris] = f;
            v->ntris++;
        }
    }
}


/******************************************************************************
Order the pointers to faces that are around a given vertex.

Entry:
        v - vertex whose face list is to be ordered
                  ******************************************************************************/

              void Polyhedron::order_Vertex_to_Tri_Ptrs(Vertex *v)
{
    int i,j;
    Triangle *f;
    Triangle *fnext;
    int nf;
    int vindex;
    int boundary;
    int count;

    nf = v->ntris;
    f = v->tris[0];

    /* go backwards (clockwise) around faces that surround a vertex */
    /* to find out if we reach a boundary */

    boundary = 0;

    for (i = 1; i <= nf; i++) {

        /* find reference to v in f */
        vindex = -1;
        for (j = 0; j < f->nverts; j++)
            if (f->verts[j] == v) {
                vindex = j;
                break;
            }

        /* error check */
        if (vindex == -1) {
            fprintf (stderr, "can't find vertex #1\n");
            exit (-1);
        }

        /* corresponding face is the previous one around v */
        fnext = other_Triangle (f->edges[vindex], f);

        /* see if we've reached a boundary, and if so then place the */
        /* current face in the first position of the vertice's face list */

        if (fnext == NULL) {
            /* find reference to f in v */
            for (j = 0; j < v->ntris; j++)
                if (v->tris[j] == f) {
                    v->tris[j] = v->tris[0];
                    v->tris[0] = f;
                    break;
                }
            boundary = 1;
            break;
        }

        f = fnext;
    }

    /* now walk around the faces in the forward direction and place */
    /* them in order */

    f = v->tris[0];
    count = 0;

    for (i = 1; i < nf; i++) {

        /* find reference to vertex in f */
        vindex = -1;
        for (j = 0; j < f->nverts; j++)
            if (f->verts[(j+1) % f->nverts] == v) {
                vindex = j;
                break;
            }

        /* error check */
        if (vindex == -1) {
            fprintf (stderr, "can't find vertex #2\n");
            exit (-1);
        }

        /* corresponding face is next one around v */
        fnext = other_Triangle (f->edges[vindex], f);

        /* break out of loop if we've reached a boundary */
        count = i;
        if (fnext == NULL) {
            break;
        }

        /* swap the next face into its proper place in the face list */
        for (j = 0; j < v->ntris; j++)
            if (v->tris[j] == fnext) {
                v->tris[j] = v->tris[i];
                v->tris[i] = fnext;
                break;
            }

        f = fnext;
    }
}


/******************************************************************************
Find the other triangle that is incident on an edge, or NULL if there is
no other.
Input: Edge *edge -- the pointer of an edge
       Triangle *tri -- pointer to current triangle
******************************************************************************/

Triangle *Polyhedron::other_Triangle(Edge *edge, Triangle *tri)
{
    /* search for any other triangle */

    for (int i = 0; i < edge->ntris; i++)
        if (edge->tris[i] != tri)
            return (edge->tris[i]);

    /* there is no such other triangle if we get here */
    return (NULL);
}



void Polyhedron::calc_Bounding_Sphere()
{
    int i;
    icVector3 min, max;

    for (i=0; i<vlist.nverts; i++) {
        if (i==0)  {
            min.set(vlist.verts[i]->x, vlist.verts[i]->y, vlist.verts[i]->z);
            max.set(vlist.verts[i]->x, vlist.verts[i]->y, vlist.verts[i]->z);
        }
        else {
            if (vlist.verts[i]->x < min.entry[0])
                min.entry[0] = vlist.verts[i]->x;
            if (vlist.verts[i]->x > max.entry[0])
                max.entry[0] = vlist.verts[i]->x;
            if (vlist.verts[i]->y < min.entry[1])
                min.entry[1] = vlist.verts[i]->y;
            if (vlist.verts[i]->y > max.entry[1])
                max.entry[1] = vlist.verts[i]->y;
            if (vlist.verts[i]->z < min.entry[2])
                min.entry[2] = vlist.verts[i]->z;
            if (vlist.verts[i]->z > max.entry[2])
                max.entry[2] = vlist.verts[i]->z;
        }
    }
    center = (min + max) * 0.5;
    radius = length(center - min);
    rot_center = center * 1.0;
}

void Polyhedron::calc_Face_Normals_and_Area()
{
    unsigned int i, j, k, l;
    icVector3 v0, v1, v2;
    Triangle *temp_t;
    double length[3];

    area = 0.0;
    for (i=0; i<tlist.ntris; i++){
        for (j=0; j<3; j++)
            length[j] = tlist.tris[i]->edges[j]->length;
        double temp_s = (length[0] + length[1] + length[2])/2.0;
        tlist.tris[i]->area = sqrt(temp_s*(temp_s-length[0])*(temp_s-length[1])*(temp_s-length[2]));

        area += tlist.tris[i]->area;
        temp_t = tlist.tris[i];
        v1.set(vlist.verts[tlist.tris[i]->verts[0]->index]->x,
               vlist.verts[tlist.tris[i]->verts[0]->index]->y,
               vlist.verts[tlist.tris[i]->verts[0]->index]->z);

        v2.set(vlist.verts[tlist.tris[i]->verts[1]->index]->x,
               vlist.verts[tlist.tris[i]->verts[1]->index]->y,
               vlist.verts[tlist.tris[i]->verts[1]->index]->z);

        v0.set(vlist.verts[tlist.tris[i]->verts[2]->index]->x,
               vlist.verts[tlist.tris[i]->verts[2]->index]->y,
               vlist.verts[tlist.tris[i]->verts[2]->index]->z);

        tlist.tris[i]->normal = cross(v0-v1, v2-v1);
        normalize(tlist.tris[i]->normal);
        //initilaize the cur_morse 12/03/2009 added by Qingqing
        tlist.tris[i]->cur_morse=false;
    }

    double signedvolume = 0.0;
    icVector3 test = center;
    for (i=0; i<tlist.ntris; i++){
        icVector3 cent(vlist.verts[tlist.tris[i]->verts[0]->index]->x,
                       vlist.verts[tlist.tris[i]->verts[0]->index]->y,
                       vlist.verts[tlist.tris[i]->verts[0]->index]->z);

        signedvolume += dot(test-cent, tlist.tris[i]->normal)*tlist.tris[i]->area;
    }
    signedvolume /= area;
    if (signedvolume<0)
        orientation = 0;
    else {
        orientation = 1;
        for (i=0; i<tlist.ntris; i++)
            tlist.tris[i]->normal *= -1.0;
    }


    /*
       Added by Guoning: compute the euler step within each triangle
    */
    for (i=0; i<tlist.ntris; i++)
    {
        tlist.tris[i]->euler_step = sqrt(tlist.tris[i]->area)/20.;
    }
}

void Polyhedron::calc_Face_Normals_and_Area_2(void)
{
    unsigned int i;
    int j, k, l;
    icVector3 v0, v1, v2;
    Triangle *face;
    int *verts;
    double temp_x[100], temp_y[100];  // assume on polygon has more than 100 vertices.
    int largest_face;

    object->area = 0.0;
    largest_face = 0;
    for (i=0; i<vlist.nverts; i++)
        vlist.verts[i]->normal.set(0.0);
    for (i=0; i<tlist.ntris; i++){
        tlist.tris[i]->area = 0.0;
        face = tlist.tris[i];
        if (face->nverts >largest_face)
            largest_face = face->nverts;
        //verts = face->verts;
        v1.set(face->verts[0]->x, face->verts[0]->y, face->verts[0]->z);
        v2.set(face->verts[1]->x, face->verts[1]->y, face->verts[1]->z);
        v0.set(face->verts[2]->x, face->verts[2]->y, face->verts[2]->z);
        tlist.tris[i]->normal = cross(v0-v1, v2-v1);
        normalize(tlist.tris[i]->normal);
        v1 -= v0;
        v2 = cross(tlist.tris[i]->normal, v1);
        normalize(v1);
        normalize(v2);
        for (j=0; j<face->nverts; j++){
            v0.set(face->verts[j]->x-face->verts[0]->x,
                   face->verts[j]->y-face->verts[0]->y,
                   face->verts[j]->z-face->verts[0]->z);
            temp_x[j] = dot(v1, v0);
            temp_y[j] = dot(v2, v0);
        }

        for (j=0; j<face->nverts; j++){
            tlist.tris[i]->area += (temp_x[j]*temp_y[(j+1) % face->nverts] - temp_x[(j+1) % face->nverts]*temp_y[j]);
            face->verts[j]->normal = face->verts[j]->normal + tlist.tris[i]->normal;
        }
        tlist.tris[i]->area /= 2;
        tlist.tris[i]->area = fabs(tlist.tris[i]->area);
        object->area += tlist.tris[i]->area;
    }
    for (i=0; i<vlist.nverts; i++)
        normalize(vlist.verts[i]->normal);

    double test_val;
    for (j=-2; j<=2; j++) {
        for (k=-2; k<=2; k++) {
            for (l=-2; l<=2; l++) {
                test_val = 0.0;
                icVector3 test((double)j, (double)k, (double)l);
                test *= radius*100;
                test += center;
                for (i=0; i<vlist.nverts; i++){
                    icVector3 cent(tlist.tris[i]->verts[0]->x,
                                   tlist.tris[i]->verts[0]->y,
                                   tlist.tris[i]->verts[0]->z);

                    test_val += dot(test-cent, tlist.tris[i]->normal)*tlist.tris[i]->area;
                }
                test_val /= object->area;
                if ((j==0) && (k==0) && (l==0)){
                    if (test_val<0) {
                        orientation = 0;
                    }
                    else {
                        orientation = 1;
                    }
                }
            }
        }
    }
    if (orientation == 1){
        for (i=0; i<tlist.ntris; i++)
            tlist.tris[i]->normal *= -1.0;
    }
}



/*
Simply average the normals of the adjacent triangles of each vertex
*/
void Polyhedron::get_vertex_normal()
{
    //int i, j;

    //for (i=0; i<vlist.nverts; i++) {
    //	vlist.verts[i]->normal = icVector3(0.0);
    //	for (j=0; j<vlist.verts[i]->ntris; j++)
    //		vlist.verts[i]->normal += vlist.verts[i]->tris[j]->normal;
    //	normalize(vlist.verts[i]->normal);
    //}
    int i, j;

    for (i=0; i<vlist.nverts; i++)
        vlist.verts[i]->normal.set(0.0);

    for (i=0; i<tlist.ntris; i++) {
        for (j=0; j<3; j++)
            tlist.tris[i]->verts[j]->normal += tlist.tris[i]->normal * tlist.tris[i]->area;
    }
    for (i=0; i<vlist.nverts; i++)
        normalize(vlist.verts[i]->normal);
}

/*
Simply average the normals of the adjacent triangles of each vertex
*/
void Polyhedron::average_Normals()
{
    int i, j;

    for (i=0; i<vlist.nverts; i++) {
        vlist.verts[i]->normal = icVector3(0.0);
        for (j=0; j<vlist.verts[i]->ntris; j++)
            vlist.verts[i]->normal += vlist.verts[i]->tris[j]->normal;
        normalize(vlist.verts[i]->normal);
    }
}


void Polyhedron::project_to_TangentPlane()
{
    int i;
    Vertex *vert;
    icVector3 N, V;



    for(i = 0; i < vlist.nverts; i++)
    {
        vert = vlist.verts[i];

        N = vert->normal;

        V = vert->g_vec - dot(vert->g_vec, N)*N;

        /****************************
          Some observations: it seems that the normalized field could get
          better ring-like structure (more smooth, fewer holes)

        However, the non-normalized vector field with certain scaling
                                   can produce good result with large tau (e.g. t=3) value as well!
                               */
                               normalize(V);    // comment out by Guoning at 2/8/2010
        //V = 1000.*V;        // modified by Guoning at 2/8/2010
        //V = 1000.*V;        // modified by Guoning at 3/9/2010 for the diesel slice

        ////Give the vector to the variable "Vector", which is on the tangent plane of the vertex
        vert->t_vec = V;
    }

    //normalize_TangentField();

    /*  for the planar slice, we may want to get rid of the bubble  */
    for (i = 0; i < vlist.nverts; i++)
    {
        vert = vlist.verts[i];

        V.set(0.);

        if (length(vert->t_vec) < 1.e-8)
        {
            for (int j=0; j<vert->nedges; j++)
            {
                Vertex *adj_v = vert->edges[j]->verts[0];

                if (adj_v == vert)
                    adj_v = vert->edges[j]->verts[1];

                V = V+adj_v->t_vec;
            }
            vert->t_vec = 1./vert->nedges* V;
        }
    }
}


//void Polyhedron::find_Opposite()
//{
//	//build 3 integer arrays
//	int *corners = new int[clist.ncorners];
//
//	int *min_v = new int [clist.ncorners];
//
//	int *max_v = new int[clist.ncorners];
//
//	int i;
//
//	//set up the three arrays
//	for(i = 0; i < clist.ncorners; i++)
//	{
//		corners[i] = i;
//		min_v[i] = min(clist.corners[clist.corners[i]->p]->v,clist.corners[clist.corners[i]->n]->v);
//		max_v[i] = max(clist.corners[clist.corners[i]->p]->v,clist.corners[clist.corners[i]->n]->v);
//	}
//
//	//sort the 3 arrays according to the max_v
//	QuickSortforCorner(corners, min_v, max_v, 0, clist.ncorners-1);
//
//	//sort the 3 arrays according to min_v, this will sort the local parts of the array only
//	//we only need to scan the whole arrays once
//	int same_max = max_v[0];
//	int counter = 0;
//	int left, right;
//
//	for(i = 1; i < clist.ncorners; i++)
//	{
//		if(max_v[i] == same_max)
//		{
//			if(counter == 0)
//				left = i-1;
//
//			counter ++;
//			continue;
//		}
//
//		right = i-1;
//		same_max = max_v[i];
//		counter = 0;
//
//		QuickSortforCorner(corners, max_v, min_v, left, right);
//
//	}
//
//	//get the opposite corners in order
//	for(i = 0; i < clist.ncorners-1; i++)
//	{
//		if(min_v[i] == min_v[i+1] && max_v[i] == max_v[i+1])
//		{
//			clist.corners[corners[i]]->o = corners[i+1];
//			clist.corners[corners[i+1]]->o = corners[i];
//			i++;
//		}
//
//		else
//		{
//			clist.corners[corners[i]]->o = -1;
//		}
//	}
//
//	delete [] corners;
//	delete [] min_v;
//	delete [] max_v;
//}
//
//
/*
Initialize some fundamental property of the mesh
*/
void Polyhedron::preprocess_Vertex(void)
{
    for(int i = 0; i < vlist.nverts; i++)
    {
        //vlist.verts[i]->x = vlist.verts[i]->x *( 0.5/radius) + 0.5;
        //vlist.verts[i]->y = vlist.verts[i]->y *( 0.5/radius) + 0.5;
        //vlist.verts[i]->z = vlist.verts[i]->z *( 0.5/radius) ;

        vlist.verts[i]->x = (vlist.verts[i]->x-center.entry[0])/radius + 0.5;
        vlist.verts[i]->y = (vlist.verts[i]->y-center.entry[1])/radius + 0.5;
        vlist.verts[i]->z = (vlist.verts[i]->z-center.entry[2])/radius ;
    }
}



/********************************************************************
Initialize the local frame for each triangle and vertex for tracing
********************************************************************/

void Polyhedron::init_Local_Frame(void)
{

    Triangle *face;
    icVector3 PP;
    int i;

    /////Build the local frame of the tangent plane on each vertex
    icVector3 T, B, N, TM; //define vectors for local frame calculation
    Vertex *temp_ver, *verts;

    icVector3 temp_vector;

    for( i = 0; i < vlist.nverts; i++)  //decrease the number of vertices for real data set visualization
    {

        verts = vlist.verts[i];

        //Calculate the local frame
        ///First, choose a vector on the tangent plane
        //here we use the projection of one edge on the tangent plane
        if(verts->edges[0]->verts[0]->index != i)
            temp_ver = verts->edges[0]->verts[0];
        else
            temp_ver = verts->edges[0]->verts[1];

        temp_vector.entry[0] = temp_ver->x - verts->x;
        temp_vector.entry[1] = temp_ver->y - verts->y;
        temp_vector.entry[2] = temp_ver->z - verts->z;

        N = verts->normal;

        TM = temp_vector - dot(temp_vector, N) * N;
        normalize(TM);

        verts->T.set(TM);    //store T;

        B = cross(N, TM);   //Get B
        normalize(B);       //Normalize B
        verts->B.set(B);    //store B

    }


    ////Build a local frame for each triangle
    for( i = 0; i < tlist.ntris; i++)
    {
        face = tlist.tris[i];
        face->index = i;

        /////using v1v0 as X axis
        PP.entry[0] = face->LX.entry[0] = face->verts[1]->x - face->verts[0]->x;
        PP.entry[1] = face->LX.entry[1] = face->verts[1]->y - face->verts[0]->y;
        PP.entry[2] = face->LX.entry[2] = face->verts[1]->z - face->verts[0]->z;

        normalize(face->LX);

        ////This is a right hand coordinate system
        face->LY = cross(face->normal,face->LX);
        normalize(face->LY);

        ////At the same time, we need to calculate the coordinates of its three vertices under local frame

        ////calculate the coordinates for verts[1]
        //PP.entry[0] = face->verts[1]->x - face->verts[0]->x;
        //PP.entry[1] = face->verts[1]->y - face->verts[0]->y;
        //PP.entry[2] = face->verts[1]->z - face->verts[0]->z;

        face->x1 = dot(PP, face->LX);

        ////Now we need to calculate the coordinates for verts[2]
        PP.entry[0] = face->verts[2]->x - face->verts[0]->x;
        PP.entry[1] = face->verts[2]->y - face->verts[0]->y;
        PP.entry[2] = face->verts[2]->z - face->verts[0]->z;

        face->x2 = dot(PP, face->LX);
        face->y2 = dot(PP, face->LY);

        ////get the orientation of the vertices of the triangle according to the sign of xy[2][1]
        ////Added at 06/10/06
        if(face->y2 < 0)
            face->orient = true;  //clockwise orientation
        else
            face->orient = false;  //counter clockwise orientation
    }

}


/***********************************************************************
Initialize the directional vectors on the vertices of each triangle
The vectors we get may be more than three because of angle deficit
***********************************************************************/
void Polyhedron::init_Trace(bool with_vec)
{
    int i, j;
    Triangle *face;
    Vertex *verts;
    Vertex *vert_prev, *vert_next;

    icVector3 temp_v;

    double vk[3];

    int vec_id;

    bool angle_deficit_yes;

    /*calculate the directional vectors under the local frame defined by the triangle*/
    //  It seems that we need to first sort the corner table before we do this calculation
    double a, b;
    for( i = 0;  i < vlist.nverts; i++)
    {
        a = dot(vlist.verts[i]->T, vlist.verts[i]->t_vec);
        b = dot(vlist.verts[i]->B, vlist.verts[i]->t_vec);

        vlist.verts[i]->t_angle = atan2(b, a);

        if(vlist.verts[i]->t_angle < 0)
            vlist.verts[i]->t_angle += 2 * M_PI;
    }


    for( i = 0; i < tlist.ntris; i++)
    {
        face = tlist.tris[i];
        vec_id = 0;

        //Extend the memory if the triangle contains vertex having angle_deficit
        angle_deficit_yes = false;
        for(j = 0; j < face->nverts; j++)
        {
            if(face->verts[j]->angle_deficit)
            {
                angle_deficit_yes = true;
                break;
            }
        }

        if(angle_deficit_yes)
        {
            face->dir_vec = new icVector2[4];
        }
        else
            face->dir_vec = new icVector2[3];

        for( j = 0; j < face->nverts; j++)
        {
            verts = face->verts[j];

            if(verts->angle_deficit) //angle deficit occurs on this vertex
            {
                ///First, we need to calculate the vector using edge verts[0] and verts[2]
                ///Get the previous corner
                ///modified 3/10/05
                vert_prev = vlist.verts[clist.corners[clist.corners[3*i + j]->p]->v];

                ////set the temparary vector
                //temp_v.entry[0] = vert_prev->x - verts->x;
                //temp_v.entry[1] = vert_prev->y - verts->y;
                //temp_v.entry[2] = vert_prev->z - verts->z;
                //
                //vk[0] = verts->x + temp_v.entry[0] * 0.1;
                //vk[1] = verts->y + temp_v.entry[1] * 0.1;
                //vk[2] = verts->z + temp_v.entry[2] * 0.1;

                vk[0] = vert_prev->x;
                vk[1] = vert_prev->y;
                vk[2] = vert_prev->z;

                face->dir_vec[vec_id] = get_Miu(verts->index, face->index, vk);

                vec_id ++;

                ///second,we need to calculate the vector used edge verts[0] and verts[1]
                ///get the next corner
                vert_next = vlist.verts[clist.corners[clist.corners[3*i + j]->n]->v];

                ////set the temparary vector
                //temp_v.entry[0] = vert_next->x - verts->x;
                //temp_v.entry[1] = vert_next->y - verts->y;
                //temp_v.entry[2] = vert_next->z - verts->z;
                //
                //vk[0] = verts->x + temp_v.entry[0] * 0.1;
                //vk[1] = verts->y + temp_v.entry[1] * 0.1;
                //vk[2] = verts->z + temp_v.entry[2] * 0.1;
                vk[0] = vert_next->x;
                vk[1] = vert_next->y;
                vk[2] = vert_next->z;

                face->dir_vec[vec_id] = get_Miu(verts->index, face->index, vk);
                vec_id ++;
            }

            else{
                ////Use the normal way to calculate the transported vector
                vk[0] = verts->x;
                vk[1] = verts->y;
                vk[2] = verts->z;

                //if (verts->index == 3729)
                //{
                //	int stop = 0;
                //}

                face->dir_vec[vec_id] = get_Miu(verts->index, face->index, vk);

                vec_id ++;
            }
        }
        face->num_dir_vecs = vec_id;
    }
}


/***********************************************************************
Initialize the directional vectors on the vertices of each triangle
The vectors we get may be more than three because of angle deficit
***********************************************************************/
void Polyhedron::recal_vec_in_triangle()
{
    int i, j;
    Triangle *face;
    Vertex *verts;
    Vertex *vert_prev, *vert_next;

    icVector3 temp_v;

    double vk[3];

    int vec_id;

    bool angle_deficit_yes;

    double a, b;
    for( i = 0;  i < vlist.nverts; i++)
    {
        a = dot(vlist.verts[i]->T, vlist.verts[i]->t_vec);
        b = dot(vlist.verts[i]->B, vlist.verts[i]->t_vec);

        vlist.verts[i]->t_angle = atan2(b, a);

        if(vlist.verts[i]->t_angle < 0)
            vlist.verts[i]->t_angle += 2 * M_PI;
    }

    /*calculate the directional vectors under the local frame defined by the triangle*/

    for( i = 0; i < tlist.ntris; i++)
    {
        face = tlist.tris[i];
        vec_id = 0;

        for( j = 0; j < face->nverts; j++)
        {
            verts = face->verts[j];

            if(verts->angle_deficit) //angle deficit occurs on this vertex
            {
                ///First, we need to calculate the vector using edge verts[0] and verts[2]
                ///Get the previous corner
                ///modified 3/10/05
                vert_prev = vlist.verts[clist.corners[clist.corners[3*i + j]->p]->v];

                ////set the temparary vector
                temp_v.entry[0] = vert_prev->x - verts->x;
                temp_v.entry[1] = vert_prev->y - verts->y;
                temp_v.entry[2] = vert_prev->z - verts->z;

                vk[0] = verts->x + temp_v.entry[0] * 0.1;
                vk[1] = verts->y + temp_v.entry[1] * 0.1;
                vk[2] = verts->z + temp_v.entry[2] * 0.1;

                face->dir_vec[vec_id] = get_Miu(verts->index, face->index, vk);

                vec_id ++;

                ///second,we need to calculate the vector used edge verts[0] and verts[1]
                ///get the next corner
                vert_next = vlist.verts[clist.corners[clist.corners[3*i + j]->n]->v];

                ////set the temparary vector
                temp_v.entry[0] = vert_next->x - verts->x;
                temp_v.entry[1] = vert_next->y - verts->y;
                temp_v.entry[2] = vert_next->z - verts->z;

                vk[0] = verts->x + temp_v.entry[0] * 0.1;
                vk[1] = verts->y + temp_v.entry[1] * 0.1;
                vk[2] = verts->z + temp_v.entry[2] * 0.1;

                face->dir_vec[vec_id] = get_Miu(verts->index, face->index, vk);
                vec_id ++;
            }

            else{
                ////Use the normal way to calculate the transported vector
                vk[0] = verts->x;
                vk[1] = verts->y;
                vk[2] = verts->z;

                face->dir_vec[vec_id] = get_Miu(verts->index, face->index, vk);

                vec_id ++;
            }
        }
        face->num_dir_vecs = vec_id;
    }
}



/***********************************************************
The core routine to implement paralle transport of vector
on vertex to any point inside a triangle
Note: we assume the origin of the local frame of the triangle
      is on the first vertex(face->verts[0])!!!!!!!
Input: p vertex
       face  current face
       vk   the point inside current face
Output: transported tensor vector
***********************************************************/
icVector2 Polyhedron::get_Miu(int vertid, int faceid, double vk[3])
{
    ////The following we temp to calculate the transported tensor vector
    ////from the specific vertex on specific triangle
    icVector2 result;
    icVector3 PK;
    double x, y;
    double Phi = 0;

    double k[3] = {0};
    Vertex *p = vlist.verts[vertid];
    Triangle *face = tlist.tris[faceid];

    x = y = 0;

    //////Note that we need to judge the point vk[3] is one of the vertices or not?


    /////Here we suppose the orientation of the triangle is counter clockwise
    ////If it is verts[0], choose edge10 as PK
    ////else if it is verts[1], choose edge21 as PK
    ////else if it is verts[2], choose edge02 as PK

    ////test shows that the orientation of the triangle maybe clockwise

    if( fabs(p->x - vk[0]) < 1e-12
        &&  fabs(p->y - vk[1]) < 1e-12
        &&  fabs(p->z - vk[2]) < 1e-12)
    {
        if( face->verts[0] == p){
            PK.entry[0] = face->verts[2]->x - p->x;
            PK.entry[1] = face->verts[2]->y - p->y;
            PK.entry[2] = face->verts[2]->z - p->z;

        }

        else if( face->verts[1] == p){
            PK.entry[0] = face->verts[0]->x - p->x;
            PK.entry[1] = face->verts[0]->y - p->y;
            PK.entry[2] = face->verts[0]->z - p->z;
        }
        else{
            PK.entry[0] = face->verts[1]->x - p->x;
            PK.entry[1] = face->verts[1]->y - p->y;
            PK.entry[2] = face->verts[1]->z - p->z;
        }

    }

    else{
        ////If it is not the vertex
        PK.entry[0] = vk[0] - p->x;
        PK.entry[1] = vk[1] - p->y;
        PK.entry[2] = vk[2] - p->z;
    }

    /////Normalize the vector PK 03/09/05
    normalize(PK);

    ////Get the two components of vector PK based on local frame of the triangle
    x = dot(PK, face->LX);
    y = dot(PK, face->LY);


    /////Get the angles for rotation

    Phi = get_Phi(p->index, faceid, vk);

    /////Rotate the vector according to local frame by angle Phi
    result.entry[0] = x * cos(Phi) - y * sin(Phi);
    result.entry[1] = x * sin(Phi) + y * cos(Phi);
    normalize(result);

    result = length(p->t_vec) * result;

    return result;
}

/*********************************************
Get the angle phi between vector on vertex p
and vector PK
*********************************************/
double Polyhedron::get_Phi(int vertid, int faceid, double vk[3])
{
    double angle, ang_pk_edge, cur_ang, PolarPK;
    int orient = 0;
    icVector3 PK;
    double dot1 = 0;

    Vertex *p = vlist.verts[vertid];
    Triangle *face = tlist.tris[faceid];

    double k[3] = {0.} ;

    int iface = face->index;
    int i;

    cur_ang = angle = ang_pk_edge = PolarPK = 0;


    if( fabs(p->x - vk[0]) < 1e-12
        && fabs (p->y - vk[1]) < 1e-12
        && fabs (p->z - vk[2]) < 1e-12 )    ////may have numerical issue
    {
        if( face->verts[0] == p){
            PK.entry[0] = face->verts[2]->x - p->x;
            PK.entry[1] = face->verts[2]->y - p->y;
            PK.entry[2] = face->verts[2]->z - p->z;
        }
        else if( face->verts[1] == p){
            PK.entry[0] = face->verts[0]->x - p->x;
            PK.entry[1] = face->verts[0]->y - p->y;
            PK.entry[2] = face->verts[0]->z - p->z;
        }
        else{
            PK.entry[0] = face->verts[1]->x - p->x;
            PK.entry[1] = face->verts[1]->y - p->y;
            PK.entry[2] = face->verts[1]->z - p->z;
        }

    }

    else{
        PK.entry[0] = vk[0] - p->x;
        PK.entry[1] = vk[1] - p->y;
        PK.entry[2] = vk[2] - p->z;
    }


    ////finding the corresponding corner for the triangle in the corner list
    for(i = 0; i < p->ncorners; i++)
    {
        if(p->corners[i]->t == iface)
        {
            break;
        }
    }

    Corner *c = p->corners[i];

    cur_ang = c->BeginAng;


    //orient = c->orient;
    if(tlist.tris[c->t]->y2 < 0)
        orient = 1;
    else
        orient = 0;

    ////calculate the angle between PK vector and the edge[0] of corner[i];
    icVector3 edge;
    ////get the edge[0] vector
    if(c->edge[0]->verts[0] == p)
    {
        edge.entry[0] = c->edge[0]->verts[1]->x - p->x;
        edge.entry[1] = c->edge[0]->verts[1]->y - p->y;
        edge.entry[2] = c->edge[0]->verts[1]->z - p->z;
    }
    else{
        edge.entry[0] = c->edge[0]->verts[0]->x - p->x;
        edge.entry[1] = c->edge[0]->verts[0]->y - p->y;
        edge.entry[2] = c->edge[0]->verts[0]->z - p->z;
    }

    normalize(edge);
    normalize(PK);

    dot1 = dot(PK, edge);
    if(fabs(dot1 - 1.0) < 1e-12) ////very important setting here!!!!!!
        dot1 = 1.;
    ang_pk_edge = acos(dot1);


    if(orient > 0)
    {
        PolarPK = c->r * ang_pk_edge + cur_ang;

        if(PolarPK >= 2 * M_PI)
            PolarPK -= 2 * M_PI;
    }

    else{
        PolarPK = cur_ang - c->r * ang_pk_edge;

        if(PolarPK < 0)
            PolarPK += 2 * M_PI;
    }

    ///////////////////////////////////////////////////////////////////////////////


    angle = p->t_angle - PolarPK;

    if(angle < -M_PI)
    {
        angle += 2 * M_PI;
    }
    if(angle > M_PI)
    {
        angle -= 2 * M_PI;
    }
    return angle;
}


/*************************************************************
Get the 2D local coordinates of a point inside a triangle
*************************************************************/
icVector2 Polyhedron::get_A_2D_Point(int faceid, double alpha[3])
{
    icVector2 result;
    Triangle *face = tlist.tris[faceid];

    result.entry[0] = alpha[1] * face->x1 + alpha[2] * face->x2;
    result.entry[1] = alpha[2] * face->y2;

    return result;
}


bool Polyhedron::is_Valid_Alpha(double alpha[])
{
    if(alpha[0] >= 0 && alpha[0] <= 1
        && alpha[1] >= 0 && alpha[1] <= 1
        && alpha[2] >= 0 && alpha[2] <= 1)
        return true;

    return false;  // not valid
}

void Polyhedron::local_To_Global(int faceid, icVector2 locpos, icVector3 &glpos)
{
    Triangle *face = tlist.tris[faceid];

    icVector3 dir3D = locpos.entry[0] * face->LX + locpos.entry[1] * face->LY;

    glpos.entry[0] = face->verts[0]->x + dir3D.entry[0];
    glpos.entry[1] = face->verts[0]->y + dir3D.entry[1];
    glpos.entry[2] = face->verts[0]->z + dir3D.entry[2];

}

/******************************************************
Get the bary centric coordinates in local frame
******************************************************/
void Polyhedron::get_2D_Barycentric_Facters(int faceid, double a, double b, double alpha[3])
{
    ////first, we need to calculate the local coordinates of the three vertices
    //// Maybe we should put them into the initialtracing routine

    Triangle *face = tlist.tris[faceid];

    alpha[2] = b/face->y2;

    if(fabs(alpha[2]) < 1e-10) alpha[2] = 0.;

    alpha[1] = (a - alpha[2]*face->x2)/ face->x1;

    if(fabs(alpha[1]) < 1e-10) alpha[1] = 0.;

    alpha[0] = 1 - alpha[1] - alpha[2];
}


/************************************************************************
Get the vector at any point inside a triangle using parallel transported
************************************************************************/
icVector2 Polyhedron::get_Vector_At_Point(int faceid, double vk[3], double alpha[3], double a, double b)
{
    icVector3 result;

    Triangle *face = tlist.tris[faceid];

    icVector2 vertsTensor[3], sum;

    Vertex *p, *pv, *pn;

    icVector2 vertsVec[2];   ////for two vectors on vertex having angle deficit

    double k[3];

    //////////////////////////////////////////////////
    ////codes for processing passing vertex 05/03/05
    if(fabs(alpha[0] - 1) < 1e-20)
    {
        p = face->verts[0];

        if(p->angle_deficit == 1 /*&& TraceRecalOn == 1*/)
        {
            pv = face->verts[2];
            k[0] = pv->x;
            k[1] = pv->y;
            k[2] = pv->z;
            vertsVec[0] = get_Miu(face->verts[0]->index, faceid, k);

            pn = face->verts[1];
            k[0] = pn->x;
            k[1] = pn->y;
            k[2] = pn->z;
            vertsVec[1] = get_Miu(face->verts[0]->index, faceid, k);

            sum = 0.5*(vertsVec[0] + vertsVec[1]);

        }
        else
            sum = get_Miu(face->verts[0]->index, faceid, vk);
    }

    else if(fabs(alpha[1] - 1) < 1e-20)
    {
        p = face->verts[1];

        if(p->angle_deficit == 1 )
        {
            pv = face->verts[0];
            k[0] = pv->x;
            k[1] = pv->y;
            k[2] = pv->z;
            vertsVec[0] = get_Miu(face->verts[1]->index, faceid, k);

            pn = face->verts[2];
            k[0] = pn->x;
            k[1] = pn->y;
            k[2] = pn->z;
            vertsVec[1] = get_Miu(face->verts[1]->index, faceid, k);

            sum = 0.5*(vertsVec[0] + vertsVec[1]);
        }
        else
            sum = get_Miu(face->verts[1]->index, faceid, vk);
    }

    else if(fabs(alpha[2] - 1) < 1e-20)
    {
        p = face->verts[2];

        if(p->angle_deficit == 1 /*&& TraceRecalOn == 1*/)
        {
            pv = face->verts[1];
            k[0] = pv->x;
            k[1] = pv->y;
            k[2] = pv->z;
            vertsVec[0] = get_Miu(face->verts[2]->index, faceid, k);

            pn = face->verts[0];
            k[0] = pn->x;
            k[1] = pn->y;
            k[2] = pn->z;
            vertsVec[1] = get_Miu(face->verts[2]->index, faceid, k);

            sum = 0.5*(vertsVec[0] + vertsVec[1]);
        }
        else
            sum = get_Miu(face->verts[2]->index, faceid, vk);
    }

    ///////////////////////////////////////////////////////////
    else{

        for(int i = 0; i < 3; i++)
            vertsTensor[i] = get_Miu(face->verts[i]->index, faceid, vk);

        /////sum the three tensors
        sum = alpha[0]*vertsTensor[0] + alpha[1]*vertsTensor[1] + alpha[2]*vertsTensor[2];
    }

    return sum;
}

void Polyhedron::initialize(){
    int i;
    icVector3 v1, v2;

    create_Pointers();
    calc_Edge_Length();

    calc_Bounding_Sphere();
    preprocess_Vertex();
    calc_Edge_Length();
    calc_Bounding_Sphere();
    calc_Face_Normals_and_Area();
    //average_Normals();
    get_vertex_normal();

    if(with_vec)
    {
        assign_color_VFmagnitude2(); /*assign color for vertex according to the vector magnitude*/

        /* normalize the vector values */
        normalize_Field();

        /* get the tangential vectors for all vertices if vector values have been read from file*/
        project_to_TangentPlane();

    }

    else
    {
        for(int i = 0; i < vlist.nverts; i++){
            //vlist.verts[i]->index = i;
            //vlist.verts[i]->nedges = 0;
            //vlist.verts[i]->ncorners = 0;
            vlist.verts[i]->visited = 0;
            vlist.verts[i]->distance = 0.;

            vlist.verts[i]->g_vec.entry[0] =
                vlist.verts[i]->g_vec.entry[1] =
                vlist.verts[i]->g_vec.entry[2] = 0.;
        }
    }

    init_Local_Frame();

    build_Corners();

    sort_Corner_on_Vertex();
    alloc_Corner_Angle();

    init_Trace(with_vec);


    cal_TexCoord();

    /*calculate the Jacobian for each triangle here*/
    cal_all_tri_Jacobians();

    mark_in_outlets();
    for (int i=0; i<tlist.ntris; i++)
    {
        /*  Modified by Guoning to avoid the region of zero vectors 03/09/2010  */
        Triangle *face = tlist.tris[i];
        bool has_zero_vert = false;
        for (int j=0; j<face->nverts; j++){
            if (length(face->verts[j]->t_vec) < 1.e-8)
            {
                has_zero_vert = true; break;
            }
        }
        if (has_zero_vert) face->has_zero_vec = true;
        else face->has_zero_vec = false;
    }

}


void
Polyhedron::mark_in_outlets()
{
    int i;

    for (i=0; i<vlist.nverts; i++)
    {
        icVector3 raw_vec = 1./1000.*vlist.verts[i]->t_vec;
        //normalize(raw_vec);
        //double re = dot (raw_vec, vlist.verts[i]->normal);

        //if (fabs(re)>1-0.1)
        //if (length(raw_vec) < 1.e-3)
        //	vlist.verts[i]->in_out_let = true;
        //else
        vlist.verts[i]->in_out_let = false;
    }
}


//void
//Polyhedron::write_file(FILE *file)
//{
//
//}



///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////



////Calculate the region normal at a specific edge
////This routine may need to be changed to 3D case
////We consider the normal in local frame ?? 12/22/05
////On 3D surfaces, an edge will have two normals for each type
void Edge::cal_normal_at_edge(Edge *cur_edge, Triangle *face, int which_face)
{
    icVector3 edge_3D, normal_3D;
    icVector2 edge_2D, normal;
    int edgecase;
    int orientation;

    ////Judge the orientation of the triangle, we need it to get the correct normal
    if(face->y2 > 0) ////It is counter clockwise orientation
    {
        orientation = 0;
    }
    else{                  ////It is clockwise orientation
        orientation = 1;
    }

    ////Get the edge vector in local frame
    edge_3D.entry[0] = cur_edge->verts[1]->x - cur_edge->verts[0]->x;
    edge_3D.entry[1] = cur_edge->verts[1]->y - cur_edge->verts[0]->y;
    edge_3D.entry[2] = cur_edge->verts[1]->z - cur_edge->verts[0]->z;
    edge_2D.entry[0] = dot(edge_3D, face->LX);
    edge_2D.entry[1] = dot(edge_3D, face->LY);

    ////
    edgecase = edge_vertex_case(cur_edge, face);
    ////calculate the outward normal for the edge
    if(orientation == 0) ////counter clockwise orientation
    {
        switch(edgecase)
        {
        case 0:   ////the same as triangle's orientation
        case 2:
        case 4:
            normal.entry[0] = edge_2D.entry[1];
            normal.entry[1] = -edge_2D.entry[0];
            break;

        case 1:   ////the inverse orientation
        case 3:
        case 5:
            normal.entry[0] = -edge_2D.entry[1];
            normal.entry[1] = edge_2D.entry[0];
            break;

        default:
            //MessageBox(NULL, "wrong edge!", "Error", MB_OK);
            fprintf(stderr, "wrong edge!\n");
        }
    }

    else{
        switch(edgecase)
        {
        case 0:   ////the same as triangle's orientation
        case 2:
        case 4:
            normal.entry[0] = -edge_2D.entry[1];
            normal.entry[1] = edge_2D.entry[0];
            break;

        case 1:   ////the inverse orientation
        case 3:
        case 5:
            normal.entry[0] = edge_2D.entry[1];
            normal.entry[1] = -edge_2D.entry[0];
            break;

        default:
            //MessageBox(NULL, "wrong edge!", "Error", MB_OK);
            fprintf(stderr, "wrong edge!\n");
        }
    }

    ////Transfer the normal back to 3D global frame
    //normal_3D = normal.entry[0] * face->LX + normal.entry[1] * face->LY;
    //normalize(normal_3D);

    //cur_edge->normal = normal_3D;

    normalize(normal);
    cur_edge->normal_2d = normal;

}



void
Edge::cal_normal_through_edge(Edge *cur_edge, Triangle *face, int which_face)
{
    /*
        compute the edge vector in 2D
    */
    icVector3 edge_3D, normal_3D;
    icVector2 edge_2D, normal;
    int edgecase;
    int orientation;

    ////Judge the orientation of the triangle, we need it to get the correct normal
    if(face->y2 > 0) ////It is counter clockwise orientation
    {
        orientation = 0;
    }
    else{                  ////It is clockwise orientation
        orientation = 1;
    }

    ////Get the edge vector in the local frame of "face"
    edge_3D.entry[0] = cur_edge->verts[1]->x - cur_edge->verts[0]->x;
    edge_3D.entry[1] = cur_edge->verts[1]->y - cur_edge->verts[0]->y;
    edge_3D.entry[2] = cur_edge->verts[1]->z - cur_edge->verts[0]->z;
    edge_2D.entry[0] = dot(edge_3D, face->LX);
    edge_2D.entry[1] = dot(edge_3D, face->LY);

    /*  rotate it 90 degree  */
    normal.entry[0] = -edge_2D.entry[1];
    normal.entry[1] = edge_2D.entry[0];

    normalize(normal);

    normal_3D = normal.entry[0]*face->LX + normal.entry[1]*face->LY;

    double test_dot = dot(edge_3D, normal_3D);

    /*
    get the triangle center
    */
    double center[3] = {0.};
    face->get_center(center);
    double mid_p[3] = {0.};
    mid_p[0] = (cur_edge->verts[0]->x + cur_edge->verts[1]->x)/2.;
    mid_p[1] = (cur_edge->verts[0]->y + cur_edge->verts[1]->y)/2.;
    mid_p[2] = (cur_edge->verts[0]->z + cur_edge->verts[1]->z)/2.;

    /*
       compute a vector from center to the edge middle point
    */
    icVector3 test_vec (mid_p[0]-center[0], mid_p[1]-center[1], mid_p[2]-center[2]);
    normalize(test_vec);

    double res = dot (test_vec, normal_3D);

    if (res < 0)
    {
        normal = -normal;
        normal_3D = -normal_3D;
    }

    cur_edge->normal_2d = normal;
}

////There are only 6 cases for each triangle
int Edge::edge_vertex_case(Edge *cur_edge, Triangle *face)
{
    if(cur_edge->verts[0] == face->verts[0] && cur_edge->verts[1] == face->verts[1])
        return 0;           ////case 0->1, the same as triangle's orientation
    else if(cur_edge->verts[0] == face->verts[1] && cur_edge->verts[1] == face->verts[0])
        return 1;           ////case 1->0, the inverse to triangle's orientation
    else if(cur_edge->verts[0] == face->verts[1] && cur_edge->verts[1] == face->verts[2])
        return 2;           ////case 1->2, the same as triangle's orientation
    else if(cur_edge->verts[0] == face->verts[2] && cur_edge->verts[1] == face->verts[1])
        return 3;           ////case 2->1, the inverse to triangle's orientation
    else if(cur_edge->verts[0] == face->verts[2] && cur_edge->verts[1] == face->verts[0])
        return 4;           ////case 2->0, the same as triangle's orientation
    else if(cur_edge->verts[0] == face->verts[0] && cur_edge->verts[1] == face->verts[2])
        return 5;           ////case 0->2, the inverse to triangle's orientation
    else
        return -1;          ////wrong
}





void Polyhedron::cal_TexCoord(void)
{
    double px, py, pz,  vx, vy, vz;

    double  dmax   = 3.13/800.;  //Hardcode here!!!

    for(int i = 0; i < vlist.nverts; i++)
    {

        vx = vlist.verts[i]->t_vec.entry[0] * dmax;
        vy = vlist.verts[i]->t_vec.entry[1] * dmax;
        vz = vlist.verts[i]->t_vec.entry[2] * dmax;

        px = vlist.verts[i]->x - vx;
        py = vlist.verts[i]->y - vy;
        pz = vlist.verts[i]->z - vz;

        vlist.verts[i]->texture_coord.entry[0] = px;
        vlist.verts[i]->texture_coord.entry[1] = py;
    }
}

/**---------------------------------------------------------------------------------------**/

/**---------------------------------------------------------------------------------------**/
/* Testing vector field */
void Polyhedron::XAxisVF(void)
{
    double nx, ny, nz;
    double vx, vy, vz;


    ////////////////////////////
    //Note: It can not judge the vertical situation!!!!!!!!
    for(int i = 0; i < vlist.nverts; i++)
    {
        nx = vlist.verts[i]->normal.entry[0];
        ny = vlist.verts[i]->normal.entry[1];
        nz = vlist.verts[i]->normal.entry[2];

        //The vector field used here is a constant vector field
        vx = (ny*ny + nz*nz);
        vy = -nx*ny;
        vz = -nx*nz;

        ///Normalize
        vlist.verts[i]->t_vec.entry[0] = vx;
        vlist.verts[i]->t_vec.entry[1] = vy;
        vlist.verts[i]->t_vec.entry[2] = vz;
        normalize(vlist.verts[i]->t_vec);

    }
    cal_TexCoord();
}


/*
Tranfer color from hsv mode to rgb mode
*/
void  HsvRgb( float hsv[3], float rgb[3] )
{
    float h, s, v;			// hue, sat, value
    float r, g, b;			// red, green, blue
    float i, f, p, q, t;		// interim values

    // guarantee valid input:
    h = hsv[0] / 60.;
    while( h >= 6. )	h -= 6.;
    while( h <  0. ) 	h += 6.;

    s = hsv[1];
    if( s < 0. )
        s = 0.;
    if( s > 1. )
        s = 1.;

    v = hsv[2];
    if( v < 0. )
        v = 0.;
    if( v > 1. )
        v = 1.;

    // if sat==0, then is a gray:
    if( s == 0.0 )
    {
        rgb[0] = rgb[1] = rgb[2] = v;
        return;
    }

    // get an rgb from the hue itself:
    i = floor( h );
    f = h - i;
    p = v * ( 1. - s );
    q = v * ( 1. - s*f );
    t = v * ( 1. - ( s * (1.-f) ) );

    switch( (int) i )
    {
    case 0:
        r = v;	g = t;	b = p;
        break;

    case 1:
        r = q;	g = v;	b = p;
        break;

    case 2:
        r = p;	g = v;	b = t;
        break;

    case 3:
        r = p;	g = q;	b = v;
        break;

    case 4:
        r = t;	g = p;	b = v;
        break;

    case 5:
        r = v;	g = p;	b = q;
        break;
    }

    rgb[0] = r;
    rgb[1] = g;
    rgb[2] = b;
}
