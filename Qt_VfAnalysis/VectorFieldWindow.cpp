#include "GL_LIB/glew.h"
#include "VectorFieldWindow.h"

// this file serves as IBFVDlg.cpp, GLView.cpp
double SCALE = 3.0;


// GLOBAL VARIABLES
const char filePath [] = "/Users/linxinw/Public/GitHub/QT-2D-Vector-Field-Analysis/Qt_VfAnalysis/Datasets/simple_fld2.ply";
int picked_node = -1;
double g_zoom_factor = 1.0;

/*------------------------------------------------------------*/
//global variables for ibfv
int     iframe = 0;
const int     Npat   = 64;        // Modified at 02/10/2010
int     alpha  = (0.06*255);  // modified for a smear out LIC texture for better visualization (TVCG revision 11/22/2010)
double  sa;
const double  tmax   = NPIX/(SCALE*NPN);
double  dmax   = 1.8/NPIX;

int Cal_Regions=0;
int ndisplay_trajs;
bool RemoveDisconnMSOn = false;
int Integrator_opt = 0;
int ndisplay_POs;
clock_t g_start, g_finish;

const double singsize = 0.009; // For visualizing the fixed points with the Morse sets together!
double trans_x=0;
double trans_y=0;

bool EnGreyTexture = false;
bool NotShiny = false;
bool DisableLighting = false;
bool DisableBackground = false;
bool EnTubeLikeVis = false;
bool ShowTDistribution = false;
bool ShowDiffECGMCG = false;
bool ShowMCGDiff = false;
bool ShowRotSumOn = false;
const int DebugOn = 1;

GLuint Textures[2];
GLubyte f_tex[NPIX][NPIX][3], b_tex[NPIX][NPIX][3], applied_tex[NPIX][NPIX][3];
GLubyte test_tex[NPIX][NPIX][3];
struct jitter_struct{
    double x;
    double y;
} jitter_para;

jitter_struct ji1[1] = {{0.0, 0.0}};
jitter_struct ji16[16] = {{0.125, 0.125}, {0.375, 0.125}, {0.625, 0.125},
                          {0.875, 0.125},{0.125, 0.375}, {0.375, 0.375}, {0.625, 0.375}, {0.875, 0.375},
                          {0.125, 0.625}, {0.375, 0.625}, {0.625, 0.625}, {0.875, 0.625},{0.125, 0.875},
                          {0.375, 0.875}, {0.625, 0.875}, {0.875, 0.875}};

// some extern varibles
extern Polyhedron *object; // declared in geometry.cpp

extern MorseDecomp *morse_decomp;
extern MorseDecomp *local_decomp;
extern TrajectoryList *separatrices;
extern PeriodicOrbitList *periodic_orbits;
extern EvenStreamlinePlace *evenplace;

extern MorseDecomp *L1_morse;
extern MorseDecomp *L2_morse;
extern ECG_Graph *ecg;
extern MCG_Graph *mcg;



VectorFieldWindow::VectorFieldWindow(QWidget *parent) : QOpenGLWidget(parent)
{
    // set surfaceFormat
//    format.setDepthBufferSize(32);
//    format.setStencilBufferSize(8);
//    //format.setAlphaBufferSize(8);
//    format.setVersion(3,2);
//    format.setProfile(QSurfaceFormat::CoreProfile);
//    format.setSwapBehavior(QSurfaceFormat::DoubleBuffer);
//    QSurfaceFormat::setDefaultFormat(format);
    qInfo() << "OpenGLWindow Default Constructor has finished";
    //qDebug() << QDir::currentPath();
}


VectorFieldWindow::~VectorFieldWindow() {
    qInfo() << "OpenGLWindow Desctuctor called";
}


void VectorFieldWindow::initializeGL(){
    this->makeCurrent();
    this->mat_ident(rotmat);
    for(int i = 0; i < 16; i++)
        this->ObjXmat[i]=0.;
    this->ObjXmat[0] = this->ObjXmat[5] = this->ObjXmat[10] = this->ObjXmat[15] = 1;

    this->init_flags();

    this->initializeOpenGLFunctions();

    // read file
    FILE *this_file = fopen(filePath, "r");
    if(this_file == 0){
        qInfo() << "Cannot open file: " << filePath;
        exit(-1);
    }
    QString str = "Reading file: " + QString(filePath);
    qInfo() << str;

    // build polyhedron
    object = new Polyhedron(this_file, 1);
    object->initialize();
    this->rot_center = object->rot_center;  ////set the rotation center!

    object->cal_TexCoord();

    morse_decomp = new MorseDecomp(); /*initialize the Morse Decomposition component*/

    L1_morse = new MorseDecomp();
    L2_morse = new MorseDecomp();

    this->InitGL(); // no drawing

    //construct the MCG here
    mcg = new MCG_Graph(); // create a new MCG
    mcg->init_MCG(); // init MCG

    //IBFV visualization initialization
    //this->makePatterns();
    //this->DrawGLScene(GL_RENDER);
    //this->ReCalTexcoord();


    qInfo() << "init gl Error: " << glGetError();
    qInfo() << "OpenGLWindow has been initialized";
}


void VectorFieldWindow::paintGL(){
    //this->DrawGLScene(GL_RENDER);
    this->draw();
    //this->update();
    qInfo() << "paintGL is Done";
}


void VectorFieldWindow::resizeGL(int w, int h){

//    qInfo() << "Resize w" << w << "  h" <<h;
//    glViewport(0, 0, w, h);
//    glMatrixMode(GL_PROJECTION);
//    glLoadIdentity();
//    glMatrixMode(GL_MODELVIEW);
//    glLoadIdentity();
}


void VectorFieldWindow::set_up_MainWindow_ptr(MainWindow *MW_ptr)
{
    this->mainWindow = MW_ptr;
}


void VectorFieldWindow::draw(){
    float r = 1.0f, g = 0.0f, b = 0.0f;
    glColor3f(r, g, b);
    glLoadIdentity();
    glTranslatef(-0.5, -0.5, 0);
    double min =10;
    double max = -10;
    glBegin(GL_TRIANGLES);
    for(int i = 0; i < object->tlist.ntris; i++){
        const Triangle * t1 = object->tlist.tris[i];
        const Vertex* v1 = t1->verts[0];
        const Vertex* v2 = t1->verts[1];
        const Vertex* v3 = t1->verts[2];
        glVertex3f(v1->x, v1->y, v1->z);
        glVertex3f(v2->x, v2->y, v2->z);
        glVertex3f(v3->x, v3->y, v3->z);
    }
    glEnd();
}


int VectorFieldWindow::InitGL( )
{
    glViewport(0, 0, (GLsizei) NPIX, (GLsizei) NPIX);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    return true;
}


void VectorFieldWindow::init_flags()
{
    this->ShowMCGOn = 0;
    this->ShowConleyCircle=0;

    this->MoveOrStop = 0;
    this->ShowFixedPtOn = 0;
    this->ShowSeparatricesOn = 0;
    this->ShowSCCsOn = 0;
    this->ShowPeriodicOrbitsOn = 0;
    this->EvenStreamlinePlacement = 0;
    this->ShowColorVFMagOn = 0;

    this->zoom_factor = 1;


    this->ShowMorseConn=0;
    this->ShowTriangleConn=0;
    this->morseC=0;
    this->triC=0;
    this->sccC=0;

    this->ShowEdgeSamplesOn = false;
    this->ShowTriMappingOn = false;
    this->selected_triangle = -1;  // 02/10/2010

    this->ShowBackward=0;//0=forward, 1=backward
    this->sampling_edge=0;

    this->ShowConnectionRegion=0;

    this->ShowBoundary=1;

    this->FlipNormalOn = false;
    this->IBFVOff = false;
}


void VectorFieldWindow::init_tex(){
    for(int i =0; i < NPIX; i++)
        for(int j =0; j < NPIX; j++){
            for(int w =0; w < 4; w++){
                f_tex[i][j][w] = 0;
                b_tex[i][j][w] = 0;
                applied_tex[i][j][w] = 0;
            }
        }
}


void VectorFieldWindow::makePatterns()
{
    this->init_tex();

    int lut[256];
    int phase[NPN][NPN];
    GLubyte pat[NPN][NPN][4];
    GLubyte spat[NPN][NPN][4];
    int i, j, k, t;
    int lut_r[256], lut_g[256], lut_b[256];

    //for (i = 0; i < 256; i++) lut[i] = i < 127 ? 0 : 255;
    for (i = 0; i < 256; i++) lut[i] = i < 127 ? 50 : 255;                 // for a brighter texture (TVCG revision 11/22/2010)

    for (i = 0; i < NPN; i++)
        for (j = 0; j < NPN; j++) phase[i][j] = rand() % 256;

    if ( EnGreyTexture )
    {
//        for (k = 0; k < Npat; k++) {
//            t = k*256/Npat;                           //t is used to control the animation of the image
//            for (i = 0; i < NPN; i++)
//                for (j = 0; j < NPN; j++) {
//                    pat[i][j][0] =
//                        pat[i][j][1] =
//                        pat[i][j][2] = lut[(t + phase[i][j]) % 255];
//                    pat[i][j][3] = alpha;

//                    spat[i][j][0] =
//                        spat[i][j][1] =
//                        spat[i][j][2] = lut[ phase[i][j] % 255];
//                    spat[i][j][3] = alpha;
//                }

//            glNewList(k + 1, GL_COMPILE);
//            glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0,
//                         GL_RGBA, GL_UNSIGNED_BYTE, pat);
//            glEndList();

//            glNewList(k + 1 + 100, GL_COMPILE);       //This is for static image
//            glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0,
//                         GL_RGBA, GL_UNSIGNED_BYTE, spat);
//            glEndList();
//        }
    }

    else
    {

        /*
           The following texture color can be adjusted to deliver the best visualization
       */
        for (i=0; i<256; i++)
        {
            lut_r[i] = 1.22*i<127?0:255;
            lut_g[i] = 1.3*i<127?0:255;
            lut_b[i] = 1.1*i<127?0:255;

        }

        for (i = 0; i < NPN; i++)
            for (j = 0; j < NPN; j++) phase[i][j] = rand() % 256;

        for (k = 0; k < Npat; k++) {
            t = k*256/Npat;                           //t is used to control the animation of the image
            for (i = 0; i < NPN; i++)
                for (j = 0; j < NPN; j++) {
                    pat[i][j][0] = lut_r[(t + phase[i][j]) % 255];
                    pat[i][j][1] = lut_g[(t + phase[i][j]) % 255];
                    pat[i][j][2] = lut_b[(t + phase[i][j]) % 255]/*lut[(t + phase[i][j]) % 255]*/;
                    pat[i][j][3] = alpha;

                    spat[i][j][0] = lut_r[ phase[i][j] % 255];
                    spat[i][j][1] = lut_g[ phase[i][j] % 255];
                    spat[i][j][2] = lut_b[ phase[i][j] % 255]/*lut[phase[i][j] % 255]*/;
                    spat[i][j][3] = alpha;
                }

//            glNewList(k + 1, GL_COMPILE);
//            glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0,
//                         GL_RGBA, GL_UNSIGNED_BYTE, pat);
//            glEndList();


//            if(k == 1){
//                output_pattern_to_file(spat, "spat");
//                exit(-5);
//            }

            glNewList(k + 1 + 100, GL_COMPILE);       //This is for static image
            glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0,
                         GL_RGBA, GL_UNSIGNED_BYTE, spat);
            glEndList();
        }
    }
    qInfo() << "makePatterns gl Error: " << glGetError();
}


int VectorFieldWindow::DrawGLScene(GLenum mode) // Here's Where We Do All The Drawing
{
    if (!IBFVOff) {
        //printf("IBFVOFF!");

        this->IBFVSEffect(mode);
    }
    else
        qInfo() << "something went wrong";
//    else
//    {
//        //printf("IBFVON!");
//        glClearColor (1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
//        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//        draw_shadedObj(mode);
//    }


//    glDisable(GL_TEXTURE_2D);


//    glEnable(GL_LIGHTING);
//    glEnable(GL_LIGHT0);

//    glEnable(GL_COLOR_MATERIAL);

//    glDepthFunc(GL_LEQUAL);
//    glEnable(GL_LINE_SMOOTH);
//    glEnable(GL_BLEND);
//    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//    glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);

//    /*display the extracted Morse sets with blending mode*/
//    glDisable(GL_LIGHTING);
//    glEnable(GL_COLOR_MATERIAL);
//    glEnable(GL_BLEND);


//    //show the rotational sum

//    glEnable(GL_BLEND);
//    if (ShowRotSumOn)
//    {
//        vis_rot_sum();
//    }
//    glDisable(GL_BLEND);

//    without_antialiasing(mode);


//    /*************************************************************
//    /*
//       Display the stored sample points along edges 02/09/2010
//    */
//    glDisable(GL_LIGHTING);
//    extern EdgeSamplePt_List *edge_samples;
//    if (ShowEdgeSamplesOn)
//    {
//        if (edge_samples != NULL)
//            //edge_samples->display(); display_sel_edges
//            edge_samples->display_sel_edges(selected_triangle,ShowBackward,sampling_edge);
//    }

//    if (ShowTriMappingOn)
//    {
//        morse_decomp->show_tri_mapping(selected_triangle);
//    }

//    /**************************************************************
//    /*
//        Test the triangle selection 02/10/2010
//    */
//    if (selected_triangle >=0 && selected_triangle < object->tlist.ntris)
//    {
//        display_sel_tri(selected_triangle);
//    }

//    glDepthFunc(GL_LESS);

//    glColor3f (1, 1, 1);


//    glDisable(GL_COLOR_MATERIAL);

    return true;
}


void VectorFieldWindow::vis_rot_sum()
{
    int i, j;
    float hsv[3] = {0, 1, 1}, rgb[3];

    float rot_rang = object->max_rot_sum - object->min_rot_sum;
    for (i=0; i<object->tlist.ntris; i++)
    {
        Triangle *t = object->tlist.tris[i];

        glBegin(GL_TRIANGLES);
        for (j=0; j<3; j++)
        {
            Vertex *v = t->verts[j];

            hsv[1] = 1. - 2.*(v->total_rot_sum - object->min_rot_sum)/rot_rang;

            if (hsv[1] < 0)
            {
                hsv[0] = 0.;
                hsv[1] = abs(hsv[1]);
            }
            else
                hsv[0] = 240.;

            HsvRgb(hsv, rgb);

            glColor4f (rgb[0], rgb[1], rgb[2], 0.5);
            glVertex3f (v->x, v->y, v->z);
        }
        glEnd();
    }
}


void VectorFieldWindow::ReCalTexcoord()
{
    double p[4], pr[4];

    double vx, vy, vz ;

    int i, j, k;

    glGetFloatv(GL_MODELVIEW_MATRIX, ObjXmat);

    for(i = 0; i < object->vlist.nverts; i++)
    {
        ////Initialize the variables nv and pr
        pr[0] = pr[1] = pr[2] = pr[3] = 0;

        object->vlist.verts[i]->g_vec = object->vlist.verts[i]->t_vec;
        normalize(object->vlist.verts[i]->g_vec);
        vx = object->vlist.verts[i]->g_vec.entry[0] * dmax * 2;
        vy = object->vlist.verts[i]->g_vec.entry[1] * dmax * 2;
        vz = object->vlist.verts[i]->g_vec.entry[2] * dmax * 2;

        ///Transform the vertex first!!! calculate the ri - vi
        p[0] = object->vlist.verts[i]->x - vx;
        p[1] = object->vlist.verts[i]->y - vy;
        p[2] = object->vlist.verts[i]->z - vz;
        p[3] = 1;


        for(k = 0; k < 4; k++)
            for(j = 0; j < 4; j++)
                pr[k] += ObjXmat[k + 4*j] * p[j];

        object->vlist.verts[i]->texture_coord.entry[0] = pr[0] ;
        object->vlist.verts[i]->texture_coord.entry[1] = pr[1] ;

        /*   compute the backward texture coordinates  */
        ///Transform the vertex first!!! calculate the ri - vi
        p[0] = object->vlist.verts[i]->x + vx;
        p[1] = object->vlist.verts[i]->y + vy;
        p[2] = object->vlist.verts[i]->z + vz;
        p[3] = 1;

        pr[0] = pr[1] = pr[2] = pr[3] = 0;

        for(k = 0; k < 4; k++)
            for(j = 0; j < 4; j++)
                pr[k] += ObjXmat[k + 4*j] * p[j];

        object->vlist.verts[i]->back_Tex.entry[0] = pr[0] ;
        object->vlist.verts[i]->back_Tex.entry[1] = pr[1] ;
    }
}


int version = 0;
void    VectorFieldWindow::IBFVSEffect(GLenum mode)
{
    qInfo() << "IBFVEFFECT Error1: " << glGetError();
    int i, j;
    Triangle *face;
    Vertex *v;

    GLfloat ambient[] = { 0.3, 0.3, 0.3, 1.0 };
    //GLfloat diffuse[] = { 1., 0.6, 0, 1.0};
    GLfloat diffuse[] = { 0.8, .8, 1., 1.0 };
    GLfloat specular[] = { 0.8, 0.8, 1.0, 1.0 };

    glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, f_tex);
    output_tex_to_file(f_tex, "f_tex");
    glClearColor (1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, f_tex);
    output_tex_to_file(f_tex, "f_tex");
    exit(0);

    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);


    glShadeModel(GL_FLAT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_TEXTURE_2D);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    ////////////////////////////////////////////////
    //glDrawBuffer(GL_BACK);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
                 GL_RGB, GL_UNSIGNED_BYTE, f_tex);

    set_view(GL_RENDER);


    glPushMatrix();
    set_scene(GL_RENDER);

    //////Texture distortion and mapping here
    for (i=0; i<object->tlist.ntris; i++) {
        face = object->tlist.tris[i];
        glBegin(GL_POLYGON);
        for (j=0; j<face->nverts; j++) {
            ///Using the storaged texture coordinates
            v = face->verts[j];
            glTexCoord2f(v->texture_coord.entry[0], v->texture_coord.entry[1]);
            glVertex3f(v->x, v->y, v->z);

        }
        glEnd();
    }

    glPopMatrix();


    ////Adding the depth judgement here!
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthFunc(GL_GREATER);

    if( MoveOrStop ==0 )   ////Static
        glCallList(iframe % Npat + 1 + 100);

//    else                   ////Moving
//        glCallList(iframe % Npat + 1);
    glBegin(GL_QUAD_STRIP);
        glTexCoord2f(0.0,  0.0);  glVertex3f(0.0, 0.0, -39.9);
        glTexCoord2f(0.0,  tmax); glVertex3f(0.0, 1.0, -39.9);
        glTexCoord2f(tmax, 0.0);  glVertex3f(1.0, 0.0, -39.9);
        glTexCoord2f(tmax, tmax); glVertex3f(1.0, 1.0, -39.9);
    glEnd();
    glDepthFunc(GL_LESS);
    glDisable(GL_BLEND);



    ////Store the vector field texture for next shading blending and next texture generating

    //glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 0, 0, NPIX, NPIX, 0);
    //glReadBuffer(GL_BACK);


    glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, f_tex);

    if(MoveOrStop == 0)
    {
//        /* Calculate backward texture */
//        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
//                     GL_RGB, GL_UNSIGNED_BYTE, b_tex);

//        glPushMatrix ();
//        set_scene(GL_RENDER);


//        //////Texture distortion and mapping here
//        for (i=0; i<object->tlist.ntris; i++) {
//            face = object->tlist.tris[i];
//            glBegin(GL_POLYGON);
//            for (j=0; j<face->nverts; j++) {
//                ///Using the storaged texture coordinates
//                v = face->verts[j];
//                glTexCoord2f(v->back_Tex.entry[0], v->back_Tex.entry[1]);
//                glVertex3f(v->x, v->y, v->z);
//            }
//            glEnd();
//        }

//        glPopMatrix();

//        iframe = iframe + 1;

//        ////Adding the depth judgement here!
//        glEnable(GL_BLEND);
//        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//        glDepthFunc(GL_GREATER);
//        //if( MoveOrStop ==0 )   ////Static
//        glCallList(iframe % Npat + 1 + 100);
//        if( MoveOrStop ==0 )   ////Static
//            glCallList(iframe % Npat + 1 + 100);
////        else                   ////Moving
////            glCallList(iframe % Npat + 1);
//        glBegin(GL_QUAD_STRIP);
//            glTexCoord2f(0.0,  0.0);  glVertex3f(0.0, 0.0, -39.9);
//            glTexCoord2f(0.0,  tmax); glVertex3f(0.0, 1.0, -39.9);
//            glTexCoord2f(tmax, 0.0);  glVertex3f(1.0, 0.0, -39.9);
//            glTexCoord2f(tmax, tmax); glVertex3f(1.0, 1.0, -39.9);
//        glEnd();
//        glDepthFunc(GL_LESS);
//        glDisable(GL_BLEND);

//        //glReadBuffer(GL_BACK);
//        glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, b_tex);

//        //blend two images

//        for(int x = 0; x < NPIX; x++)
//        {
//            for(int y = 0; y < NPIX; y++)
//            {
//                /**/
//                int temp_color = (int)(f_tex[x][y][0] + b_tex[x][y][0])/2.;
//                applied_tex[x][y][0] = (int)(f_tex[x][y][0] + b_tex[x][y][0])/2;
//                applied_tex[x][y][1] = (int)(f_tex[x][y][1] + b_tex[x][y][1])/2;
//                applied_tex[x][y][2] = (int)(f_tex[x][y][2] + b_tex[x][y][2])/2;
//            }
//        }
    }


    ///////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    ////Blending shading here 12/1

//    /*  First, draw the background  */
//    if (!DisableLighting && !DisableBackground)
//    {
//        glClearColor (0.0, 0.0, 0.0, 1.0);  // background for rendering color coding and lighting
//        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
//        glDisable(GL_DEPTH_TEST);
//        glDisable(GL_LIGHTING);
//        glDisable(GL_TEXTURE_2D);
//        glShadeModel(GL_SMOOTH);

//        glMatrixMode(GL_PROJECTION);
//        glPushMatrix();
//        glLoadIdentity();
//        gluOrtho2D(-1., 1., -1., 1.);

//        glMatrixMode(GL_MODELVIEW);
//        glPushMatrix();
//        glLoadIdentity();

//        glBegin(GL_QUADS);
//        //red color
//        //glColor3f(0.4,0.4,0.3);
//        glColor3f(0.1,0.1,0.3);
//        glVertex2f(-1.0,-1.0);
//        glVertex2f(1.0,-1.0);
//        //blue color
//        glColor3f(1.,1.,1.);
//        glVertex2f(1.0, 1.0);
//        glVertex2f(-1.0, 1.0);
//        glEnd();

//        glMatrixMode(GL_PROJECTION);
//        glPopMatrix();
//        glMatrixMode(GL_MODELVIEW);
//        glPopMatrix();
//    }
//    ////////////////////////////////////////////

    glEnable(GL_DEPTH_TEST);

    set_scene(GL_RENDER);

    glShadeModel(GL_SMOOTH);

    if (!DisableLighting)
    {
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        //glEnable(GL_LIGHT1);
        glEnable(GL_NORMALIZE);
    }

    else
    {
        glDisable(GL_LIGHTING);
    }
    //glEnable(GL_RESCALE_NORMAL);


    glEnable(GL_TEXTURE_2D);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    //glDisable(GL_COLOR_MATERIAL);

    if (MoveOrStop == 0)
    {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
                     GL_RGB, GL_UNSIGNED_BYTE, applied_tex);
    }
//    else
//        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
//                     GL_RGB, GL_UNSIGNED_BYTE, f_tex);


    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT_AND_BACK/*GL_FRONT*/, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT_AND_BACK/*GL_FRONT*/, GL_SPECULAR, specular);
    //glMaterialf(GL_FRONT, GL_SHININESS, (GLfloat)shiny);
    glMateriali(GL_FRONT_AND_BACK/*GL_FRONT*/, GL_SHININESS, 80);

//    for (i=0; i<object->tlist.ntris; i++) {
//        if (mode == GL_SELECT)
//            glLoadName(i+1);

//        face = object->tlist.tris[i];

//        /*
//            Do not display the texture if it is part of the in/outlets of the data
//        */
//        if (face->verts[0]->in_out_let || face->verts[1]->in_out_let || face->verts[2]->in_out_let)
//            continue;

//        glBegin(GL_POLYGON);
//        for (j=0; j<face->nverts; j++) {
//            v = face->verts[j];
//            glNormal3d(v->normal.entry[0],
//                       v->normal.entry[1],
//                       v->normal.entry[2]);
//            glTexCoord2f(v->texture_coord.entry[0], v->texture_coord.entry[1]);
//            glVertex3f(v->x, v->y, v->z);
//        }
//        glEnd();

//    }

    glDisable(GL_TEXTURE_2D);
    qInfo() << "IBFVEFFECT Error2: " << glGetError();
    version++;
    iframe = iframe + 1;
}


void    VectorFieldWindow::draw_shadedObj(GLenum mode)
{
    int i, j;
    Triangle *face;

    glClearColor (1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    GLfloat ambient[] = { 0.5, 0.5, 0.5, 1.0 };
    GLfloat diffuse[] = { 0.8, 0.8, 0.8, 1.0 };
    GLfloat specular[] = { 0.9, 0.9, 0.9, 1.0 };
    int shiny = 100;

    set_view(GL_RENDER);

    set_scene(GL_RENDER);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT, GL_SHININESS, (GLfloat)shiny);


    for (i=0; i<object->tlist.ntris; i++) {
        if (mode == GL_SELECT)
            glLoadName(i+1);

        face = object->tlist.tris[i];

        glBegin(GL_POLYGON);
        for (j=0; j<face->nverts; j++) {
            glNormal3d(face->verts[j]->normal.entry[0],
                       face->verts[j]->normal.entry[1],
                       face->verts[j]->normal.entry[2]);
            glVertex3f(face->verts[j]->x, face->verts[j]->y, face->verts[j]->z);
        }
        glEnd();
    }
}


void    VectorFieldWindow::without_antialiasing(GLenum mode)
{
    if(ShowConnectionRegion==1)
    {
        display_MCG_connections();
    }

    if(ShowSCCsOn == 1)
    {
        display_SCCs(mode);
    }



    if(ShowColorVFMagOn == 1)
    {
        display_color_VFMag();
    }

    if (ShowTDistribution)
        display_tau_distribution();

    if (ShowDiffECGMCG)
        display_diff_ECG_MCG();

    if (ShowMCGDiff)
        display_diff_MCGs();

    glDisable(GL_BLEND);

    ////show separatrices
    if(ShowSeparatricesOn == 1)
        display_separatrices(mode);


    //////Display all the detected limit cycles
    if(ShowPeriodicOrbitsOn == 1)
        display_periodicorbits(mode);


    //glDisable(GL_COLOR_MATERIAL);
    if(EvenStreamlinePlacement == 1)
    {
        display_even_streamlines(mode); /*display the calculated evenly placed streamlines*/
    }

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHT2);

    /* Display the legends for the objects in the vector field */
    if(ShowFixedPtOn == 1)
        display_FixedPtsIcon(mode);
}


void    VectorFieldWindow::set_view(GLenum mode)
{
    icVector3 up, ray, view;
    GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
    GLfloat light_diffuse0[] = { 0.6, 0.6, 0.6, 1.0 };
    GLfloat light_specular0[] = { 1.0, 0.2, 0.0, 1.0 };
    GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat light_diffuse1[] = { 0.6, 0.6, 0.6, 1.0 };
    GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat light_ambient2[] = { 0.3, 0.3, 0.3, 1.0 };
    GLfloat light_diffuse2[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_specular2[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_position[] = { 0.0, 0.0, 0.0, 1.0 };

    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);


    glMatrixMode(GL_PROJECTION);

    if (mode == GL_RENDER)
        glLoadIdentity();

    glOrtho(0., 1., 0., 1., -1000.0, 4000.0);
    //gluPerspective(45.0, 1.0, 0.1, 40.0);


    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    ////////////////////////////////////////////////////////////////////////

    GLfloat fLightPos[4]   = { -1000.0f, 0., 50.0f, 1.0f };  // Point source
    GLfloat fNoLight[] = { 0.0f, 0.0f, 0.0f, 0.0f };
    GLfloat fLowLight[] = { 0.3f, 0.3f, 0.3f, 1.0f };
    GLfloat fBrightLight[] = { 0.6f, 0.6f, 0.6f, 1.0f };

    if (!NotShiny)
    {
        glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);
    }
    else
    {
        glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SINGLE_COLOR/*GL_SEPARATE_SPECULAR_COLOR*/);
    }

    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, fNoLight);
    glLightfv(GL_LIGHT0, GL_AMBIENT, fLowLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, fBrightLight);
    glLightfv(GL_LIGHT0, GL_SPECULAR, fBrightLight);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    // Position light before any other transformations
    glLightfv(GL_LIGHT0, GL_POSITION, fLightPos);

    GLfloat fLightPos2[4]   = { 1000.0f, -10000.0f, -5000.0f, 1.0f };  // Point source
    glLightfv(GL_LIGHT1, GL_AMBIENT, fLowLight);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, fBrightLight);
    glLightfv(GL_LIGHT1, GL_SPECULAR, fBrightLight);
    glLightfv(GL_LIGHT1, GL_POSITION, fLightPos2);
}


void    VectorFieldWindow::set_scene(GLenum mode)
{
    //glTranslatef(trans_x, trans_y, 0);
    //glTranslatef(0.0, 0.0, -10.0);
    glTranslatef(this->rot_center.entry[0] ,
                 this->rot_center.entry[1],
                 this->rot_center.entry[2]);


    //multmatrix( this->rotmat );

    glScalef(this->zoom_factor, this->zoom_factor, this->zoom_factor);

    glScalef(0.7, 0.7, 0.7);

    glTranslatef(-this->rot_center.entry[0] ,
                 -this->rot_center.entry[1] ,
                 -this->rot_center.entry[2] );
}

void VectorFieldWindow::set_ColorByType(unsigned char type)
{
    if(type == SOURCE){
        glColor4f(0.,1., 0.,1);
    }

    else if(type == SINK){
        glColor4f(1.,0., 0.,1);
    }

    else if(type == SADDLE){
        glColor4f(0.,0.,1., 1);
    }

    else if(type == CWCENTER){
        glColor4f(1.,0., 1.,1);
    }

    else if(type == CCWCENTER){
        //glColor4f(0.3, 1., 1.,1);
        glColor4f(1.,0., 1.,1);
    }

    else if(type == AFOCUS){
        //glColor4f(1.,0.5, 0, 1);
        glColor4f(1.,0., 0.,1);
    }

    else if(type == RFOCUS){
        //glColor4f(0, 1, 0.7, 1);
        glColor4f(0.,1., 0.,1);
    }
}

void    VectorFieldWindow::draw_ASolidSphere(double x, double y, double z, double r)
{
    glPushMatrix();
    glTranslatef(x, y, z);
    //glutSolidSphere(r, 80, 80 );

    // a replacement of glutSolidSphere
    GLUquadric* quadric = gluNewQuadric();

    glTranslated(x, y, z);
    //glColor3f(R, G, B);
    gluSphere(quadric, r, 16, 16);
    glPopMatrix();
    gluDeleteQuadric(quadric);
}


void    VectorFieldWindow::multmatrix(const Matrix m)
{
    int i,j, index = 0;

    GLfloat mat[16];

    for ( i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            mat[index++] = m[i][j];

    glMultMatrixf (mat);
}


void    VectorFieldWindow::mat_ident(Matrix m)
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


void    VectorFieldWindow::display_MCG_connections()
{
    int i, j, k;

    int node1, node2;
    Vertex *v;
    icVector3 z(0, 0, -1.);
    float rot_mat[4][4];
    GLUquadricObj *q;
    q = gluNewQuadric();

    for (i=0; i<object->tlist.ntris; i++)
    {
        object->tlist.tris[i]->visited = false;
    }

    for(i=0; i<mcg->elist->nedges; i++)
    {
        if(mcg->elist->edges[i]->cancel)
            continue;

        node1 = mcg->elist->edges[i]->node_index1;
        node2 = mcg->elist->edges[i]->node_index2;

        for(j=0; j< mcg->elist->edges[i]->ntris; j++)
        {

            float rgb[3] = {0.};
            glColor4f(0.8, 0.8, 0.8, 0.4);


            Triangle *face = object->tlist.tris[mcg->elist->edges[i]->triangles[j]];

            if (face->visited)
                continue;

            face->visited = true;

            if (morse_decomp->scclist->scccomponents[morse_decomp->dg->nlist->dirnodes[face->index]->sscomp_index]->node_index >=0)
                continue;

            glBegin(GL_TRIANGLES);
            for(k=0; k<3; k++)
            {
                Vertex* v = face->verts[k];
                glVertex3f(v->x+0.001*v->normal.entry[0], v->y+0.001*v->normal.entry[1], v->z+0.001*v->normal.entry[2]);

            }
            glEnd();

            double center[3] = {0.};
            Vertex* v0=face->verts[0];
            Vertex* v1=face->verts[1];
            Vertex* v2=face->verts[2];

            center[0] = (v0->x+v1->x+v2->x)/3.;
            center[1] = (v0->y+v1->y+v2->y)/3.;
            center[2] = (v0->z+v1->z+v2->z)/3.;


            if(mcg->nlist->mnodes[node1]->type == 0)
                glColor4f(0, 1., 0, 0.6);
            else if(mcg->nlist->mnodes[node2]->type == 1)
                glColor4f(1., 0, 0, 0.6);
            else
                glColor4f(0.0, 0.0, 1., 0.6);


            extern double g_zoom_factor;
            extern void
            cal_rot_mat(icVector3 vec1, icVector3 vec2, float rotmat[4][4]);
            extern void
            transpose_mat(float mat[4][4]);
            cal_rot_mat(z, face->normal, rot_mat);
            extern void
            get_rot_between_two_vecs(double vec1[3], double vec2[3], float rot_mat[4][4]);

            //get_rot_between_two_vecs(z.entry, face->normal.entry, rot_mat);
            transpose_mat(rot_mat);

            glPushMatrix();
            glTranslatef(center[0]+0.002*face->normal.entry[0],
                         center[1]+0.002*face->normal.entry[1],
                         center[2]+0.002*face->normal.entry[2]);
            // we need to define a rotation here
            multmatrix(rot_mat);
            //gluDisk (q, 0.0, 0.003/g_zoom_factor, 20, 1); // 0.003 is for the diesel slices
            gluDisk (q, 0.0, 0.004/g_zoom_factor, 20, 1);
            glPopMatrix();

        }
    }
}


void    VectorFieldWindow::display_SCCs(GLenum mode)
{
    int i, j, k;
    Triangle *face;
    Vertex *vert;
    icVector3 outn;
    int colorcount = 0;

    if(morse_decomp == NULL || morse_decomp->scclist == NULL)
        return;


    for(i = 0; i < morse_decomp->scclist->nsccs; i++)
    {
        if(morse_decomp->scclist->scccomponents[i] == NULL)
            return;

        if(morse_decomp->scclist->scccomponents[i]->nnodes <= 2
            && morse_decomp->scclist->scccomponents[i]->nfixedpoints <= 0)
            continue;

        else if(morse_decomp->scclist->scccomponents[i]->nnodes > 2
                 && !morse_decomp->scclist->scccomponents[i]->valid)
            continue;

        else
        {

            if (morse_decomp->scclist->scccomponents[i]->node_index < 0)
                continue;

            float alpha_channel = 0.4;
            float color_rgb[3];

            switch(mcg->nlist->mnodes[morse_decomp->scclist->scccomponents[i]->node_index]->type)
            {
            case 0:glColor4f(0, 1, 0, 0.4);//0-source_like, 1-sink_like, 2-saddle_like
                color_rgb[0]=0; color_rgb[1]=1; color_rgb[2]=0;
                break;
            case 1:glColor4f(1, 0, 0, 0.4);
                color_rgb[0]=1; color_rgb[1]=0; color_rgb[2]=0;
                break;
            case 2:glColor4f(0, 0, 1, 0.4);
                color_rgb[0]=0; color_rgb[1]=0; color_rgb[2]=1;
                break;

            case 3: glColor4f (0, 0, 1, 0.4); // trivial index
                color_rgb[0]=0; color_rgb[1]=0; color_rgb[2]=1;
                break;

            }

            if (mcg->nlist->mnodes[morse_decomp->scclist->scccomponents[i]->node_index]->conley[0]==0
                &&mcg->nlist->mnodes[morse_decomp->scclist->scccomponents[i]->node_index]->conley[1]==0
                && mcg->nlist->mnodes[morse_decomp->scclist->scccomponents[i]->node_index]->conley[2]==0)
            {
                //glColor4f (1, 0, 1, 0.4); // trivial index
            }


            glBegin(GL_TRIANGLES);
            for(j = 0; j < morse_decomp->scclist->scccomponents[i]->nnodes; j++)
            {
                face = object->tlist.tris[morse_decomp->scclist->scccomponents[i]->nodes[j]];
                /////////////////////////////////////////////////////////////////
                //check if it is in the refined morse set, 12/03/2009 by Qingqing
                /////////////////////////////////////////////////////////////////
                if(face->cur_morse)continue;

                for(k = 0; k < face->nverts; k++)
                {
                    vert = face->verts[k];

                    if (FlipNormalOn)
                        outn = -vert->normal;
                    else
                        outn = vert->normal;
                    glNormal3dv(vert->normal.entry);
                    glVertex3d(vert->x + 0.0016*outn.entry[0],
                               vert->y + 0.0016*outn.entry[1],
                               vert->z + 0.0016*outn.entry[2]);
                }
            }
            glEnd();

            //draw the boundary
            if(ShowBoundary==1)
            {
                glDepthFunc(GL_LEQUAL);
                glColor4f(0,0,0,0.6);
                glBegin(GL_LINES);
                glLineWidth(10.);
                int morse_index=morse_decomp->scclist->scccomponents[i]->node_index;
                int EN=mcg->nlist->mnodes[morse_index]->boundaryN;
                for(j=0;j<EN;j++)
                {
                    Edge *e=object->elist.edges[mcg->nlist->mnodes[morse_index]->boundary[j]];
                    Vertex *v1=e->verts[0];
                    Vertex *v2=e->verts[1];
                    glVertex3f(v1->x+0.0016*v1->normal.entry[0],v1->y+0.0016*v1->normal.entry[1],v1->z+0.0016*v1->normal.entry[2]);
                    glVertex3f(v2->x+0.0016*v2->normal.entry[0],v2->y+0.0016*v2->normal.entry[1],v2->z+0.0016*v2->normal.entry[2]);

                }
                glEnd();

                glDepthFunc(GL_LESS);
            }
        }
    }


    /*Display the highlighted SCC*/
    if( picked_node >= 0 )
    {
        int *tris = morse_decomp->scclist->scccomponents[mcg->nlist->mnodes[picked_node]->scc_index]->nodes;
        int ntris = morse_decomp->scclist->scccomponents[mcg->nlist->mnodes[picked_node]->scc_index]->nnodes;
        for(i = 0; i < ntris; i++)
        {
            glColor3f(1, 1, 1);

            if(tris[i] < 0)
                continue;

            face = object->tlist.tris[tris[i]];

            glBegin(GL_LINE_LOOP);
            for(j = 0; j < 3; j++)
            {
                //glVertex2f(Object.vlist[face->verts[j]]->x, Object.vlist[face->verts[j]]->y);
                vert = face->verts[j];

                if (FlipNormalOn)
                    outn = -vert->normal;
                else
                    outn = vert->normal;
                glVertex3d(vert->x + 0.005*outn.entry[0],
                           vert->y + 0.005*outn.entry[1],
                           vert->z + 0.005*outn.entry[2]);
            }
            glEnd();
        }
    }

    ///////////////////////////////////////////////////
    //local morse set to black, 12/03/2009 by Qingqing
    ////////////////////////////////////////////////////
    if(local_decomp == NULL || local_decomp->scclist == NULL)
        return;

    //if(

    for(i = 0; i < local_decomp->scclist->nsccs; i++)
    {
        if(local_decomp->scclist->scccomponents[i] == NULL)
            return;

        if(local_decomp->scclist->scccomponents[i]->nnodes <= 2
            && local_decomp->scclist->scccomponents[i]->nfixedpoints <= 0)
            continue;

        else if(local_decomp->scclist->scccomponents[i]->nnodes > 2
                 && !local_decomp->scclist->scccomponents[i]->valid)
            continue;

        else
        {

            //glColor4f(0, 0, 0, 1);//black
            int gSCC=local_decomp->scclist->scccomponents[i]->global_SCC;
            int mcgID=morse_decomp->scclist->scccomponents[gSCC]->node_index;
            switch(mcg->nlist->mnodes[mcgID]->type)
            {
            case 0:glColor4f(0, 1, 0, 0.4);//0-source_like, 1-sink_like, 2-saddle_like
                break;
            case 1:glColor4f(1, 0, 0, 0.4);
                break;
            case 2:glColor4f(0, 0, 1, 0.4);
                break;

            }

            glBegin(GL_TRIANGLES);
            for(j = 0; j < local_decomp->scclist->scccomponents[i]->nnodes; j++)
            {
                int nid=local_decomp->dg->nlist->dirnodes[local_decomp->scclist->scccomponents[i]->nodes[j]]->global_index;
                face = object->tlist.tris[nid];

                for(k = 0; k < face->nverts; k++)
                {
                    vert = face->verts[k];
                    outn = vert->normal;
                    glNormal3dv(vert->normal.entry);
                    glVertex3d(vert->x + 0.001*outn.entry[0],
                               vert->y + 0.001*outn.entry[1],
                               vert->z + 0.001*outn.entry[2]);
                }
            }
            glEnd();
        }
    }
}


void    VectorFieldWindow::display_color_VFMag()
{
    int i, j;
    Triangle *face;
    Vertex *v;
    icVector3 outn;

    for(i=0; i<object->tlist.ntris; i++)
    {
        face = object->tlist.tris[i];

        glBegin(GL_TRIANGLES);
        for(j = 0; j < face->nverts; j++)
        {
            v = face->verts[j];
            //glColor4f(v->mag_rgb[0], v->mag_rgb[1], v->mag_rgb[2], 0.3);  // for other data, we use 0.3
            glColor4f(v->mag_rgb[0], v->mag_rgb[1], v->mag_rgb[2], 0.55);  // for earthquake data

            if (FlipNormalOn)
                outn = -v->normal;
            else
                outn = v->normal;
            glVertex3d(v->x + 0.001*outn.entry[0],
                       v->y + 0.001*outn.entry[1],
                       v->z + 0.001*outn.entry[2]);
        }
        glEnd();
    }
}

void    VectorFieldWindow::display_separatrices(GLenum mode)
{
    //Draw the separatrices begin from each saddle and long their major direction

    int i;

    glLineWidth(2.);

    if(separatrices == NULL)
        return;

    for(i = 0; i < separatrices->ntrajs; i+=4)
    {
        display_single_traj(separatrices->trajs[i], 1);

        display_single_traj(separatrices->trajs[i+1], 2);

        display_single_traj(separatrices->trajs[i+2], 3);

        display_single_traj(separatrices->trajs[i+3], 4);

    }
}


void    VectorFieldWindow::display_periodicorbits(GLenum mode)
{
    int i, j;
    icVector3 outn;

    if(periodic_orbits == NULL)
        return;

    /***********************************************************************************/
    /*   Display the trajectory as tube-like    */
    GLUquadricObj *quadratic = NULL;				// Storage For Our Quadratic Objects
    if (quadratic == NULL)
        quadratic=gluNewQuadric();			// Create A Pointer To The Quadric Object ( NEW )
    gluQuadricNormals(quadratic, GLU_SMOOTH);	// Create Smooth Normals ( NEW )
    gluQuadricTexture(quadratic, GL_TRUE);
    /////////////////////////////////////////////////////////////////////////////////////


    glDepthFunc(GL_LEQUAL);
    glLineWidth(3.0);
    for(i = 0; i < ndisplay_POs/*periodic_orbits->nporbits*/; i++)
    {
        if(periodic_orbits->polist == NULL || periodic_orbits->polist[i] == NULL)
            return;

        if(periodic_orbits->polist[i]->type == 0)
            glColor3f(0, 1., 0);
        else
            glColor3f(1, 0., 0);

        if(periodic_orbits->polist[i]->traj == NULL)
            return;

        for(j = 0; j < periodic_orbits->polist[i]->traj->nlinesegs; j++)
        {
            if(periodic_orbits->polist[i]->traj->linesegs == NULL )//||
                //periodic_orbits->polist[i]->traj->linesegs[j] == NULL)
                return;

            if (FlipNormalOn)
                outn = -object->tlist.tris[periodic_orbits->polist[i]->traj->linesegs[j].Triangle_ID]->normal;

            else
                outn = object->tlist.tris[periodic_orbits->polist[i]->traj->linesegs[j].Triangle_ID]->normal;

            glBegin(GL_LINES);
            glVertex3f(periodic_orbits->polist[i]->traj->linesegs[j].gstart.entry[0]+ 0.002*outn.entry[0],
                       periodic_orbits->polist[i]->traj->linesegs[j].gstart.entry[1]+ 0.002*outn.entry[1],
                       periodic_orbits->polist[i]->traj->linesegs[j].gstart.entry[2]+ 0.002*outn.entry[2]);
            glVertex3f(periodic_orbits->polist[i]->traj->linesegs[j].gend.entry[0]+ 0.002*outn.entry[0],
                       periodic_orbits->polist[i]->traj->linesegs[j].gend.entry[1]+ 0.002*outn.entry[1],
                       periodic_orbits->polist[i]->traj->linesegs[j].gend.entry[2]+ 0.002*outn.entry[2]);
            glEnd();
            /*
                If we display the separatrices as tube-like
            */
            if (EnTubeLikeVis)
            {
                glEnable(GL_LIGHTING);
                glEnable(GL_LIGHT0);
                glEnable(GL_LIGHT1);
                icVector3 Z (periodic_orbits->polist[i]->traj->linesegs[j].gend.entry[0]-periodic_orbits->polist[i]->traj->linesegs[j].gstart.entry[0],
                            periodic_orbits->polist[i]->traj->linesegs[j].gend.entry[1]-periodic_orbits->polist[i]->traj->linesegs[j].gstart.entry[1],
                            periodic_orbits->polist[i]->traj->linesegs[j].gend.entry[2]-periodic_orbits->polist[i]->traj->linesegs[j].gstart.entry[2]);
                double len = length (Z);

                normalize (Z);

                //That will give you the Z axis of the rotation matrix you need to set.
                //Then, make another vector which is

                icVector3 Y (Z.entry[1], -Z.entry[2], 0);
                normalize (Y);

                icVector3 X = cross (Y, Z);
                normalize (X);

                icVector3 Y2 = cross (X, Z);
                normalize (Y2);

                float rotMat[4][4];

                matrix_ident (rotMat);

                rotMat[0][0] = X.entry[0];
                rotMat[0][1] = X.entry[1];
                rotMat[0][2] = X.entry[2];

                rotMat[1][0] = Y2.entry[0];
                rotMat[1][1] = Y2.entry[1];
                rotMat[1][2] = Y2.entry[2];

                rotMat[2][0] = Z.entry[0];
                rotMat[2][1] = Z.entry[1];
                rotMat[2][2] = Z.entry[2];

                glMatrixMode(GL_MODELVIEW);
                glPushMatrix();

                glTranslatef (periodic_orbits->polist[i]->traj->linesegs[j].gstart.entry[0],
                             periodic_orbits->polist[i]->traj->linesegs[j].gstart.entry[1],
                             periodic_orbits->polist[i]->traj->linesegs[j].gstart.entry[2]);

                glMultMatrixf ((GLfloat *) rotMat);

                gluCylinder (quadratic, 0.006*object->radius/zoom_factor, 0.006*object->radius/zoom_factor, len, 10, 10);
                glMatrixMode(GL_MODELVIEW);
                glPopMatrix();
            }
        }
    }

    glDepthFunc(GL_LESS);
}


void    VectorFieldWindow::display_even_streamlines(GLenum mode)
{
    int i;

    for(i = 0; i < ndisplay_trajs/*evenplace->evenstreamlines->ntrajs*/; i++)
    {

        if(evenplace->evenstreamlines->trajs[i] != NULL)
            display_single_traj(evenplace->evenstreamlines->trajs[i], 5);
    }
}


void    VectorFieldWindow::display_FixedPtsIcon(GLenum mode)
{
    int i;

    ////set the surface property
    GLfloat ambient[] = { 0.8, 0.8, 0.8, 1.0 };
    GLfloat diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat specular[] = { 1.0, 1.0, 1.0, 1.0 };
    int shiny = 100;

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT, GL_SHININESS, (GLfloat)shiny);

    glEnable(GL_COLOR_MATERIAL);

    for(i = 0; i < object->slist.nsingularities; i++)
    {
        if(mode == GL_SELECT)
            glLoadName(NAMEOFSINGULARITY+i);

        set_ColorByType(object->slist.slist[i]->type);

        draw_ASolidSphere(object->slist.slist[i]->gpos.entry[0],
                          object->slist.slist[i]->gpos.entry[1],
                          object->slist.slist[i]->gpos.entry[2], singsize/zoom_factor);
    }
}


void    VectorFieldWindow::display_single_traj(Trajectory *traj, int sep_id)
{
    icVector3 outn;

    if(traj == NULL)
        return;

    /***********************************************************************************/
    /*   Display the trajectory as tube-like    */
    GLUquadricObj *quadratic = NULL;				// Storage For Our Quadratic Objects
    if (quadratic == NULL)
        quadratic=gluNewQuadric();			// Create A Pointer To The Quadric Object ( NEW )
    gluQuadricNormals(quadratic, GLU_SMOOTH);	// Create Smooth Normals ( NEW )
    gluQuadricTexture(quadratic, GL_TRUE);
    ////set the surface property
    GLfloat ambient[] = { 0.8, 0.8, 0.8, 1.0 };
    GLfloat diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat specular[] = { 1.0, 1.0, 1.0, 1.0 };
    int shiny = 100;

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT, GL_SHININESS, (GLfloat)shiny);

    glEnable(GL_COLOR_MATERIAL);
    /////////////////////////////////////////////////////////////////////////////////////

    glDepthFunc(GL_LEQUAL);

    if(sep_id == 1 || sep_id == 3)
        glColor3f(1, 0, 0);

    else if(sep_id == 2 || sep_id == 4)
        glColor3f(0, 1, 0);

    else
        glColor3f(0, 0, 0);

    for(int i = 0; i < traj->nlinesegs; i++)
    {
        if(traj->linesegs == NULL )//|| traj->linesegs[i] == NULL)
            return;

        if (FlipNormalOn)
            outn = -object->tlist.tris[traj->linesegs[i].Triangle_ID]->normal;
        else
            outn = object->tlist.tris[traj->linesegs[i].Triangle_ID]->normal;

        glBegin(GL_LINES);
        glVertex3f(traj->linesegs[i].gstart.entry[0]+ 0.002* outn.entry[0],\
                                                                                 traj->linesegs[i].gstart.entry[1]+ 0.002* outn.entry[1],\
                                                                                                             traj->linesegs[i].gstart.entry[2]+ 0.002* outn.entry[2]);
        glVertex3f(traj->linesegs[i].gend.entry[0]+ 0.002* outn.entry[0],\
                                                                               traj->linesegs[i].gend.entry[1]+ 0.002* outn.entry[1],\
                                                                                                           traj->linesegs[i].gend.entry[2]+ 0.002* outn.entry[2]);
        glEnd();

        /*
                If we display the separatrices as tube-like
            */
        if (EnTubeLikeVis)
        {
            glEnable(GL_LIGHTING);
            glEnable(GL_LIGHT0);
            glEnable(GL_LIGHT1);
            icVector3 Z (traj->linesegs[i].gend.entry[0]-traj->linesegs[i].gstart.entry[0],
                        traj->linesegs[i].gend.entry[1]-traj->linesegs[i].gstart.entry[1],
                        traj->linesegs[i].gend.entry[2]-traj->linesegs[i].gstart.entry[2]);
            double len = length (Z);

            normalize (Z);

            //That will give you the Z axis of the rotation matrix you need to set.
            //Then, make another vector which is

            icVector3 Y (Z.entry[1], -Z.entry[2], 0);
            normalize (Y);

            icVector3 X = cross (Y, Z);
            normalize (X);

            icVector3 Y2 = cross (X, Z);
            normalize (Y2);

            float rotMat[4][4];

            matrix_ident (rotMat);

            rotMat[0][0] = X.entry[0];
            rotMat[0][1] = X.entry[1];
            rotMat[0][2] = X.entry[2];

            rotMat[1][0] = Y2.entry[0];
            rotMat[1][1] = Y2.entry[1];
            rotMat[1][2] = Y2.entry[2];

            rotMat[2][0] = Z.entry[0];
            rotMat[2][1] = Z.entry[1];
            rotMat[2][2] = Z.entry[2];

            glMatrixMode(GL_MODELVIEW);
            glPushMatrix();

            glTranslatef (traj->linesegs[i].gstart.entry[0],
                         traj->linesegs[i].gstart.entry[1],
                         traj->linesegs[i].gstart.entry[2]);

            glMultMatrixf ((GLfloat *) rotMat);

            gluCylinder (quadratic, 0.005*object->radius/zoom_factor, 0.005*object->radius/zoom_factor, len, 10, 10);
            glMatrixMode(GL_MODELVIEW);
            glPopMatrix();
        }
    }

    glDepthFunc(GL_LESS);
}


void VectorFieldWindow::display_tau_distribution()
{
    int i, j;
    Triangle *face;
    Vertex *v;
    icVector3 outn;

    float hsv[3], rgb[3];

    double tau_diff = morse_decomp->tau_max - morse_decomp->tau_min;

    for(i=0; i<object->tlist.ntris; i++)
    {
        face = object->tlist.tris[i];

        if (tau_diff == 0)
        {
            glColor4f (0, 0, 0.5, 0.4);
        }

        else
        {
            hsv[0] = 256.-256.*(face->used_tau-morse_decomp->tau_min)/tau_diff;
            hsv[1] = 0.5;
            hsv[2] = 1;
            HsvRgb(hsv, rgb);
            glColor4f (rgb[0], rgb[1], rgb[2], 0.4);
        }

        glBegin(GL_TRIANGLES);
        for(j = 0; j < face->nverts; j++)
        {
            v = face->verts[j];
            outn = v->normal;
            glVertex3d(v->x + 0.001*outn.entry[0],
                       v->y + 0.001*outn.entry[1],
                       v->z + 0.001*outn.entry[2]);
        }
        glEnd();
    }
}


void VectorFieldWindow::display_diff_ECG_MCG()
{
    int i, j;
    Triangle *face;
    Vertex *v;
    icVector3 outn;

    double tau_diff = morse_decomp->tau_max - morse_decomp->tau_min;

    for(i=0; i<object->tlist.ntris; i++)
    {
        face = object->tlist.tris[i];

        //if (face->diff_ECG_MCG)
        if (face->exclude)
        {
            glColor4f(0.8, 0.8, 0.8, 0.3);
        }
        else
            glColor4f(0.2,0.2,0.2,0.5);

        glBegin(GL_TRIANGLES);
        for(j = 0; j < face->nverts; j++)
        {
            v = face->verts[j];
            outn = v->normal;
            glVertex3d(v->x + 0.001*outn.entry[0],
                       v->y + 0.001*outn.entry[1],
                       v->z + 0.001*outn.entry[2]);
        }
        glEnd();
    }
}


void VectorFieldWindow::display_diff_MCGs()
{
    int i, j;
    Triangle *face;
    Vertex *v;
    icVector3 outn;

    double tau_diff = morse_decomp->tau_max - morse_decomp->tau_min;

    for(i=0; i<object->tlist.ntris; i++)
    {
        face = object->tlist.tris[i];

        if (face->counter%2==1)
        {
            glColor4f(0.8, 0.0, 0.2, 0.4);
        }
        else
            glColor4f(0.2,0.2,0.2,0.0);

        glBegin(GL_TRIANGLES);
        for(j = 0; j < face->nverts; j++)
        {
            v = face->verts[j];
            outn = v->normal;
            glVertex3d(v->x + 0.001*outn.entry[0],
                       v->y + 0.001*outn.entry[1],
                       v->z + 0.001*outn.entry[2]);
        }
        glEnd();
    }
}
