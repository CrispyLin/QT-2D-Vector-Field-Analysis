// GlView.cpp : implementation file
//Created and Modified by Guoning Chen
//

#include <math.h>

#include "GlView.h"
#include "StreamlineCalculate/EvenStreamlines.h"
#include "Graph2D.h"

extern std::vector<Contour_Graph> isolines;

/*------------------------------------------------------------*/
//some constants for ibfv, some may be useless in future!!!06/23/05
//#define M_PI 3.1415926      ////Move to the datastructure.h file 06/23/05
#define	NPN 64
#define NMESH  100
#define DM  ((double) (1.0/(NMESH-1.0)))
#define NPIX  /*512*/800
double SCALE = 3.0;
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
//some global variables for ibfv
int     iframe = 0;
int     Npat   = 64;        // Modified at 02/10/2010
int     alpha  = (0.06*255);  // modified for a smear out LIC texture for better visualization (TVCG revision 11/22/2010)
double  sa;
double  tmax   = NPIX/(SCALE*NPN);
//double  dmax   = 0.9*SCALE/NPIX;
double  dmax   = 1.8/NPIX;
int npn = NPN;
int npix = /*512*/800;


/*------------------------------------------------------------*/

//get
GLuint Textures[2];
GLubyte f_tex[NPIX][NPIX][3], b_tex[NPIX][NPIX][3], applied_tex[NPIX][NPIX][3];

HDC  g_hDC;		        // GDI Device Context
HGLRC	*g_hglRC;		// Rendering Context

struct jitter_struct{
    double x;
    double y;
} jitter_para;

jitter_struct ji1[1] = {{0.0, 0.0}};
jitter_struct ji16[16] = {{0.125, 0.125}, {0.375, 0.125}, {0.625, 0.125},
                          {0.875, 0.125},{0.125, 0.375}, {0.375, 0.375}, {0.625, 0.375}, {0.875, 0.375},
                          {0.125, 0.625}, {0.375, 0.625}, {0.625, 0.625}, {0.875, 0.625},{0.125, 0.875},
                          {0.375, 0.875}, {0.625, 0.875}, {0.875, 0.875}};

//int MoveOrStop = 0;

extern Polyhedron *object;
extern MorseDecomp *morse_decomp;
extern MorseDecomp *local_decomp;
extern TrajectoryList *separatrices;
extern PeriodicOrbitList *periodic_orbits;
extern	EvenStreamlinePlace *evenplace;
extern MCG_Graph *mcg;
extern int picked_node;


double singsize = 0.009; // For visualizing the fixed points with the Morse sets together!
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




// CGlView

CGlView::CGlView()
{
}

CGlView::~CGlView()
{
}



// CGlView message handlers
int CGlView::OnCreate()
{

    m_hDC = ::GetDC(this->m_hWnd);

    if(!SetPixelformat(m_hDC))
    {
        ::MessageBox(::GetFocus(),"SetPixelformat Failed!","Error",MB_OK);
        return -1;
    }

    m_hglRC = wglCreateContext(m_hDC);

    g_hglRC = &m_hglRC;

    int i= wglMakeCurrent(m_hDC,m_hglRC);

    InitGL();

    return 0;
}


void  CGlView::get_rotation(double p[3], icVector3 rot_axis, double rotang, double rotmatrix[16])
{
    wglMakeCurrent(m_hDC, m_hglRC);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glTranslatef(p[0], p[1], p[2]);
    glRotatef(rotang, rot_axis.entry[0], rot_axis.entry[1], rot_axis.entry[2]);
    glTranslatef(-p[0], -p[1], -p[2]);

    glGetDoublev(GL_MODELVIEW_MATRIX, (double *)rotmatrix);
    glPopMatrix();
}


void  CGlView::get_rotation(double rotang, icVector3 rot_axis, double rotmatrix[16])
{
    wglMakeCurrent(m_hDC, m_hglRC);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glRotatef(rotang, rot_axis.entry[0], rot_axis.entry[1], rot_axis.entry[2]);
    glGetDoublev(GL_MODELVIEW_MATRIX, (double *)rotmatrix);
    glPopMatrix();
}


void CGlView::OnDestroy()
{
    CWnd::OnDestroy();

    // TODO: Add your message handler code here
    wglMakeCurrent(NULL,NULL);
    wglDeleteContext(m_hglRC);
}

BOOL CGlView::OnEraseBkgnd(CDC* pDC)
{
    // TODO: Add your message handler code here and/or call default

    return CWnd::OnEraseBkgnd(pDC);
}

/*--------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------*/
// Other extended operations for OpenGl window setting

BOOL CGlView::SetPixelformat(HDC hdc)
{

    PIXELFORMATDESCRIPTOR *ppfd;
    int pixelformat;

    PIXELFORMATDESCRIPTOR pfd = {
        sizeof(PIXELFORMATDESCRIPTOR),  //  size of this pfd
        1,                     // version number
        PFD_DRAW_TO_WINDOW |   // support window
            PFD_SUPPORT_OPENGL |   // support OpenGL
            PFD_GENERIC_FORMAT |
            PFD_DOUBLEBUFFER,      // double buffered
        PFD_TYPE_RGBA,         // RGBA type
        32,                    // 24-bit color depth
        0, 0, 0, 0, 0, 0,      // color bits ignored
        8,                     // no alpha buffer
        0,                     // shift bit ignored
        8,                     // no accumulation buffer
        0, 0, 0, 0,            // accum bits ignored
        64,                    // 32-bit z-buffer
        8,                     // no stencil buffer
        8,                     // no auxiliary buffer
        PFD_MAIN_PLANE,        // main layer
        0,                     // reserved
        0, 0, 0                // layer masks ignored
    };


    ppfd = &pfd;


    if ( (pixelformat = ChoosePixelFormat(hdc, ppfd)) == 0 )
    {
        ::MessageBox(NULL, "ChoosePixelFormat failed", "Error", MB_OK);
        return FALSE;
    }

    if (SetPixelFormat(hdc, pixelformat, ppfd) == FALSE)
    {
        ::MessageBox(NULL, "SetPixelFormat failed", "Error", MB_OK);
        return FALSE;
    }

    return TRUE;

}


int CGlView::InitGL(GLvoid)								// All Setup For OpenGL Goes Here
{
    //glViewport(0, 0, (GLsizei) NPIX, (GLsizei) NPIX);
    //glMatrixMode(GL_PROJECTION);
    //glLoadIdentity();
    //glTranslatef(-1.0, -1.0, 0.0);
    //glScalef(2.0, 2.0, 1.0);
    //glTexParameteri(GL_TEXTURE_2D,
    //				GL_TEXTURE_WRAP_S, GL_REPEAT);
    //glTexParameteri(GL_TEXTURE_2D,
    //				GL_TEXTURE_WRAP_T, GL_REPEAT);
    //glTexParameteri(GL_TEXTURE_2D,
    //				GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    //glTexParameteri(GL_TEXTURE_2D,
    //				GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    //glTexEnvf(GL_TEXTURE_ENV,
    //				GL_TEXTURE_ENV_MODE, GL_REPLACE);
    //glEnable(GL_TEXTURE_2D);
    //glShadeModel(GL_FLAT);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //glClear(GL_COLOR_BUFFER_BIT);

    glViewport(0, 0, (GLsizei) NPIX, (GLsizei) NPIX);


    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    //glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    //if (NotShiny)
    //	set_light();

    return TRUE;										// Initialization Went OK
}
/*--------------------------------------------------------------------------------*/



void CGlView::init_flags()
{
    MoveOrStop = 0;
    ShowFixedPtOn = 0;
    ShowSeparatricesOn = 0;
    ShowSCCsOn = 0;
    ShowPeriodicOrbitsOn = 0;
    EvenStreamlinePlacement = 0;
    ShowColorVFMagOn = 0;

    zoom_factor = 1;


    ShowMorseConn=0;
    ShowTriangleConn=0;
    morseC=0;
    triC=0;
    sccC=0;

    ShowEdgeSamplesOn=false;
    ShowTriMappingOn = false;
    selected_triangle = -1;  // 02/10/2010

    ShowBackward=0;//0=forward, 1=backward
    sampling_edge=0;

    ShowConnectionRegion=0;

    ShowBoundary=1;

    FlipNormalOn = false;
    IBFVOff = false;
}

void
display_sel_tri(int tri)
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




int CGlView::DrawGLScene(GLenum mode)					// Here's Where We Do All The Drawing
{

    //if (!NotShiny)

    if (!IBFVOff) {
        printf("IBFVOFF!");
        IBFVSEffect(mode);
    }
    else
    {
        printf("IBFVON!");
        glClearColor (1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        draw_shadedObj(mode);
    }
    //else
    //	IBFVSEffect_noshiny(mode);


    ////The following display other attribute on the surfaces
    ////such as streamlines, critical points ...
    //glDisable(GL_TEXTURE_2D);

    //   glDisable(GL_LIGHTING);


    glDisable(GL_TEXTURE_2D);

    //glEnable(GL_COLOR_MATERIAL);
    //   glEnable(GL_LINE_SMOOTH);
    //   glEnable(GL_BLEND);
    //   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //   glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);

    //	glEnable(GL_LIGHTING);
    //glEnable(GL_LIGHT0);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    glEnable(GL_COLOR_MATERIAL);

    glDepthFunc(GL_LEQUAL);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);

    /*display the extracted Morse sets with blending mode*/
    glDisable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_BLEND);

    /*
            This calling has been commented out by Guoning 03/01/2010
            What is the purpose of this calling?
        */
    ////////if(ShowMorseConn==1)
    ////////{
    ////////	display_morse_connection(mode);
    ////////}


    ////////if(ShowTriangleConn==1)
    ////////{
    ////////	display_triangle_connection(mode);
    ////////}

    //if(ShowConnectionRegion==1)
    //{
    //	display_MCG_connections();
    //}

    //if(ShowSCCsOn == 1)
    //{
    //	display_SCCs(mode);
    //}

    //

    //if(ShowColorVFMagOn == 1)
    //{
    //	display_color_VFMag();
    //}
    //	glDisable(GL_BLEND);

    ////////we may first mark those triangles where we add elements at
    ////if(DrawTrajectoryOn == 1)   ////Draw a single trajectory or whole??
    ////    display_single_traj();
    //
    ////glEnable(GL_POLYGON_OFFSET_FILL);
    ////glPolygonOffset (1., 1.);

    //////show separatrices
    //if(ShowSeparatricesOn == 1)
    //	display_separatrices(mode);
    //

    ////////Display all the detected limit cycles
    //if(ShowPeriodicOrbitsOn == 1)
    //	display_periodicorbits(mode);
    ////	DisplayLimitCycleLegends(mode);  //display periodic orbit handler to edit it

    ////glDisable(GL_BLEND);

    ////glEnable(GL_LIGHTING);

    ////glDisable(GL_COLOR_MATERIAL);
    //if(EvenStreamlinePlacement == 1)
    //{
    //	display_even_streamlines(mode); /*display the calculated evenly placed streamlines*/
    //}

    //glEnable(GL_LIGHTING);
    //glEnable(GL_LIGHT1);
    //glEnable(GL_LIGHT2);

    ///* Display the legends for the objects in the vector field */
    //if(ShowFixedPtOn == 1)
    //	display_FixedPtsIcon(mode);

    // show the rotational sum

    //glEnable(GL_BLEND);
    if (ShowRotSumOn)
    {
        vis_rot_sum();
    }
    //glDisable(GL_BLEND);

    without_antialiasing(mode);
    //with_antialiasing(mode);


    ///////*************************************************************
    ///////*
    //////   Display the stored sample points along edges 02/09/2010
    //////*/
    //////extern EdgeSamplePt_List *edge_samples;
    //////if (ShowEdgeSamplesOn)  // you can use different name of the variable
    //////{
    ////// if (edge_samples != NULL)
    //////  edge_samples->display();
    //////}




    /*************************************************************
    /*
       Display the stored sample points along edges 02/09/2010
    */
    glDisable(GL_LIGHTING);
    extern EdgeSamplePt_List *edge_samples;
    if (ShowEdgeSamplesOn)
    {
        if (edge_samples != NULL)
            //edge_samples->display(); display_sel_edges
            edge_samples->display_sel_edges(selected_triangle,ShowBackward,sampling_edge);
    }

    if (ShowTriMappingOn)
    {
        morse_decomp->show_tri_mapping(selected_triangle);
    }

    /**************************************************************
    /*
        Test the triangle selection 02/10/2010
    */
    if (selected_triangle >=0 && selected_triangle < object->tlist.ntris)
    {
        display_sel_tri(selected_triangle);
    }

    glDepthFunc(GL_LESS);




    glColor3f (1, 1, 1);


    glDisable(GL_COLOR_MATERIAL);


    //////draw_shadedObj(mode);


    SwapBuffers(m_hDC);
    return TRUE;										// Keep Going
}

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
////create pattern for texture mapping to create flow effect visualization
////the size of the pattern is 64x64
void CGlView::makePatterns(void)
{
    int lut[256];
    int phase[NPN][NPN];
    GLubyte pat[NPN][NPN][4];
    GLubyte spat[NPN][NPN][4];
    int i, j, k, t;
    int lut_r[256], lut_g[256], lut_b[256];

    //for (i = 0; i < 256; i++) lut[i] = i < 127 ? 0 : 255;
    for (i = 0; i < 256; i++) lut[i] = i < 127 ? 50 : 255;                 // for a brighter texture (TVCG revision 11/22/2010)
    //for (i = 0; i < 256; i++) lut[i] = 255.0 * pow(belta, (double)i/255); // exponential decay noise
    //for (i = 0; i < 256; i++) lut[i] = 255-i;                             // saw-tooth noise
    //for (i = 0; i < 256; i++) lut[i] = 255.0*cos(2*M_PI*i/255.0);         // cosine noise

    for (i = 0; i < NPN; i++)
        for (j = 0; j < NPN; j++) phase[i][j] = rand() % 256;

    if ( EnGreyTexture)
    {
        for (k = 0; k < Npat; k++) {
            t = k*256/Npat;                           //t is used to control the animation of the image
            for (i = 0; i < NPN; i++)
                for (j = 0; j < NPN; j++) {
                    pat[i][j][0] =
                        pat[i][j][1] =
                        pat[i][j][2] = lut[(t + phase[i][j]) % 255];
                    pat[i][j][3] = alpha;

                    spat[i][j][0] =
                        spat[i][j][1] =
                        spat[i][j][2] = lut[ phase[i][j] % 255];
                    spat[i][j][3] = alpha;
                }
            glNewList(k + 1, GL_COMPILE);
            glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0,
                         GL_RGBA, GL_UNSIGNED_BYTE, pat);
            glEndList();

            glNewList(k + 1 + 100, GL_COMPILE);       //This is for static image
            glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0,
                         GL_RGBA, GL_UNSIGNED_BYTE, spat);
            glEndList();
        }
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

            //lut_r[i] = 1.22*i<127?70:255; // for a brighter texture (TVCG revision 11/22/2010)
            // lut_g[i] = 1.3*i<127?70:255;
            //lut_b[i] = 1.1*i<127?70:255;
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
            glNewList(k + 1, GL_COMPILE);
            glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0,
                         GL_RGBA, GL_UNSIGNED_BYTE, pat);
            glEndList();

            glNewList(k + 1 + 100, GL_COMPILE);       //This is for static image
            glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0,
                         GL_RGBA, GL_UNSIGNED_BYTE, spat);
            glEndList();
        }
    }
}

void CGlView::getDP(double x, double y, double *px, double *py)
{
    double dx, dy, vx, vy, r;

    dx = x - 0.5;
    dy = y - 0.5;
    r  = dx*dx + dy*dy;
    if (r < 0.0001) r = 0.0001;
    vx = sa*dx/r + 0.02;
    vy = sa*dy/r;
    r  = vx*vx + vy*vy;
    if (r > dmax*dmax) {
        r  = sqrt(r);
        vx *= dmax/r;
        vy *= dmax/r;
    }
    *px = x + vx;
    *py = y + vy;
}

void CGlView::IBFVSEffect_noshiny(GLenum mode)
{

    int i;
    unsigned char j;
    Triangle *face;
    Vertex *v;

    GLfloat ambient[] = { 0.5, 0.5, 0.5, 1.0 };
    //GLfloat diffuse[] = { 1., 1.0, 1., 1.0 };
    GLfloat diffuse[] = { 1.0, 0.6, 0.0, 1.0 };
    GLfloat specular[] = { 1.0, 1.0, 1.0, 1.0 };
    int shiny = 100;

    if (EnGreyTexture)
    {
        diffuse[0] = diffuse[1] = diffuse[2] = 1.;
    }
    else
    {
        diffuse[0] = 1.;
        diffuse[1] = 0.6;
        diffuse[2] = 0.;
    }

    glClearColor (1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glShadeModel(GL_FLAT);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);
    glDisable(GL_LIGHT1);
    glDisable(GL_LIGHT2);
    glEnable(GL_DEPTH_TEST);

    glEnable(GL_TEXTURE_2D);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

    ////////////////////////////////////////////////
    glDrawBuffer(GL_BACK);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
                 GL_RGB, GL_UNSIGNED_BYTE, f_tex);

    set_view(GL_RENDER);

    glPushMatrix ();
    set_scene(GL_RENDER);

    ///* Calculate forward texture */
    //glDrawBuffer(GL_BACK);


    //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 512, 512, 0,
    //	GL_RGB, GL_UNSIGNED_BYTE, f_tex);

    //////Texture distortion and mapping here
    for (i=0; i<object->tlist.ntris; i++) {
        face = object->tlist.tris[i];

        glBegin(GL_POLYGON);
        for (j=0; j<face->nverts; j++) {
            v = face->verts[j];
            ///Using the storaged texture coordinates
            glTexCoord2f(v->texture_coord.entry[0], v->texture_coord.entry[1]);
                //glTexCoord2f(v->back_Tex.entry[0], v->back_Tex.entry[1]);
            glVertex3f(v->x, v->y, v->z);
        }
        glEnd();
    }

    glPopMatrix();

    iframe = iframe + 1;

    ////Adding the depth judgement here!
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthFunc(GL_GREATER);
    if( MoveOrStop ==0 )   ////Static
        glCallList(iframe % Npat + 1 + 100);
    else                   ////Moving
        glCallList(iframe % Npat + 1);
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
    glReadBuffer(GL_BACK);
    glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, f_tex);
    ///////////////////////////////////////////////////////////////////

    /*    for backward frame  */
    if(MoveOrStop == 0)
    {
        glClearColor (1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        /* Calculate backward texture */
        glDrawBuffer(GL_BACK);
        glEnable(GL_TEXTURE_2D);
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
                     GL_RGB, GL_UNSIGNED_BYTE, b_tex);

        glPushMatrix ();
        set_scene(GL_RENDER);


        //////Texture distortion and mapping here
        for (i=0; i<object->tlist.ntris; i++) {
            face = object->tlist.tris[i];
            glBegin(GL_POLYGON);
            for (j=0; j<face->nverts; j++) {
                v = face->verts[j];
                ///Using the storaged texture coordinates
                glTexCoord2f(v->back_Tex.entry[0], v->back_Tex.entry[1]);
                //glTexCoord2f(v->texture_coord.entry[0], v->texture_coord.entry[1]);
                glVertex3f(v->x, v->y, v->z);
            }
            glEnd();
        }

        glPopMatrix();

        //iframe = iframe - 1;

        ////Adding the depth judgement here!
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDepthFunc(GL_GREATER);
        glCallList(iframe % Npat + 1 + 100);
        glBegin(GL_QUAD_STRIP);
        glTexCoord2f(0.0,  0.0);  glVertex3f(0.0, 0.0, -39.9);
        glTexCoord2f(0.0,  tmax); glVertex3f(0.0, 1.0, -39.9);
        glTexCoord2f(tmax, 0.0);  glVertex3f(1.0, 0.0, -39.9);
        glTexCoord2f(tmax, tmax); glVertex3f(1.0, 1.0, -39.9);
        glEnd();
        glDepthFunc(GL_LESS);
        //glDisable(GL_BLEND);

        glReadBuffer(GL_BACK);
        glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, b_tex);

        //blend two images

        for(int x = 0; x < NPIX; x++)
        {
            for(int y = 0; y < NPIX; y++)
            {
                /**/
                int temp_color = (int)(f_tex[x][y][0] + b_tex[x][y][0])/2.;
                applied_tex[x][y][0] = (int)(f_tex[x][y][0] + b_tex[x][y][0])/2;
                applied_tex[x][y][1] = (int)(f_tex[x][y][1] + b_tex[x][y][1])/2;
                applied_tex[x][y][2] = (int)(f_tex[x][y][2] + b_tex[x][y][2])/2;
            }
        }
    }

    ////////////////////////////////////////////////////////////////////
    ////Blending shading here 12/1

    set_scene(GL_RENDER);

    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    //glEnable(GL_LIGHT1);
    //glEnable(GL_LIGHT2);

    //glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
    //glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
    //glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
    //glMaterialf(GL_FRONT, GL_SHININESS, (GLfloat)shiny);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT_AND_BACK/*GL_FRONT*/, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT_AND_BACK/*GL_FRONT*/, GL_SPECULAR, specular);
    //glMaterialf(GL_FRONT, GL_SHININESS, (GLfloat)shiny);
    glMateriali(GL_FRONT_AND_BACK/*GL_FRONT*/, GL_SHININESS, 80);

    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

    if (MoveOrStop == 0)
    {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
                     GL_RGB, GL_UNSIGNED_BYTE, applied_tex);
    }
    else
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
                     GL_RGB, GL_UNSIGNED_BYTE, f_tex);

    for (i=0; i<object->tlist.ntris; i++) {
        if (mode == GL_SELECT)
            glLoadName(i+1);

        face = object->tlist.tris[i];

        glBegin(GL_POLYGON);
        for (j=0; j<face->nverts; j++) {
            v = face->verts[j];
            glNormal3dv(v->normal.entry);
            glTexCoord2f(v->texture_coord.entry[0], v->texture_coord.entry[1]);
            glVertex3f(v->x, v->y, v->z);
        }
        glEnd();

    }

    glDisable(GL_TEXTURE_2D);

}

//
void CGlView::IBFVSEffect(GLenum mode)
{

    int i, j;
    Triangle *face;
    Vertex *v;

    GLfloat ambient[] = { 0.3, 0.3, 0.3, 1.0 };
    //GLfloat diffuse[] = { 1., 0.6, 0, 1.0};
    GLfloat diffuse[] = { 0.8, .8, 1., 1.0 };
    GLfloat specular[] = { 0.8, 0.8, 1.0, 1.0 };
    int shiny = 100;

    glClearColor (1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glShadeModel(GL_FLAT);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);

    glEnable(GL_TEXTURE_2D);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

    ////////////////////////////////////////////////
    glDrawBuffer(GL_BACK);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
                 GL_RGB, GL_UNSIGNED_BYTE, f_tex);

    set_view(GL_RENDER);

    glPushMatrix ();
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

    iframe = iframe + 1;

    ////Adding the depth judgement here!
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthFunc(GL_GREATER);
    if( MoveOrStop ==0 )   ////Static
        glCallList(iframe % Npat + 1 + 100);
    else                   ////Moving
        glCallList(iframe % Npat + 1);
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
    glReadBuffer(GL_BACK);
    glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, f_tex);

    if(MoveOrStop == 0)
    {
        /* Calculate backward texture */
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
                     GL_RGB, GL_UNSIGNED_BYTE, b_tex);

        glPushMatrix ();
        set_scene(GL_RENDER);


        //////Texture distortion and mapping here
        for (i=0; i<object->tlist.ntris; i++) {
            face = object->tlist.tris[i];
            glBegin(GL_POLYGON);
            for (j=0; j<face->nverts; j++) {
                ///Using the storaged texture coordinates
                v = face->verts[j];
                glTexCoord2f(v->back_Tex.entry[0], v->back_Tex.entry[1]);
                glVertex3f(v->x, v->y, v->z);
            }
            glEnd();
        }

        glPopMatrix();

        iframe = iframe + 1;

        ////Adding the depth judgement here!
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDepthFunc(GL_GREATER);
        //if( MoveOrStop ==0 )   ////Static
        glCallList(iframe % Npat + 1 + 100);
        //else                   ////Moving
        //	glCallList(iframe % Npat + 1);
        glBegin(GL_QUAD_STRIP);
        glTexCoord2f(0.0,  0.0);  glVertex3f(0.0, 0.0, -39.9);
        glTexCoord2f(0.0,  tmax); glVertex3f(0.0, 1.0, -39.9);
        glTexCoord2f(tmax, 0.0);  glVertex3f(1.0, 0.0, -39.9);
        glTexCoord2f(tmax, tmax); glVertex3f(1.0, 1.0, -39.9);
        glEnd();
        glDepthFunc(GL_LESS);
        glDisable(GL_BLEND);

        glReadBuffer(GL_BACK);
        glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, b_tex);

        //blend two images

        for(int x = 0; x < NPIX; x++)
        {
            for(int y = 0; y < NPIX; y++)
            {
                /**/
                int temp_color = (int)(f_tex[x][y][0] + b_tex[x][y][0])/2.;
                applied_tex[x][y][0] = (int)(f_tex[x][y][0] + b_tex[x][y][0])/2;
                applied_tex[x][y][1] = (int)(f_tex[x][y][1] + b_tex[x][y][1])/2;
                applied_tex[x][y][2] = (int)(f_tex[x][y][2] + b_tex[x][y][2])/2;
            }
        }
    }


    ///////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    ////Blending shading here 12/1

    /*  First, draw the background  */
    if (!DisableLighting && !DisableBackground)
    {
        glClearColor (0.0, 0.0, 0.0, 1.0);  // background for rendering color coding and lighting
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);
        glDisable(GL_TEXTURE_2D);
        glShadeModel(GL_SMOOTH);

        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        gluOrtho2D(-1., 1., -1., 1.);

        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();

        glBegin(GL_QUADS);
        //red color
        //glColor3f(0.4,0.4,0.3);
        glColor3f(0.1,0.1,0.3);
        glVertex2f(-1.0,-1.0);
        glVertex2f(1.0,-1.0);
        //blue color
        glColor3f(1.,1.,1.);
        glVertex2f(1.0, 1.0);
        glVertex2f(-1.0, 1.0);
        glEnd();

        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();
    }
    ////////////////////////////////////////////

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
    else
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
                     GL_RGB, GL_UNSIGNED_BYTE, f_tex);


    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT_AND_BACK/*GL_FRONT*/, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT_AND_BACK/*GL_FRONT*/, GL_SPECULAR, specular);
    //glMaterialf(GL_FRONT, GL_SHININESS, (GLfloat)shiny);
    glMateriali(GL_FRONT_AND_BACK/*GL_FRONT*/, GL_SHININESS, 80);

    for (i=0; i<object->tlist.ntris; i++) {
        if (mode == GL_SELECT)
            glLoadName(i+1);

        face = object->tlist.tris[i];

        /*
            Do not display the texture if it is part of the in/outlets of the data
        */
        if (face->verts[0]->in_out_let || face->verts[1]->in_out_let || face->verts[2]->in_out_let)
            continue;

        glBegin(GL_POLYGON);
        for (j=0; j<face->nverts; j++) {
            v = face->verts[j];
            glNormal3d(v->normal.entry[0],
                       v->normal.entry[1],
                       v->normal.entry[2]);
            glTexCoord2f(v->texture_coord.entry[0], v->texture_coord.entry[1]);
            glVertex3f(v->x, v->y, v->z);
        }
        glEnd();

    }

    glDisable(GL_TEXTURE_2D);

}


////View port setting routines
void CGlView::set_view(GLenum mode)
{
    ////////icVector3 up, ray, view;
    ////////GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
    ////////GLfloat light_diffuse0[] = { 0.5, 0.5, 0.5, 1.0 };
    ////////GLfloat light_specular0[] = { 1.0, 0.2, 0.0, 1.0 };
    ////////GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
    ////////GLfloat light_diffuse1[] = { 0.6, 0.6, 0.6, 1.0 };
    ////////GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };
    ////////GLfloat light_ambient2[] = { 1.0, 1.0, 1.0, 1.0 };
    ////////GLfloat light_diffuse2[] = { 1.0, 1.0, 1.0, 1.0 };
    ////////GLfloat light_specular2[] = { 1.0, 1.0, 1.0, 1.0 };
    ////////GLfloat light_position[] = { 0.0, 0.0, 0.0, 1.0 };

    ////////glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
    ////////glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
    ////////glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);
    ////////glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
    ////////glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
    ////////glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


    //////   glMatrixMode(GL_PROJECTION);
    //////
    //////if (mode == GL_RENDER)
    //////	glLoadIdentity();

    //////glOrtho(0., 1., 0., 1., -400.0, 4000.0);


    //////glMatrixMode(GL_MODELVIEW);
    //////   glLoadIdentity();
    ////////light_position[0] = 5.5;
    ////////light_position[1] = 0.0;
    ////////light_position[2] = 4.0;
    ////////glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    //////////light_position[0] = -0.1;
    ////////light_position[0] = -1.1;
    ////////light_position[1] = 0.3;
    ////////light_position[2] =4.0;
    ////////glLightfv(GL_LIGHT2, GL_POSITION, light_position);

    /*   new set_view  03/01/2010  */

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

    // Light and material Data (added by Guoning 11/08/2009)
    //GLfloat fLightPos[4]   = { -1000.0f, 10000.0f, 50.0f, 1.0f };  // Point source
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

    //else
    //	set_light();

    ///////////////////////////////////////////////////////////////////////


}


////View port setting routines
void CGlView::set_light()
{
    icVector3 up, ray, view;
    GLfloat light_ambient0[] = { 0.5, 0.5, 0.5, 1.0 };
    GLfloat light_diffuse0[] = { 0.5, 0.5, 0.5, 1.0 };
    //GLfloat light_specular0[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_specular0[] = { 1.0, 0.5, 0., 1.0 };
    GLfloat light_ambient1[] = { 0.5, 0.5, 0.5, 1.0 };
    GLfloat light_diffuse1[] = { 0.6, 0.6, 0.6, 1.0 };
    //GLfloat light_specular1[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_specular1[] = { 1.0, 0.5, 0., 1.0 };
    GLfloat light_ambient2[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_diffuse2[] = { 1.0, 1.0, 1.0, 1.0 };
    //GLfloat light_specular2[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_specular2[] = { 1.0, 0.5, 0., 1.0 };
    GLfloat light_position[] = { 0.0, 3.0, -3.0, 1.0 };

    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);
    glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
    glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);
    glLightfv(GL_LIGHT2, GL_DIFFUSE, light_diffuse2);
    glLightfv(GL_LIGHT2, GL_SPECULAR, light_specular2);
    glLightfv(GL_LIGHT2, GL_SPECULAR, light_specular2);

    //glMatrixMode(GL_MODELVIEW);

    light_position[0] = 15.5;
    light_position[1] = 40.0;
    light_position[2] = 40.0;
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    light_position[0] = 15.1;
    light_position[1] = 40.3;
    light_position[2] =40.0;
    glLightfv(GL_LIGHT1, GL_POSITION, light_position);

    light_position[0] = -20.1;
    light_position[1] = 0.0;
    light_position[2] =70.0;
    glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}


void CGlView::set_scene(GLenum mode)
{
    glTranslatef(trans_x, trans_y, 0);
    glTranslatef(0.0, 0.0, -3.0);
    glTranslatef(rot_center.entry[0] ,
                 rot_center.entry[1],
                 rot_center.entry[2]);

    //glTranslatef(0.5 ,
    //	0.5,
    //	rot_center.entry[2]);

    multmatrix( rotmat );

    //glRotatef(180, 1, 0, 0); /*show the back of the object here 10/26/2007*/


    glScalef(zoom_factor, zoom_factor, zoom_factor);

    //glScalef(1.5, 1.5, 1.5);
    glScalef(0.7, 0.7, 0.7);
    //glScalef(0.51, 0.51, 0.51);

    glTranslatef(-rot_center.entry[0] ,
                 -rot_center.entry[1] ,
                 -rot_center.entry[2] );


}


/////For recalculate all the texture coordinates
void CGlView::ReCalTexcoord(void)
{
    double p[4], pr[4];

    double vx, vy, vz ;

    int i, j, k;

    wglMakeCurrent(m_hDC, m_hglRC);
    glGetFloatv(GL_MODELVIEW_MATRIX, ObjXmat);


    for(i = 0; i < object->vlist.nverts; i++)
    {
        ////Initialize the variables nv and pr
        pr[0] = pr[1] = pr[2] = pr[3] = 0;


        //vx = object->vlist.verts[i]->t_vec.entry[0] * dmax;
        //vy = object->vlist.verts[i]->t_vec.entry[1] * dmax;
        //vz = object->vlist.verts[i]->t_vec.entry[2] * dmax;

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

void CGlView::getInitTex(void)
{
    glDrawBuffer(GL_BACK);

    glCallList(0);

    //   glBegin(GL_QUAD_STRIP);
    //	glTexCoord2f(0.0,  0.0);  glVertex2f(0.0, 0.0);
    //	glTexCoord2f(0.0,  tmax); glVertex2f(0.0, 1.0);
    //	glTexCoord2f(tmax, 0.0);  glVertex2f(1.0, 0.0);
    //	glTexCoord2f(tmax, tmax); glVertex2f(1.0, 1.0);
    //glEnd();

    glBegin(GL_QUAD_STRIP);
    glTexCoord2f(0.0,  0.0);  glVertex3f(0.0, 0.0, -39.9);
    glTexCoord2f(0.0,  tmax); glVertex3f(0.0, 1.0, -39.9);
    glTexCoord2f(tmax, 0.0);  glVertex3f(1.0, 0.0, -39.9);
    glTexCoord2f(tmax, tmax); glVertex3f(1.0, 1.0, -39.9);
    glEnd();

    glReadBuffer(GL_BACK);
    glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, f_tex);

    glReadBuffer(GL_BACK);
    glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, b_tex);
}



/******************************************************************
Draw a solid sphere for marking the captured fixed points
******************************************************************/
void CGlView::display_FixedPtsIcon(GLenum mode)
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


////Select the color for the visual icons according to the type of the singularity
/*current color scheme  08/25/05
source: pure green (0, 1, 0)
repeller: light green (0, 1, 0.5)
sink:   pure red (1, 0, 0)
attractor:  orange (1, 0.5, 0)
saddle: pure blue (0, 0, 1)
center: light red (1, 0, 1)
*/

void CGlView::set_ColorByType(unsigned char type)
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

void CGlView::draw_ASolidSphere(double x, double y, double z, double r)
{
    glPushMatrix();
    glTranslatef(x, y, z);
    glutSolidSphere(r, 80, 80 );
    glPopMatrix();
}



/*
Draw all the separatrices in the field
*/
void CGlView::display_separatrices(GLenum mode)
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

    //glLineWidth(1.);
}

int ndisplay_trajs;
/*
Display the obtained evenly placed streamlines
*/
void CGlView::display_even_streamlines(GLenum mode)
{
    int i;

    for(i = 0; i < ndisplay_trajs/*evenplace->evenstreamlines->ntrajs*/; i++)
    {
        /*Testing code*/
        //FILE *fp = fopen("display_test", "w");
        //fprintf(fp, "display the %d th streamline\n", i);
        //fclose(fp);

        if(evenplace->evenstreamlines->trajs[i] != NULL)
            display_single_traj(evenplace->evenstreamlines->trajs[i], 5);
    }
}


void CGlView::draw_shadedObj(GLenum mode)
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
    //glEnable(GL_LIGHT1);
    //glEnable(GL_LIGHT2);
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

void matrix_ident( float m[4][4])
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

/******************************************************************
Draw specific trajectory according to its index
It seems that we should move the streamline along the outward normal direction
******************************************************************/

void CGlView::display_single_traj(Trajectory *traj, int sep_id)
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
    //glLineWidth(1.);

    //glBegin(GL_LINES);
    //for(int i = 0; i < traj->nlinesegs; i++)
    //{
    //	outn = object->tlist.tris[traj->linesegs[i].Triangle_ID]->normal;

    //	glVertex3f(traj->linesegs[i].gstart.entry[0]+ 0.002* outn.entry[0],\
    //		traj->linesegs[i].gstart.entry[1]+ 0.002* outn.entry[1],\
    //		traj->linesegs[i].gstart.entry[2]+ 0.002* outn.entry[2]);
    //}
    //glEnd();

    //}

    //else
    //{
    //	////set the surface property
    //	GLfloat ambient[] = { 0.2, 0.2, 0.2, 1.0 };
    //	GLfloat diffuse[] = { 0.2, 0.2, 1.0, 1.0 };

    //	if(MoveOrStop == 0)
    //	{
    //		if(sep_id == 1 || sep_id == 3)
    //		{
    //			diffuse[0] = 1.0;
    //			diffuse[1] = 0.0;
    //			diffuse[2] = 0.0;
    //		}
    //		else if(sep_id == 2 || sep_id == 4)
    //		{
    //			diffuse[0] = 0.0;
    //			diffuse[1] = 1.0;
    //			diffuse[2] = 0.0;
    //		}
    //	}


    //	GLfloat specular[] = { 1.0, 1.0, 1.0, 1.0 };
    //	int shiny = 100;

    //	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
    //	glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
    //	glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
    //	glMaterialf(GL_FRONT, GL_SHININESS, (GLfloat)shiny);

    //	glDisable(GL_COLOR_MATERIAL);

    //	icVector3 dir;
    //	for(int i = 0; i < traj->nlinesegs; i++)
    //	{
    //   dir.entry[0] = traj->linesegs[i].gend.entry[0]-traj->linesegs[i].gstart.entry[0];
    //dir.entry[1] = traj->linesegs[i].gend.entry[1]-traj->linesegs[i].gstart.entry[1];
    //dir.entry[2] = traj->linesegs[i].gend.entry[2]-traj->linesegs[i].gstart.entry[2];

    //		/* If we want to do animation, we have to the color mapping here 01/30/07 */
    //		if(MoveOrStop == 1)
    //		{
    //			if(sep_id == 1 || sep_id == 3)
    //			{
    //				if(color_pattern[(i-movecolor+length_color_pattern)%length_color_pattern])
    //				{
    //					diffuse[0] = 1.0;
    //					diffuse[1] = 0.0;
    //					diffuse[2] = 0.0;
    //					glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
    //				}
    //				else
    //				{
    //					diffuse[0] = 1.0;
    //					diffuse[1] = 1.0;
    //					diffuse[2] = 1.0;
    //					glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
    //				}
    //			}
    //
    //			else if(sep_id == 2 || sep_id == 4)
    //			{
    //				if(color_pattern[(i+movecolor)%length_color_pattern])
    //				{
    //					diffuse[0] = 0.0;
    //					diffuse[1] = 1.0;
    //					diffuse[2] = 0.0;
    //					glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
    //				}
    //				else
    //				{
    //					diffuse[0] = 0.9;
    //					diffuse[1] = 0.9;
    //					diffuse[2] = 0.9;
    //					glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
    //				}
    //			}
    //		}

    //		int num_leaves = 15;
    //		if((sep_id == 1 || sep_id == 3) && i>=numlinesegs-num_leaves)
    //		{
    //			dir.entry[0] = trajectories[index][numlinesegs-5].gend[0]
    //			-trajectories[index][i].gstart[0];
    //			dir.entry[1] = trajectories[index][numlinesegs-5].gend[1]
    //			-trajectories[index][i].gstart[1];
    //			dir.entry[2] = trajectories[index][numlinesegs-5].gend[2]
    //			-trajectories[index][i].gstart[2];

    //			DrawACylider(singsize+0.003, 0, length(dir), dir,
    //				trajectories[index][i].gstart);
    //			return;
    //		}

    //		else if((sep_id == 2 || sep_id == 4) && i==0)
    //		{
    //			if(i+10 < numlinesegs)
    //			{
    //				dir.entry[0] = trajectories[index][i+10].gend[0]
    //				-trajectories[index][i].gstart[0];
    //				dir.entry[1] = trajectories[index][i+10].gend[1]
    //				-trajectories[index][i].gstart[1];
    //				dir.entry[2] = trajectories[index][i+10].gend[2]
    //				-trajectories[index][i].gstart[2];
    //
    //				i+=10;

    //				DrawACylider(0, singsize+0.003, length(dir), dir,
    //					trajectories[index][i].gstart);
    //			}
    //			else
    //				DrawACylider(0, singsize+0.003, length(dir), dir,
    //					trajectories[index][i].gstart);
    //		}
    //		else
    //			DrawACylider(singsize-0.002, singsize-0.002, length(dir), dir,
    //				trajectories[index][i].gstart);

    //	}
    //}

}


void CGlView::get_Color(int num, float rgb[3])
{
    float hsv[3] = {0, 1., 1.};

    hsv[0] = 360 * (double)num/20.;

    HsvRgb(hsv,rgb);
}


void CGlView::get_Color_frac(int num, int frac, float rgb[3])
{
    float hsv[3] = {0, 1., 1.};

    hsv[0] = 360 * (double)num/frac;

    HsvRgb(hsv,rgb);
}


/*
Display the extracted strongly connected components
*/
void CGlView::display_SCCs(GLenum mode)
{
    int i, j, k;
    Triangle *face;
    Vertex *vert;
    icVector3 outn;
    int colorcount = 0;

    if(morse_decomp == NULL || morse_decomp->scclist == NULL)
        return;

    //if(

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

            //float rgb[3] = {0.};
            ////get_Color(i%20, rgb);
            //get_Color_frac(colorcount%8, 8, rgb);
            //colorcount++;
            //glColor4f(rgb[0], rgb[1], rgb[2], 0.4);

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

                //case 3: glColor4f (1, 0, 1, 0.4); // trivial index
                //	break;

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

                //	if(face->index==6338)glColor4f(0,0,0,0.6);
                //	else
                //		if(face->index==6397 || face->index==6935 ||face->index==7104 ||face->index==7103 ||face->index==6936||face->index==7105)
                //			glColor4f(1,1,0,0.6);
                //		else
                //			if(face->index==6335||face->index==6298||face->index==6375||face->index==6297||face->index==6376||face->index==6334||face->index==6336||face->index==6259)
                //				glColor4f(1,0,1,0.6);
                //			else
                //			{
                //					switch(mcg->nlist->mnodes[morse_decomp->scclist->scccomponents[i]->node_index]->type)
                //{
                //case 0:glColor4f(0, 1, 0, 0.4);//0-source_like, 1-sink_like, 2-saddle_like
                //	break;
                //case 1:glColor4f(1, 0, 0, 0.4);
                //	break;
                //case 2:glColor4f(0, 0, 1, 0.4);
                //	break;

                //}


                //			}


                ///// This is only for the hierarchy MCGs for cooling jacket data 07/14/2010

                //if (face->exclude)
                //	alpha_channel = 0.4;
                //else
                //	alpha_channel = 0.17;

                //glColor4f(color_rgb[0], color_rgb[1], color_rgb[2], alpha_channel);


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
    if(/*ShowMCGOn == 1 && */picked_node >= 0)
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


    /*Display the highlighted SCC*/
    //if(/*ShowMCGOn == 1 && */picked_node >= 0)
    //{
    //	int *tris = morse_decomp->scclist->scccomponents[mcg->nlist->mnodes[picked_node]->scc_index]->nodes;
    //	int ntris = morse_decomp->scclist->scccomponents[mcg->nlist->mnodes[picked_node]->scc_index]->nnodes;
    //	for(i = 0; i < ntris; i++)
    //	{
    //		glColor3f(1, 1, 1);

    //		if(tris[i] < 0)
    //			continue;

    //		face = object->tlist.tris[tris[i]];

    //		glBegin(GL_LINE_LOOP);
    //		for(j = 0; j < 3; j++)
    //		{
    //			//glVertex2f(Object.vlist[face->verts[j]]->x, Object.vlist[face->verts[j]]->y);
    //				vert = face->verts[j];
    //				outn = vert->normal;
    //				glVertex3d(vert->x + 0.005*outn.entry[0],
    //					vert->y + 0.005*outn.entry[1],
    //					vert->z + 0.005*outn.entry[2]);
    //		}
    //		glEnd();
    //	}
    //}







}


void CGlView::display_color_VFMag()
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


    /*
       Also display the vector field magnitude iso-contour
       For the earthquake data, we have to use the flipped normal for some reason
    */

    //for (i=0; i<isolines.size(); i++)
    //{
    //	Contour_Graph &cg = isolines[i];
    //	//cg.display();
    //	glColor3f (0, 0, 0);
    //	icVector3 norm;
    //	glBegin(GL_LINES);
    //	for (int j=0; j<cg.elist.size(); j++)
    //	{
    //		ContourEdge &e=cg.elist[j];
    //		ContourVertex &v1=cg.vlist[e.verts[0]];
    //		ContourVertex &v2=cg.vlist[e.verts[1]];
    //		norm = object->tlist.tris[e.at_triangle]->normal;

    //		icVector3 offset = v1.pos-0.002*norm;
    //		glVertex3dv(offset.entry);
    //		offset = v2.pos-0.002*norm;
    //		glVertex3dv(offset.entry);
    //	}
    //	glEnd();
    //}

}


void CGlView::display_sepandatt_pts(GLenum mode)
{
    int j, k;
    Triangle *face;
    Edge *cur_e;
    glPointSize(3.);
    glBegin(GL_POINTS);
    icVector3 outn;

    for(j = 0; j < object->tlist.ntris; j++)
    {
        face = object->tlist.tris[j];
        outn = face->normal;

        for(int k = 0; k < 3; k++)
        {
            cur_e = face->edges[k];
            glColor4f(0, 0, 0, 1);   //sep points

            if(cur_e->find_sep)
                glVertex3f(cur_e->sep->x+0.002*outn.entry[0],
                           cur_e->sep->y+0.002*outn.entry[1],
                           cur_e->sep->z+0.002*outn.entry[2]);

            glColor4f(0.8, 1, 0, 1); //att points
            if(cur_e->find_attp)
                glVertex3f(cur_e->attp->x+0.002*outn.entry[0],
                           cur_e->attp->y+0.002*outn.entry[1],
                           cur_e->attp->z+0.002*outn.entry[2]);
        }

    }

    glEnd();
    glPointSize(1.);
}

int ndisplay_POs;

/******************************************************************************
This routine is used to display the closed streamline of all the limit cycles
******************************************************************************/
void CGlView::display_periodicorbits(GLenum mode)
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
    //glLineWidth(1.);
}

/*-------------------------------------------------------------*/
//////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
void CGlView::multmatrix(const Matrix m)
{
    int i,j, index = 0;

    GLfloat mat[16];

    for ( i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            mat[index++] = m[i][j];

    glMultMatrixf (mat);
}

void CGlView::mat_ident(Matrix m)
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
void CGlView::display_morse_connection(GLenum mode)
{
    int sccc=mcg->nlist->mnodes[morseC]->scc_index;

    for(int i=0;i<morse_decomp->scclist->scccomponents[sccc]->nnodes;i++)
    {
        int node=morse_decomp->scclist->scccomponents[sccc]->nodes[i];
        triC=node;
        display_triangle_connection(mode);
    }
}

void CGlView::display_triangle_connection(GLenum mode)
{
    FILE* conn=fopen("connection.txt","w");
    int tri=triC;

    glBegin(GL_TRIANGLES);
    //draw the triangle
    glColor4f(0,0,0,1);
    for(int i=0;i<object->tlist.tris[tri]->nverts;i++)
    {
        Vertex *v=object->tlist.tris[tri]->verts[i];
        icVector3 outn=v->normal;
        glNormal3dv(v->normal.entry);
        glVertex3d(v->x+0.003*outn.entry[0],v->y+0.003*outn.entry[1],v->z+0.003*outn.entry[2]);
    }
    glEnd();


    fprintf(conn,"SCC=%d, The triangles connecting to triangle %d:\n",mcg->nlist->mnodes[morseC]->scc_index,tri);
    fprintf(conn,"from\tto\n");
    //draw other triangles connecting to the triangle by directed edges
    glBegin(GL_TRIANGLES);
    for(int i=0;i<morse_decomp->dg->nlist->dirnodes[tri]->nedges;i++)
    {
        int node;
        int e=morse_decomp->dg->nlist->dirnodes[tri]->edges[i];
        int n1=morse_decomp->dg->elist->edges[e]->node_index1;
        int n2=morse_decomp->dg->elist->edges[e]->node_index2;
        fprintf(conn,"%d\t%d\n",n1,n2);
        if(n1==tri)//out
        {
            glColor4f(1,1,0,0.4);//yellow
            node=n2;
        }
        else
        {
            glColor4f(1,0,1,0.4);//purple
            node=n1;
        }
        for(int j=0;j<object->tlist.tris[node]->nverts;j++)
        {
            Vertex *v=object->tlist.tris[node]->verts[j];
            icVector3 outn=v->normal;
            glNormal3dv(v->normal.entry);
            glVertex3d(v->x+0.002*outn.entry[0],v->y+0.002*outn.entry[1],v->z+0.002*outn.entry[2]);
        }
    }
    glEnd();
    //print out the triangles around the triangle
    fprintf(conn,"The triangles around %d:\n",tri);

    for(int i=0;i<3;i++)
    {
        Edge *e=object->tlist.tris[tri]->edges[i];
        if(e->tris[0]->index==tri)
            fprintf(conn,"%d\n",e->tris[1]->index);
        else
            fprintf(conn,"%d\n",e->tris[0]->index);
    }
    fclose(conn);
}

void display_Obj(GLenum mode)
{
    int i, j;
    Triangle *face;
    int *verts;

    for (i=0; i<object->tlist.ntris; i++) {
        if (mode == GL_SELECT)
            glLoadName(i+1);

        face = object->tlist.tris[i];

        glBegin(GL_POLYGON);
        for (j=0; j<face->nverts; j++) {
            glVertex3f(face->verts[j]->x, face->verts[j]->y, face->verts[j]->z);
        }
        glEnd();
    }

}


void CGlView::HitProcess(double ss, double st)
{
    GLuint  selectBuffer[128];
    GLint	vp[4] = {0, 0 ,1, 1};
    int hits;

    //glGetIntegerv(GL_VIEWPORT, vp);

    ////Build the selection buffer here
    glSelectBuffer(128, selectBuffer);

    //glMatrixMode(GL_MODELVIEW);
    //glPushMatrix();


    (void)glRenderMode(GL_SELECT);
    glInitNames();
    glPushName(0);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    //gluPickMatrix(ss, st, 0.005, 0.005, vp );  ////set smaller pick window for triangle selection
    gluPickMatrix(ss, st, 1.e-6, 1.e-6, vp );  ////set smaller pick window for triangle selection

    glMatrixMode( GL_MODELVIEW );

    glPushMatrix ();
    set_view(GL_SELECT);
    set_scene(GL_SELECT);

    display_Obj(GL_SELECT);
    glPopMatrix();


    ////If one of the element has been selected for being edited
    hits = glRenderMode(GL_RENDER);

    if(hits > 0)
    {
        selected_triangle = SelecteObj(hits, selectBuffer) - 1;

    }

    else
        selected_triangle = -1;

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();

    glMatrixMode( GL_MODELVIEW );
}


int CGlView::SelecteObj(int hits, GLuint buffer[])
{
    unsigned int i, j;
    GLuint names, *ptr;
    double smallest_depth=1.0e+30, current_depth;
    int seed_id=-1;
    unsigned char need_to_update;

    printf("hits = %d\n", hits);
    ptr = (GLuint *) buffer;
    for (i = 0; i < hits; i++) {  /* for each hit  */
        need_to_update = 0;
        names = *ptr;
        ptr++;

        current_depth = (double) *ptr/0x7fffffff;
        if (current_depth < smallest_depth) {
            smallest_depth = current_depth;
            need_to_update = 1;
        }
        ptr++;
        current_depth = (double) *ptr/0x7fffffff;
        if (current_depth < smallest_depth) {
            smallest_depth = current_depth;
            need_to_update = 1;
        }
        ptr++;
        for (j = 0; j < names; j++) {  /* for each name */
            if (need_to_update == 1)
                seed_id = *ptr;
            ptr++;
        }
    }

    return seed_id;
}
void CGlView::display_MCG_connections(void)
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
            /**/
            //float rgb[3] = {0.};
            //glColor4f(0.3, 0.3, 0.3, 0.8);

            //if(mcgnodes[node1].type == 0)
            //	glColor4f(0, 0.35, 0, 0.8);
            //else if(mcgnodes[node2].type == 1)
            //	glColor4f(0.35, 0, 0, 0.8);
            //else
            //	glColor4f(0.0, 0.0, 0.4, 0.8);

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

            //DrawSolidCircle_size(center[0], center[1], 0.007);
            //draw_ASolidSphere(center[0],center[1],center[2],0.005);

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

            //draw the edges
            //for(int e=0;e<morse_decomp->dg->nlist->dirnodes[face->index]->nedges;e++)
            //{
            //	int E=morse_decomp->dg->nlist->dirnodes[face->index]->edges[e];
            //	int dirN1=morse_decomp->dg->elist->edges[E]->node_index1;
            //	int dirN2=morse_decomp->dg->elist->edges[E]->node_index2;

            //	Triangle *t1=object->tlist.tris[dirN1];
            //	Triangle *t2=object->tlist.tris[dirN2];

            //	//Find the centers of the two triangles
            //	double c1[3],c2[3];

            //	Vertex* v10=t1->verts[0];
            //	Vertex* v11=t1->verts[1];
            //	Vertex* v12=t1->verts[2];
            //	Vertex* v20=t2->verts[0];
            //	Vertex* v21=t2->verts[1];
            //	Vertex* v22=t2->verts[2];

            //	c1[0] = (v10->x+v11->x+v12->x)/3.;
            //	c1[1] = (v10->y+v11->y+v12->y)/3.;
            //	c1[2] = (v10->z+v11->z+v12->z)/3.;
            //	c2[0] = (v20->x+v21->x+v22->x)/3.;
            //	c2[1] = (v20->y+v21->y+v22->y)/3.;
            //	c2[2] = (v20->z+v21->z+v22->z)/3.;

            //	glBegin(GL_LINES);
            //	glVertex3d(c1[0]+0.001*t1->normal.entry[0],c1[1]+0.001*t1->normal.entry[1],c1[2]+0.001*t1->normal.entry[2]);
            //	glVertex3d(c2[0]+0.001*t2->normal.entry[0],c2[1]+0.001*t2->normal.entry[1],c2[2]+0.001*t2->normal.entry[2]);
            //	glEnd();
            //}
        }
    }
}


void
CGlView::without_antialiasing(GLenum mode)
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

    //////we may first mark those triangles where we add elements at
    //if(DrawTrajectoryOn == 1)   ////Draw a single trajectory or whole??
    //    display_single_traj();

    //glEnable(GL_POLYGON_OFFSET_FILL);
    //glPolygonOffset (1., 1.);

    ////show separatrices
    if(ShowSeparatricesOn == 1)
        display_separatrices(mode);


    //////Display all the detected limit cycles
    if(ShowPeriodicOrbitsOn == 1)
        display_periodicorbits(mode);
    //	DisplayLimitCycleLegends(mode);  //display periodic orbit handler to edit it

    //glDisable(GL_BLEND);

    //glEnable(GL_LIGHTING);

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

void
CGlView::with_antialiasing(GLenum mode)
{
    glClear(GL_ACCUM_BUFFER_BIT);
    for(int i = 0; i < 16; i++)
    {
        glPushMatrix ();
        glTranslatef (ji16[i].x*1.0/(npix/*512*/*zoom_factor), ji16[i].y*1.0/(npix/*512*/*zoom_factor), 0.0);

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
        glDisable(GL_BLEND);

        //////we may first mark those triangles where we add elements at
        //if(DrawTrajectoryOn == 1)   ////Draw a single trajectory or whole??
        //    display_single_traj();

        //glEnable(GL_POLYGON_OFFSET_FILL);
        //glPolygonOffset (1., 1.);

        ////show separatrices
        if(ShowSeparatricesOn == 1)
            display_separatrices(mode);


        //////Display all the detected limit cycles
        if(ShowPeriodicOrbitsOn == 1)
            display_periodicorbits(mode);
        //	DisplayLimitCycleLegends(mode);  //display periodic orbit handler to edit it

        //glDisable(GL_BLEND);

        //glEnable(GL_LIGHTING);

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

        glPopMatrix ();
        glAccum(GL_ACCUM, 1.0/16);
    }
    glAccum (GL_RETURN, 1.0);
    glReadBuffer(GL_BACK);
    glLineWidth(1.);
}



void
CGlView::display_tau_distribution()
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


void
CGlView::display_diff_ECG_MCG()
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

void
CGlView::display_diff_MCGs()
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


void
CGlView::vis_rot_sum()
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
