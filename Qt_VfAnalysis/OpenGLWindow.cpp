#include "OpenGLWindow.h"
// this file serves as IBFVDlg.cpp, GLView.cpp
#include <QDir>


// GLOBAL VARIABLES
const char filePath [] = "../Qt_VfAnalysis/datasets/simple_fld2.ply";

const int win_height=800;
const int win_width=800;

/*------------------------------------------------------------*/
//some global variables for ibfv
double SCALE = 3.0;
int     iframe = 0;
int     Npat   = 64;        // Modified at 02/10/2010
int     alpha  = (0.06*255);  // modified for a smear out LIC texture for better visualization (TVCG revision 11/22/2010)
double  sa;
double  tmax   = NPIX/(SCALE*NPN);
double  dmax   = 1.8/NPIX;
int npn = NPN;
int npix = 800;

int Cal_Regions=0;
int ndisplay_trajs;
bool RemoveDisconnMSOn = false;
int Integrator_opt = 0;
int ndisplay_POs;
clock_t g_start, g_finish;


// some extern varibles
extern Polyhedron *object; // declared in geometry.cpp

extern MorseDecomp *morse_decomp;
extern MorseDecomp *L1_morse;
extern MorseDecomp *L2_morse;
extern ECG_Graph *ecg;

OpenGLWindow::OpenGLWindow(QWidget *parent) : QOpenGLWidget(parent)
{
    qInfo() << "OpenGLWindow Default Constructor called";
    qDebug() << QDir::currentPath();

}


OpenGLWindow::~OpenGLWindow() {
    qInfo() << "OpenGLWindow Desctuctor called";
}


void OpenGLWindow::build_Polyhedron()
{
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

    object->cal_TexCoord();

    morse_decomp = new MorseDecomp(); /*initialize the Morse Decomposition component*/

    L1_morse = new MorseDecomp();
    L2_morse = new MorseDecomp();

    //mcg = new MCG_Graph();
    //mcg->init_MCG();
}


void OpenGLWindow::initializeGL(){
    this->build_Polyhedron();
    this->initializeOpenGLFunctions();
    this->InitGL();

    glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    qInfo() << "OpenGLWindow has been initialized";
}


void OpenGLWindow::paintGL(){
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    this->draw();
}


void OpenGLWindow::resizeGL(int w, int h){

    qInfo() << "Resize w" << w << "  h" <<h;
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}


void OpenGLWindow::draw(){
    float r = 0.0f, g = 0.0f, b = 0.0f;
    glColor3f(r, g, b);
    glBegin(GL_LINE_STRIP);
        glVertex3f(0, 0, 0);
        glVertex3f(0.5, 0.5, 0);
        glVertex3f(0, 0.5, 0.5);
        glVertex3f(0, 0, 0);
    glEnd();

}


void OpenGLWindow::set_up_MainWindow_ptr(MainWindow *MW_ptr)
{
    this->mainWindow = MW_ptr;
}


int OpenGLWindow::InitGL( )
{
    glViewport(0, 0, (GLsizei) NPIX, (GLsizei) NPIX);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    return TRUE;
}


void OpenGLWindow::init_flags()
{
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

    this->ShowEdgeSamplesOn=false;
    this->ShowTriMappingOn = false;
    this->selected_triangle = -1;  // 02/10/2010

    this->ShowBackward=0;//0=forward, 1=backward
    this->sampling_edge=0;

    this->ShowConnectionRegion=0;

    this->ShowBoundary=1;

    this->FlipNormalOn = false;
    this->IBFVOff = false;
}
