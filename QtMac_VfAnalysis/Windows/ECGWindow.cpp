#include "Windows/ECGWindow.h"
#include <GLUT/glut.h>
#include <QPainter>

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


ECGWindow::ECGWindow(QWidget *parent) : QOpenGLWidget(parent)
{
    qInfo() << "ECG Window Default Constructor has finished";
}


ECGWindow::~ECGWindow()
{
    qInfo() << "ECG Window Default Destructor has finished";
}

void ECGWindow::draw_nodes(GLenum mode)
{
    int i;
    if( !ShowECGOn )
        return;

    for(i = 0; i < ecg->cur_ecgnode_index; i++)
    {
        if(ecg->nlist == nullptr || ecg->nlist->enodes[i] == nullptr)
            return;

        if(ecg->nlist->enodes[i]->cancelled)
            continue;

        if(mode == GL_SELECT)
            glLoadName(i+1);

        this->set_color(ecg->nlist->enodes[i]->type);

        this->draw_solid_circle(ecg->nlist->enodes[i]->pos.entry[0],
            ecg->nlist->enodes[i]->pos.entry[1], 0.04);

        if (ecg->nlist->enodes[i]->periodicorbitID>=0)
        {
            draw_circle(ecg->nlist->enodes[i]->pos.entry[0],
            ecg->nlist->enodes[i]->pos.entry[1], 0.06, ecg->nlist->enodes[i]->type);
        }
    }
}


void ECGWindow::draw_solid_circle(double cx, double cy, double R)
{
    int i;
    double theta, deta ;
    deta = 2 * M_PI/80.;
    double x, y;
    theta = 0.;
    glBegin(GL_POLYGON);
    for(i = 0; i < 80; i++, theta += deta)
    {
        x =  cx + R * cos(theta);
        y =  cy + R * sin(theta);

        glVertex2f(x, y);
    }
    glEnd();
}

void ECGWindow::set_color(int nodetype)
{
    if(nodetype == 0) // source-like
        glColor3f(0, 1, 0);
    else if (nodetype == 1) // sink-like
        glColor3f(1, 0, 0);
    else // saddle-like
        glColor3f(0, 0, 1);
}

void ECGWindow::draw_edges()
{
    int i;
    int node1, node2;
    icVector2 direct;
    double head[2];

    glLineWidth(1.5);

    if( !ShowECGOn )
        return;

    for(i = 0; i < ecg->cur_ecgedge_index; i++)
    {
        node1 = ecg->elist->edges[i]->node_index1;
        node2 = ecg->elist->edges[i]->node_index2;

        if(	ecg->nlist->enodes[node1]->cancelled || ecg->nlist->enodes[node2]->cancelled )
            continue;


        glBegin(GL_LINES);
        set_color(ecg->nlist->enodes[node1]->type);
        glVertex2f(ecg->nlist->enodes[node1]->pos.entry[0],
            ecg->nlist->enodes[node1]->pos.entry[1]);

        set_color(ecg->nlist->enodes[node2]->type);
        glVertex2f(ecg->nlist->enodes[node2]->pos.entry[0],
            ecg->nlist->enodes[node2]->pos.entry[1] + 0.02);
        glEnd();

        ////Draw the wings of the arrow
        direct.entry[0] = ecg->nlist->enodes[node2]->pos.entry[0] -
            ecg->nlist->enodes[node1]->pos.entry[0];
        direct.entry[1] = ecg->nlist->enodes[node2]->pos.entry[1] -
            ecg->nlist->enodes[node1]->pos.entry[1];
        normalize(direct);

        head[0] = ecg->nlist->enodes[node2]->pos.entry[0] ;
        head[1] = ecg->nlist->enodes[node2]->pos.entry[1] + 0.02;
        this->draw_wings(head, direct);
    }
}

void ECGWindow::draw_wings(double head[], icVector2 direct)
{
    glPushMatrix();
    glTranslatef(head[0], head[1], 0);
    glRotatef(atan2(direct.entry[1], direct.entry[0])*360/(2*M_PI), 0, 0, 1);
    glScalef(0.33, 0.33, 1);

    ////Draw the wings of the arrow
    glBegin(GL_TRIANGLES);
    glVertex2f(0, 0);
    glVertex2f(-0.35, 0.12);
    glVertex2f(-0.35, -0.12);
    glEnd();

    glPopMatrix();
}

void ECGWindow::display_label(int x, int y, char* string)
{
    int len, i;
    glRasterPos2f(x, y);
    len = (int) strlen(string);
    for (i = 0; i < len; i++)
    {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, string[i]);
    }

//    QPainter painter(this);
//    painter.setPen(Qt::black);
//    painter.setFont(QFont("Arial", 10));
//    painter.drawText(x, y, QString(string));

}

void ECGWindow::draw_labels()
{
    int i;
    int posx, posy;
    char strings[5];
    //int r_counter, a_counter, s_counter;

    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    glMatrixMode(GL_PROJECTION);		// Select The Projection Matrix
    glLoadIdentity();					// Reset The Projection Matrix


    // Calculate The Aspect Ratio Of The Window
    glOrtho(viewport[0],viewport[2],viewport[1],viewport[3],0.0f,50.0f);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    if( ! ShowECGOn )
        return;

    for(i = 0; i < ecg->cur_ecgnode_index; i++)
    {
        if(ecg->nlist->enodes[i]->cancelled)
            continue;

        glPushMatrix();

        posx = (int)((ecg->nlist->enodes[i]->pos.entry[0]/4.)*(float)viewport[2])-35;
        posy = (int)((ecg->nlist->enodes[i]->pos.entry[1]/2.)*(float)viewport[3])-5;

        glColor3f(0,0,0);
        if(ecg->nlist->enodes[i]->type == 0) // a repeller
        {
            strings[0] = 'R';
        }

        else if(ecg->nlist->enodes[i]->type == 1) // an attractor
        {
            strings[0] = 'A';
        }

        else   //saddle
        {
            strings[0] = 'S';
        }

        if(ecg->nlist->enodes[i]->labelindex+1 >= 100)
        {
            int hundred = floor((ecg->nlist->enodes[i]->labelindex+1)/100.);
            int ten = floor((ecg->nlist->enodes[i]->labelindex + 1 - hundred * 100)/10.);
            int number = (ecg->nlist->enodes[i]->labelindex+1) % 10;
            strings[1] = '000'+hundred;
            strings[2] = '000'+ten;
            strings[3] = '000'+number;
            strings[4] = '\0';
        }

        else if(ecg->nlist->enodes[i]->labelindex+1 >= 10)
        {
            int ten = floor((ecg->nlist->enodes[i]->labelindex+1)/10.);
            int number = (ecg->nlist->enodes[i]->labelindex+1) % 10;
            strings[1] = '000'+ten;
            strings[2] = '000'+number;
            strings[3] = '\0';
        }

        else{
            strings[1] = '001'+ecg->nlist->enodes[i]->labelindex;
            strings[2] = '\0';
        }

        this->display_label(posx+5, posy, strings);

        glPopMatrix();
    }

    ////Reset the viewport back to original setting
    glMatrixMode(GL_PROJECTION);		// Select The Projection Matrix
    glLoadIdentity();					// Reset The Projection Matrix

    // Calculate The Aspect Ratio Of The Window
    glOrtho(0.0f,4.0f,0.0f,2.0f,0.0f,50.0f);
}

void ECGWindow::DrawGLScene(GLenum mode)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer

    glDisable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);

    draw_edges();          //Draw edges first

    draw_nodes(mode);          //Draw nodes

    draw_labels();
}

void ECGWindow::draw_circle(double cx, double cy, double radius, int type)
{
    int i;
    double R = radius;
    double theta, deta ;
    deta = 2 * M_PI/49.;
    double x, y;
    theta = 0.;

    if(type>=0)
    this->set_color(type);
    glBegin(GL_LINE_LOOP);
    for(i = 0; i < 50; i++, theta += deta)
    {
        x =  cx + R * cos(theta);
        y =  cy + R * sin(theta);

        glVertex2f(x, y);
    }
    glEnd();
}

void ECGWindow::update_scene()
{
    this->update();
}

void ECGWindow::initializeGL()
{
    this->initializeOpenGLFunctions();
    qInfo() << "ECG Window has been initialized";
}


void ECGWindow::paintGL()
{
    glClearColor (1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    this->DrawGLScene(GL_RENDER);
    this->update();
}

void ECGWindow::resizeGL(int w, int h)
{
    qInfo() << "ECG Window has been resized";
}
