#include "MCGWindow.h"
#include <GLUT/glut.h>

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

extern bool showRealIndex;

MCGWindow::MCGWindow(QWidget *parent) : QOpenGLWidget(parent)
{
    qInfo() << "MCG Window Default Constructor has finished";
}


MCGWindow::~MCGWindow()
{
    qInfo() << "MCG Window Default Destructor has finished";
}

void MCGWindow::draw_nodes(GLenum mode)
{
    int i;
    if( !this->ShowMCGOn || mcg == nullptr )
        return;

    for(i = 0; i < mcg->cur_mcgnode_index; i++)
    {
        if(mcg->nlist == NULL || mcg->nlist->mnodes[i] == NULL)
            return;

        if(mcg->nlist->mnodes[i]->cancelled)
            continue;

        if(mode == GL_SELECT)
            glLoadName(i+1);

        set_color(mcg->nlist->mnodes[i]->type);

        draw_solid_circle(mcg->nlist->mnodes[i]->pos.entry[0],
            mcg->nlist->mnodes[i]->pos.entry[1], 0.04);

        //draw the new conley index out circle
        if(ShowConleyCircle==true)
        {
            double incR=0.025;
            double R=0.041; // only for planar case example 03/07/2010 original value = 0.05

            if(mcg->nlist->mnodes[i]->type==0||mcg->nlist->mnodes[i]->type==1)//source ||sink
            {
                for(int a=0;a<mcg->nlist->mnodes[i]->conley[1];a++)
                    draw_circle(mcg->nlist->mnodes[i]->pos.entry[0],mcg->nlist->mnodes[i]->pos.entry[1], R+incR*(a+1),-1);
            }
            if(mcg->nlist->mnodes[i]->type==2)//saddle
            {
                if(mcg->nlist->mnodes[i]->conley[0]==0 && mcg->nlist->mnodes[i]->conley[2]==1 )
                    draw_circle(mcg->nlist->mnodes[i]->pos.entry[0],mcg->nlist->mnodes[i]->pos.entry[1], R+incR,0);
                else if(mcg->nlist->mnodes[i]->conley[0]==1 && mcg->nlist->mnodes[i]->conley[2]==0 )
                    draw_circle(mcg->nlist->mnodes[i]->pos.entry[0],mcg->nlist->mnodes[i]->pos.entry[1], R+incR,1);
                if(mcg->nlist->mnodes[i]->conley[1]!=0)
                {
                    for(int a=0;a<mcg->nlist->mnodes[i]->conley[1]/*min(2,mcg->nlist->mnodes[i]->conley[1])*/;a++)
                    draw_circle(mcg->nlist->mnodes[i]->pos.entry[0],mcg->nlist->mnodes[i]->pos.entry[1], R+incR*(a+1),-1);
                }
            }
        }
    }
}

void MCGWindow::draw_solid_circle(double cx, double cy, double R)
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

void MCGWindow::set_color(int nodetype)
{
    if(nodetype == 0)
        glColor3f(0, 1, 0);
    else if (nodetype == 1)
        glColor3f(1, 0, 0);
    else
        glColor3f(0, 0, 1);
}

void MCGWindow::draw_edges()
{
    int i;
    int node1, node2;
    icVector2 direct;
    double head[2];

    glLineWidth(1.5);

    if( !ShowMCGOn || mcg == nullptr)
        return;

    for(i = 0; i < mcg->cur_mcgedge_index; i++)
    {
        node1 = mcg->elist->edges[i]->node_index1;
        node2 = mcg->elist->edges[i]->node_index2;

        if(mcg->elist->edges[i]->cancel)
            continue;

        if(	mcg->nlist->mnodes[node1]->cancelled  || mcg->nlist->mnodes[node2]->cancelled )
            continue;


        glBegin(GL_LINES);
        set_color(mcg->nlist->mnodes[node1]->type);
        glVertex2f(mcg->nlist->mnodes[node1]->pos.entry[0],
            mcg->nlist->mnodes[node1]->pos.entry[1]);

        set_color(mcg->nlist->mnodes[node2]->type);
        glVertex2f(mcg->nlist->mnodes[node2]->pos.entry[0],
            mcg->nlist->mnodes[node2]->pos.entry[1] + 0.02);
        glEnd();

        ////Draw the wings of the arrow
        direct.entry[0] = mcg->nlist->mnodes[node2]->pos.entry[0] - mcg->nlist->mnodes[node1]->pos.entry[0];
        direct.entry[1] = mcg->nlist->mnodes[node2]->pos.entry[1] - mcg->nlist->mnodes[node1]->pos.entry[1];
        normalize(direct);

        head[0] = mcg->nlist->mnodes[node2]->pos.entry[0];
        head[1] = mcg->nlist->mnodes[node2]->pos.entry[1] + 0.02;
        draw_wings(head, direct);
    }
}

void MCGWindow::draw_wings(double head[], icVector2 direct)
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

void MCGWindow::display_label(int x, int y, char *string)
{
    int len, i;
    glRasterPos2f(x, y);
    len = (int) strlen(string);
    for (i = 0; i < len; i++)
    {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, string[i]);
    }
}

void MCGWindow::draw_labels()
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

    if( !ShowMCGOn || mcg == nullptr )
        return;

    for(i = 0; i < mcg->cur_mcgnode_index; i++)
    {
        if(mcg->nlist->mnodes[i]->cancelled)
            continue;

        glPushMatrix();

        posx = (int)((mcg->nlist->mnodes[i]->pos.entry[0]/4.)*(float)viewport[2])-35;
        posy = (int)((mcg->nlist->mnodes[i]->pos.entry[1]/2.)*(float)viewport[3])-5;

        glColor3f(0,0,0);
        if(mcg->nlist->mnodes[i]->type == 0) // a repeller
        {
            strings[0] = 'R';
        }

        else if(mcg->nlist->mnodes[i]->type == 1) // an attractor
        {
            strings[0] = 'A';
        }

        else   //saddle
        {
            strings[0] = 'S';
        }

        if(mcg->nlist->mnodes[i]->labelindex+1 >= 100)
        {
            int hundred = floor((mcg->nlist->mnodes[i]->labelindex+1)/100.);
            int ten = floor((mcg->nlist->mnodes[i]->labelindex + 1 - hundred * 100)/10.);
            int number = (mcg->nlist->mnodes[i]->labelindex+1) % 10;
            strings[1] = '000'+hundred;
            strings[2] = '000'+ten;
            strings[3] = '000'+number;
            strings[4] = '\0';
        }

        else if(mcg->nlist->mnodes[i]->labelindex+1 >= 10)
        {
            int ten = floor((mcg->nlist->mnodes[i]->labelindex+1)/10.);
            int number = (mcg->nlist->mnodes[i]->labelindex+1) % 10;
            strings[1] = '000'+ten;
            strings[2] = '000'+number;
            strings[3] = '\0';
        }

        else{
            strings[1] = '001'+mcg->nlist->mnodes[i]->labelindex;
            strings[2] = '\0';
        }

        if(showRealIndex)
        {
            int hundred = floor((mcg->nlist->mnodes[i]->node_index)/100.);
            int ten = floor((mcg->nlist->mnodes[i]->node_index  - hundred * 100)/10.);
            int number = (mcg->nlist->mnodes[i]->node_index) % 10;
            strings[0] = '000'+hundred;
            strings[1] = '000'+ten;
            strings[2] = '000'+number;
            strings[3] = '\0';
        }

        display_label(posx, posy, strings);

        glPopMatrix();
    }
    ////Reset the viewport back to original setting
    glMatrixMode(GL_PROJECTION);		// Select The Projection Matrix
    glLoadIdentity();					// Reset The Projection Matrix

    // Calculate The Aspect Ratio Of The Window
    glOrtho(0.0f,4.0f,0.0f,2.0f,0.0f,50.0f);
}

extern int picked_node;
void MCGWindow::DrawGLScene(GLenum mode)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer

    glDisable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);

    if(picked_node >= 0)
        draw_highlights();     //Draw highlighted objects if exist
    draw_edges();   //Draw edges first
    draw_nodes(mode);          //Draw nodes
    draw_labels();
}


void MCGWindow::draw_highlights()
{
    glColor3f(1, 0.6, 0.3);   ////set yellow as the highlighted color
    glLineWidth(3);

    draw_circle(mcg->nlist->mnodes[picked_node]->pos.entry[0],
            mcg->nlist->mnodes[picked_node]->pos.entry[1],0.04,-1);
}

void MCGWindow::draw_circle(double cx, double cy, double radius, int type)
{
    int i;
    double R = radius;
    double theta, deta ;
    deta = 2 * M_PI/49.;
    double x, y;
    theta = 0.;

    if(type>=0)
    set_color(type);
    glBegin(GL_LINE_LOOP);
    for(i = 0; i < 50; i++, theta += deta)
    {
        x =  cx + R * cos(theta);
        y =  cy + R * sin(theta);

        glVertex2f(x, y);
    }
    glEnd();
}

void MCGWindow::update_scene()
{
    this->update();
}

void MCGWindow::HitProcessforGraph(double ss, double tt)
{
    GLuint  selectBuffer[SELECTBUFFERSIZE];
    GLint	vp[4] = {0, 0 ,4, 2};
    int hits;

    ////Build the selection buffer here
    glSelectBuffer(SELECTBUFFERSIZE, selectBuffer);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    glRenderMode(GL_SELECT);
    glInitNames();
    glPushName(0);

    // 0.01 is better than 0.005
    gluPickMatrix(ss, tt, 0.01, 0.01, vp );  ////set a larger pick window for element selection
    glOrtho(0, 4,  0, 2,  0, 50);

    draw_nodes(GL_SELECT);

    hits = glRenderMode(GL_RENDER);

    if(hits>0)
    {
        picked_node = selectBuffer[3]-1;
    }
    else{
        picked_node = -1;
    }

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();

    glMatrixMode( GL_MODELVIEW );
}


void MCGWindow::initializeGL()
{
    this->initializeOpenGLFunctions();

    qInfo() << "MCG Window has been initialized";
}


void MCGWindow::paintGL()
{
    glClearColor (1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    this->DrawGLScene(GL_RENDER);
    this->update();
}


void MCGWindow::resizeGL(int w, int h)
{
    qInfo() << "MCG Window has been resized";
}

void MCGWindow::mouseMoveEvent(QMouseEvent *event)
{
    double s, t;

    int firstwin_leftbottom_x = 0;
    int firstwin_leftbottom_y = 150;
    int firstwin_rightx = 400;
    int firstwin_bottomy = 150;

    QPointF p = event->pos();

    ScreenToSecondWin(p.x(), p.y(), firstwin_leftbottom_x, firstwin_leftbottom_y,
                    firstwin_rightx, firstwin_bottomy, 0, 0, 4, 2, s, t);

    if( (s<0) || (s > 4) || (t < 0) || (t > 2)){
        return;
    }
    this->HitProcessforGraph(s, t);
}
