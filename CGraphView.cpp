/*
This file contains the implementation of the memeber functions of the class CGraphView

Created and modified by Guoning Chen
        copyright @2007
    */

#include "CGraphView.h"
#include "Analysis/MorseDecomp.h"

extern MCG_Graph *mcg;
extern ECG_Graph *ecg;

int picked_node = -1;                  ////The picked up node when mouse moves to it

#define SELECTBUFFERSIZE 128


CGraphView::CGraphView(void)
{

}

CGraphView::~CGraphView(void)
{

}


int CGraphView::OnCreate()
{

    m_hDC = ::GetDC(this->m_hWnd);

    if(!SetPixelformat(m_hDC))
    {
        ::MessageBox(::GetFocus(),"SetPixelformat Failed!","Error",MB_OK);
        return -1;
    }

    m_hglRC = wglCreateContext(m_hDC);
    int i= wglMakeCurrent(m_hDC,m_hglRC);

    InitGL();

    return 0;
}

/*--------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------*/
// Other extended operations for OpenGl window setting

BOOL CGraphView::SetPixelformat(HDC hdc)
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

/* --------------------------------------------------
*  Initialize the opengl environment
---------------------------------------------------*/
int CGraphView::InitGL(GLvoid)								// All Setup For OpenGL Goes Here
{
    glMatrixMode(GL_PROJECTION);		// Select The Projection Matrix
    glLoadIdentity();					// Reset The Projection Matrix

    // Calculate The Aspect Ratio Of The Window
    glOrtho(0.0f,4.0f,0.0f,2.0f,0.0f,50.0f);

    glMatrixMode(GL_MODELVIEW);							// Select The Modelview Matrix
    glLoadIdentity();									// Reset The Modelview Matrix

    glClearColor(1.0f, 1.0f, 1.0f, 1.f);				// Black Background

    return TRUE;										// Initialization Went OK
}



void CGraphView::init_flags()
{
    ShowMCGOn = 0;
    ShowConleyCircle=0;
}

void  CGraphView::get_rotation(double p[3], icVector3 rot_axis, double rotang, double rotmatrix[16])
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

    /*Write down the pure rotation*/
    //FILE *fp = fopen("rot_opengl.txt", "w");
    //for(int i=0; i<16; i++)
    //		fprintf(fp, "%f\n", rotmatrix[i]);
    //fclose(fp);
}


void  CGraphView::get_rotation(double rotang, icVector3 rot_axis, double rotmatrix[16])
{
    wglMakeCurrent(m_hDC, m_hglRC);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glRotatef(rotang, rot_axis.entry[0], rot_axis.entry[1], rot_axis.entry[2]);
    glGetDoublev(GL_MODELVIEW_MATRIX, (double *)rotmatrix);
    glPopMatrix();
}

void CGraphView::rightMultiply16(double old_rotmat[16], double rot_mat[16])
{
    wglMakeCurrent(m_hDC, m_hglRC);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMultMatrixd ((double*)old_rotmat);
    glMultMatrixd ((double*)rot_mat);
    glGetDoublev(GL_MODELVIEW_MATRIX, (double *)old_rotmat);
    glPopMatrix();
}


//
/* --------------------------------------------------
*  Draw the nodes of the graph
---------------------------------------------------*/
void CGraphView::draw_nodes(GLenum mode)
{
    int i;

    if(ShowMCGOn == 0)
    {
        if(ecg == NULL)
            return;

        for(i = 0; i < ecg->cur_ecgnode_index/*ecg->nlist->nenodes*/; i++)
        {
            if(ecg->nlist == NULL || ecg->nlist->enodes[i] == NULL)
                return;

            if(ecg->nlist->enodes[i]->cancelled)
                continue;

            if(mode == GL_SELECT)
                glLoadName(i+1);

            set_color(ecg->nlist->enodes[i]->type);

            //if(mcg->nlist->mnodes[i]->LimitCycleID >= 0)
            //{
            //	//DrawSolidSquare(graphnodes[i].pos_x, graphnodes[i].pos_y);
            //	int type = limitcycles[graphnodes[i].LimitCycleID].type;
            //	DrawLimitCircle(graphnodes[i].pos_x, graphnodes[i].pos_y, type);
            //	DrawSolidCircle(graphnodes[i].pos_x, graphnodes[i].pos_y, 0.04);
            //}

            //else{

            draw_solid_circle(ecg->nlist->enodes[i]->pos.entry[0],
                              ecg->nlist->enodes[i]->pos.entry[1], 0.04);
            //}

            if (ecg->nlist->enodes[i]->periodicorbitID>=0)
            {
                draw_circle(ecg->nlist->enodes[i]->pos.entry[0],
                            ecg->nlist->enodes[i]->pos.entry[1], 0.06, ecg->nlist->enodes[i]->type);
            }
        }
    }
    else{
        if(mcg == NULL)
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
            if(ShowConleyCircle==1)
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
}



void CGraphView::draw_solid_circle(double cx, double cy, double R)
{
    int i;
    //double R = 0.06;
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


void CGraphView::set_color(int nodetype)
{
    if(nodetype == 0)
        glColor3f(0, 1, 0);
    else if (nodetype == 1)
        glColor3f(1, 0, 0);
    else
        glColor3f(0, 0, 1);
}



void CGraphView::draw_edges()
{
    int i;
    int node1, node2;
    icVector2 direct;
    double head[2];

    glLineWidth(1.5);
    if(ShowMCGOn == 0)
    {
        if(ecg == NULL)
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
            draw_wings(head, direct);
        }
    }

    else
    {
        if(mcg == NULL)
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
}


void CGraphView::draw_wings(double head[2], icVector2 direct)
{
    glPushMatrix();
    glTranslatef(head[0], head[1], 0);
    glRotatef(atan2(direct.entry[1], direct.entry[0])*360/(2*M_PI), 0, 0, 1);
    glScalef(0.33, 0.33, 1);

    ////Draw the wings of the arrow
    //glBegin(GL_LINES);
    //glVertex2f(0, 0);
    //glVertex2f(-0.2, 0.16);

    //glVertex2f(0, 0);
    //glVertex2f(-0.2, -0.16);
    //glEnd();
    glBegin(GL_TRIANGLES);
    glVertex2f(0, 0);
    glVertex2f(-0.35, 0.12);
    glVertex2f(-0.35, -0.12);
    glEnd();

    glPopMatrix();
}

////display the labels
void CGraphView::display_label(int x, int y, char *string)
{
    int len, i;
    glRasterPos2f(x, y);
    len = (int) strlen(string);
    for (i = 0; i < len; i++)
    {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, string[i]);
    }
}

extern bool showrealindex;

void CGraphView::draw_labels()
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

    if(ShowMCGOn == 0)
    {
        if(ecg == NULL)
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

            display_label(posx+5, posy, strings);

            glPopMatrix();
        }
    }

    else
    {
        if(mcg == NULL)
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

            if(showrealindex)
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
    }
    ////Reset the viewport back to original setting
    glMatrixMode(GL_PROJECTION);		// Select The Projection Matrix
    glLoadIdentity();					// Reset The Projection Matrix

    // Calculate The Aspect Ratio Of The Window
    glOrtho(0.0f,4.0f,0.0f,2.0f,0.0f,50.0f);
}


void CGraphView::DrawGLScene(GLenum mode)
{
    ////Draw the conley relation graph according to the graphnodes and graphedges structure
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer

    glDisable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);

    if(/*num_related_edges > 0 && */picked_node >= 0)
        draw_highlights();     //Draw highlighted objects if exist

    draw_edges();          //Draw edges first

    draw_nodes(mode);          //Draw nodes

    draw_labels();

    SwapBuffers(m_hDC);
}


/*--------------------------------------------------------------------------------*/
/*Mouse selection routine
/*--------------------------------------------------------------------------------*/
void
CGraphView::HitProcessforGraph(double ss, double tt)
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

    gluPickMatrix(ss, tt, 0.005, 0.005, vp );  ////set a larger pick window for element selection
    glOrtho(0, 4,  0, 2,  0, 50);

    draw_nodes(GL_SELECT);
    //picked_node = -1;

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



void
CGraphView::draw_highlights()
{
    int i;
    int node1, node2;

    glColor3f(1, 0.6, 0.3);   ////set yellow as the highlighted color
    glLineWidth(3);

    if(ShowMCGOn == 0)
    {
        ////1. Draw the highlighted edges first
        //glBegin(GL_LINES);
        //for(i = 0; i < num_related_edges; i++)
        //{
        //	node1 = graphedges[related_edges[i]].node_index1;
        //	node2 = graphedges[related_edges[i]].node_index2;

        //	////10/28/05 temporary
        //	if(node1 < 0 || node2 < 0)
        //		return;

        //	if(graphnodes[node1].cancelled == 1 || graphnodes[node2].cancelled == 1)
        //		continue;

        //	////Line1
        //	glVertex2f(graphnodes[node1].pos_x, graphnodes[node1].pos_y-0.01);
        //	glVertex2f(graphnodes[node2].pos_x, graphnodes[node2].pos_y-0.01);
        //
        //	////Line2
        //	glVertex2f(graphnodes[node1].pos_x, graphnodes[node1].pos_y+0.01);
        //	glVertex2f(graphnodes[node2].pos_x, graphnodes[node2].pos_y+0.01);
        //}
        //glEnd();

        //////2. Draw the highlighted nodes
        //DrawCircle(graphnodes[picked_node].pos_x, graphnodes[picked_node].pos_y);

        //for(i = 0; i < num_related_edges; i++)
        //{
        //	node1 = graphedges[related_edges[i]].node_index1;

        //	if(graphnodes[node1].cancelled == 1)
        //		continue;

        //	if(node1 == picked_node)
        //		node1 = graphedges[related_edges[i]].node_index2;

        //	DrawCircle(graphnodes[node1].pos_x, graphnodes[node1].pos_y);
        //}
    }

    else
    {
        ////2. Draw the highlighted nodes
        draw_circle(mcg->nlist->mnodes[picked_node]->pos.entry[0],
                    mcg->nlist->mnodes[picked_node]->pos.entry[1],0.04,-1);

        //for(i = 0; i < num_related_edges; i++)
        //{
        //	node1 = mcgedges[related_edges[i]].node_index1;

        //	if(mcgnodes[node1].cancelled == 1)
        //		continue;

        //	if(node1 == picked_node)
        //		node1 = mcgedges[related_edges[i]].node_index2;

        //	DrawCircle(mcgnodes[node1].pos_x, mcgnodes[node1].pos_y);
        //}
    }
}



/* --------------------------------------------------
*  Draw a hollow circle for highlighting the nodes
*  that selected by user in the graph
---------------------------------------------------*/
void
CGraphView::draw_circle(double cx, double cy, double radius,int type)
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
