#include "ConleyIndex.h"
#include "Analysis/MorseDecomp.h"

extern MorseDecomp *morse_decomp;
extern MorseDecomp *local_decomp;
extern MCG_Graph *mcg;
extern Polyhedron *object;
extern int refined_morse;
//ConleyIndex *con;

ConleyIndex::ConleyIndex(void)
{
}

ConleyIndex::~ConleyIndex(void)
{
}

void ConleyIndex::build_morse_nlist(int scc,int local)
{
    if(morse_nlist.size())morse_nlist.clear();
    if(local)
    {
        for(int nodeID=0;nodeID<local_decomp->scclist->scccomponents[scc]->nnodes;nodeID++)
            morse_nlist.push_back(local_decomp->dg->nlist->dirnodes[local_decomp->scclist->scccomponents[scc]->nodes[nodeID]]->global_index);
    }
    else
    {
        for(int nodeID=0;nodeID<morse_decomp->scclist->scccomponents[scc]->nnodes;nodeID++)
            morse_nlist.push_back(morse_decomp->scclist->scccomponents[scc]->nodes[nodeID]);
    }
}

void ConleyIndex::compute_conley_index(int local)
{
    //initialize
    if(pro.size())pro.clear();

    if(local)
    {
        isLocal=true;
        for(int i=0;i<local_decomp->scclist->nsccs;i++)

            //for(int i=0;i<3;i++)
            GetConley(i,local);
        //copy temp to the scc
        for(int i=0;i<pro.size();i++)
        {
            local_decomp->scclist->scccomponents[i]->XM=pro[i].XM;
            local_decomp->scclist->scccomponents[i]->XL=pro[i].XL;
            local_decomp->scclist->scccomponents[i]->classification=pro[i].classification;
            local_decomp->scclist->scccomponents[i]->Conley0=pro[i].conleyB0;
            local_decomp->scclist->scccomponents[i]->Conley1=pro[i].conleyB1;
            local_decomp->scclist->scccomponents[i]->Conley2=pro[i].conleyB2;
        }
    }
    else
    {
        //initialize boundary
        for(int i=0;i<object->elist.nedges;i++)
            object->elist.edges[i]->boundary=false;

        isLocal=false;

        for(int i=0;i<morse_decomp->scclist->nsccs;i++)
        {
            GetConley(i,local);

        }
        //copy temp to the scc
        for(int i=0;i<pro.size();i++)
        {
            morse_decomp->scclist->scccomponents[i]->XM=pro[i].XM;
            morse_decomp->scclist->scccomponents[i]->XL=pro[i].XL;
            morse_decomp->scclist->scccomponents[i]->classification=pro[i].classification;
            morse_decomp->scclist->scccomponents[i]->Conley0=pro[i].conleyB0;
            morse_decomp->scclist->scccomponents[i]->Conley1=pro[i].conleyB1;
            morse_decomp->scclist->scccomponents[i]->Conley2=pro[i].conleyB2;

        }
    }
    //output pro
    /*FILE *xm=fopen("conley.txt","w");
    fprintf(xm,"scc\tXM\tXL\tTYPE\tB0\tB1\tB2\n");
    for(int i=0;i<pro.size();i++)
    {
        fprintf(xm,"%d\t%d\t%d\t",i,pro[i].XM,pro[i].XL);
        switch(pro[i].classification)
        {
            case 0: fprintf(xm,"trivial\t");break;
            case 1: fprintf(xm,"source\t");break;
            case 2: fprintf(xm,"sink\t");break;
            case 3: fprintf(xm,"saddle\t");break;
        }
        fprintf(xm,"%d\t%d\t%d\t%d\n",pro[i].conleyB0,pro[i].conleyB1,pro[i].conleyB2,pro[i].size);
    }
    fclose(xm);*/

    //output morse
    /*FILE *xm=fopen("conley.txt","w");
    fprintf(xm,"scc\tmorse\tXM\tXL\tTYPE\tB0\tB1\tB2\tnodes#\n");
    for(int i=0;i<morse_decomp->scclist->nsccs;i++)
    {
        SCComponent *s=morse_decomp->scclist->scccomponents[i];
        fprintf(xm,"%d\t%d\t%d\t%d\t",i,s->node_index,s->XM,s->XL);
        switch(s->classification)
        {
            case 0: fprintf(xm,"trivial\t");break;
            case 1: fprintf(xm,"source\t");break;
            case 2: fprintf(xm,"sink\t");break;
            case 3: fprintf(xm,"saddle\t");break;
        }
        fprintf(xm,"%d\t%d\t%d\t%d\n",s->Conley0,s->Conley1,s->Conley2,s->nnodes);
    }
    fclose(xm);*/

}

void ConleyIndex::GetConley(int scc,int local)
{
    //initial compute a new morse set
    SCCProperties tempPro;
    pro.push_back(tempPro);

    //initialize
    build_morse_nlist(scc,local);
    pro[scc].size=morse_nlist.size();
    cur_scc=scc;

    //compute X(M),X(L),Conley Index
    pro[cur_scc].XM=computeXM();
    computeXL();
    computeBetti();
}


int ConleyIndex::computeXM(void)
{
    //set edge and vertex unvisited, must be the whole mesh, since two morse sets can be neighbour
    for(int i=0;i<object->vlist.nverts;i++)
        object->vlist.verts[i]->visited=false;
    for(int i=0;i<object->elist.nedges;i++)
        object->elist.edges[i]->visited=false;

    //count face#, vertex#,and edge#
    int tN=0,vN=0,eN=0;
    tN=morse_nlist.size();
    for(int n=0;n<morse_nlist.size();n++)
    {
        Triangle *t=object->tlist.tris[morse_nlist[n]];

        //count vertex#
        for(int i=0;i<3;i++)
        {
            Vertex *v=t->verts[i];
            if(v->visited==false)
            {
                vN++;
                v->visited=true;
            }
        }

        //count edge#
        for(int i=0;i<3;i++)
        {
            Edge *e=t->edges[i];
            if(e->visited==false)
            {
                eN++;
                e->visited=true;
            }
        }
    }
    return (vN+tN-eN);
}


void ConleyIndex::computeXL(void)
{
    //get the edgelist
    BoundaryEdgeExtraction();
    //get the exit_elist
    get_exit_elist();
    ////get the exit_vlist
    get_exit_vlist();

    int n_all_e=boundary_edgelist.size();
    int n_e_exit=exit_elist.size();
    int n_v_exit=exit_vlist.size();

    //fill in XL
    //pro[cur_scc].XL=n_v_exit-n_e_exit;
    XL=n_v_exit-n_e_exit;


    //fill in type
    if(n_e_exit==0)
        //pro[cur_scc].classification=1;//sink
        classification=1;//sink
    else
    {
        if(n_all_e==n_e_exit)
            //pro[cur_scc].classification=2;//source
            classification=0;//source
        else
            //pro[cur_scc].classification=3;//saddle, need to be reassigned
            classification=2;//saddle, need to be reassigned
    }
}

void ConleyIndex::BoundaryEdgeExtraction(void)
{
    Triangle *t;
    Edge *e;
    Triangle *t1,*t2;
    int scc1, scc2;

    //initialize the boundary_edgelist
    if(boundary_edgelist.size())boundary_edgelist.clear();
    if(boundary_tlist.size())boundary_tlist.clear();

    //set all edge.visited to 0
    for(int i=0;i<object->elist.nedges;i++)
    {
        object->elist.edges[i]->visited=0;
    }
    //start
    //if(isLocal)
    //{
    //	for(int i=0;i<morse_nlist.size();i++)
    //	{
    //		t=object->tlist.tris[morse_nlist[i]];
    //		//record edge
    //		for(int j=0;j<3;j++)
    //		{
    //			e=t->edges[j];
    //			if(e->visited==false)
    //			{
    //				e->visited=true;
    //				bool t0in,t1in;
    //				t0in=IsInScc(e->tris[0]->index);
    //				t1in=IsInScc(e->tris[1]->index);
    //
    //				if(t0in && t1in)
    //					continue;
    //				boundary_edgelist.push_back(e->index);
    //				e->boundary=true;
    //				if(t0in)
    //					boundary_tlist.push_back(e->tris[0]->index);
    //				else
    //					boundary_tlist.push_back(e->tris[1]->index);
    //			}
    //		}
    //	}
    //}
    //else
    {
        for(int i=0;i<morse_nlist.size();i++)
        {
            t=object->tlist.tris[morse_nlist[i]];
            //record edge
            for(int j=0;j<3;j++)
            {
                e=t->edges[j];
                if(e->visited==false)
                {
                    //get the two scc of the two triangles
                    //if(e->tris[1]<0)
                    if(e->tris[1]==nullptr)
                        scc2=-1;
                    else
                        scc2=morse_decomp->dg->nlist->dirnodes[e->tris[1]->index]->sscomp_index;
                    scc1=morse_decomp->dg->nlist->dirnodes[e->tris[0]->index]->sscomp_index;
                    if(scc1!=scc2)
                    {
                        boundary_edgelist.push_back(e->index);
                        e->boundary=true;
                        if(scc1==cur_scc)
                            boundary_tlist.push_back(e->tris[0]->index);
                        else
                            boundary_tlist.push_back(e->tris[1]->index);
                    }
                    e->visited=true;
                }
            }
        }
    }
}


void ConleyIndex::get_exit_elist(void)
{
    if(exit_elist.size())exit_elist.clear();
    for(int eN=0;eN<boundary_edgelist.size();eN++)
    {
        //int type=GetDir_Tau(eN);
        int type;
        /*if(isLocal)
            type=GetDir_Geo(eN);
        else*/
        type=GetDir_Tau(eN);
        if(type==EXIT)
            exit_elist.push_back(boundary_edgelist[eN]);
    }
}
void ConleyIndex::get_exit_vlist(void)
{
    if(exit_vlist.size())exit_vlist.clear();
    for(int i=0;i<object->vlist.nverts;i++)
        object->vlist.verts[i]->visited=false;
    for(int i=0;i<exit_elist.size();i++)
    {
        Edge *e;
        e=object->elist.edges[exit_elist[i]];
        int v1,v2;
        v1=e->verts[0]->index;
        v2=e->verts[1]->index;
        if(object->vlist.verts[v1]->visited==false)
        {
            exit_vlist.push_back(v1);
            object->vlist.verts[v1]->visited=true;
        }
        if(object->vlist.verts[v2]->visited==false)
        {
            exit_vlist.push_back(v2);
            object->vlist.verts[v2]->visited=true;
        }
    }
}
int ConleyIndex::GetDir_Geo(int edge)
{
    //initialize
    Edge *e;
    Vertex  *v1,*v2,*v3;
    Triangle *t;
    int scc1, scc2;
    double evx,evy,epx,epy,ev2x,ev2y,dot1,dot2,v1x,v1y,v2x,v2y;
    double sqrtv1,sqrtv2,sqrtep,sqrtev,sqrtev2;

    //////////////////////
    //get the right v1,v2
    //////////////////////
    e=object->elist.edges[boundary_edgelist[edge]];

    v1=object->vlist.verts[e->verts[0]->index];
    v2=object->vlist.verts[e->verts[1]->index];

    //getting v3
    //get another vertex in this face
    if(IsInBoundary(e->tris[0]->index))
        t=object->tlist.tris[e->tris[0]->index];
    else
        t=object->tlist.tris[e->tris[1]->index];
    for(int findV=0;findV<3;findV++)
    {
        if(t->verts[findV]->index!=v1->index && t->verts[findV]->index!=v2->index)
        {
            v3=object->vlist.verts[t->verts[findV]->index];
            break;
        }
    }
    ///////////////////////
    //compute edge status
    ///////////////////////
    //compute ev
    evx=v2->x	-	v1->x;
    evy=v2->y	-	v1->y;

    sqrtev=sqrt(evx*evx+evy*evy);
    evx=evx/sqrtev;
    evy=evy/sqrtev;

    //compute ep
    epx=-evy;
    epy=evx;

    //compute ev2
    ev2x=v3->x	-	v1->x;
    ev2y=v3->y	-	v1->y;

    sqrtev2=sqrt(ev2x*ev2x+ev2y*ev2y);
    ev2x=ev2x/sqrtev2;
    ev2y=ev2y/sqrtev2;



    //compute ep*ev2
    if(epx*ev2x+epy*ev2y>0) /*> ???????*/
    {
        epx=-epx;
        epy=-epy;
    }


    //compute dot1, dot2
    v1x=v1->g_vec.entry[0];
    v1y=v1->g_vec.entry[1];
    v2x=v2->g_vec.entry[0];
    v2y=v2->g_vec.entry[1];


    sqrtv1=sqrt(v1x*v1x+v1y*v1y);
    v1x=v1x/sqrtv1;
    v1y=v1y/sqrtv1;

    sqrtv2=sqrt(v2x*v2x+v2y*v2y);
    v2x=v2x/sqrtv2;
    v2y=v2y/sqrtv2;


    dot1=v1x*epx   +   v1y*epy;//dot1=v1.x*ep.x+v1.y*ep.y
    dot2=v2x*epx   +   v2y*epy;//dot2=v2.x*ep.x+v2.y*ep.y

    //get current edge characteristic
    if(dot1*dot2<0)
        return MIX;
    else
    {
        if(dot1>0)
            return EXIT;
        else
            return ENTRANCE;
    }
}

int ConleyIndex::GetDir_Tau(int edge)
{
    Edge *e;
    int scc1, scc2;
    Triangle *t;
    int scc_edge;
    int test_node;
    int test_scc;
    int test_morse;

    e=object->elist.edges[boundary_edgelist[edge]];
    //on the boundary
    if(e->tris[1]== nullptr)
        return GetDir_Geo(edge);

    //get outer node
    //which one is in the scc
    bool isIn0,isIn1;
    isIn0=IsInBoundary(e->tris[0]->index);
    isIn1=IsInBoundary(e->tris[1]->index);

    if(isIn0)
        t=object->tlist.tris[e->tris[1]->index];
    else
        t=object->tlist.tris[e->tris[0]->index];

    //get at least one edge of sccnode in scc
    int eN;
    int nE;
    //if(isLocal)
    //{
    //	if(IsInLocal(t->index))//in local just do the normal job
    //	{
    //		return GetDir_Geo(edge);

    //		nE=local_decomp->dg->nlist->dirnodes[t->local_index]->nedges;
    //		for(eN=0;eN<nE;eN++)
    //		{
    //			scc_edge=local_decomp->dg->nlist->dirnodes[t->local_index]->edges[eN];
    //			test_node=local_decomp->dg->elist->edges[scc_edge]->node_index1;
    //			test_scc=local_decomp->dg->nlist->dirnodes[test_node]->sscomp_index;

    //			if(test_scc==cur_scc)break;
    //		}
    //		if(eN==local_decomp->dg->nlist->dirnodes[t->local_index]->nedges)return ENTRANCE;
    //		else return EXIT;
    //	}
    //	else//outside the local,compute global
    //	{
    //		//return GetDir_Geo(edge);
    //		nE=morse_decomp->dg->nlist->dirnodes[t->index]->nedges;
    //		for(eN=0;eN<nE;eN++)
    //		{
    //			scc_edge=morse_decomp->dg->nlist->dirnodes[t->index]->edges[eN];
    //			test_node=morse_decomp->dg->elist->edges[scc_edge]->node_index1;
    //			test_scc=morse_decomp->dg->nlist->dirnodes[test_node]->sscomp_index;
    //			test_morse=morse_decomp->scclist->scccomponents[test_scc]->node_index;

    //			if(test_morse==refined_morse)break;
    //		}
    //	if(eN==morse_decomp->dg->nlist->dirnodes[t->index]->nedges)return ENTRANCE;
    //	else return EXIT;
    //	}

    //}
    //else
    {
        nE=morse_decomp->dg->nlist->dirnodes[t->index]->nedges;

        int inE=0;
        for(eN=0;eN<nE;eN++)
        {
            scc_edge=morse_decomp->dg->nlist->dirnodes[t->index]->edges[eN];
            if(morse_decomp->dg->elist->edges[scc_edge]->node_index1==t->index)
            {
                test_node=morse_decomp->dg->elist->edges[scc_edge]->node_index2;
                test_scc=morse_decomp->dg->nlist->dirnodes[test_node]->sscomp_index;
                if(test_scc==cur_scc)
                    inE++;
            }

            else
            {
                test_node=morse_decomp->dg->elist->edges[scc_edge]->node_index1;
                test_scc=morse_decomp->dg->nlist->dirnodes[test_node]->sscomp_index;
                if(test_scc==cur_scc)
                {

                    return EXIT;
                }

            }
        }



        if(inE==0)return GetDir_Geo(edge);
        else return ENTRANCE;

    }

}
void ConleyIndex::computeBetti(void)
{
    //if(pro[cur_scc].classification==3)//saddle or trivial
    //{
    //	pro[cur_scc].conleyB0=0;
    //	pro[cur_scc].conleyB2=0;
    //}
    //if(pro[cur_scc].classification==1)//source
    //{
    //	pro[cur_scc].conleyB0=0;
    //	pro[cur_scc].conleyB2=1;
    //}
    //if(pro[cur_scc].classification==2)//sink
    //{
    //	pro[cur_scc].conleyB0=1;
    //	pro[cur_scc].conleyB2=0;
    //}
    //pro[cur_scc].conleyB1=pro[cur_scc].conleyB0+pro[cur_scc].conleyB2+pro[cur_scc].XL-pro[cur_scc].XM;
    //
    //if(pro[cur_scc].classification==3 && pro[cur_scc].conleyB1==0)//trivial
    //	pro[cur_scc].classification=0;


    if(classification==2)//saddle or trivial
    {
        conleyB0=0;
        conleyB2=0;
    }
    if(classification==0)//source
    {
        conleyB0=0;
        conleyB2=1;
    }
    if(classification==1)//sink
    {
        conleyB0=1;
        conleyB2=0;
    }
    conleyB1=conleyB0+conleyB2+XL-XM;

    if (conleyB1<0) conleyB1 = 0;  // added by Guoning on 07/15/2010 to correct the bug

    //if(pro[cur_scc].classification==3 && pro[cur_scc].conleyB1==0)//trivial
    //	pro[cur_scc].classification=0;
}

bool ConleyIndex::IsInScc(int node)
{
    for(int i=0;i<morse_nlist.size();i++)
    {
        if(node==morse_nlist[i])
            return true;
    }
    return false;
}

bool ConleyIndex::IsInBoundary(int node)
{
    for(int i=0;i<boundary_tlist.size();i++)
    {
        if(node==boundary_tlist[i])
            return true;
    }
    return false;
}

bool ConleyIndex::IsInLocal(int node)
{

    for(int i=0;i<local_decomp->dg->nlist->ndirnodes;i++)
    {
        if(node==local_decomp->dg->nlist->dirnodes[i]->global_index)
            return true;
    }
    return false;
}

//for auto refine
//int ConleyIndex::select_morse(double min_priority)
//{
//	double c1,c2,c3;
//	c1=0.2;
//	c2=0.5;
//	c3=0.3;
//
//	for(int i=0;i<mcg->cur_mcgnode_index;i++)
//	{
//
//
//		//pri[i].priority=	(pri[i].tris_num)/(pri[i].B1+1)	;
//
//		pri[i].priority=	c1*pri[i].B1	+	c2*pri[i].tris_num	+c3*pri[i].variance_vector;
//
//		//pri[i].priority=(pri[i].variance_vector+1)*(pri[i].B1+1)*pri[i].tris_num;
//		if(pri[i].priority<min_priority)pri[i].out=true;
//	}
//	double max_pri=-1;
//	int first_mID;
//
//	for(int i=0;i<mcg->cur_mcgnode_index;i++)
//	{
//		bool discard=false;
//		int scc_id=0;
//		//for(int j=0;j<Object.nfaces;j++)
//		for(int sccN=0;sccN<morse_decomp->scclist->nsccs;sccN++)
//		{
//			//if morseID of the node == the current morse set, set discard
//			if(morse_decomp->scclist->scccomponents[sccN]->node_index==i)
//			{
//				scc_id = sccN;
//				for(int node=0;node<morse_decomp->scclist->scccomponents[sccN]->nnodes;node++)
//				{
//					if(object->tlist.tris[morse_decomp->scclist->scccomponents[sccN]->nodes[node]]->exclude)
//					{discard=true;
//					break;}
//				}
//				break;
//			}
//
//
//
//			//if(scclist.scccomponents[sccnodes[j].sscomp_index].node_index==i && Object.flist[j]->discard==true)
//			//{
//			//	//pri[i].out==true;//need to discard it
//			//	//funcP<<"morseID:"<<i<<" tri#:"<<j<<endl;
//			//	discard=true;
//			//	//break;
//			//}
//		}
//
//		//mcgN<<"scclist.scccomponents[scc_id].num_nodes:"<<scclist.scccomponents[scc_id].num_nodes<<" pri["<<i<<"].tris_num:"<<pri[i].tris_num<<endl;
//		if(discard)continue;
//		if(pri[i].out)continue;
//		if(morse_decomp->scclist->scccomponents[scc_id]->nnodes<5) continue;
//
//		if(pri[i].priority>max_pri)
//		{
//			max_pri=pri[i].priority;
//			first_mID=i;
//		}
//		//funcP<<"first_mID:"<<first_mID<<endl;
//		//funcP<<"MorseID:"<<i<<" priority: "<<pri[i].priority<<endl;
//	}
//
//	if(max_pri==-1)return -1;
//
//	return first_mID;
//}

double ConleyIndex::compute_variance_vector(void)
{
    Vertex  *v[3];
    double vx[3], vy[3],vz[3];
    double cur_vec;
    double max_vec=0, min_vec=10;
    double vec_var;
    double vecx, vecy,vecz;
    for(int i=0;i<morse_nlist.size();i++)
    {

        for(int j=0;j<3;j++)
        {
            v[j]=object->vlist.verts[object->tlist.tris[morse_nlist[i]]->verts[j]->index];
            vx[j]=1000*v[j]->g_vec.entry[0];
            vy[j]=1000*v[j]->g_vec.entry[1];
            vz[j]=1000*v[j]->g_vec.entry[2];
        }
        vecx=(vx[0]+vx[1]+vx[2])/3;
        vecy=(vy[0]+vy[1]+vy[2])/3;
        vecz=(vz[0]+vz[1]+vz[2])/3;
        cur_vec=sqrt(vecx*vecx+vecy*vecy+vecz*vecz);
        if(cur_vec>max_vec)max_vec=cur_vec;
        if(cur_vec<min_vec)min_vec=cur_vec;
    }
    vec_var=max_vec/min_vec;

    return vec_var;
}



int ConleyIndex::GetMCGType(int scc)
{
    //initialize
    build_morse_nlist(scc,0);
    cur_scc=scc;

    //compute X(M),X(L),Conley Index
    XM=computeXM();
    computeXL();
    computeBetti();

    return classification;
}

int* ConleyIndex::GetMCGBoundary(int & bn)
{

    bn=boundary_edgelist.size();
    int *t=new int[bn];


    for(int i=0;i<bn;i++)
        t[i]=boundary_edgelist[i];
    return t;
}

double ConleyIndex::GetMCGPriority(void)
{
    //initialize
    //weight
    /*double c1,c2,c3;
    c1=0.2;
    c2=0.5;
    c3=0.3;*/
    double priority=0;

    //tris_num
    int tris_num=morse_nlist.size();
    //variance_vector
    //double variance_vector=compute_variance_vector();

    //priority
    //1.
    if (conleyB0==0 && conleyB1==0 && conleyB2==0)
        priority = tris_num*100;
    else
        priority=tris_num*(conleyB1+1);
    //2.
    //priority=c1*conleyB1+c2*tris_num+c3*variance_vector;
    //3.
    //priority=(variance_vector+1)*(ConleyB1+1)*tris_num;
    //4.





    return priority;
}

void ConleyIndex::GetMCGConley(int  conley[])
{
    conley[0]=conleyB0;
    conley[1]=conleyB1;
    conley[2]=conleyB2;
}
