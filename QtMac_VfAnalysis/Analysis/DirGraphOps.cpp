/*
This file contains routines for Directed graph operations, such DFS, finding SCCs

Created and Modified by Guoning Chen
        copyright @2007
    */

#include "VField.h"

extern const int DebugOn;

extern Polyhedron *object;

void DirGraph::DFS(int inverse)
{
    int i;
    int time = 0;
    num_sccomps = 0;

    //initialization
    for(i = 0; i < nlist->ndirnodes; i++)
    {
        nlist->dirnodes[i]->visited = 0;
        nlist->dirnodes[i]->sscomp_index = 0;
        nlist->dirnodes[i]->parent = i;
        nlist->dirnodes[i]->levels = 0;
    }

    //
    for(i = 0; i < nlist->ndirnodes; i++)
    {
        if(inverse == 0)
        {
            if(nlist->dirnodes[i]->visited == 0)
            {
                //DFS_visit(i, time, inverse);
                DFS_visit_nonrecursive(i, num_sccomps, inverse);
            }
        }

        else
        {
            if(nlist->dirnodes[cur_nodes_order[i]]->visited == 0)
            {
                //DFS_visit(cur_nodes_order[i], time, inverse);
                DFS_visit_nonrecursive(cur_nodes_order[i], num_sccomps, inverse);
                num_sccomps ++;           //for each subtree, increase the index
            }
        }
    }
}


void DirGraph::DFS_visit(int u, int &time, int inverse)
{
    int i;
    Graph_Edge *cur_e;
    int adj_node;

    nlist->dirnodes[u]->visited = 1;
    time ++;

    //nodes[u].levels = time;

    for(i = 0; i < nlist->dirnodes[u]->nedges; i++)
    {
        cur_e = elist->edges[nlist->dirnodes[u]->edges[i]];

        adj_node = cur_e->node_index2;

        if(adj_node == u)
            continue;

        if(nlist->dirnodes[adj_node]->visited == 0)
        {
            nlist->dirnodes[adj_node]->sscomp_index = num_sccomps;  //mark 'adj_node' as the same sub-tree (the same scc)
                //as u
            nlist->dirnodes[adj_node]->parent = u;                  //set the parent
            DFS_visit(adj_node, time, inverse); //recursive to find out all the nodes in the subtree
        }
    }


    nlist->dirnodes[u]->visited = 2;
    time ++;
    nlist->dirnodes[u]->levels = time;
    nlist->dirnodes[u]->sscomp_index = num_sccomps;

    ///Add the node into the list instead of sorting it later
    if(cur_sccnode_index >= 0)
    {
        cur_nodes_order[cur_sccnode_index] = u;
        cur_sccnode_index --;
    }
}



/*
The non-recursive DFS
*/
/*
Here we need a stack to keep track of the searching
*/



void DirGraph::add_to_dfsstack(int elem)
{
    //if(top_dfsstack == curMaxDFSStack)
    //{
    //	/*extend it*/
    //	int* temp_stack;
    //	temp_stack=(int*)malloc(sizeof(int)*top_dfsstack);
    //
    //	for(int i=0;i<top_dfsstack;i++)
    //		temp_stack[i]=dfs_stack[i];

    //	//if(dfs_stack != NULL)
    //	free(dfs_stack);

    //	//free(dfs_stack);

    //	dfs_stack=(int*)malloc(sizeof(int)*(curMaxDFSStack+500));

    //	//dfs_stack = (int*)realloc(dfs_stack, sizeof(int)*(curMaxDFSStack+500));

    //	if(dfs_stack == NULL)
    //	{
    //		char rout[256], var[256];
    //		sprintf(rout, "%s", "DirGraph::add_to_dfsstack");
    //		sprintf(var, "%s", "dfs_stack");

    //		write_mem_error(rout, var, 0);
    //		curMaxDFSStack = 0;
    //		exit(-1);
    //	}

    //	//curMaxDFSStack+=800;
    //	curMaxDFSStack+=500;

    //	for(int i=0;i<top_dfsstack;i++)
    //		dfs_stack[i]=temp_stack[i];
    //	free(temp_stack);
    //}


    /*top_dfsstack++;
    dfs_stack[top_dfsstack] = elem;*/
    dfs_stack.push(elem);
}



/*non-recursive DFS*/
void DirGraph::DFS_visit_nonrecursive(int u, int scc_index, int inverse)
{
    int i;
    Graph_Edge *cur_e;
    int cur_node, adj_node;
    int nchildren = 0;

    //nodes[u].visited = 1;

    /*initialize the stack*/
    /*if(dfs_stack != NULL)
        free(dfs_stack);*/
    while(!dfs_stack.empty())
        dfs_stack.pop();

    /*curMaxDFSStack = nlist->ndirnodes*2;
    dfs_stack = (int*)malloc(sizeof(int)*curMaxDFSStack);
    if(dfs_stack == NULL)
    {
        char rout[256], var[256];
        sprintf(rout, "%s", "DirGraph::DFS_visit_nonrecursive");
        sprintf(var, "%s", "dfs_stack");

    write_mem_error(rout, var, 0);
    curMaxDFSStack = 0;
    exit(-1);
}*/


    //dfs_stack[0] = u;

    dfs_stack.push(u);

//top_dfsstack = 0;

//nodes[u].levels = time;


//while(top_dfsstack >= 0)/*keep searching until the stack is empty*/
while(!dfs_stack.empty())
{
    /*pop up the last node*/
    //cur_node = dfs_stack[top_dfsstack];
    cur_node=dfs_stack.top();

    if(nlist->dirnodes[cur_node]->visited == 2)
    {
        /*remove it from the stack*/
        //top_dfsstack --;
        dfs_stack.pop();
        continue;
    }

    nlist->dirnodes[cur_node]->visited = 1; /*mark it as 'visited' but not 'finished' yet*/

    nchildren = 0;
    //for(i = nodes[cur_node].nedges-1; i > 0; i--)
    for(i = 0; i < nlist->dirnodes[cur_node]->nedges; i++)
    {
        cur_e = elist->edges[nlist->dirnodes[cur_node]->edges[i]];

        adj_node = cur_e->node_index2;

        if(adj_node == cur_node)/*consider outgoing edges only*/
            continue;

        if(nlist->dirnodes[adj_node]->visited == 0)
        {
            //nodes[adj_node].sscomp_index = scc_index;  //mark 'adj_node' as the same sub-tree (the same scc)
            nlist->dirnodes[adj_node]->parent = cur_node;         //set the parent

            /*add to the stack*/
            add_to_dfsstack(adj_node);
            nchildren ++;
        }
    }

    if(nchildren == 0)/*this is a leaf or no more 'unvisited' children*/
    {
        nlist->dirnodes[cur_node]->visited = 2;  /*the searching is 'finished' for this node*/
        nlist->dirnodes[cur_node]->sscomp_index = scc_index;

        /*put it to the cur_nodes_order array*/
        if(cur_sccnode_index >= 0)
        {
            cur_nodes_order[cur_sccnode_index] = cur_node;
            cur_sccnode_index --;
        }

        //else{
        //	if(inverse == 0)
        //	{
        //		FILE *fp= fopen("nodesearch_error.txt", "a");
        //		fprintf(fp, "cur_nodes_order = %d \n", cur_sccnode_index);
        //		fprintf(fp, "current node: %d\n", cur_node);
        //		fprintf(fp, "the root is: %d\n", u);
        //		fclose(fp);
        //	}
        //}

        /*remove it from the stack*/
        //top_dfsstack --;
        dfs_stack.pop();
    }
}

//free(dfs_stack);
//dfs_stack = NULL;
while(!dfs_stack.empty())
    dfs_stack.pop();
}



void DirGraph::find_SCCS()
{
    /*first allocate the space for cur_nodes_order*/
    cur_nodes_order = new int[nlist->ndirnodes];

    int i;
    for(i = 0; i < nlist->ndirnodes; i++)
        cur_nodes_order[i] = i;

    cur_sccnode_index = nlist->ndirnodes-1;

    DFS(0);

    reverse_edges();

    DFS(1);

    reverse_edges();

    delete [] cur_nodes_order;
}

void DirGraph::reverse_edges()
{
    int i, temp;

    for(i = 0; i < elist->nedges; i++)
    {
        temp = elist->edges[i]->node_index1;
        elist->edges[i]->node_index1 = elist->edges[i]->node_index2;
        elist->edges[i]->node_index2 = temp;
    }
}



/*-----------------------------------------------------------
* Add the edge "edgeindex" to the edge list of "node"
-----------------------------------------------------------*/
void DirGraph::add_edge_to_node(int node, int edgeindex)
{
    //we need to first extend the space of the edge list
    if(nlist->dirnodes[node]->nedges == 0)
    {
        if(nlist->dirnodes[node]->edges != NULL)
            free(nlist->dirnodes[node]->edges);
        nlist->dirnodes[node]->edges = (int*)malloc(sizeof(int)*1);

        if(nlist->dirnodes[node]->edges == NULL)
        {
            char rout[256], var[256];
            sprintf(rout, "%s", "DirGraph::add_to_edge_to_node");
            sprintf(var, "%s", "nlist->dirnodes[]->edges");

            //write_mem_error(rout, var, 1);
            exit(-1);
        }
    }

    else
    {
        int *temp = nlist->dirnodes[node]->edges;

        nlist->dirnodes[node]->edges = (int*)malloc(sizeof(int)*(nlist->dirnodes[node]->nedges+1));
        if(nlist->dirnodes[node]->edges == NULL)
        {
            char rout[256], var[256];
            sprintf(rout, "%s", "DirGraph::add_to_edge_to_node");
            sprintf(var, "%s", "nlist->dirnodes[]->edges");

            //write_mem_error(rout, var, 1);
            exit(-1);
        }

        for(int i = 0; i < nlist->dirnodes[node]->nedges; i++)
        {
            nlist->dirnodes[node]->edges[i] = temp[i];
        }

        free(temp);
    }

    nlist->dirnodes[node]->edges[nlist->dirnodes[node]->nedges] = edgeindex;
    nlist->dirnodes[node]->nedges ++;
}


bool DirGraph::is_repeated_edge(int from, int to)
{
    int i;

    for(i = 0; i < elist->nedges; i++)
    {
        if(elist->edges[i]->node_index1 == from && elist->edges[i]->node_index2 == to)
            return true;

        //if(elist->edges[i]->node_index1 != from)
        //	continue;

        //if(elist->edges[i]->node_index2 == to)
        //	return true;

    }

    return false;
}


/*   The updated for the repeated_edge pending  */
bool DirGraph::is_repeated_edge_2(int from, int to)
{
    int i;
    Graph_Edge *e;

    for(i = 0; i < nlist->dirnodes[from]->nedges; i++)
    {
        e = elist->edges[nlist->dirnodes[from]->edges[i]];

        if(e->node_index1 != from)
            continue;

        if(e->node_index2 == to)
            return true;

    }
    return false;
}


void DirGraph::add_to_edgelist(int node1, int node2, int &cur_index)
{

    if(elist->nedges >= elist->curMaxNumGedges)
    {
        int oldnum = elist->curMaxNumGedges;
        if(!elist->extend(object->tlist.ntris/*10(int)elist->nedges/2*/))
        {
            char rout[256], var[256];
            sprintf(rout, "%s", "DirGraph::add_to_edgelist");
            sprintf(var, "%s", "elist->edges");

            //write_mem_error(rout, var, 1);
            exit(-1);
        }

        /*
        we need to add new space for the new allocate pointers
        */
        for(int i = oldnum; i < elist->curMaxNumGedges; i++)
        {
            elist->edges[i] = new Graph_Edge();
        }
    }


    elist->edges[elist->nedges]->node_index1 = node1;
    elist->edges[elist->nedges]->node_index2 = node2;
    elist->edges[elist->nedges]->edge_index = elist->nedges;
    elist->nedges++;

    cur_index = elist->nedges;

}


