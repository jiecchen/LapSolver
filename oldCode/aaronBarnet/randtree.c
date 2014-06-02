// create a random tree from a graph (probability of picking edge inv. proportional to weight) 
// (not really random spanning tree: try cyclepoprandtree)
// if SIM_COMPLETE_G==true: randly select pairs of verts for edges
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

#include "components.h"
#include "weightedselection.h"
#include "graph.h"
//#include "priorityq.h"

const char* DOC_STRING = "run like ./randtree GRAPH_IN.BINIJV TREE_OUT.PARRAY\n";
const bool SIM_COMPLETE_G = false;

typedef struct List_s
{
    int u;     // keep u smaller than v
    int v;
    struct List_s* next;
} List;

typedef struct
{
    int n;     // num nodes in complete graph
    List* list;
} SimCompleteData;

SimCompleteData* newSimCompleteData(int n)
{
    SimCompleteData* scd = malloc(sizeof(SimCompleteData));
    scd->n = n;
    scd->list = NULL;
    return scd;
}

void freeSimCompleteData(SimCompleteData* scd)
{
    List* curList = scd->list;
    while (curList)
    {
	List* tmp = curList->next;
	free(curList);
	curList = tmp;
    }
}

bool scdFind(List* list, SimCompleteData* scd)
{
    List* curList = scd->list;
    while (curList)
    {
	if (list->u == curList->u
	    && list->v == curList->v)
	    return true;
	curList = curList->next;
    }
    return false;
}

void scdSelect(SimCompleteData* scd, int* u, int* v)
{
    while (1)
    {
	int upick = (double)rand()/(double)RAND_MAX * scd->n;
	int vpick = (double)rand()/(double)RAND_MAX * scd->n;
	if (upick == vpick)
	    continue;
	if (upick < vpick)
	{
	    *u = upick;
	    *v = vpick;
	}
	else
	{
	    *u = vpick;
	    *v = upick;
	}
	break;
    }
}

/*     List* newList = malloc(sizeof(List)); */
/*     while (1) */
/*     { */
/* 	int u = (double)rand()/(double)RAND_MAX * scd->n; */
/* 	int v = (double)rand()/(double)RAND_MAX * scd->n; */
/* 	if (u == v) */
/* 	    continue; */
/* 	if (u < v) */
/* 	{ */
/* 	    newList->u = u; */
/* 	    newList->v = v; */
/* 	} */
/* 	else */
/* 	{ */
/* 	    newList->u = v; */
/* 	    newList->v = u; */
/* 	} */
/* 	if (scdFind(newList, scd)) */
/* 	    continue; */
/* 	break; */
/*     } */
/*     *u = newList->u; */
/*     *v = newList->v; */
/*     newList->next = scd->list; */
/*     scd->list = newList; */
/* } */

// return a an array of the number 0 -> n-1 in a random order
// freshly malloc'd so must be freed!
/* int* randPerm(int n) */
/* { */
/*     PriorityQ* pq = newPriorityQ(n); */
/*     int* out = malloc(sizeof(int) * n); */
/*     int ctr = 0; */
/*     for (ctr = 0; ctr < n; ctr++) */
/*     { */
/* 	pqInsert(pq, ctr, (double)rand()); */
/*     } */
/*     for (ctr = 0; ctr < n; ctr++) */
/*     { */
/* 	int num = pqTop(pq); */
/* 	pqPop(pq); */
/* 	out[ctr] = num; */
/*     } */
    
/*     freePriorityQ(pq); */

/*     return out; */
/* } */



/* typedef struct */
/* { */
/*     int nnz; */
/*     int* i; */
/*     int* j; */
/* } Edges; */

/* Edges* newEdges(int nnz) */
/* { */
/*     Edges* out = malloc(sizeof(Edges)); */
/*     out->i = (int*)malloc(sizeof(int) * nnz); */
/*     out->j = (int*)malloc(sizeof(int) * nnz); */
/*     return out; */
/* } */

/* void freeEdges(Edges* in) */
/* { */
/*     free(in->i); */
/*     free(in->j); */
/*     free(in); */
/* } */

/* // could just use ijv used to make graph */
/* Edges* enumerateEdges(myGraph* graph) */
/* { */
    

/* } */

int main(int argc, char* argv[])
{
    srand(time(NULL));
    if (argc != 3)
	{
	    fprintf(stderr, "%s", DOC_STRING);
	    exit(1);
	}

    FILE* fpIn = NULL;
    ijvType* ijv = NULL;
    //myGraph* graph = NULL;
    if ((fpIn = fopen(argv[1], "rt")) == NULL)
        {
	    fprintf(stderr, "Error: could not read %s\n", argv[1]);
	    exit(1);
	}
    ijv = binReadIJV(fpIn);
    fclose(fpIn);
    //graph = ijv2graph(ijv);
    //pArray* parray = newPArray(graph->n);
    Components* comps = initGraphComponents(ijv->n);
    //int* shuffledEdges = randPerm(ijv->nnz);   // order we wil go through edges
    // create weightelements for random selection
    WeightedElements* we = newWeightedElements(ijv->nnz);   // the elements are the indexes to the edges in ijv
    int ctr = 0;
    for (ctr = 0; ctr < ijv->nnz; ctr++)
    {
	weInsert(we, ctr, 1.0/ijv->v[ctr]);
    }

    ijvType* treeOut = makeIJV(ijv->n - 1);
    treeOut->n = ijv->n;
    int treeEdgeCtr = 0;
    
    SimCompleteData* scd = NULL;
    if (SIM_COMPLETE_G)
    {
	scd = newSimCompleteData(ijv->n);
    }

    

    // randomly select edges until we have a spanning tree (only one comp)
    int numComponents = ijv->n;
    while (numComponents > 1)
    {
	int u = -1;
	int v = -1;
	double weight = 0.0;
	if (SIM_COMPLETE_G)
	{
	    scdSelect(scd, &u, &v);
	    weight = 0.0;
	}
	else
	{
	    int randEdgeIdx = weSelect(we);
	    u = ijv->i[randEdgeIdx];
	    v = ijv->j[randEdgeIdx];
	    weight = ijv->v[randEdgeIdx];
	}

	int uComp = comps->vertComps[u];
	int vComp = comps->vertComps[v];
	if (uComp == vComp)
	{
	    continue;
	}
	// let u be the smaller component 
	if (comps->compSizes[vComp] < comps->compSizes[uComp])
	{
	    int tmp = u;
	    u = v;
	    v = tmp;
	    tmp = uComp;
	    uComp = vComp;
	    vComp = tmp;
	}

	// make uComp part of vComp
	CompCell* curUcell = comps->comps[uComp];
	while (1)
	{
	    comps->vertComps[curUcell->vert] = vComp;
	    if (curUcell->next)
		curUcell = curUcell->next;
	    else
		break;
	}
	curUcell->next = comps->comps[vComp];
	comps->comps[vComp] = comps->comps[uComp];
	comps->comps[uComp] = NULL;
	comps->compSizes[vComp] += comps->compSizes[uComp];
	comps->compSizes[uComp] = 0;

	// add edge to tree
	treeOut->i[treeEdgeCtr] = u;
	treeOut->j[treeEdgeCtr] = v;
	treeOut->v[treeEdgeCtr++] = weight;
	numComponents--;
    }

    //go through all the edges always taking an edge that does not create a cycle
/*     for (ctr = 0; ctr < ijv->nnz; ctr++) */
/*     { */
/* 	int u = ijv->i[shuffledEdges[ctr]]; */
/* 	int v = ijv->j[shuffledEdges[ctr]]; */
/* 	int uComp = comps->vertComps[u]; */
/* 	int vComp = comps->vertComps[v]; */
/* 	if (uComp == vComp) */
/* 	    continue; */
/* 	// let u be the smaller component  */
/* 	if (comps->compSizes[vComp] < comps->compSizes[uComp]) */
/* 	{ */
/* 	    int tmp = u; */
/* 	    u = v; */
/* 	    v = tmp; */
/* 	    tmp = uComp; */
/* 	    uComp = vComp; */
/* 	    vComp = tmp; */
/* 	} */

/* 	// make uComp part of vComp */
/* 	Cell* curUcell = comps->comps[uComp]; */
/* 	while (1) */
/* 	{ */
/* 	    comps->vertComps[curUcell->vert] = vComp; */
/* 	    if (curUcell->next) */
/* 		curUcell = curUcell->next; */
/* 	    else */
/* 		break; */
/* 	} */
/* 	curUcell->next = comps->comps[vComp]; */
/* 	comps->comps[vComp] = comps->comps[uComp]; */
/* 	comps->comps[uComp] = NULL; */
/* 	comps->compSizes[vComp] += comps->compSizes[uComp]; */
/* 	comps->compSizes[uComp] = 0; */

/* 	// add edge to tree */
/* 	treeOut->i[treeEdgeCtr] = u; */
/* 	treeOut->j[treeEdgeCtr] = v; */
/* 	treeOut->v[treeEdgeCtr++] = ijv->v[shuffledEdges[ctr]]; */
/*     } */

    myGraph* treeOutGraph = ijv2graph(treeOut);
    pArray* parray = graph2pArray(treeOutGraph);
    binWritePArray(parray, argv[2]);
//    binWriteIJV(treeOut, argv[2]);
    
    //free(shuffledEdges);
    freeWeightedElements(we);
    freeComponents(comps);
    //freePArray(parray);
    //freeGraph(graph);
    freeGraph(treeOutGraph);
    freePArray(parray);
    freeIJV(ijv);
    freeIJV(treeOut);
    if (scd)
	freeSimCompleteData(scd);

    printf("exiting at %lf seconds\n", (double)clock()/(double)CLOCKS_PER_SEC);
    
    return 0;
}
