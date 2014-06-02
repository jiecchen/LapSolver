// calculate the probability that an edge in a graph will by chosen by randtree for spanning tree
// outputs in clear text to stdout prob of each edge in binijv (matched by line)
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

#include "components.h"
//#include "weightedselection.h"
#include "graph.h"
//#include "priorityq.h"

const char* DOC_STRING = "run like ./randedgeprob GRAPH_IN.BINIJV";
const int NUM_ITERATIONS = 1000; // how many samples we take

// return numbers 0 ... n-1 in random order
int* randPerm(int n)
{
    int* out = malloc(sizeof(int) * n);
    int ctr = 0;
    for (ctr = 0; ctr < n; ctr++)
	out[ctr] = ctr;
    for (ctr = 0; ctr < n; ctr++)
    {
	int tmp = out[ctr];
	int numLeft = n - ctr;
	int randPick = (double)rand()/(double)RAND_MAX * numLeft;
	out[ctr] = out[ctr + randPick];
	out[ctr + randPick] = tmp;
    }
    return out;
}

int main(int argc, char* argv[])
{
    srand(time(NULL));
    if (argc != 2)
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
    int numRuns = 0;
    int* countedEdges = malloc(sizeof(int) * ijv->nnz);  // where we keep track of which edges get used
    int ctr = 0;
    for (ctr = 0; ctr < ijv->nnz; ctr++)
	countedEdges[ctr] = 0;
    
    while (numRuns++ < NUM_ITERATIONS)
    {
	Components* comps = initGraphComponents(ijv->n);
	int* shuffledEdges = randPerm(ijv->nnz);   // order we wil go through edges
	//ijvType* treeOut = makeIJV(ijv->n - 1);
	//treeOut->n = ijv->n;
	//int treeEdgeCtr = 0;
	
	//go through all the edges always taking an edge that does not create a cycle
	for (ctr = 0; ctr < ijv->nnz; ctr++)
	{
	    int u = ijv->i[shuffledEdges[ctr]];
	    int v = ijv->j[shuffledEdges[ctr]];
	    int uComp = comps->vertComps[u];
	    int vComp = comps->vertComps[v];
	    if (uComp == vComp)
		continue;
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
	    
	    // note we used the edge
	    countedEdges[shuffledEdges[ctr]]++;
	}
	free(shuffledEdges);
	freeComponents(comps);
    }

    // output results
    for (ctr = 0; ctr < ijv->nnz; ctr++)
	printf("%f\n", (double)countedEdges[ctr]/(double)numRuns);

    freeIJV(ijv);
    free(countedEdges);
    //printf("exiting at %lf seconds\n", (double)clock()/(double)CLOCKS_PER_SEC);
    
    return 0;
}
