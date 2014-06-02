// take a tree and replace some of the edge randomly. keep result if stretch is improved
// rinse and repeat.

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include "graph.h"
#include "weightedselection.h"
#include "components.h"
#include "calcstretchctree.h"

const char* DOC_STRING = "run like ./randimprove GRAPH_IN.BINIJV TREE_IN.PARRAY TREE_OUT.PARRAY\n";

const int NUM_ITERATIONS = 100;
const double REPLACE_FRACTION = .01;

// stretch is expected to already be known for tree at beginning.
// at end stretch will be updated
Tree* randImproveTree(Tree* tree, myGraph* graph, ijvType* ijv, double* stretch)
{
    ijvType* newTree = makeIJV(tree->n - 1);
    
    Components* comps = initGraphComponents(tree->n);
    int numComponents = ijv->n;
    bool* selected = malloc(sizeof(bool)* tree->n);
    int ctr = 0;
    for (ctr = 0; ctr < tree->n; ctr++)
	selected[ctr] = false;
    int newTreeEdgeCtr = 0;
    // first select which edges to keep
    
    while (newTreeEdgeCtr < (double)tree->n * REPLACE_FRACTION)
    {
	int randNode;
	while (true)
	{
	    randNode = (double)rand()/((double)RAND_MAX + .1)*tree->n;
	    if (!selected[randNode] && randNode != tree->root->v)
		break;
	}
	selected[randNode] = true;
	int u = tree->verts[randNode]->v;
	int v = tree->verts[randNode]->p->v;
	int uComp = comps->vertComps[u];
	int vComp = comps->vertComps[v];

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
	newTree->i[newTreeEdgeCtr] = u;
	newTree->j[newTreeEdgeCtr] = v;
	// avoiding using weights in ijv
	newTreeEdgeCtr++;
	numComponents--;
    }
    free(selected);
    
    //printf("after selecting edges to keep\n");

    // now add all allowable edges to weightedselection
    WeightedElements* we = newWeightedElements(ijv->nnz); // won't nec. fill it up completely!
    for (ctr = 0; ctr < ijv->nnz; ctr++)
    {
	int u = ijv->i[ctr];
	int v = ijv->j[ctr];
	if (comps->vertComps[u] != comps->vertComps[v])
	{
	    weInsert(we, ctr, 1.0/ijv->v[ctr]);
	}
    }
    we->n = we->cur;

    //printf("after filling up we\n");

    while (numComponents > 1)
    {
	int randEdgeIdx = weSelect(we);
	int u = ijv->i[randEdgeIdx];
	int v = ijv->j[randEdgeIdx];
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
	newTree->i[newTreeEdgeCtr] = u;
	newTree->j[newTreeEdgeCtr] = v;
	// avoid tree->v
	newTreeEdgeCtr++;
	numComponents--;
    }
    
    //printf("after finished picking new tree edges\n");
    newTree->nnz = newTreeEdgeCtr;
    newTree->n = graph->n;

    myGraph* newTreeGraph = ijv2graph(newTree);
    pArray* parray = graph2pArray(newTreeGraph);
    Tree* newTreeTree = pArray2tree(parray, graph);
    
    double newStretch = calcStretch(newTreeTree, graph);
    Tree* treeToReturn = NULL;
    if (newStretch < *stretch)
    {
	*stretch = newStretch;
	treeToReturn = newTreeTree;
    }
    else
    {
	freeTree(newTreeTree);
	treeToReturn = tree;
    }

    freeIJV(newTree);
    freeGraph(newTreeGraph);
    freeComponents(comps);
    freePArray(parray);
    //printf("returning from improvetree\n");
    return treeToReturn;
}


int main(int argc, char* argv[])
{
    srand(time(NULL));

    Tree* tree = NULL;
    Tree* improvedTree = NULL;
    
    if (argc != 4)
    {
	fprintf(stderr, "%s", DOC_STRING);
	exit(1);
    }
    
    FILE* fpIn = NULL;
    ijvType* ijv = NULL;
    myGraph* graph = NULL;
    if ((fpIn = fopen(argv[1], "rt")) == NULL)
    {
	fprintf(stderr, "Error: could not read %s\n", argv[1]);
	exit(1);
    }

    ijv = binReadIJV(fpIn);
    fclose(fpIn);
    graph = ijv2graph(ijv);
    pArray* parray = binReadPArray(argv[2]);
    tree = pArray2tree(parray, graph);
    
    //printf("before improvement\n");

    int iterations = 0;
    double stretch = calcStretch(tree, graph);
    improvedTree = randImproveTree(tree, graph, ijv, &stretch);
    iterations++;
    while (iterations++ < NUM_ITERATIONS)
    {
	Tree* tmp = randImproveTree(improvedTree, graph, ijv, &stretch);
	if (tmp != improvedTree && improvedTree != tree)
	    freeTree(improvedTree);
	improvedTree = tmp;
    }

    pArray* parrayOut = tree2pArray(improvedTree);
    binWritePArray(parrayOut, argv[3]);

    //printf("after output of parray\n");
//    binWritePArray(parray, const char* fileName); // vert ints must be zero based (this function will increment them)
    
    // free mem
    freePArray(parray);
    freePArray(parrayOut);
    freeIJV(ijv);
    freeGraph(graph);
    //printf("about to free trees\n");
    if (tree != improvedTree)
	freeTree(improvedTree);
    freeTree(tree);
    
    printf("exiting at %lf seconds\n", (double)clock()/(double)CLOCKS_PER_SEC);

    return 0;
}
