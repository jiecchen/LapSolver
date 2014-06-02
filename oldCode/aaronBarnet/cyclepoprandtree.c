// This program takes in a myGraph graph and outputs a random spanning tree
// By random, I mean that each possible tree will occur with probability proportional to the product of all the weights of the used edges
// Based on Wilson's paper

// pick root node for tree (with prob proportional to degree of node)
// then randomly generate tree by taking random node and creating random path until joined up with tree
// anytime a cycle appears, just go back to right before the area that was partof cycle

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "graph.h"
#include "cyclepoprandtreemodule.h"

const char* DOC_STRING = "run like ./cyclepoprandtree GRAPH_IN.BINIJV TREE_OUT.PARRAY\n";

/* // pick a node at random such that the prob of any node being picked is proportional to its degree */
/* // implement as linear time function */
/* int pickRoot(ijvType* ijv) */
/* { */
/*     // we will pick a random edge (with prob prop. to weight) */
/*     double totalWeights = 0.0; */
/*     int ctr = 0; */
/*     for (ctr = 0; ctr < ijv->nnz; ctr++) */
/*     { */
/* 	totalWeights += ijv->v[ctr]; */
/*     } */
/*     double randWeight = (double)rand()/(double)(RAND_MAX) * totalWeights; */
/*     double weightCtr = 0.0; */
/*     int randEdge = 0; */
/*     for (ctr = 0; ctr < ijv->nnz; ctr++) */
/*     { */
/* 	weightCtr += ijv->v[ctr]; */
/* 	if (weightCtr >= randWeight) */
/* 	{ */
/* 	    randEdge = ctr; */
/* 	    break; */
/* 	} */
/*     } */
    
/*     // now just pick either i or j of randEdge; */
/*     int randVert = 0; */
/*     if ((double)rand()/(double)RAND_MAX < .5) */
/* 	randVert = ijv->i[randEdge]; */
/*     else */
/* 	randVert = ijv->j[randEdge]; */

/*     return randVert; */
/* } */

/* // vert is the node for which you want a nbr, */
/* // avoid is nbr that's not aloud (stopped using avoid) */
/* int pickNeighbor(int vert, myGraph* graph, int avoid) */
/* { */

/*     double totalWeights = 0.0; */
/*     int ctr; */
/*     for (ctr = 0; ctr < graph->deg[vert]; ctr++) */
/* 	totalWeights += graph->wts[vert][ctr]; */
/*     int randNbr = 0; */
/*     //while (randNbr == avoid) */
/*     //{ */
/* 	double randWeight = (double)rand()/(double)(RAND_MAX) * totalWeights; */
/* 	double weightCtr = 0.0; */
/* 	for (ctr = 0; ctr < graph->deg[vert]; ctr++) */
/* 	{ */
/* 	    weightCtr += graph->wts[vert][ctr]; */
/* 	    if (weightCtr >= randWeight) */
/* 	    { */
/* 		randNbr = graph->nbrs[vert][ctr]; */
/* 		break; */
/* 	    } */
	    
/* 	} */
/*     //} */
/*     return randNbr; */
/* } */

/* pArray* randRootTree(int root, myGraph* graph) */
/* { */
/*     pArray* parray = newPArray(graph->n); */
/*     bool* inTree = malloc(sizeof(bool) * graph->n); */
/*     bool* inPath = malloc(sizeof(bool) * graph->n); */
/*     int* stack = malloc(sizeof(int) * graph->n); */
/*     int stackTop = 0;   // always where to put next entry on stack (i.e. it == stack size) */

/*     int ctr = 0; */
/*     for (ctr = 0; ctr < graph->n; ctr++) */
/*     { */
/* 	inTree[ctr] = false; */
/* 	inPath[ctr] = false; */
/*     } */
/*     inTree[root] = true; */
/*     parray->array[root] = root; */
    
/*     int vertItr = 0; */
/*     for (vertItr = 0; vertItr < graph->n; vertItr++) */
/*     { */
/* 	int vert = vertItr; */
/* 	//printf("new path\n"); */
/* 	if (inTree[vert]) */
/* 	    continue; */

/* 	int ctr = 0; */
/* 	for (ctr = 0; ctr < graph->n; ctr++) */
/* 	{ */
/* 	    if (inPath[ctr]) */
/* 		printf("bad inPath\n"); */
/* 	} */

/* 	stack[stackTop++] = vert; */
/* 	inPath[vert] = true; */
/* 	// keep building up path until it reaches tree */
/* 	while (1) */
/* 	{ */
/* /\* 	    printf("vert: %d\n", vert); *\/ */
/* /\* 	    printf("stack: "); *\/ */
/* /\* 	    for (ctr = 0; ctr < stackTop; ctr++) *\/ */
/* /\* 	    { *\/ */
/* /\* 		printf("%d ", stack[ctr]); *\/ */
/* /\* 	    } *\/ */
/* /\* 	    printf("\n"); *\/ */
/* 	    // don't pick previous node in path */
/* 	    int nbr = 0; */
/* 	    if (stackTop > 1) */
/* 		nbr = pickNeighbor(vert, graph, stack[stackTop-2]); */
/* 	    else */
/* 		nbr = pickNeighbor(vert, graph, -1); */
		
/* 	    //printf("nbr: %d\n", nbr); */
/* 	   /\*  if (nbr >= graph->n) *\/ */
/* /\* 		fprintf(stderr, "node out of range\n"); *\/ */
/* /\* 	    if (stackTop >= graph->n) *\/ */
/* /\* 		fprintf(stderr, "stackTop out of range\n"); *\/ */
/* 	    if (inPath[nbr])      // we have cycle to pop */
/* 	    { */
/* 		//printf("found cycle\n"); */
/* 		// keeping popping off from stack until the first instance of nbr is top of stack */
/* 		/\* int stackItr = 0; *\/ */
/* /\* 		for (stackItr = 0; stackItr < stackTop; stackItr++) *\/ */
/* /\* 		{ *\/ */
/* /\* 		    if (stack[stackItr] == nbr) *\/ */
/* /\* 			printf("found nbr in stack\n"); *\/ */
/* /\* 		} *\/ */
/* /\* 		fprintf(stderr, "before while loop\n"); *\/ */
/* 		while (stack[stackTop - 1] != nbr) */
/* 		{ */
/* 		    int nodeToDelete = stack[--stackTop]; */
/* 		    inPath[nodeToDelete] = false; */
/* 		} */
/* 		vert = nbr; // (the old nbr); */
/* 		//fprintf(stderr, "after while loop\n"); */
/* 	    } */
/* 	    else if (inTree[nbr]) // we are done with path */
/* 	    { */
/* 		//printf("found tree\n"); */
/* 		int stackItr = 0; */
/* 		for (stackItr = 0; stackItr < stackTop - 1; stackItr++) */
/* 		{ */
/* 		    int nodeInTree = stack[stackItr]; */
/* 		    inPath[nodeInTree] = false; */
/* 		    inTree[nodeInTree] = true; */
/* 		    parray->array[nodeInTree] = stack[stackItr + 1]; */
/* 		} */
/* 		inPath[stack[stackTop - 1]] = false; */
/* 		inTree[stack[stackTop - 1]] = true; */
/* 		parray->array[stack[stackTop - 1]] = nbr; */
/* 		stackTop = 0; */
/* 		break; */
/* 	    } */
/* 	    else                  // nothing special */
/* 	    { */
/* 		//printf("nothing special\n"); */
/* 		stack[stackTop++] = nbr; */
/* 		inPath[nbr] = true; */
/* 		vert = nbr; */
/* 	    } */
/* 	} */

/*     } */
    
/*     free(inTree); */
/*     free(inPath); */
/*     free(stack); */

/*     return parray; */
/* } */


int main(int argc, char* argv[])
{
    srand(time(NULL));
    if (argc != 3)
    {
	fprintf(stderr, "%s", DOC_STRING);
	exit(1);
    }
    
    // get graph
    FILE* fpIn = NULL;
    ijvType* ijv = NULL;
    myGraph* graph = NULL;
    if ((fpIn = fopen(argv[1], "rt")) == NULL)
    {
	fprintf(stderr, "Error: could not read %s\n", argv[1]);
	exit(1);
    }
    ijv = binReadIJVcosts(fpIn);  // in this case we want costs not lenghts!
    fclose(fpIn);
    graph = ijv2graph(ijv);

    // choose root
    int root = pickRoot(ijv);
    
    // get tree
    pArray* parray = randRootTree(root, graph);
	
    // output tree
/*     printf("parray: "); */
/*     int ctr = 0; */
/*     for (ctr = 0; ctr < graph->n; ctr++) */
/* 	printf("%d ", parray->array[ctr]); */
    binWritePArray(parray, argv[2]);

    freeGraph(graph);
    freeIJV(ijv);
    freePArray(parray);

    return 0;
}
