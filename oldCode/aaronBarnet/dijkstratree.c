// takes in a graph and outputs a tree (parray)
// take random node and run dijkstra. Tree is the paths
// dijkstra chooses.
#include <stdio.h>
#include <stdlib.h>
#include "priorityq.h"
#include "graph.h"

const char* DOC_STRING = "run like ./dijkstratree GRAPH.BINIJV TREEOUT.PARRAY\n";
enum COLORS {WHITE, GREY, BLACK}; // = {unvisited, in Q, processed}


int main(int argc, char* argv[])
{
    if (argc != 3)
    {
	fprintf(stderr, "%s", DOC_STRING);
	exit(1);
    }
    
    myGraph* graph = getGraph(argv[1]);
    if (graph->n < 1)
    {
	fprintf(stderr, "ERROR: graph has no vertices.\n");
	exit(2);
    }

    // data structs for dijkstra
    pArray* parray = newPArray(graph->n);
    double* dists = malloc(sizeof(double) * graph->n);
    int* colorMap = malloc(sizeof(int) * graph->n);
    PriorityQ* pq = newPriorityQ(graph->n);

    int ctr = 0;
    for (ctr = 0; ctr < graph->n; ctr++)
	colorMap[ctr] = WHITE;

    // we will let vert 0 be our start node
    parray->array[0] = 0;
    dists[0] = 0.0;
    colorMap[0] = GREY;
    pqInsert(pq, 0, 0.0);

    while (!pqIsEmpty(pq))
    {
	int vert = pqTop(pq);
	pqPop(pq);
	colorMap[vert] = BLACK;

	int nbrItr = 0;
	for (nbrItr = 0; nbrItr < graph->deg[vert]; nbrItr++)
	{
	    int nbrVert = graph->nbrs[vert][nbrItr];
	    if (colorMap[nbrVert] == WHITE)
	    {
		parray->array[nbrVert] = vert;
		colorMap[nbrVert] = GREY;
		dists[nbrVert] = dists[vert] + graph->wts[vert][nbrItr];
		pqInsert(pq, nbrVert, dists[nbrVert]);
	    }
	    else if (colorMap[nbrVert] == GREY)
	    {
		double newDist = dists[vert] + graph->wts[vert][nbrItr];
		if (newDist < dists[nbrVert])
		{
		    dists[nbrVert] = newDist;
		    parray->array[nbrVert] = vert;
		    pqDecreaseKey(pq, nbrVert, newDist);
		}		
	    }
/* 	    else if (colorMap[nbrVert] == BLACK) */
/* 	    { */
/* 		// pass */
/* 	    } */
/* 	    else */
/* 	    { */
/* 	    } */
	}

    }

    for (ctr = 0; ctr < parray->n; ctr++)
    {
	if (colorMap[ctr] != BLACK)
	    printf("node %d was not processed\n", ctr);
    }
    printf("parray->n == %d\n", parray->n);

    // output
    binWritePArray(parray, argv[2]);

    // mem free
    free(dists);
    free(colorMap);
    freePriorityQ(pq);
    freePArray(parray);
    freeGraph(graph);

    return 0;
}
