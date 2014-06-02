// cluster the nodes in small dijkstra trees, then run dijkstra to combine these small trees into one big tree (try to get low stretch)

#include <stdio.h>
#include <stdlib.h>
#include "graph.h"
#include "priorityq.h"


const char* DOC_STRING = "run like ./doubledijkstra GRAPH_IN.BINIJV TREEOUT.PARRAY\n";

const int CLUST_SIZE = 20; // how big each mini dijkstra tre should be

enum COLORS { WHITE = -3, GREY = -2, BLACK = -1}; // BLACK (meaning processed) can also be the comp number which changes


// figure out what order to go through the nodes in (run dijkstra twice)
// for now assume we alaways start with node 0...
// take start node and find furthest node away.
// Now set ordered nodes so in order of distance from "furthest node"
// inputs:
//  - graph: cf graph.h
//  - pq: cf priorityq.h
//  - colorMap: will be initialized (will overwrite current values)
//  - dists: will overwrite current values
// output: orderedNodes is modified to be output
void orderNodes(myGraph* graph, PriorityQ* pq, int* colorMap, double* dists, int* orderedNodes)
{
    int startNode = 0; // temp
    
    int ctr = 0;
    for (ctr = 0; ctr < graph->n; ctr++)
	colorMap[ctr] = WHITE;
    
    // dijkstra it up
    dists[startNode] = 0.0;
    colorMap[startNode] = GREY;
    pqInsert(pq, startNode, 0.0);

    int vert = 0;
    while (!pqIsEmpty(pq))
    {
	vert = pqTop(pq);
	pqPop(pq);
	colorMap[vert] = BLACK;

	int nbrItr = 0;
	for (nbrItr = 0; nbrItr < graph->deg[vert]; nbrItr++)
	{
	    int nbrVert = graph->nbrs[vert][nbrItr];
	    if (colorMap[nbrVert] == WHITE)
	    {
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
		    pqDecreaseKey(pq, nbrVert, newDist);
		}		
	    }
	}

    }

    // vert should now be furthest node from startNode
    // rerun dijkstra and this time keep track of orderedNodes
    ctr = 0;
    for (ctr = 0; ctr < graph->n; ctr++)
	colorMap[ctr] = WHITE;
    
    dists[vert] = 0.0;
    colorMap[vert] = GREY;
    pqInsert(pq, vert, 0.0);
    
    ctr = 0;
    while (!pqIsEmpty(pq))
    {
	vert = pqTop(pq);
	pqPop(pq);
	colorMap[vert] = BLACK;
	orderedNodes[ctr++] = vert;

	int nbrItr = 0;
	for (nbrItr = 0; nbrItr < graph->deg[vert]; nbrItr++)
	{
	    int nbrVert = graph->nbrs[vert][nbrItr];
	    if (colorMap[nbrVert] == WHITE)
	    {
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
		    pqDecreaseKey(pq, nbrVert, newDist);
		}		
	    }
	}

    }
}

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
    pArray* parray = newPArray(graph->n); // will have multiple trees in same array
    double* parrayDists = malloc(sizeof(double) * graph->n); // keep track of dists for ijv output
    int* components = malloc(sizeof(int) * graph->n);  // each tree gets its own number (# > 0 for comps; -2 for untouched nodes, -1 for nodes in q)
    double* dists = malloc(sizeof(double) * graph->n);
    //int* colorMap = malloc(sizeof(int) * graph->n);
    PriorityQ* pq = newPriorityQ(graph->n);
    int* orderedNodes = malloc(sizeof(int) * graph->n); // tells what order to go through nodes in
    int orderedNodesItr = 0;

    int* finalNodeOrder = malloc(sizeof(int) * graph->n); // record of order all nodes were processed in
    int finalNodeOrderItr = 0;
	    
    orderNodes(graph, pq, components, dists, orderedNodes);

    int ctr = 0;
    for (ctr = 0; ctr < graph->n; ctr++)
	components[ctr] = WHITE;

    int comp = 0; // keeps track of what comp we are on
    while (1)
    {
	// find the first node in orderedNodes not assigned a comp
	while (orderedNodesItr < graph->n)
	{
	    if (components[orderedNodes[orderedNodesItr]] < 0)
		break;
	    else
		orderedNodesItr++;
	}
    
	// stop loop when we run out of nodes in orderedNodes
	if (orderedNodesItr >= graph->n)
	    break;

	int startNode = orderedNodes[orderedNodesItr];
	parray->array[startNode] = startNode;
	dists[startNode] = 0.0;
	components[startNode] = GREY;
	pqInsert(pq, startNode, 0.0);


	ctr = 0;
	while (!pqIsEmpty(pq) && ctr++ < CLUST_SIZE)
	{
	    int vert = pqTop(pq);
	    pqPop(pq);
	    components[vert] = comp;
	    finalNodeOrder[finalNodeOrderItr++] = vert;

	    int nbrItr = 0;
	    for (nbrItr = 0; nbrItr < graph->deg[vert]; nbrItr++)
	    {
		int nbrVert = graph->nbrs[vert][nbrItr];
		if (components[nbrVert] == WHITE)
		{
		    parray->array[nbrVert] = vert;
		    parrayDists[nbrVert] = graph->wts[vert][nbrItr];
		    components[nbrVert] = GREY;
		    dists[nbrVert] = dists[vert] + graph->wts[vert][nbrItr];
		    pqInsert(pq, nbrVert, dists[nbrVert]);
		}
		else if (components[nbrVert] == GREY)
		{
		    double newDist = dists[vert] + graph->wts[vert][nbrItr];
		    if (newDist < dists[nbrVert])
		    {
			dists[nbrVert] = newDist;
			parray->array[nbrVert] = vert;
			parrayDists[nbrVert] = graph->wts[vert][nbrItr];
			pqDecreaseKey(pq, nbrVert, newDist);
		    }		
		}
	    }

	}
    
	// clean up Q (empty it out and turn all the grey nodes white)
	while (!pqIsEmpty(pq))
	{
	    int vert = pqTop(pq);
	    pqPop(pq);
	    components[vert] = WHITE;
	}

	// increment comp
	//printf("incrementing com\n");
	comp++;
    }

    printf("comp: %d\n", comp);
    // now we have to join everything together
    

    // make a second graph: each node represents a components, edges are connections between components
    // (when there are multiple edge between two components takes the shorter one
    // since don't know how many edges we will need: just malloc enough for every edge in first graph (overkill)
    myGraph* h = newGraph(comp, graph->nnz); // extreme overestimate of number of edges
    int** h_is = malloc(sizeof(int*) * h->n);  // these correspond to the actaul verts in graph for each edge in h
    int** h_js = malloc(sizeof(int*) * h->n);
    int* h_isBlock = malloc(sizeof(int) * h->nnz);
    int* h_jsBlock = malloc(sizeof(int) * h->nnz);
    finalNodeOrderItr = 0;
    int totalEdgeItr = 0; // where we are in nbrsBlock and wtsBlock 
    for (ctr = 0; ctr < h->n; ctr++)
    {
	h->deg[ctr] = 0;
	h->nbrs[ctr] = &(h->nbrsBlock[totalEdgeItr]);
	h->wts[ctr] = &(h->wtsBlock[totalEdgeItr]);
	h_is[ctr] = &(h_isBlock[totalEdgeItr]);
	h_js[ctr] = &(h_jsBlock[totalEdgeItr]);

	while (finalNodeOrderItr < graph->n   // avoid buffer overrun
	       && components[finalNodeOrder[finalNodeOrderItr]] == ctr)
	{
//	    if (ctr == 0)
//		printf("dealing with comp 0 in h\n");

	    int node = finalNodeOrder[finalNodeOrderItr];
	    int nbrItr = 0;
	    for (nbrItr = 0; nbrItr < graph->deg[node]; nbrItr++)
	    {
		int nbr = graph->nbrs[node][nbrItr];
		if (components[nbr] != ctr)
		{
		    //    if (ctr== 0)
		    //	printf("found nbr for comp 0\n");

		    h->deg[ctr]++;
		    h->nbrsBlock[totalEdgeItr] = components[nbr];
		    h->wtsBlock[totalEdgeItr] = graph->wts[node][nbrItr];
		    h_isBlock[totalEdgeItr] = node;
		    h_jsBlock[totalEdgeItr++] = nbr;
		    if (totalEdgeItr > h->nnz)
			printf ("totalEdgeItr too big\n");
		}
	    }
	    finalNodeOrderItr++;
	}	
    }
    // now we have multigraph that might have multiple edges between same nodes
    // we will run dijkstra starting with comp 0
    int* colorMap = malloc(sizeof(int) * h->n);
    for (ctr = 0; ctr < h->n; ctr++)
	colorMap[ctr] = WHITE;           // components is now the colormap for graph h
    int* hi = malloc(sizeof(int) * h->n);    // these are the actual edges (in graph) that join the different components
    int* hj = malloc(sizeof(int) * h->n);    // each comp has an ede hi[comp] -> hj[comp] (note that the root will have a self-loop)
    pArray* hpArray = newPArray(h->n); // will have multiple trees in same array
    double* hpArrayDists = malloc(sizeof(double) * h->n); // so we can back out the last edge in dists (in case we used an overlylong in in disjkstras)
    for (ctr = 0; ctr < h->n; ctr++)
	hpArray->array[ctr] = -1;

    int startNode = 0;
    hpArray->array[startNode] = startNode;
    dists[startNode] = 0.0;
    colorMap[startNode] = GREY;
    pqInsert(pq, startNode, 0.0);

    while (!pqIsEmpty(pq))
    {
	int vert = pqTop(pq);
	// check associated edge in graph
	/* if (hpArray->array[vert] == vert) */
/* 	{ */
/* 	    printf("vert is parent\n"); */
/* 	} */
/* 	else */
/* 	{ */
/* 	    int i = hi[vert]; */
/* 	    int j = hj[vert]; */
/* 	    if (components[j] != vert || components[i] != hpArray->array[vert]) */
/* 	    { */
/* 		printf("problem\n"); */
/* 	    } */
/* 	    else */
/* 	    { */
/* 		printf("ok\n"); */
/* 	    } */
/* 	} */
	pqPop(pq);
	colorMap[vert] = BLACK;
	
	int nbrItr = 0;
	//printf("h->deg[%d]: %d\n", vert, h->deg[vert]);
	for (nbrItr = 0; nbrItr < h->deg[vert]; nbrItr++)
	{
	    //printf("found child\n");
	    int nbrVert = h->nbrs[vert][nbrItr];
	    //printf("determined  nbrVert: %d\n", nbrVert);
	    if (colorMap[nbrVert] == WHITE)
	    {
		hpArray->array[nbrVert] = vert;
		hpArrayDists[nbrVert] = h->wts[vert][nbrItr];
		hi[nbrVert] = h_is[vert][nbrItr];
		hj[nbrVert] = h_js[vert][nbrItr];
		colorMap[nbrVert] = GREY;
		dists[nbrVert] = dists[vert] + h->wts[vert][nbrItr];
		pqInsert(pq, nbrVert, dists[nbrVert]);
	    }
	    // case where we are already using an overly long edge to go from vert to nbrVert
	    else if (hpArray->array[nbrVert] == vert)
	    {
		// do we have a shorter edge now?
		if (h->wts[vert][nbrItr] < hpArrayDists[nbrVert])
		{
		    dists[nbrVert] -= (hpArrayDists[nbrVert] - h->wts[vert][nbrItr]);
		    pqDecreaseKey(pq, nbrVert, dists[nbrVert]);
		    hpArrayDists[nbrVert] = h->wts[vert][nbrItr];
		    hi[nbrVert] = h_is[vert][nbrItr];
		    hj[nbrVert] = h_js[vert][nbrItr];
		}
		
	    }
	    else if (colorMap[nbrVert] == GREY)
	    {
		double newDist = dists[vert] + h->wts[vert][nbrItr];
		if (newDist < dists[nbrVert])
		{
		    dists[nbrVert] = newDist;
		    hpArray->array[nbrVert] = vert;
		    hpArrayDists[nbrVert] = h->wts[vert][nbrItr];
		    hi[nbrVert] = h_is[vert][nbrItr];
		    hj[nbrVert] = h_js[vert][nbrItr];
		    pqDecreaseKey(pq, nbrVert, newDist);
		}		
	    }
	}
	    
    }
    
    // now we just output all edges as one giant ijv list (avoid combining everything into one parray tree (would require tree rehanging)
    FILE* fp = NULL;
    if ((fp = fopen(argv[2], "wt")) == NULL)
    {
	fprintf(stderr, "ERROR: could not open output file\n");
	exit(1);
    }
    // expectation that everything will be incremented by one for matlab compatibility
    // also we have been dealing in lengths but files should have costs (reciprical of length)
    int num_verts = graph->n;
    int num_edges = num_verts - 1;
    fwrite(&num_verts, sizeof(int), 1, fp);
    fwrite(&num_edges, sizeof(int), 1, fp);
    // write is
    for (ctr = 0; ctr < graph->n; ctr++)
    {
	if (ctr != parray->array[ctr])
	    {
		int tmp = ctr + 1;
		fwrite(&tmp, sizeof(int), 1, fp);
	    }
    }
    for (ctr = 0; ctr < h->n; ctr++)
    {
	if (hpArray->array[ctr] != ctr)
	{
	    int tmp =hi[ctr] + 1; 
	    fwrite(&tmp, sizeof(int), 1, fp);
	}
    }
    // write js
    for (ctr = 0; ctr < graph->n; ctr++)
    {
	if (ctr != parray->array[ctr])
	{
	    int tmp = parray->array[ctr] + 1; 
	    fwrite(&tmp, sizeof(int), 1, fp);
	}
    }
    for (ctr = 0; ctr < h->n; ctr++)
    {
	if (hpArray->array[ctr] != ctr)
	{
	    int tmp = hj[ctr] + 1;
	    fwrite(&tmp, sizeof(int), 1, fp);
	}
    }
    // write vs
    for (ctr = 0; ctr < graph->n; ctr++)
    {
	if (ctr != parray->array[ctr])
	{
	    double tmp = 1.0/parrayDists[ctr]; 
	    fwrite(&tmp, sizeof(double), 1, fp);
	}
    }
    for (ctr = 0; ctr < h->n; ctr++)
    {
	if (hpArray->array[ctr] != ctr)
	{
	    double tmp = 1.0/hpArrayDists[ctr];
	    fwrite(&tmp, sizeof(double), 1, fp);
	}
    }

    fclose(fp);

    free(hi);
    free(hj);
    free(h_is);
    free(h_js);
    free(h_isBlock);
    free(h_jsBlock);
    

    // output
    //binWritePArray(parray, argv[2]);

    // mem free
    free(dists);
    free(colorMap);
    free(components);
    free(orderedNodes);
    freePriorityQ(pq);
    freePArray(parray);
    free(parrayDists);
    freePArray(hpArray);
    free(hpArrayDists);
    freeGraph(graph);
    freeGraph(h);

    return 0;
}
