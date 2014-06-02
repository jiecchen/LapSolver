// figure out prob of an edge in a graph being in a truly random spanning tree
// output prob of each edge to line of stdout in ascii
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "cyclepoprandtreemodule.h"
#include "graph.h"

const char* DOC_STRING = "run like cyclepoprandedgeprob GRAPH_IN.BINGIJV\n";
const int NUM_ITERATIONS = 1000;  // how many samples to take

// stuff for heapsorting the edges of graph
// data structures
// edge sorted by u and then v nodes in ascending order (u must be smaller than v)
typedef struct
{
    int u;        // the first vert of edge
    int v;        // the second vert of edge
    int key;      // this is the edge number (according to binijv)
    int heapIdx;
} entry;

typedef struct
{
    int capacity; // how many entries there (not mutable!)
    int nextOpen; // the lowest idx free slot
    entry** heap; // actual 'q'
    entry* values;// raw data that gets orginized in the heap
} PriorityQ;

// malloc and return a fresh data structure
PriorityQ* newPriorityQ(int cap)
{
    PriorityQ* newQ = malloc(sizeof(PriorityQ));
    newQ->capacity = cap;
    newQ->nextOpen = 0;
    newQ->heap = malloc(sizeof(entry*)*cap);
    newQ->values = malloc(sizeof(entry)*cap);
    return newQ;
};

// no safety checks against double frees etc.
void freePriorityQ(PriorityQ* pq)
{
    free(pq->heap);
    free(pq->values);
    free(pq);
    pq = NULL;
}

int pqIsEmpty(PriorityQ* pq)
{
    if (pq->nextOpen == 0)
        return 1;
    else
        return 0;
}

// no check to avoid buffer overrun
void pqInsert(PriorityQ* pq, int key, int u, int v)
{
    int curSlot = pq->nextOpen++;
    pq->values[key].u = u;
    pq->values[key].v = v;
    pq->values[key].key = key;
    pq->values[key].heapIdx = curSlot;

    pq->heap[curSlot] = &(pq->values[key]);


    // make sure 3/2 = 1 not 2
    while (curSlot > 0 && 
	   (pq->heap[(curSlot-1)/2]->u > pq->heap[curSlot]->u
	    || (pq->heap[(curSlot-1)/2]->u == pq->heap[curSlot]->u && pq->heap[(curSlot-1)/2]->v > pq->heap[curSlot]->v)))
    {
        // swap with parent
        int parentIdx = (curSlot-1)/2;
        entry* parentEntry = pq->heap[parentIdx];
        pq->heap[parentIdx] = pq->heap[curSlot];
        pq->heap[parentIdx]->heapIdx = parentIdx;

        parentEntry->heapIdx = curSlot;
        pq->heap[curSlot] = parentEntry;

        curSlot = parentIdx;
    }
}

// return key of min entry (don't pop, no check for empty!)
int pqTop(PriorityQ* pq)
{
    return pq->heap[0]->key;
}


// pop the top entry of heap (don't return it)
void pqPop(PriorityQ* pq)
{
    pq->heap[0]->key = -1;
    int curSlot = 0; // refers to gap we need to fill

    // fill in gap left at root with children
    while ((curSlot*2) + 1 < pq->nextOpen)
    {
        int leftChildIdx = curSlot*2 + 1;
        int rightChildIdx = curSlot*2 + 2;
        // left child is only option  or left child is smaller than right
        if (rightChildIdx >= pq->nextOpen
            || (pq->heap[leftChildIdx]->u < pq->heap[rightChildIdx]->u
		|| (pq->heap[leftChildIdx]->u == pq->heap[rightChildIdx]->u && pq->heap[leftChildIdx]->v < pq->heap[rightChildIdx]->v)))
        {
            entry* leftChild = pq->heap[leftChildIdx];
            //pq->heap[leftChildIdx] = pq->heap[curSlot];
            //pq->heap[leftChildIdx]->heapIdx = leftChildIdx;
            pq->heap[curSlot] = leftChild;
            pq->heap[curSlot]->heapIdx = curSlot;
            curSlot = leftChildIdx;
        }
        // right child is smaller
        else
        {
            entry* rightChild = pq->heap[rightChildIdx];
            //pq->heap[rightChildIdx] = pq->heap[curSlot];
            //pq->heap[rigthChildIdx]->heapIdx = rightChildIdx;
            pq->heap[curSlot] = rightChild;
            pq->heap[curSlot]->heapIdx = curSlot;
            curSlot = rightChildIdx;
        }    
    }

    // there may be gap at bottom of tree that we should fill witht the highest index node
    pq->nextOpen--;
    if (curSlot < pq->nextOpen)
    {
        pq->heap[curSlot] = pq->heap[pq->nextOpen];
        // now restore heap property (just like insert from now on)
        // curSlot is idx of node that might need to be moved up
        while (curSlot > 0 && 
	       (pq->heap[(curSlot-1)/2]->u > pq->heap[curSlot]->u
		|| (pq->heap[(curSlot-1)/2]->u == pq->heap[curSlot]->u && pq->heap[(curSlot-1)/2]->v > pq->heap[curSlot]->v)))
        {
            // swap with parent
            int parentIdx = (curSlot-1)/2;
            entry* parentEntry = pq->heap[parentIdx];
            pq->heap[parentIdx] = pq->heap[curSlot];
            pq->heap[parentIdx]->heapIdx = parentIdx;
            
            parentEntry->heapIdx = curSlot;
            pq->heap[curSlot] = parentEntry;
            
            curSlot = parentIdx;
        }
    }
}

// struct to store the sorted edges (as elements of array)
typedef struct
{
    int u;
    int v;
    int ijvIdx; // where  edge is in ijv;
} Edge;
  
// find the ijvIdx for the given edge in sorted Edge array using binary search
int findEdgeIdx(Edge* sortedEdges, int u, int v, int numEdges)
{
    int lowBound = 0;
    int highBound = numEdges - 1;
    
    while (lowBound != highBound)
    {
	int center = (highBound + lowBound)/2;
	
	// look up
	if (sortedEdges[center].u < u 
	    || (sortedEdges[center].u == u && sortedEdges[center].v < v))
	{
	    lowBound = center + 1;
	}
	// look down
	else if (sortedEdges[center].u > u 
	    || (sortedEdges[center].u == u && sortedEdges[center].v > v))
	{
	    highBound = center - 1;
	}
	// we found it
	else
	{
	    highBound = center;
	    lowBound = center;
	}

    }
    return sortedEdges[lowBound].ijvIdx;
    
}

void checkSortedEdges(Edge* sortedEdges, int n)
{
    int ctr = 0;
    for (ctr = 0; ctr < n - 1; ctr++)
    {
	printf("u: %d   v: %d\n", sortedEdges[ctr].u, sortedEdges[ctr].v);
	if (sortedEdges[ctr].u < sortedEdges[ctr+1].u)
	    continue;
	if (sortedEdges[ctr].u == sortedEdges[ctr+1].u
	    && sortedEdges[ctr].v == sortedEdges[ctr+1].v)
	    continue;
	fprintf(stderr, "error in sorted edges!\n");
    }


}

int main(int argc, char* argv[])
{
    srand(time(NULL));
    if (argc != 2)
    {
	fprintf(stderr, "%s", DOC_STRING);
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

    Edge* sortedEdges = malloc(sizeof(Edge) * ijv->nnz);
    int* countedEdges = malloc(sizeof(int) * ijv->nnz);

    PriorityQ* pq = newPriorityQ(ijv->nnz);
    int ctr = 0;
    for (ctr = 0; ctr < ijv->nnz; ctr++)
    {
	if (ijv->i[ctr] > ijv->j[ctr])
	{
	    int tmp = ijv->i[ctr];
	    ijv->i[ctr] = ijv->j[ctr];
	    ijv->j[ctr] = tmp;
	}
	pqInsert(pq, ctr, ijv->i[ctr], ijv->j[ctr]);
    }
    
    ctr = 0;
    while (!pqIsEmpty(pq))
    {
	int edgeIdx = pqTop(pq);
	pqPop(pq);
	sortedEdges[ctr].u = ijv->i[edgeIdx];
	sortedEdges[ctr].v = ijv->j[edgeIdx];
	sortedEdges[ctr].ijvIdx = edgeIdx;
	ctr++;
    }

    freePriorityQ(pq);


    //checkSortedEdges(sortedEdges, ijv->nnz);
    
    //printf("before while loop\n");

    int numIterations = 0;
    while (numIterations++ < NUM_ITERATIONS)
    {
	// choose root
	int root = pickRoot(ijv);
    
	// get tree
	pArray* parray = randRootTree(root, graph);
	printf("found tree\n");
	int ctr = 0;
	for (ctr = 0; ctr < ijv->n; ctr++)
	{
	    int u = ctr;
	    int v = parray->array[ctr];
	    if (u == v)
		continue;
	    if (v < u)
	    {
		int tmp = u;
		u = v;
		v = tmp;
	    }
	    int edgeIjvIdx = findEdgeIdx(sortedEdges, u, v, ijv->nnz);
	    countedEdges[edgeIjvIdx]++;
	}
	freePArray(parray);
    }

    // output
    ctr = 0;
    for (ctr = 0; ctr < ijv->nnz; ctr++)
    {
	printf("%lf\n", (double)countedEdges[ctr]/(double)NUM_ITERATIONS);
    }
    freeIJV(ijv);
    freeGraph(graph);
    free(sortedEdges);
    free(countedEdges);

    return 0;
}
