// this is a specialized heap based priority Q designed for use in dijkstras
// you must know in advance how big the capacity will need to be
// moroever, the 'key's (vertices) must go from 0 -> capicity-1
// (perhaps 'key' & value are swapped
#include <stdio.h>
#include <stdlib.h>
#include "priorityq.h"


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
void pqInsert(PriorityQ* pq, int key, double value)
{
    int curSlot = pq->nextOpen++;
    pq->values[key].value = value;
    pq->values[key].key = key;
    pq->values[key].heapIdx = curSlot;

    pq->heap[curSlot] = &(pq->values[key]);


    // make sure 3/2 = 1 not 2
    while (curSlot > 0 && pq->heap[(curSlot-1)/2]->value > pq->heap[curSlot]->value)
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
	    || pq->heap[leftChildIdx]->value < pq->heap[rightChildIdx]->value)
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
	while (curSlot > 0 && pq->heap[(curSlot-1)/2]->value > pq->heap[curSlot]->value)
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


void pqDecreaseKey(PriorityQ* pq, int key, double newValue)
{
    entry* modifiedEntry = &(pq->values[key]);
    modifiedEntry->value = newValue;
    int curSlot = modifiedEntry->heapIdx;
    while (curSlot > 0 && pq->heap[(curSlot-1)/2]->value > pq->heap[curSlot]->value)
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

// for testing only
/* int main(int argc, char* argv[]) */
/* { */
/*     PriorityQ* pq = newPriorityQ(100); */
/*     int ctr = 0; */
/*     for (ctr = 0; ctr< 100; ctr++) */
/*     { */
/* 	pqInsert(pq, ctr, (double)(-1 * ctr)); */
/* 	int top = pqTop(pq); */
/* 	printf(" %d ", top); */
/*     } */


/*     printf("\n\n\n"); */

/*     pqDecreaseKey(pq, 50, -10000.0); */
/*     pqDecreaseKey(pq, 30, -10001.0); */
    

/*     while (!pqIsEmpty(pq)) */
/*     { */
/* 	int top = pqTop(pq); */
/* 	pqPop(pq); */
/* 	printf(" %d ", top); */
/*     } */

/*     freePriorityQ(pq); */

/*     return 0; */
/* } */
