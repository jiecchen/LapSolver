#include "weightedselection.h"
#include <stdlib.h>
#include <stdio.h>

WeightedElements* newWeightedElements(int n)
{
    WeightedElements* out = malloc(sizeof(WeightedElements));
    out->n = n;
    out->cur = 0;
    out->totalWeight = 0.0;
    out->elements = malloc(sizeof(int) * n);
    out->cumulativeWeights = malloc(sizeof(double) * n);
    out->avoid = malloc(sizeof(bool) * n);
    out->numAvoids = 0;
    out->weightAvoids = 0.0;
    int ctr = 0;
    for (ctr = 0; ctr < n; ctr++)
    {
	out->avoid[ctr] = false;
    }
    return out;
}

void freeWeightedElements(WeightedElements* in)
{
    free(in->elements);
    free(in->cumulativeWeights);
    free(in->avoid);
    free(in);
}

// remove all avoid elements from we
void weCleanUp(WeightedElements* we)
{
    //printf("running cleanup\n");
    int insertIdx = 0;
    double subtractWeight = 0.0;
    int ctr = 0;
    for (ctr = 0; ctr < we->n; ctr++)
    {
	if (we->avoid[ctr])
	{
	    if (ctr == 0)
		subtractWeight += we->cumulativeWeights[ctr];
	    else
		subtractWeight += we->cumulativeWeights[ctr] - we->cumulativeWeights[ctr-1];
	    continue;
	}
	
	we->elements[insertIdx] = we->elements[ctr];
	we->cumulativeWeights[insertIdx] = we->cumulativeWeights[ctr] - subtractWeight;
	we->avoid[insertIdx] = false;
	insertIdx++;
    }
    we->n -= we->numAvoids;
    we->numAvoids = 0;
    we->weightAvoids = 0.0;
    we->totalWeight -= subtractWeight;
    we->cur = we->n;
}

// all elements must be sorted from low to high
void weAvoid(WeightedElements* we, int element)
{
    int highBound = we->n - 1;
    int lowBound = 0;
    while (highBound != lowBound)
    {
	//printf("highBound: %d   lowBound: %d\n", highBound, lowBound);
	int center = (highBound + lowBound)/2;
	// move up
	if (we->elements[center] < element)
	{
	    lowBound = center + 1;
	}
	// move down
	else if (we->elements[center] > element)
	{
	    highBound = center - 1;
	}
	// center is what we want
	else
	{
	    highBound = center;
	    lowBound = center;
	}
    }
    if (we->elements[lowBound] != element)
    {
	fprintf(stderr, "ERROR: can not find element: %d.\n", element);
	exit(1);
    }
    else
    {
	we->avoid[lowBound] = true;
	we->numAvoids++;
	if (lowBound == 0)
	    we->weightAvoids += we->cumulativeWeights[lowBound];
	else
	    we->weightAvoids += we->cumulativeWeights[lowBound] - we->cumulativeWeights[lowBound - 1];	
    }
    if (we->weightAvoids > we->totalWeight / 2)
	weCleanUp(we);
}

// must insert elements from low to high to use weAvoid
void weInsert(WeightedElements* we, int i, double w)
{
    if (we->cur >= we->n)
    {
	fprintf(stderr, "WeightedElements is full!\n");
	exit(1);
    }
    
    we->elements[we->cur] = i;
    if (we->cur == 0)
	we->cumulativeWeights[we->cur] = w;
    else
	we->cumulativeWeights[we->cur] = w + we->cumulativeWeights[we->cur - 1];
    we->totalWeight = we->cumulativeWeights[we->cur];
    we->cur++;
}

// assume we never want to get the same value more than once!!!
int weSelect(WeightedElements* we)
{
    int selectionIdx = -1;
    int lowBound = 0;
    while (selectionIdx < 0 || we->avoid[selectionIdx])
    {
	double randWeight = (double)rand()/(double)(RAND_MAX) * (double)we->totalWeight;
	// binary search for corresponding element
	//int curPos = we->n / 2;
	int highBound = we->n - 1;
	lowBound = 0;
	while (highBound != lowBound)
	{
	    //printf("highBound: %d   lowBound: %d\n", highBound, lowBound);
	    int center = (highBound + lowBound)/2;
	    // move up
	    if (we->cumulativeWeights[center] < randWeight)
	    {
		lowBound = center + 1;
	    }
	    // move down
	    else if (center > 0 && we->cumulativeWeights[center - 1] > randWeight)
	    {
		highBound = center - 1;
	    }
	    // center is what we want
	    else
	    {
		highBound = center;
		lowBound = center;
	    }
	}
	selectionIdx = lowBound;
    }
    we->avoid[selectionIdx] = true;
    we->numAvoids++;
    if (lowBound == 0)  // WHAT IS THIS?
	we->weightAvoids += we->cumulativeWeights[selectionIdx];
    else
	we->weightAvoids += we->cumulativeWeights[selectionIdx] - we->cumulativeWeights[selectionIdx - 1];	
    int out = we->elements[selectionIdx];
    if (we->weightAvoids > we->totalWeight / 2)
	weCleanUp(we);
    return out;
}

