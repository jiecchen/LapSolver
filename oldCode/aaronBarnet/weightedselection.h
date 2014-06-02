#ifndef WEIGHTEDSELECTION_H
#define WEIGHTEDSELECTION_H
// module to allow random selection of a set of ints according to a set of probabilistic weights


#include <stdbool.h>

// struct for selecting weighted random elements
typedef struct
{
    int n;                     // number of elements
    int cur;                   // the next free slot to insert something
                               // once cur == n, everything is ready
    int numAvoids;
    double weightAvoids;       // how much weight are we avoiding
    double totalWeight;        // total weight of all elements
    bool* avoid;               // what we should never select
    int* elements;             // what we are actually selecting
    double* cumulativeWeights; // sums of all weights
} WeightedElements;

WeightedElements* newWeightedElements(int n);

void freeWeightedElements(WeightedElements* in);

// remove all avoid elements from we
void weCleanUp(WeightedElements* we);

// all elements must be sorted from low to high (necesary to make locating elements efficient)
void weAvoid(WeightedElements* we, int element);

// must insert elements from low to high to use weAvoid
void weInsert(WeightedElements* we, int i, double w);

// assume we never want to get the same value more than once!!!
int weSelect(WeightedElements* we);



#endif // WEIGHTEDSELECTION_H
