#pragma once

#include <vector>
#include <algorithm>
#include "structures/graph.h"

using namespace std;

class PartialCholeskyOrder
{
public:
	PartialCholeskyOrder(const Graph &g);					// The constructor computes 
	~PartialCholeskyOrder();

	vector <int> DegreeTwoOrdering();						// The ordering in which the degree two vertices should be processed
	vector <int> DegreeTwoChainOrdering();					// The ordering in which the degree two vertices present in chains should be processed
	vector <int> CheckIfCycle();							// Checks if the graphs is a cycle. Treat this case accordingly.

	vector <int> DegreeOneOrdering(); 						// Gets the ordering for degree one vertices

	vector <int> ConstructFinalOrdering(vector <int> v);	// Turns the given vertices into a permutation

private:
	int n;
	Graph graph;

	int *permutation;
	int removal_count;

	int *removed;
    int *updated_deg;

    int *un_removable;
};