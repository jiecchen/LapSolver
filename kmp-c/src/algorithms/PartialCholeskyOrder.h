#pragma once

#include <vector>
#include "structures/graph.h"

using namespace std;

class PartialCholeskyOrder
{
public:
	PartialCholeskyOrder(const Graph &g);			// The constructor computes 
	~PartialCholeskyOrder();

	//EnhancedArray OrderDegreeTwo();				// Gets the ordering for degree two vertices
	//EnhancedArray OrderDegreeTwoChains();			// Orders the vertices in degree two chains
	//EnhancedArray CheckIfCycle();					// Check if the graph is composed of a cycle

	//EnhancedArray OrderDegreeOne();				// Get the Cholesky order for degree one vertices
	//EnhancedArray OrderBFS();						// Find the order using a BFS search

	//EnhancedArray ConstructPermutation();			// The answer generated should be a permutation

private:
	struct EnhancedArray {
		int *v;
		int size;
	} answer;

	int n;
	Graph graph;
};