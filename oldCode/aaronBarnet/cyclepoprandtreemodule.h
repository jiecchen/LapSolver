#ifndef CYCLEPOPRANDTREEMODULE_H
#define CYCLEPOPRANDTREEMODULE_H
// basic functions that form basis of cycle pop random walk

#include "graph.h"

int pickRoot(ijvType* ijv); // pick a node to be tree root (with prob prop. to node weight)

// pick neighboring node with prob prop. to reciprical of edge length
// vert is the node for which you want a nbr
int pickNeighbor(int vert, myGraph* graph, int avoid);

// return random  spanning tree for graph given root node
pArray* randRootTree(int root, myGraph* graph);

#endif //CYCLEPOPRANDTREEMODULE_H
