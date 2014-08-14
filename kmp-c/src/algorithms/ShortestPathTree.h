#pragma once

#include "structures/Graph.h"

class ShortestPathTree {
public:
    ShortestPathTree(const Graph &g, int source);
    ~ShortestPathTree();

    const double * getDistances() const { return dist; }
    const int * getParentArray() const { return parent; }
    const double * getWeights() const { return weight; }

private:
    double *dist;
    int *parent;
    double *weight;
};
