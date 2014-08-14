#pragma once

#include "../structures/graph.h"

struct ShortestPathTree {
    ShortestPathTree(const Graph &g, int source);
    ~ShortestPathTree();

    double *dist;
    int *parent;
    int *parentIndex;
    double *weight;
};
