#pragma once

#include "structures/Graph.h"

struct StarDecompositionTree {
    StarDecompositionTree (const Graph& g_, int source);
    ~StarDecompositionTree ();

    void BuildTree(int source, int level);

    const Graph* g;
    int n;
    int* currentLevel;
    int* parentIndex;
};