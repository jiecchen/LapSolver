#include "structures/Graph.h"
#include "StarDecompositionTree.h"
#include "algorithms/ShortestPathTree.h"

StarDecompositionTree::StarDecompositionTree (const Graph& g_, int source) {
    g = &g_;
    n = g->nv();
    parentIndex = new int[n];
    currentLevel = new int[n];
}

StarDecompositionTree::~StarDecompositionTree () {
    delete[] currentLevel;
    delete[] parentIndex;
}

void StarDecompositionTree::BuildTree(int source, int level) {
    ShortestPathTree* spt = new ShortestPathTree(*g, source);
    parentIndex[0:n] = spt->parentIndex[0:n];
}