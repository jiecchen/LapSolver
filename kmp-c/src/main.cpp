/* Copyright 2014 Yale Institute for Network Science
 * Author: Alex Reinking
 *
 * Example benchmark for vector arithmetic
 */
#include <cstdlib>
#include "util/aligned.h"
#include "structures/graph.h"
#include "algorithms/ShortestPathTree.h"

int main(int argc, char const *argv[])
{
    aligned_vector<int> u    = { 0, 1, 2, 2, 3 };
    aligned_vector<int> v    = { 1, 2, 3, 4, 4 };
    aligned_vector<double> w = { 1, 1, 1, 1, 5 };

    EdgeList edges(u,v,w);
    Graph g(edges);

    g.debugPrint();

    printf("\nAdjacency List:\n");
    for (int v = 0; v < g.nv; ++v) {
        int deg = g.getDegree(v);
        auto nbrs = g.getNeighbors(v);
        printf("%d(%d):", v, deg);
        for (int i = 0; i < deg; ++i) {
            printf(" %d", nbrs[i]);
        }
        printf("\n");
    }

    ShortestPathTree spt(g, 4);
    printf("Dijkstra:");
    for (int i = 0; i < g.nv; ++i)
        printf(" %d", spt.getParentArray()[i]);
    printf("\n");

    return 0;
}
