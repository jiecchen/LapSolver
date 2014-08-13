/* Copyright 2014 Yale Institute for Network Science
 * Author: Alex Reinking
 *
 * Example benchmark for vector arithmetic
 */
#include <sys/time.h>
#include <mkl.h>
#include <vector>
#include <numeric>
#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <algorithm>
#include "util/aligned.h"
#include "structures/graph.h"

int main(int argc, char const *argv[])
{
    aligned_vector<int> u = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1 };
    aligned_vector<int> v = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 5, 8 };
    aligned_vector<double> w = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

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

    return 0;
}
