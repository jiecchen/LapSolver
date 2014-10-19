/* Copyright 2014 Yale Institute for Network Science
 * Author: Alex Reinking
 */
#include <cstdlib>
#include <algorithm>
#include <structures/GraphLoader.h>
#include <algorithms/ShortestPathTree.h>
#include <lsst/StarDecompositionTree.h>

int main(int argc, char **argv)
{
    atexit(mkl_free_buffers);
    Graph g = GraphLoader::fromStdin();

    StarDecompositionTree spt(g, 0);

    printf("%d %d\n", g.nv(), g.nv() - 1);
    for (int i = 0; i < g.nv(); ++i) {
        if (i != 0) {
            int parent = g.neighbor(i, spt.parentIndex[i]);
            int u = std::min(i, parent);
            int v = std::max(i, parent);
            double w = 1.0; //spt.weight[i];
            printf("%d %d %lf\n", u, v, w);
        }
    }

    return 0;
}
