/* Copyright 2014 Yale Institute for Network Science
 * Author: Alex Reinking
 *
 * Example benchmark for vector arithmetic
 */
#include <cstdlib>
#include <string>
#include <structures/GraphLoader.h>
#include <algorithms/ShortestPathTree.h>

int main(int argc, char const *argv[])
{
    atexit(MKL_Free_Buffers);
    Graph g = GraphLoader::fromStdin();

    ShortestPathTree spt(g, 0);
    for (int i = 0; i < g.nv; ++i)
        printf("%d,%d ", spt.parent[i], i);
    printf("\n");

    return 0;
}
