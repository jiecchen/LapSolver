/* Copyright 2014 Yale Institute for Network Science
 * Author: Alex Reinking
 *
 * Example benchmark for vector arithmetic
 */
#include <cstdlib>
#include "structures/GraphFactory.h"
#include "algorithms/ShortestPathTree.h"

int main(int argc, char const *argv[])
{
    Graph g;
    try
    {
        if (argc < 2)
            g = GraphFactory::fromFile("inG.ijv");
        else
            g = GraphFactory::fromFile(argv[1]);
    }
    catch (int e)
    {
        if (e > 0)
            fprintf(stderr, "Failed to parse input file (line %d)\n", e + 1);
        else
            fprintf(stderr, "Failed to open file\n");
        return 0;
    }

    ShortestPathTree spt(g, 0);
    for (int i = 0; i < g.nv; ++i)
        printf("%d ", spt.getParentArray()[i]);
    printf("\n");

    return 0;
}
