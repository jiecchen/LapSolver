/* Copyright 2014 Yale Institute for Network Science
 * Author: Alex Reinking
 *
 * Example benchmark for vector arithmetic
 */
#include <cstdlib>
#include <string>
#include "structures/GraphFactory.h"
#include "algorithms/ShortestPathTree.h"

int main(int argc, char const *argv[])
{
    Graph g;
    std::string inFileName;
    try
    {
        if (argc < 2)
            inFileName = "inG.ijv";
        else
            inFileName = argv[1];
        g = GraphFactory::fromFile(inFileName);
    }
    catch (int e)
    {
        if (e > 0)
            fprintf(stderr, "Parse error in '%s' on line %d\n", inFileName.c_str(), e + 1);
        else
            fprintf(stderr, "Failed to open file '%s'\n", inFileName.c_str());
        return 0;
    }

    ShortestPathTree spt(g, 0);
    for (int i = 0; i < g.nv; ++i)
        printf("%d ", spt.parent[i]);
    printf("\n");

    return 0;
}
