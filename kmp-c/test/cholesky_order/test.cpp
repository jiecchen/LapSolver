// Test for Cholesky elimination order

#include <cstdlib>
#include <string>
#include <structures/Graph.h>
#include <structures/GraphLoader.h>
#include <algorithms/PartialCholeskyOrder.h>

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
        g = GraphLoader::fromFile(inFileName);
    }
    catch (int e)
    {
        if (e > 0)
            fprintf(stderr, "Parse error in '%s' on line %d\n", inFileName.c_str(), e + 1);
        else
            fprintf(stderr, "Failed to open file '%s'\n", inFileName.c_str());
        return 0;
    }

    PartialCholeskyOrder gvr(g);
    printf("Eliminate %d vertices\n", gvr.removal_count);

    printf("Eliminate:\n");
    for (int i = 0; i < gvr.removal_count; i++) {
        printf("%d\n", gvr.permutation[i]);
    }

    printf("Keep:\n");
    for (int i = gvr.removal_count; i < g.nv; i++) {
        printf("%d\n", gvr.permutation[i]);
    }

    return 0;
}
