// Test for Cholesky elimination order

#include <cstdlib>
#include <string>
#include <structures/Graph.h>
#include <structures/GraphLoader.h>
#include <algorithms/PartialCholeskyOrder.h>

int main(int argc, char const *argv[])
{
    atexit(mkl_free_buffers);
    Graph g = GraphLoader::fromStdin();

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
