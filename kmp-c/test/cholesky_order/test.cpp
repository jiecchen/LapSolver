// Test for Cholesky elimination order

#include <cstdlib>
#include <string>
#include <memory>
#include <structures/Graph.h>
#include <structures/GraphLoader.h>
#include <algorithms/PartialCholeskyOrder.h>
#include <util/Benchmark.h>

int main(int argc, char *argv[])
{
    atexit(mkl_free_buffers);
    Graph g = GraphLoader::fromStdin();

    std::shared_ptr<PartialCholeskyOrder> gvr;

    auto bench = make_benchmark(argc, argv, [&] () {
        gvr.reset(new PartialCholeskyOrder(g));
    });
    printf("Eliminate %d vertices\n", gvr->removal_count);

    printf("Eliminate:\n");
    for (int i = 0; i < gvr->removal_count; i++) {
        printf("%d\n", gvr->permutation[i]);
    }

    printf("Keep:\n");
    for (int i = gvr->removal_count; i < g.nv(); i++) {
        printf("%d\n", gvr->permutation[i]);
    }

    return 0;
}
