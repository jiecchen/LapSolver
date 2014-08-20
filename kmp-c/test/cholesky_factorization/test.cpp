// Test for Cholesky elimination order

#include <cstdlib>
#include <iostream>
#include <string>
#include <memory>
#include <structures/Graph.h>
#include <structures/GraphLoader.h>
#include <algorithms/PartialCholeskyFactorization.h>
#include <util/Benchmark.h>

using namespace std;

int main(int argc, char *argv[])
{
    atexit(mkl_free_buffers);
    Graph g = GraphLoader::fromStdin();

    aligned_vector <double> diag_values;
    for (int i = 0; i < g.nv(); i++) {
        double x;
        scanf("%lf", &x);
        diag_values.push_back(x);
    }

    int num_steps;
    scanf("%d", &num_steps);

    std::shared_ptr<PartialCholeskyFactorization> ldl;

    auto bench = make_benchmark(argc, argv, [&] () {
        ldl.reset(new PartialCholeskyFactorization(g, diag_values, num_steps));
    });

    // printf("Eliminate %d vertices\n", gvr->removal_count);

    // printf("Eliminate:\n");
    // for (int i = 0; i < gvr->removal_count; i++) {
    //     printf("%d\n", gvr->permutation[i]);
    // }

    // printf("Keep:\n");
    // for (int i = gvr->removal_count; i < g.nv(); i++) {
    //     printf("%d\n", gvr->permutation[i]);
    // }

    return 0;
}
