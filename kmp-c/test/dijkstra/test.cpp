/* Copyright 2014 Yale Institute for Network Science
 * Author: Alex Reinking
 *
 * Example benchmark for vector arithmetic
 */
#include <cstdlib>
#include <string>
#include <functional>
#include <memory>
#include <structures/GraphLoader.h>
#include <algorithms/ShortestPathTree.h>
#include <util/Benchmark.h>

int main(int argc, char **argv)
{
    atexit(mkl_free_buffers);
    Graph g = GraphLoader::fromStdin();

    std::shared_ptr<ShortestPathTree> spt;

    auto bench = make_benchmark(argc, argv, [&] () {
        spt.reset(new ShortestPathTree(g, 0));
    });

    for (int i = 0; i < g.nv(); ++i)
        printf("%d,%d ", spt->parent[i], i);
    printf("\n");

    return 0;
}
