/* Copyright 2014 Yale Institute for Network Science
 * Author: Cyril Zhang
 */
#include <cstdlib>
#include <string>
#include <functional>
#include <memory>
#include <structures/GraphLoader.h>
#include <algorithms/ShortestPathTree.h>
#include <util/Benchmark.h>
#include <algorithms/UnionFind.h>
#include <lsst/StarDecompositionTree.h>

int main(int argc, char **argv)
{
    atexit(mkl_free_buffers);
    Graph g = GraphLoader::fromStdin();

    std::shared_ptr<StarDecompositionTree> sdt;
    int source = 0;

    auto bench = make_benchmark(argc, argv, [&] () {
        sdt.reset(new StarDecompositionTree(g, source));
    });

    UnionFind uf(g.nv());

    // check that it's actually a spanning tree: add edges to union-find
    for (int i = 0; i < g.nv(); ++i) {
        if (i == source) continue;
        int j = g.neighbor(i, sdt->parentIndex[i]);
        uf.link(i, j);
    }

    // check that everyone's in the same component
    int c0 = uf.find_set(0);
    bool bad = false;
    for (int i = 1; i < g.nv(); ++i) {
        if (uf.find_set(i) != c0) {
            bad = true;
            printf ("SHIT\n");
            break;
        }
    }

    if (!bad) printf("OK\n");

    return 0;
}
