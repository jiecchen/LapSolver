// Test for TreeChildren

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <structures/TreeChildren.h>
#include <util/Benchmark.h>

using namespace std;

int main(int argc, char *argv[])
{
    // read a parent array from stdin
    const int MAX_N = 1000;
    int nv, parent[MAX_N], root;
    cin >> nv;
    for (int i = 0; i <= nv; i++) {
        cin >> parent[i];
        if (parent[i] == i) root = i;
    }

    printf("Root: %d\n", root);

    // convert to child array
    std::shared_ptr<TreeChildren> tChildren;
    auto bench = make_benchmark(argc, argv, [&] () {
        tChildren.reset(new TreeChildren(nv, parent));
    });

    for (int i = 0; i < nv; i++) {
        printf("%d:", i);
        for (int j = tChildren->offset[i]; j < tChildren->offset[i+1]; j++) {
            printf(" %d", tChildren->child[j]);
        }
        printf("\n");
    }

    return 0;
}
