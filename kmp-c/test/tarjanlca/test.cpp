// Test for Tarjan LCA

#include <cstdlib>
#include <iostream>
#include <string>
#include <algorithm>
#include <algorithms/TarjanLCA.h>
#include <structures/TreeChildren.h>
#include <structures/GraphLoader.h>

using namespace std;

int main(int argc, char *argv[])
{
    atexit(mkl_free_buffers);

    // read graph
    Graph g = GraphLoader::fromStdin();

    // read spanning tree
    int n = g.nv();
    vector<int> parent(n);
    int root;
    for (int i = 0; i < n; i++) {
        cin >> parent[i];
        if (parent[i] == i) root = i;
    }

    // get children
    TreeChildren tChildren(n, parent.data());

    // run lca
    TarjanLCA lca(g, root, parent.data(), tChildren);

    // canonize lca results
    for (int i = 0; i < lca.ne; i++) {
        if (lca.u[i] > lca.v[i]) swap(lca.u[i], lca.v[i]);
    }

    vector<int> ord(lca.ne);
    for (int i = 0; i < lca.ne; i++) ord[i] = i;
    sort(ord.begin(), ord.end(), [&] (int a, int b) {
        if (lca.u[a] != lca.u[b]) return lca.u[a] < lca.u[b];
        return lca.v[a] < lca.v[b];
    });

    // read off sorted lca results
    for (int i = 0; i < lca.ne; i++) {
        printf("lca(%d, %d) = %d\n", lca.u[ord[i]], lca.v[ord[i]], lca.lca[ord[i]]);
    }

    return 0;
}
