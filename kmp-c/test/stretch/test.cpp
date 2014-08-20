// Test for Tarjan LCA

#include <cstdlib>
#include <iostream>
#include <string>
#include <algorithm>
#include <algorithms/TarjanLCA.h>
#include <algorithms/TreeDFS.h>
#include <algorithms/Stretch.h>
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

    // precompute children and dfs order
    TreeChildren tChildren(n, parent.data());
    vector<int> dfsOrder = TreeDFS(n, root, tChildren);

    // crudely get parent indices
    vector<int> parentIndex(n);
    for (int u = 0; u < n; u++) {
        if (u == root) continue;
        for (int i = 0; i < g.degree(u); i++) {
            if(g.neighbor(u, i) == parent[u]) parentIndex[u] = i;
        }
    }

    // run lca then stretch
    TarjanLCA lca(g, root, parent.data(), tChildren, dfsOrder.data());
    Stretch st(g, root, parentIndex.data(), dfsOrder.data(), lca);

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
        printf("stretch(%d, %d) = %.6lf\n", lca.u[ord[i]], lca.v[ord[i]], st.stretch[ord[i]]);
    }

    return 0;
}
