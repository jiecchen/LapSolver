#include "structures/Graph.h"
#include "structures/TreeChildren.h"
#include "TarjanLCA.h"
#include "Stretch.h"

Stretch::Stretch (const Graph& g, int tRoot, const int* tParentIndex, int* dfsOrder, const TarjanLCA& lca) {
    int n = g.nv();

    // compute resistances to root top-down
    double* resistanceToRoot = new double[n];
    resistanceToRoot[tRoot] = 0;
    for (int i = 1; i < n; i++) {
        int u = dfsOrder[i];
        int parent = g.neighbor(u, tParentIndex[u]);
        resistanceToRoot[u] = resistanceToRoot[parent] + 1/g.weight(u, tParentIndex[u]);
    }

    // get stretches
    total = 0;
    stretch = new double[lca.ne];
    for (int i = 0; i < lca.ne; i++) {
        stretch[i] = resistanceToRoot[lca.u[i]] + resistanceToRoot[lca.v[i]] - 2*resistanceToRoot[lca.lca[i]];
        total += stretch[i];
    }
}

Stretch::~Stretch () {
    delete[] stretch;
}
