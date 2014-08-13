#include "UnionFind.h"

UnionFind::UnionFind(int n) {
    parent = new int[n];
    for (int i = 0; i < n; i++) {
        parent[i] = i;
    }
}

UnionFind::~UnionFind() {
    delete[] parent;
}

void UnionFind::link(int u, int v) {
    int pu = find_set(u), pv = find_set(v);
    if (pu != pv) parent[pu] = pv;
}

int UnionFind::find_set(int u) {
    if (parent[u] == u) return u;
    return parent[u] = find_set(parent[u]);
}
