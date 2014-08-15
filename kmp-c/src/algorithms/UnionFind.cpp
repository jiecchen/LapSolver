#include "UnionFind.h"

UnionFind::UnionFind(int n)
{
    to_compress = new int[n];
    parent = new int[n];
    for (int i = 0; i < n; i++)
        parent[i] = i;
}

UnionFind::~UnionFind()
{
    delete[] parent;
}

void UnionFind::link(int u, int v)
{
    int pu = find_set(u), pv = find_set(v);
    if (pu != pv) parent[pu] = pv;
}

int UnionFind::find_set(int u)
{
    int pos = 0;

    // follow pointers until self-loop
    while (parent[u] != u)
    {
        to_compress[pos++] = u;
        u = parent[u];
    }

    // Path compression
    parent[to_compress[0:pos]] = u;
    return u;
}
