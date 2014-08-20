#include <iostream>

#include "structures/Graph.h"
#include "structures/TreeChildren.h"
#include "TarjanLCA.h"
#include "UnionFind.h"
#include "TreeDFS.h"

using namespace std;

TarjanLCA::TarjanLCA (const Graph& g, int tRoot, const int* tParent, const TreeChildren& tChildren, int* dfsOrder) {
    int n = g.nv();

    // run LCA
    int* ancestor = new int[n];
    bool* black = new bool[n];
    int* childrenLeft = new int[n];

    UnionFind components(n);

    int edgePos = 0;
    ne = g.ne() - n + 1;
    u = new int[ne];
    v = new int[ne];
    lca = new int[ne];

    black[0:n] = false; 
    childrenLeft[0:n] = tChildren.offset[1:n] - tChildren.offset[0:n];
    for (int i = 0; i < n; i++) ancestor[i] = i;

    for (int i = n-1; i >= 0; i--) {
        int curNode = dfsOrder[i];
        int parent = tParent[curNode];

        if (curNode != tRoot) {
            childrenLeft[parent]--;
        }

        // done with all children
        if (!childrenLeft[curNode]) {
            black[curNode] = true;

            for (int i = 0; i < g.degree(curNode); i++) {
                int queryNode = g.neighbor(curNode, i);
                if (black[queryNode] && tParent[queryNode] != curNode && tParent[curNode] != queryNode) {
                    u[edgePos] = curNode;
                    v[edgePos] = queryNode;
                    lca[edgePos] = ancestor[components.find_set(queryNode)];
                    edgePos++;
                }
            }

            components.link(parent, curNode);
            ancestor[components.find_set(parent)] = parent;
        }
    }

    delete[] ancestor;
    delete[] black;
    delete[] childrenLeft;
}

TarjanLCA::~TarjanLCA () {
    delete[] u;
    delete[] v;
    delete[] lca;
}
