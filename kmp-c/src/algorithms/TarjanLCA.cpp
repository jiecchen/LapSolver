#include <iostream>

#include "TarjanLCA.h"
#include "UnionFind.h"
#include "structures/Graph.h"
#include "structures/TreeChildren.h"

using namespace std;

TarjanLCA::TarjanLCA (const Graph& g, int tRoot, const int* tParent, const TreeChildren& tChildren) {
    // get dfs order
    int n = g.nv();
    int* dfsStack = new int[n]; dfsStack[0] = tRoot;
    int dfsStackPos = 1;
    int* dfsOrder = new int[n];
    int dfsOrderPos = 0;

    while (dfsStackPos > 0) {
        int curNode = dfsStack[--dfsStackPos];
        dfsOrder[dfsOrderPos++] = curNode;
        for (int i = tChildren.offset[curNode]; i < tChildren.offset[curNode+1]; i++) {
            dfsStack[dfsStackPos++] = tChildren.child[i];
        }
    }

    delete[] dfsStack;

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

        // push children
        for (int i = tChildren.offset[curNode]; i < tChildren.offset[curNode+1]; i++) {
            dfsStack[dfsStackPos++] = tChildren.child[i];
        }
    }

    delete[] ancestor;
    delete[] black;
    delete[] childrenLeft;
    delete[] dfsOrder;
}

TarjanLCA::~TarjanLCA () {
    delete[] u;
    delete[] v;
    delete[] lca;
}
