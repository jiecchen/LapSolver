#include <vector>
#include "structures/TreeChildren.h"

using namespace std;

vector<int> TreeDFS(int nv, int root, const TreeChildren& children) {
    int* dfsStack = new int[nv]; dfsStack[0] = root;
    int dfsStackPos = 1;
    vector<int> dfsOrder;
    dfsOrder.reserve(nv);

    while (dfsStackPos > 0) {
        int curNode = dfsStack[--dfsStackPos];
        dfsOrder.push_back(curNode);
        for (int i = children.offset[curNode]; i < children.offset[curNode+1]; i++) {
            dfsStack[dfsStackPos++] = children.child[i];
        }
    }

    delete[] dfsStack;

    return dfsOrder;
}
