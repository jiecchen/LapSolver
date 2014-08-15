#include "TreeChildren.h"

TreeChildren::TreeChildren(int nv, const int* parent) {
    // offset = cumulative sums of child counts
    offset = new int[nv+1];
    offset[0:nv+1] = 0;

    for (int i = 0; i < nv; i++) {
        if(parent[i] != i) offset[parent[i]+1]++;
    }
    for (int i = 1; i <= nv; i++) {
        offset[i] += offset[i-1];
    }

    // fill child array with these offsets
    child = new int[nv];
    int* childIndex = new int[nv];
    childIndex[0:nv] = 0;

    for (int i = 0; i < nv; i++) {
        int p = parent[i];
        if (p != i) {
            child[offset[p] + childIndex[p]++] = i;
        }
    }

    delete[] childIndex;
}

TreeChildren::~TreeChildren() {
    delete[] child;
    delete[] offset;
}