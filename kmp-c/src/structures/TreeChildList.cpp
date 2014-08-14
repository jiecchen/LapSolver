#include <cstring>

#include "TreeChildList.h"

TreeChildList::TreeChildList(int nv, const int* parent) {
    // offset = cumulative sums of child counts
    offset = new int[nv+1];
    memset(offset, 0, (nv+1)*sizeof(*offset));
    for (int i = 0; i < nv; i++) {
        if(parent[i] != i) offset[parent[i]+1]++;
    }
    for (int i = 1; i <= nv; i++) {
        offset[i] += offset[i-1];
    }

    // fill child array with these offsets
    child = new int[nv];
    int* childIndex = new int[nv];
    std::memset(childIndex, 0, sizeof(*childIndex) * nv);
    for (int i = 0; i < nv; i++) {
        int p = parent[i];
        if (p != i) {
            child[offset[p] + childIndex[p]++] = i;
        }
    }

    delete[] childIndex;
}

TreeChildList::~TreeChildList() {
    delete[] child;
    delete[] offset;
}