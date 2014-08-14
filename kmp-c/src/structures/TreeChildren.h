// TreeChildren: a structure to compute and store child arrays from a parent array.
// child[offset[v], offset[v]+1, ..., offset[v+1]-1] are the children of vertex v.

#pragma once

struct TreeChildren {
    TreeChildren(int nv, const int* parent);
    ~TreeChildren();

    int* child;
    int* offset;
};
