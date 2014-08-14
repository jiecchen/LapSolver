// TreeChildList: a structure to compute and store child arrays from a parent array.
// child[offset[v], offset[v]+1, ..., offset[v+1]-1] are the children of vertex v.

struct TreeChildList {
    TreeChildList(int nv, const int* parent);
    ~TreeChildList();

    int* child;
    int* offset;
};
