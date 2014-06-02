#ifndef COMPONENTS_H
#define COMPONENTS_H
// module used by randtree and randtreeimprove to keep track of node components
// linked list for specifying components
typedef struct CompCell_s
{
    int vert;
    struct CompCell_s* next;
} CompCell;

typedef struct Components_s
{
    CompCell* rawCells;   // where memory is actually allocated (1st rawCEll is for 1st vertex etc)
    CompCell** comps;     // linked list for 5th comp in 5th in array
    int* compSizes;   // number of verts in each comp
    int* vertComps;   // tells what component each vertex is in
} Components;


// just mallocs
Components* newComponents(int n);
Components* initGraphComponents(int n);
void freeComponents(Components* in);

#endif // COMPONENTS_H
