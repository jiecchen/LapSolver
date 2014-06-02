#include "components.h"
#include <stdlib.h>

// just mallocs
Components* newComponents(int n)
{
    Components* out = malloc(sizeof(Components));
    out->rawCells = (CompCell*)malloc(sizeof(CompCell) * n);
    out->comps = (CompCell**)malloc(sizeof(CompCell*) * n);
    out->compSizes = (int*)malloc(sizeof(int) * n);
    out->vertComps = (int*)malloc(sizeof(int) * n);
    return out;
}


void freeComponents(Components* in)
{
    free(in->rawCells);
    free(in->comps);
    free(in->compSizes);
    free(in->vertComps);
    free(in);
}


// given the size of a graph (verts) initialize a components struct so that each vert gets its own comp
Components* initGraphComponents(int n)
{
    Components* out = newComponents(n);
    int ctr = 0;
    for (ctr = 0; ctr < n; ctr++)
    {
	out->rawCells[ctr].vert = ctr;
	out->rawCells[ctr].next = NULL;
	out->comps[ctr] = &(out->rawCells[ctr]);
	out->vertComps[ctr] = ctr;
	out->compSizes[ctr] = 1;
    }
    return out;
}
