#pragma once

#include "structures/Graph.h"
#include "structures/TreeChildren.h"

struct TarjanLCA {
    TarjanLCA (const Graph& g, int tRoot, const int* tParent, const TreeChildren& tChildren);
    ~TarjanLCA ();

    int ne;
    int* u;
    int* v;
    int* lca;
};