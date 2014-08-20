#pragma once

#include "structures/Graph.h"
#include "structures/TreeChildren.h"
#include "TarjanLCA.h"

struct Stretch {
    Stretch (const Graph& g, int tRoot, const int* tParentIndex, int* dfsOrder, const TarjanLCA& lca);
    ~Stretch ();

    double total;
    double* stretch;
};