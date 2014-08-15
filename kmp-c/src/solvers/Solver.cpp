#include "Solver.h"

Solver::Solver(const Graph &g)
    : g(g),
      dim(g.nv)
{

}

Solver::Solver(Graph &&g)
    : g(g),
      dim(g.nv)
{

}