#pragma once
#include <util/aligned.h>
#include "structures/Graph.h"
#include "VectorOperator.h"

class LaplacianOperator : public VectorOperator
{
public:
    using VectorOperator::apply;
    explicit LaplacianOperator(const Graph &g)
        : g(g), dim(g.nv) {}
    explicit LaplacianOperator(Graph &&g)
        : g(g), dim(g.nv) {}
    ~LaplacianOperator() {}

    virtual void apply(double *x, double *b) {

    }

    virtual int getDimension() const
    {
        return dim;
    }

protected:
    Graph g;
    const int dim;
    aligned_vector<double> d;
};