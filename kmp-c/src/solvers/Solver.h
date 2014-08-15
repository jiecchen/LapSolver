#pragma once
#include <util/aligned.h>
#include <structures/graph.h>
#include "LinearOperator.h"

class Solver : public LinearOperator
{
public:
	explicit Solver(const Graph &g);
	explicit Solver(Graph &&g);
	~Solver() {}

	virtual void apply(double *) = 0;
    virtual int getDimension() const { return dim; }

protected:
	Graph g;
	const int dim;
	aligned_vector<double> d;
};