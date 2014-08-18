#pragma once

#import <utility>

#import "util/aligned.h"
#import "structures/Graph.h"

class PartialCholeskyFactorization {
public:
	PartialCholeskyFactorization(const Graph &g, const aligned_vector<double> &diag_values, const int num_steps);
	~PartialCholeskyFactorization();

	void DoFactorization();

	void FactorizeDegreeOne(int u);
	void FactorizeDegreeTwoChains(int start, int stop, int outer_start, int outer_stop);

	EdgeList* SanitizeEdgeList(EdgeList edges);

	bool IsNeighbor(int u, int v);

	void AddToL(int u, int v, double weight);
	void AddToD(int u, int v, double weight);

	void InitL();
	void InitD();

	static aligned_vector<double> ApplyLInv(EdgeList L, aligned_vector<double> x);
	static aligned_vector<double> ApplyLTransInv(EdgeList L, aligned_vector<double> x);

	Graph GetReducedGraph(EdgeList D, int shift);

	Graph graph;
	aligned_vector<double> diag_values;
	int num_steps;
	int n;

	int *updated_degree;

	EdgeList *L;
	EdgeList *D;
};