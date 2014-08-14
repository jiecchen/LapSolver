#include "PartialCholeskyOrder.h"

PartialCholeskyOrder::PartialCholeskyOrder(const Graph &g) : graph(g) {
	n = g.nv;

	permutation = new int[n];
	removed = new int[n];
	un_removable = new int[n];
	updated_deg = new int[n];

	vector <int> deg1ordering = DegreeOneOrdering();
	vector <int> deg2ordering = DegreeTwoOrdering();

	vector <int> all_ordering; all_ordering.reserve(n);
	all_ordering.insert(all_ordering.end(), deg1ordering.begin(), deg1ordering.end());
	all_ordering.insert(all_ordering.end(), deg2ordering.begin(), deg2ordering.end());

	ConstructFinalOrdering(all_ordering);
}

PartialCholeskyOrder::~PartialCholeskyOrder() {
	delete[] permutation;
	delete[] removed;
	delete[] updated_deg;
	delete[] un_removable;
}

vector <int> PartialCholeskyOrder::DegreeTwoOrdering() {
	vector <int> ordering = CheckIfCycle();
	if (ordering.size() > 0)
		return ordering;

	return DegreeTwoChainOrdering();
}

vector <int> PartialCholeskyOrder::DegreeTwoChainOrdering() {
	for (int i = 0; i < n; i++)
		if (updated_deg[i] > 2)
			un_removable[i] = 1;

	vector <int> chains; chains.reserve(n);

	for (int i = 0; i < n; i++)
		if (un_removable[i] == 1) {
			for (int j = 0; j < graph.getDegree(i); j++) {
				int u = graph.getNeighbors(i)[j];

				if (updated_deg[u] == 2 && removed[u] == 0) {
					// This is a degree two chain
					while (updated_deg[u] == 2) {
						removed[u] = 1;
						chains.push_back(u);

						int v = u;
						for (int k = 0; k < graph.getDegree(u); k++) {
							if (updated_deg[graph.getNeighbors(u)[k]] == 2 && removed[graph.getNeighbors(u)[k]] == 0)
								v = graph.getNeighbors(u)[k];
						}

						if (v == u) {
							// Try to exit the chain 
							for (int k = 0; k < graph.getDegree(u); k++) {
								if (updated_deg[graph.getNeighbors(u)[k]] > 2 && graph.getNeighbors(u)[k] != i)
									v = graph.getNeighbors(u)[k];
							}

							if (v == u) {
								// The chain is connected to the same vertex
								v = i;
							}
						}

						u = v;
					}

					if (u == i) {
						// If the chain is connected to the same outer vertex at both ends, then its end should be counted only once
						un_removable[*chains.rbegin()] = 1;
						chains.pop_back();
					}
				}
			}
		}

	return chains;
}

vector <int> PartialCholeskyOrder::CheckIfCycle() {
	vector <int> v; v.reserve(n);

	int max_deg = 0;
	for (int i = 0; i < n; i++)
		if (updated_deg[i] > max_deg)
			max_deg = updated_deg[i];

	if (max_deg > 2 || max_deg == 0)
		return v;

	// Add all but two vertices from the cycle (so that Cholesky gives a valid output)
	int last = 0;
	for (int i = 0; i < n; i++)
		if (updated_deg[i] == 2)
			last = i;

	for (int i = last - 1; i >= 0; i--)
		if (updated_deg[i] == 2) {
			last = i;
			break;
		}

	for (int i = 0; i < last; i++)
		if (updated_deg[i] == 2)
			v.push_back(i);

	return v;
}

vector <int> PartialCholeskyOrder::DegreeOneOrdering() {
	vector <int> ordering; ordering.reserve(n);
	for (int i = 0; i < n; i++)
		if (graph.getDegree(i) == 1) {
			ordering.push_back(i);
			removed[i] = 1;
		}

	vector <int>::iterator it = ordering.begin();
	while (it != ordering.end()) {
		int u = *it;

		for (int i = 0; i < graph.getDegree(u); i++) {
			int v = graph.getNeighbors(u)[v];
			updated_deg[v]--;

			if (removed[v] == 0 && updated_deg[v] == 1) {
				ordering.push_back(v);
				removed[v] = 1;
			}
		}

		++it;
	}

	if (ordering.size() == n) {
		// In case the graph is a tree, remove all vertices but one
		removed[ordering[n - 1]] = 0;
		updated_deg[ordering[n - 1]] = 0;
		ordering.pop_back();
	}

	return ordering;
}

vector <int> PartialCholeskyOrder::ConstructFinalOrdering(vector <int> ordering) {
	int *use = new int[n];

	int index = 0;
	for (vector <int>::iterator it = ordering.begin(); it != ordering.end(); ++it) {
		use[*it] = 1;
		permutation[index++] = *it;
	}
	removal_count = index;

	for (int i = 0; i < n; i++)
		if (use[i] == 0)
			permutation[index++] = i;
}

