#include "PartialCholeskyFactorization.h"


PartialCholeskyFactorization::PartialCholeskyFactorization(const Graph &g, const aligned_vector<double> &diag_values, const int num_steps) 
	: graph(g), 
	  diag_values(diag_values), 
	  num_steps(num_steps) {
	
	n = graph.nv();
	updated_degree = new int[n];

	InitL();
	InitD();

	DoFactorization();
}

PartialCholeskyFactorization::~PartialCholeskyFactorization() {
	delete[] updated_degree;
	delete L;
	delete D;
}

void PartialCholeskyFactorization::DoFactorization() {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < graph.degree(i); j++)
			diag_values[i] = diag_values[i] + graph.weight(i, j);
		updated_degree[i] = graph.degree(i);

		AddToD(i, i, diag_values[i]);
	}

	for (int i = 0; i < num_steps; i++)
		if (updated_degree[i] == 1)
			FactorizeDegreeOne(i);
		else {
			// Will eliminate chains starting at 'start' and ending at 'stop'
			int start = i;
			int stop = i;
			while (stop < num_steps - 1 && IsNeighbor(stop, stop + 1))
				stop++;	

			// outer_start and outer_stop represent the higher degree vertices connected to the chain
            int outer_start = 0, outer_stop = 0;

            // pick outer_start and outer_stop
            if (start == stop) {
            	for (int j = 0; j < graph.degree(start); j++)
            		if (graph.neighbor(start, j) >= num_steps) {
            			outer_start = graph.neighbor(start, j);

            			for (int k = j + 1; k < graph.degree(start); k++) 
            				if (graph.neighbor(start, k) >= num_steps) {
            					outer_stop = graph.neighbor(start, k);
            					break;
            				}

            			break;
            		}
            } else {
            	for (int j = 0; j < graph.degree(start); j++)
            		if (graph.neighbor(start, j) >= num_steps) {
            			outer_start = graph.neighbor(start, j);
            			break;
            		}

            	for (int j = 0; j < graph.degree(stop); j++) 
            		if (graph.neighbor(stop, j) >= num_steps) {
            			outer_stop = graph.neighbor(stop, j);
            			break;
            		}
            }

            FactorizeDegreeTwoChains(start, stop, outer_start, outer_stop);

            i = stop;
		}
}

void PartialCholeskyFactorization::FactorizeDegreeTwoChains(int start, int stop, int outer_start, int outer_stop) {
	// outer_value is used to compute the L values of type (outer_start, u)
	double outer_value = -1;
	for (int i = 0; i < graph.degree(start); i++)
		if (graph.neighbor(start, i) == outer_start)
			outer_value = -1.0 * graph.weight(start, i);

	// new_value_in_lap is used to mimic the update of the laplacian matrix
	double new_value_in_lap = outer_value;

	// go through the chain and update values
	for (int u = start; u <= stop; u++) {
		outer_value = outer_value / diag_values[u];

		diag_values[outer_start] = diag_values[outer_start] - new_value_in_lap * outer_value;

		AddToL(outer_start, u, outer_value);
		AddToD(outer_start, outer_start, -new_value_in_lap * outer_value);

		for (int i = 0; i < graph.degree(u); i++) {
			int v = graph.neighbor(u, i);
			double weight = graph.weight(u, i);

			if (updated_degree[v] > 1 && v > u && v != outer_start) {
				diag_values[v] = diag_values[v] - weight * weight / diag_values[u];
				AddToL(v, u, -weight / diag_values[u]);
				AddToD(v, v, -weight * weight / diag_values[u]);

				outer_value = outer_value * weight;
				new_value_in_lap = new_value_in_lap * weight / diag_values[u];
			}
		}
	}

	AddToD(outer_stop, outer_start, new_value_in_lap);
	AddToD(outer_start, outer_stop, new_value_in_lap);
}

void PartialCholeskyFactorization::FactorizeDegreeOne(int u) {
	for (int i = 0; i < graph.degree(u); i++) {
		int v = graph.neighbor(u, i);
		double weight = graph.weight(u, i);

		if (u < v) {
			AddToL(v, u, -weight / diag_values[u]);
			AddToD(v, v, -weight * weight / diag_values[u]);

			diag_values[v] = diag_values[v] - weight * weight / diag_values[u];

			updated_degree[v]--;
		}
	}
}

EdgeList* PartialCholeskyFactorization::SanitizeEdgeList(EdgeList edges) {
	int size = edges.ne;
	aligned_vector <int> u(size);
	aligned_vector <int> v(size);
	aligned_vector <double> w(size);

	// move edges' elements to the three vectors, and get rid of (0,0,0) triplets
	for (int i = 0; i < size; i++)
		if (!(edges.u[i] == 0 && edges.v[i] == 0 && edges.w[i] == 0)) {
			u.push_back(edges.u[i]);
			v.push_back(edges.v[i]);
			w.push_back(edges.w[i]);
		}
	size = u.size();

	// sort the remaining values
	int index[size];
	for (int i = 0; i < size; i++)
		index[i] = i;

	std::sort(index, index + size, [&] (int lhs, int rhs) {
		if (edges.u[lhs] != edges.u[rhs])
			return edges.u[lhs] < edges.u[rhs];
		return edges.v[lhs] < edges.v[rhs];
	});

	aligned_vector <int> fin_u(size);
	aligned_vector <int> fin_v(size);
	aligned_vector <double> fin_w(size);

	int last = 0;
	for (int i = 0; i < size; i++)
		if (i > 0 && u[i] == *fin_u.rbegin() && v[i] == *fin_v.rbegin()) {
			*fin_u.rbegin() += w[index[i]];
		}
		else {
			fin_u.push_back(u[index[i]]);
			fin_v.push_back(v[index[i]]);
			fin_w.push_back(w[index[i]]);
		}

	fin_u.resize(fin_u.size());
	fin_v.resize(fin_v.size());
	fin_w.resize(fin_w.size());

	return new EdgeList(std::move(fin_u), std::move(fin_v), std::move(fin_w));
}

bool PartialCholeskyFactorization::IsNeighbor(int u, int v) {
	for (int i = 0; i < graph.degree(u); i++)
		if (graph.neighbor(u, i) == v)
			return true;
	return false;
}

void PartialCholeskyFactorization::InitL() {
	int cnt = n;
	for (int i = 0; i < num_steps; i++)
		cnt = cnt + graph.degree(i) + 1;

	L = new EdgeList(cnt);
	for (int i = 0; i < n; i++)
		AddToL(i, i, 1);
}

void PartialCholeskyFactorization::InitD() {
	int cnt = n + 2 * num_steps;
	for (int i = 0; i < num_steps; i++)
		cnt = cnt + 2 * graph.degree(i);

	D = new EdgeList(cnt);

	for (int i = num_steps; i < n; i++)
		for (int j = 0; j < graph.degree(i); j++)
			if (graph.neighbor(i, j) > i) {
				AddToD(i, graph.neighbor(i, j), -graph.weight(i, j));
				AddToD(graph.neighbor(i, j), i, -graph.weight(i, j));
			}
}

void PartialCholeskyFactorization::AddToL(int u, int v, double weight) {
	L->u.push_back(u);
	L->v.push_back(v);
	L->w.push_back(weight);
}

void PartialCholeskyFactorization::AddToD(int u, int v, double weight) {
	D->u.push_back(u);
	D->v.push_back(v);
	D->w.push_back(weight);
}

static aligned_vector<double> ApplyLInv(EdgeList L, aligned_vector<double> x) {
	int size = L.ne;
	int index[size];
	for (int i = 0; i < size; i++)
		index[i] = i;

	// sort L
	std::sort(index, index + size, [&] (int lhs, int rhs) {
		if (L.v[lhs] != L.v[rhs])
			return L.v[lhs] < L.v[rhs];
		return L.u[lhs] < L.u[rhs];
	});

	for (int i = 0; i < size; i++)
		x[L.u[index[i]]] -= L.w[index[i]] * x[L.v[index[i]]];

	return x;
}

static aligned_vector<double> ApplyLTransInv(EdgeList L, aligned_vector<double> x) {
	int size = L.ne;
	int index[size];
	for (int i = 0; i < size; i++)
		index[i] = i;

	// sort L
	std::sort(index, index + size, [&] (int lhs, int rhs) {
		if (L.v[lhs] != L.v[rhs])
			return L.v[lhs] > L.v[rhs];
		return L.u[lhs] > L.u[rhs];
	});

	for (int i = 0; i < size; i++)
		x[L.v[index[i]]] -= L.w[index[i]] * x[L.u[index[i]]];

	return x;
}

Graph PartialCholeskyFactorization::GetReducedGraph(EdgeList D, int shift) {
	EdgeList* sanitized_D = SanitizeEdgeList(D);

	aligned_vector <int> edges_to_add(sanitized_D->ne);

	for (int i = 0; i < sanitized_D->ne; i++) {
		if (sanitized_D->u[i] >= sanitized_D->v[i]) continue;
		if (sanitized_D->u[i] >= shift && sanitized_D->v[i] >= shift) {
			edges_to_add.push_back(i);
		}
	}

	EdgeList reduced_sparsifier_edges = EdgeList(edges_to_add.size());

	int position = 0;
	for (int i : edges_to_add) {
		reduced_sparsifier_edges.u[position] = sanitized_D->u[i] - shift;
		reduced_sparsifier_edges.v[position] = sanitized_D->v[i] - shift;
		reduced_sparsifier_edges.w[position] = -sanitized_D->w[i];
		position++;
	}

	return Graph(std::move(reduced_sparsifier_edges));
}


