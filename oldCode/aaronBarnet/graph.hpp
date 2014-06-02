#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <cstdio>

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/pending/mutable_queue.hpp>
using namespace boost;

// boost typedef shortcuts
typedef adjacency_list <listS, vecS, undirectedS, property <vertex_distance_t, float>,
			property <edge_weight_t, float> > graph_t;

typedef property_map<graph_t, edge_weight_t>::type weight_map_t;
typedef graph_traits< graph_t >::vertex_descriptor vertex_descriptor;
typedef graph_traits< graph_t >::edge_descriptor edge_descriptor;
typedef std::pair<int, int> edge_pair_t;
typedef graph_traits<graph_t>::out_edge_iterator out_edge_it;
typedef graph_traits<graph_t>::adjacency_iterator adjacency_it;
typedef property_map<graph_t, vertex_distance_t>::type distance_map_t;
typedef property_map<graph_t, edge_weight_t>::type weight_map_t;

enum colors { WHITE, GREY, BLACK };   // white = unvisited; grey = in Q; black = already processed

const float SMIDGEN = .0000000000000001;
const float BIG_FLOAT = 3.3e38;


// class to read in graph from binary file
template <class G>
class IJVclass
{
public:
    IJVclass(const char* fileName) : n(0), nnz(0), i(NULL), j(NULL), v(NULL)
	{

	    FILE* fpr;
	    if ((fpr = fopen(fileName, "rb")) == NULL)
		exit(1);

	    if (!fread(&(n), sizeof(int), 1, fpr))
		exit(1);//cError ("error reading n from file\n");
    
	    if (!fread(&(nnz), sizeof(int), 1, fpr))
		exit(1);//cError ("error reading nnz from file\n");
    
	    i = new int[nnz];
	    j = new int[nnz];
	    v = new double[nnz];

	    if (!fread(i, sizeof(int), nnz, fpr))
		exit(1); //cError ("error reading i from file\n");
    
	    if (!fread(j, sizeof(int), nnz, fpr))
		exit(1); //cError ("error reading j from file\n");
    
	    if (!fread(v, sizeof(double), nnz, fpr))
		exit(1); // cError ("error reading v from file\n");
    
	    for (int x = 0; x < nnz; x++)
	    {
		i[x] -= 1;
		j[x] -= 1;
	    }

	    fclose(fpr);
	}
    ~IJVclass()
	{
	    if (i != NULL)
		delete [] i;
	    if (j != NULL)
		delete [] j;
	    if (v != NULL)
		delete [] v;
	}

    G* getBoostGraph()
	{ 
	    G* g_ptr = new G(n);
	    G& g = *g_ptr;
	    //   EdgeWeight weightmap = get(edge_weight, g);

	    for (int x = 0; x < nnz; x++)
	    {
		edge_descriptor e;
		bool inserted;
	
		tie(e, inserted) = add_edge(i[x], j[x], v[x], g);
//	weightmap[e] = v[x];
	    }
    
	    return g_ptr;
	}

    int n;
    int nnz;
    int *i;
    int *j;
    double *v;
};


// from before IJVclass was templatized
// class IJVclass
// {
// public:
//     IJVclass(const char* fileName);
//     ~IJVclass()
// 	{
// 	    if (i != NULL)
// 		delete [] i;
// 	    if (j != NULL)
// 		delete [] j;
// 	    if (v != NULL)
// 		delete [] v;
// 	}

//     graph_t* getBoostGraph();
    
//     int n;
//     int nnz;
//     int *i;
//     int *j;
//     double *v;
// };

// // function to write out graph to file
// template <typename G, typename G_EDGE_IT>
// void binWriteGraph(G& graph, const char* file_name)
// {
    
//     int n_verts = num_vertices(graph);
//     int n_edges = num_edges(graph);
//     int* i = new int[n_edges];
//     int* j = new int[n_edges];
//     float* v = new float[n_edges];
//     int edge_ctr = 0;
//     property_map<G, edge_weight_t> weightmap = get(edge_weight, graph);
//     //graph_traits<G>::edge_iterator e, end;
//     G_EDGE_IT e, end;
//     for (tie(e, end) = edges(graph); e != end; e++)
//     {
// 	//graph_traits<G>::vertex_descriptor unode, vnode;
// 	unode = source(*e, graph);
// 	vnode = target(*e, graph);
// 	float weight = weightmap[*e];
// 	i[edge_ctr] = unode;
// 	j[edge_ctr] = vnode;
// 	v[edge_ctr] = weight;
// 	edge_ctr++;
//     }

//     FILE* fp;
//     if ((fp = fopen(file_name, "wt")) == NULL)
// 	exit(1);

//     fwrite(&n_verts, sizeof(int), 1, fp);
//     fwrite(&n_edges, sizeof(int), 1, fp);
//     fwrite(i, sizeof(int), n_edges, fp);
//     fwrite(j, sizeof(int), n_edges, fp);
//     fwrite(v, sizeof(float), n_edges, fp);
//     fclose(fp);

//     delete [] i;
//     delete [] j;
//     delete [] v;
// }


// used by graphclusterwpr
template <typename T>
void binWriteArray(T* array, int n, FILE* fp) //const char* const fileName)
{
    //FILE* fp;
    //if ((fp = fopen(fileName, "wt")) == NULL)
//	exit(1);
    
    fwrite(&n, sizeof(int), 1, fp);
    fwrite(array, sizeof(T), n, fp);
    
    //  fclose(fp);
}



// used by graphclusterpr
template <typename T>
void binWriteArray(T* array, int n, const char* const fileName)
{
    FILE* fp;
    if ((fp = fopen(fileName, "wt")) == NULL)
	exit(1);
    
    fwrite(&n, sizeof(int), 1, fp);
    fwrite(array, sizeof(T), n, fp);
    
    fclose(fp);
}


#endif // GRAPH_HPP
