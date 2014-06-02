// graphcluster.cpp
// Aaron Barnet 2007-07-09
// Clustering alg based on the ACL PageRank algorithm
// This version uses weights!
// See char* DOC_STRING for usage

#include "graph.hpp"
#include "dijkstra.hpp"
#include "weightedPageRank.hpp"
#include "profiler.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

using namespace boost;

const char DOC_STRING[] = "use: graphclusterpr CLUSTER_SIZE ALPHA EPSILON_DIVIDER MAX_CONDUCTANCE MAX_RECLUSTERS IJV_FILE OUT_FILE\n";

enum ordering_e { DISTANCE, WEIGHT, DEGREE };

class WeightCompare
{
 public:
    WeightCompare(distance_map_t& weightmap): weights(weightmap) {}
    bool operator()(const vertex_descriptor& x, const vertex_descriptor& y) const
	{
	    return weights[x] > weights[y];
	}
 private:
    distance_map_t& weights;
	};

const bool VERBOSE = true;
const enum ordering_e ORDERING = DISTANCE;

// see DOC_STRING for expected argv parameters
int main(int argc, char** argv)
{
    if (argc < 8)
    {
	printf(DOC_STRING);
	exit(1);
    }
    
    // user specified parameters
    int CLUSTER_SIZE = -1;
    float ALPHA = -1.0;
    float EPSILON_DIVIDER = -1.0;
    float MAX_CONDUCTANCE = -1.0;
    int MAX_RECLUSTERS = -1;
    
    const    char* IN_FILE = argv[6];
    const    char* OUT_FILE = argv[7];

    if (1 != sscanf(argv[1], "%d", &CLUSTER_SIZE))
	exit(1);
    if (1 != sscanf(argv[2], "%f", &ALPHA))
	exit(1);
    if (1 != sscanf(argv[3], "%f", &EPSILON_DIVIDER))
	exit(1);
    if (1 != sscanf(argv[4], "%f", &MAX_CONDUCTANCE))
	exit(1);
    if (1 != sscanf(argv[5], "%d", &MAX_RECLUSTERS))
	exit(1);
    

    Profiler profiler(VERBOSE);

    // get graph from file
    IJVclass<graph_t>* ijv = new IJVclass<graph_t>(IN_FILE);
    profiler.log("done reading file into memory.");
    graph_t* g_ptr = ijv->getBoostGraph();
    graph_t& g = *g_ptr;
    int num_verts = num_vertices(g);
    delete ijv;
    distance_map_t distancemap = get(vertex_distance, g);
    weight_map_t weightmap = get(edge_weight, g);        
    profiler.log("finished creating graph.");
    profiler.start_alg();
    
    std::vector<vertex_descriptor> ordered_nodes(num_verts);
    if (ORDERING == DISTANCE)
    {
	// set up Q
	dist_comp_class mycomp(distancemap);
	identity_property_map ident;
	MutableQueue Q(num_verts, mycomp, ident);
	profiler.log("created Q.");
    
	// find random vertex; then find furthest node from that
	srand(0);//time(NULL));
	int random = (int)floor((double)num_verts * (double)rand()/((double)(RAND_MAX) + 1.0));
	vertex_descriptor rand_v = random;
        
	// dijkstra it to find furthest node
	dijkstra(rand_v, g_ptr, &Q, &ordered_nodes);   
	vertex_descriptor start_node = ordered_nodes.back();
	profiler.log("ran dijkstra and found furthest node.");
    
	// Find out what order to cluster nodes in
	dijkstra(start_node, g_ptr, &Q, &ordered_nodes);
	profiler.log("ran dijkstra again.");	
    }
    // we order the nodes by their weighted degree (decreasing)
    else if (ORDERING == WEIGHT)
    {
	graph_traits < graph_t >::vertex_iterator i, end;
	for (tie(i, end) = vertices(g); i != end; ++i)
	{
	    vertex_descriptor v = *i;

	    float total_weight = 0;
	    out_edge_it edge_e, edge_end;
	    for (tie(edge_e, edge_end) = out_edges(v, g); edge_e != edge_end; ++edge_e)
	    {
		total_weight += weightmap[*edge_e];
	    }
	    distancemap[v] = total_weight;	    
	}

	for (int x = 0; x < num_verts; x++)
	{
	    ordered_nodes[x] = x;
	}

	std::sort(ordered_nodes.begin(), ordered_nodes.end(), WeightCompare(distancemap));
    }
    else if (ORDERING == DEGREE)
    {	
	graph_traits < graph_t >::vertex_iterator i, end;
	for (tie(i, end) = vertices(g); i != end; ++i)
	{
	    vertex_descriptor v = *i;
	    distancemap[v] = out_degree(v, g);
	}
	
	for (int x = 0; x < num_verts; x++)
	{
	    ordered_nodes[x] = x;
	}
	
	std::sort(ordered_nodes.begin(), ordered_nodes.end(), WeightCompare(distancemap));
    }
    else
    {
	exit(45);
    }
    // Now run weighted page rank
    WeightedPageRank wpr(g, profiler, ordered_nodes, CLUSTER_SIZE, ALPHA, EPSILON_DIVIDER, MAX_CONDUCTANCE, MAX_RECLUSTERS);
    profiler.log("done with weighted page rank.");
    profiler.end_alg();

    // Write out our clusters
    FILE* fp_out;
    if ((fp_out = fopen(OUT_FILE, "wt")) == NULL)
	exit(1);
    binWriteArray(wpr.get_clusters(), num_verts, fp_out);
    profiler.write_clust_info(fp_out);
    fclose(fp_out);

    // Time to free memory!
    delete g_ptr;

    profiler.log("done outputting to file.");
    profiler.cluster_stats();
    profiler.log_alg();

    return 0;
}
