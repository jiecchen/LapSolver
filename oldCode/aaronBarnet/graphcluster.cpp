// graphcluster.cpp
// Aaron Barnet 2007-07-02
// Implementation of Dan Spielman's graph clustering heuristic:
// 1) pick any node and find furthest node away (call it start_node) (Dijkstra)
// 2) While there are nodes without cluster:
//       Find unclustered node closest to start_node and start building list of closest
//         nodes to it using modified dijkstra (when find second path to a node the new dist is 1/(1/current_dist + 1/new_dist)).
//         Never let distance become less than edge weight you are looking at. Otherwise same as Dijkstra.
//       Go through this list of nodes and find group of elements indexed 0, 1,..., n that maximized EDGES_IN_CLUSTER/NODES_IN_CLUSTER
//       Make this group a cluster

// See char* DOC_STRING for usage

#include "graph.hpp"
#include "dijkstra.hpp"
//#include <boost/config.hpp>
  //#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
   //#include <algorithm>
//#include <boost/graph/graph_traits.hpp>
//#include <boost/graph/adjacency_list.hpp>
   //#include <boost/pending/indirect_cmp.hpp>
//#include <boost/pending/mutable_queue.hpp>
   //#include <boost/graph/dijkstra_shortest_paths.hpp>

using namespace boost;

// // boost typedef shortcuts
// typedef adjacency_list <listS, vecS, undirectedS, property <vertex_distance_t, float>,
// 			property <edge_weight_t, float> > graph_t;

// typedef graph_traits< graph_t >::vertex_descriptor vertex_descriptor;
// typedef graph_traits< graph_t >::edge_descriptor edge_descriptor;
// //typedef std::pair<int, int> Edge;

// typedef property_map<graph_t, vertex_distance_t>::type distance_map_t;
// typedef property_map<graph_t, edge_weight_t>::type weight_map_t;


// // functor: used by priority queue to organize heap by distance map
// class mycomp_class
// {
// public:
//     mycomp_class(distance_map_t distance_map) : d_map(distance_map)
// 	{

// 	}
    
//     bool operator()(const vertex_descriptor& x,const vertex_descriptor& y) const
// 	{
// 	    return d_map[x] < d_map[y];
// 	}

// private:
//     distance_map_t d_map;
// };
 
// // the priority queue type
// typedef mutable_queue<vertex_descriptor, std::vector<vertex_descriptor>, mycomp_class , identity_property_map> MutableQueue;


// class NodeWithDist
// {
// public:
//     //NodeWithDist(float d, vertex_descriptor v) : dist(d), node(v)
// //	{}
//     float dist;
//     vertex_descriptor node;
// };

// bool operator<(const NodeWithDist& x, const NodeWithDist& y)
// {
//     return x.dist < y.dist;
// }

const char DOC_STRING[] = "use: graphcluster CLUSTER_SIZE IJV_FILE OUT_FILE\n";
//enum colors { WHITE, GREY, BLACK };   // white = unvisited grey= in Q, black = done with

// WARNING: ordered_nodes must already be big enough to hold all vertex descriptors in graph!
// custom version of dijkstra's shortest paths algorithm for boost graphs
// inputs:
//    - s: vertex descriptor that specifies starting point for alg
//    - g: pointer to boost graph (adjacency list) that alg is to work on
//    - Q: pointer to (empty) priority Q. Must order nodes by distance (must already have access to distance map) and support: push, update, pop, and top
//    - ordered_nodes: pointer to vector of nodes. Nodes will be placed in vector according to distance from s
// returns: None
// void dijkstra(vertex_descriptor s, graph_t* g, MutableQueue* Q, std::vector<vertex_descriptor>* ordered_nodes)
// {
//      // g must have a weightmap an distance map already specified. Note that distance map will be overwritten
//      weight_map_t weightmap = get(edge_weight, *g);
//      distance_map_t distancemap = get(vertex_distance, *g);
     
//      // white: node has not been seen yet
//      // grey: node in Q
//      // black: node has been processed
//      std::vector<int> color_map(num_vertices(*g), WHITE);
     
//      distancemap[s] = 0;
//      color_map[s] = GREY;
//      Q->push(s);

//      int node_ctr = 0;

     
//      while (!Q->empty())
//      {
// 	 vertex_descriptor processing = Q->top();
// 	 Q->pop();
// 	 color_map[processing] = BLACK;
// 	 float processing_dist = distancemap[processing];

// 	 (*ordered_nodes)[node_ctr++] = processing;
	 

// 	 graph_traits < graph_t >::adjacency_iterator neighbor_i, neighbor_end;
// 	 for (tie(neighbor_i, neighbor_end) = adjacent_vertices(processing, *g); neighbor_i != neighbor_end; ++neighbor_i)
// 	 {
	     
// 	     vertex_descriptor v = *neighbor_i;
// 	     edge_descriptor e;
// 	     bool b;
// 	     tie(e, b) = edge(processing, v, *g);


// 	     if (color_map[v] == WHITE)
// 	     {
// 		 distancemap[v] =  processing_dist + weightmap[e];
// 		 color_map[v] = GREY;
// 		 Q->push(v);
// 	     }
// 	     else if (color_map[v] == GREY)
// 	     {
// 		 if (processing_dist + weightmap[e] < distancemap[v])
// 		 {
// 		     distancemap[v] = processing_dist + weightmap[e];
// 		     Q->update(v);
// 		 }
// 	     }
// 	 }
//      }
// }


// see DOC_STRING for expected argv parameters
int main(int argc, char** argv)
{
    if (argc < 4)
    {
	printf(DOC_STRING);
	exit(1);
    }
    
    int CLUSTER_SIZE(0);
    if (1 != sscanf(argv[1], "%d", &CLUSTER_SIZE))
	exit(1);


    // read file into new graph. File format is first line: 'NUM_NODES NUM_EDGES \n'.
    // For each edge there is line 'i j v \n'
    FILE *fp;
    
    if ((fp = fopen(argv[2], "rt")) == NULL)
	exit(1);
    
    int n, nnz;
    
    fscanf(fp, "%d %d\n", &n, &nnz);
    graph_t* g_ptr = new graph_t(n);
    graph_t& g = *g_ptr;
    weight_map_t weightmap = get(edge_weight, g);
    distance_map_t distancemap = get(vertex_distance, g);

    
    for (int x = 0; x < nnz; x++)
    {
	int i, j;
	float v;
	fscanf(fp, "%d %d %f\n", &i, &j, &v);
	
	edge_descriptor e;
	bool inserted;
	
	tie(e, inserted) = add_edge(i-1, j-1, g);
	
	weightmap[e] = v;
    }
    fclose(fp);

    
    clock_t clock_ticks = clock();
    double seconds = (double)clock_ticks / (double)CLOCKS_PER_SEC;
    printf("%f seconds elapsed: done reading file.\n", seconds);
    double time_started_alg = seconds;

    // set up Q
    dist_comp_class* mycomp = new dist_comp_class(distancemap);
    identity_property_map ident;
    MutableQueue* Q = new MutableQueue(num_vertices(g), *mycomp, ident);
    
    clock_ticks = clock();
    seconds = (double)clock_ticks / (double)CLOCKS_PER_SEC;
    printf("%f seconds elapsed: created Q.\n", seconds);
    

    // find random vertex; then find furthest node from that
    srand(time(NULL));
    int random = (int)floor((double)num_vertices(g) * (double)rand()/((double)(RAND_MAX) + 1.0));
    vertex_descriptor rand_v = random;
    
    //printf("randomly selected node %d \n", rand_v);


    // dijkstra it to find furthest node
    std::vector<vertex_descriptor> ordered_nodes(num_vertices(g));
    dijkstra(rand_v, g_ptr, Q, &ordered_nodes);
    
    graph_traits < graph_t >::vertex_iterator i, end;
    
    // I am under assumption we only deal with connected graphs!!
    
    vertex_descriptor start_node = ordered_nodes.back();
    
    clock_ticks = clock();
    seconds = (double)clock_ticks / (double)CLOCKS_PER_SEC;
    printf("%f seconds elapsed: ran dijkstra and found furthest node.\n", seconds);
    

    // Find out what order to cluster nodes in
    //printf("will start at node %d \n", start_node);
    dijkstra(start_node, g_ptr, Q, &ordered_nodes);

    clock_ticks = clock();
    seconds = (double)clock_ticks / (double)CLOCKS_PER_SEC;
    printf("%f seconds elapsed: ran dijkstra again.\n", seconds);
    

    std::vector<int> cluster_map(num_vertices(g), -1); // property map: graph_node -> cluster number
    std::vector<int> colors(num_vertices(g), WHITE);   // white = not seen ; black = already processed ; grey = in priority queue

    
    clock_ticks = clock();
    seconds = (double)clock_ticks / (double)CLOCKS_PER_SEC;
    printf("%f seconds elapsed: allocated cluster_map and colors vectors.\n", seconds);


    //////// TESTING QUEUE : Sort everything and test that way
//     for (tie(i, end) = vertices(g); i != end; ++i)
//     {
// 	Q->push(*i);
//     }
//     while (!Q->empty())
//     {
// 	Q->top();
// 	Q->pop();
//     }
//     clock_ticks = clock();
//     seconds = (double)clock_ticks / (double)CLOCKS_PER_SEC;
//     printf("%f seconds elapsed: sorted using queue.\n", seconds);
///////


    // Do actual clustering
    int max_nodes_consider = (int)((4.0/3.0)*(double)CLUSTER_SIZE);    
    int cluster_num = 0;
    
    int cur_node_idx = 0;     // num_vertices() returns unsigned int annoyingly
 
    std::vector<vertex_descriptor> cluster_consider(max_nodes_consider);  // bug, why are these so big 
    std::vector<float> cluster_score(max_nodes_consider);                // num_inside_edges/num_nodes in cluster

//    vertex_descriptor cluster_consider[max_nodes_consider];
//    float cluster_score[max_nodes_consider];
    int num_inside_edges = 0;
    int num_nodes_considering = 0;


    while (cur_node_idx < num_vertices(g))
    {
	num_inside_edges = 0;
	num_nodes_considering = 0;
	
	vertex_descriptor cur_node = ordered_nodes[cur_node_idx];

	colors[cur_node] = GREY;
	distancemap[cur_node] = 0;
	Q->push(cur_node);
   
	while (!Q->empty() && num_nodes_considering < max_nodes_consider)
	{
	   
	    vertex_descriptor processing = Q->top();
	    Q->pop();

	    num_nodes_considering++;
	    cluster_consider[num_nodes_considering - 1] = processing;
	    float processing_dist = distancemap[processing];
	    colors[processing] = BLACK;
	    //if (debug)
	    //	printf("processing vertex %d with final distance of %f\n", processing, processing_dist);
	    graph_traits < graph_t >::adjacency_iterator neighbor_i, neighbor_end;
	    for (tie(neighbor_i, neighbor_end) = adjacent_vertices(processing, g); neighbor_i != neighbor_end; ++neighbor_i)
	    {


		vertex_descriptor v = *neighbor_i;
		//if (colors[v] != black && debug)
		//    printf("   Examining neighbor %d\n", v);
		edge_descriptor e;
		bool b;
		tie(e, b) = edge(processing, v, g);


		if (cluster_map[v] >= 0)
		{
		    // ignore nodes that are already part of clusters
		}
		else if (colors[v] == BLACK)
		{
		    num_inside_edges++;
		}
		else if (colors[v] == WHITE)
		{
		    float new_path_dist = processing_dist + weightmap[e];
		    distancemap[v] = new_path_dist;
		    colors[v] = GREY;
		    Q->push(v);
		}
		else if (colors[v] == GREY)
		{
		    // we calculate new distance like resisters in parallel (i.e. the more paths to a node the closer it is, not just absolute dist)
		    float curDistance = distancemap[v];
		    float new_path_dist = processing_dist + weightmap[e];
		    float new_dist = 1.0/(1.0/new_path_dist + 1.0/curDistance);
		    if (new_dist > weightmap[e])
		    {
			distancemap[v] = new_dist;
			Q->update(v);
		    }
		    else if (weightmap[e] < curDistance)
		    {
			distancemap[v] = weightmap[e];
			Q->update(v);
		    }
		}
	    }

	    cluster_score[num_nodes_considering - 1] = (float)num_inside_edges/(float)num_nodes_considering;
	}
    
	float best_score = 0;
	int where_best = 0;

	// we select size of cluster to maximize the ratio between the number of edges inside the cluter and the number of nodes in cluster
	for (int x = (int)((2.0/3.0)*CLUSTER_SIZE); x < num_nodes_considering; x++)
	{
	    if (cluster_score[x] > best_score)
	    {
		best_score = cluster_score[x];
		where_best = x;
	    }
	}
	// in case we are forced to have a cluster smaller than the normal min
	if (where_best == 0)
	    where_best = num_nodes_considering - 1;

	// now assign the new cluster
	for (int x = 0; x <= where_best; x++)
	{
	    cluster_map[cluster_consider[x]] = cluster_num;
	}

	cluster_num++;
    
	// clean up
	for (int x = 0; x < num_nodes_considering; x++)
	{
	    colors[cluster_consider[x]] = WHITE;
	}

	while (!Q->empty())
	{
	    vertex_descriptor v = Q->top();
	    Q->pop();
	    colors[v] = WHITE;
	}
	
	// find the clost node to start_node that is not in a cluster
	while (cur_node_idx < num_vertices(g) && cluster_map[ordered_nodes[cur_node_idx]] >= 0)
	{
	    //printf("in curnode loop\n");
	    cur_node_idx++;
	}
    }

    clock_ticks = clock();
    seconds = (double)clock_ticks / (double)CLOCKS_PER_SEC;
    printf("%f seconds elapsed: done with alg.\n", seconds);
    double time_ended_alg = seconds;


    // output file in simple format: one line for each node (in order of numbering). Each line has string giving its cluster number.
    FILE* cluster_out;
    if ((cluster_out = fopen(argv[3], "wt")) == NULL)
	exit(1);
   
 
    for (tie(i, end) = vertices(g); i != end; ++i)
    {
	vertex_descriptor v = *i;
	fprintf(cluster_out, "%d\n", cluster_map[v]);
	//printf("vertex, %d, of %d is in clust %d\n", v , *(end-1), cluster_map[v]);
    }

    fclose(cluster_out);


    // Time to free memory!
    delete Q;
    delete mycomp;
    delete g_ptr;

    clock_t end_clustering = clock();
    seconds = (double)end_clustering / (double)CLOCKS_PER_SEC;
    printf("%f seconds elapsed: done outputting to file.\n", seconds);


    printf("***Whole algorithm took %f seconds.***\n", time_ended_alg - time_started_alg);

    return 0;
}
