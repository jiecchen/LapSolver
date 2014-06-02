#include "dijkstra.hpp"


using namespace boost;



// WARNING: ordered_nodes must already be big enough to hold all vertex descriptors in graph!
// custom version of dijkstra's shortest paths algorithm for boost graphs
// inputs:
//    - s: vertex descriptor that specifies starting point for alg
//    - g: pointer to boost graph (adjacency list) that alg is to work on
//    - Q: pointer to (empty) priority Q. Must order nodes by distance (must already have access to distance map) and support: push, update, pop, and top
//    - ordered_nodes: pointer to vector of nodes. Nodes will be placed in vector according to distance from s
// returns: None
void dijkstra(vertex_descriptor s, graph_t* g, MutableQueue* Q, std::vector<vertex_descriptor>* ordered_nodes)
{
     // g must have a weightmap an distance map already specified. Note that distance map will be overwritten
     weight_map_t weightmap = get(edge_weight, *g);
     distance_map_t distancemap = get(vertex_distance, *g);
     
     // white: node has not been seen yet
     // grey: node in Q
     // black: node has been processed
     std::vector<int> color_map(num_vertices(*g), WHITE);
     
     distancemap[s] = 0;
     color_map[s] = GREY;
     Q->push(s);

     int node_ctr = 0;

     
     while (!Q->empty())
     {
	 vertex_descriptor processing = Q->top();
	 Q->pop();
	 color_map[processing] = BLACK;
	 float processing_dist = distancemap[processing];

	 (*ordered_nodes)[node_ctr++] = processing;
	 

	 graph_traits < graph_t >::adjacency_iterator neighbor_i, neighbor_end;
	 for (tie(neighbor_i, neighbor_end) = adjacent_vertices(processing, *g); neighbor_i != neighbor_end; ++neighbor_i)
	 {
	     
	     vertex_descriptor v = *neighbor_i;
	     edge_descriptor e;
	     bool b;
	     tie(e, b) = edge(processing, v, *g);


	     if (color_map[v] == WHITE)
	     {
		 distancemap[v] =  processing_dist + weightmap[e];
		 color_map[v] = GREY;
		 Q->push(v);
	     }
	     else if (color_map[v] == GREY)
	     {
		 if (processing_dist + weightmap[e] < distancemap[v])
		 {
		     distancemap[v] = processing_dist + weightmap[e];
		     Q->update(v);
		 }
	     }
	 }
     }
}

