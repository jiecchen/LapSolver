#ifndef DIJKSTRA_HPP
#define DIJKSTRA_HPP

// #include <boost/config.hpp>
// #include <boost/graph/graph_traits.hpp>
// #include <boost/graph/adjacency_list.hpp>
// #include <boost/pending/mutable_queue.hpp>

#include "graph.hpp"

using namespace boost;

// functor: used by priority queue to organize heap by distance map
class dist_comp_class
{
public:
    dist_comp_class(distance_map_t distance_map) : d_map(distance_map) {}
    
    bool operator()(const vertex_descriptor& x,const vertex_descriptor& y) const
	{
	    return d_map[x] < d_map[y];
	}
private:
    distance_map_t d_map;
};

typedef mutable_queue<vertex_descriptor, std::vector<vertex_descriptor>, dist_comp_class , identity_property_map> MutableQueue;

void dijkstra(vertex_descriptor s, graph_t* g, MutableQueue* Q, std::vector<vertex_descriptor>* ordered_nodes);

#endif // DIJKSTRA_HPP
