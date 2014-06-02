#include "graph.hpp"

const bool AVOID_SINGLETONS = false;

// functor: used by priority queue to organize heap by distance map
class p_comp_class
{
public:
    p_comp_class(std::vector<float>* prop_map) : p_map(prop_map) {}
    
    bool operator()(const vertex_descriptor& x,const vertex_descriptor& y) const
	{
	    return (*p_map)[x] > (*p_map)[y];
	}
private:
     std::vector<float>* p_map;
};
 
// the priority queue type
typedef mutable_queue<vertex_descriptor, std::vector<vertex_descriptor>, p_comp_class , identity_property_map> PropertyMaxQueue;

class Profiler;

class WeightedPageRank
{
public:
    WeightedPageRank(graph_t &graph, Profiler& pfiler, std::vector<vertex_descriptor>& on, int c, float a, float e, float mc, int mr) :
	g(graph), profiler(pfiler), ordered_nodes(on), cluster_size(c), alpha(a), epsilon_divider(e), max_conductance(mc), max_reclusters(mr)
	{
	    // allocate our memory
	    num_verts = num_vertices(g);
	    p.resize(num_verts, 0.0); //= new std::vector<float>(num_verts, 0.0);
	    r.resize(num_verts, 0.0); // = new std::vector<float>(num_verts, 0.0);
	    color_map.resize(num_verts, WHITE);// = new std::vector<int>(num_verts, WHITE);
	    weighted_degree.resize(num_verts, 0);// = new std::vector<float>(num_verts, 0);
	    cluster_map = new int[num_verts];
	    for (int x = 0; x < num_verts; x++)
	    {
		cluster_map[x] = -1;
	    }
	    
	    p_comp = new p_comp_class(&p);
	    cluster_consider = new PropertyMaxQueue(num_verts, *p_comp, ident);
	    
	    degreemap = get(vertex_distance, g);
	    weightmap = get(edge_weight, g);

	    // calculate rest of params
	    epsilon = 1.0 / (epsilon_divider * (float)cluster_size);
	    max_nodes_consider = (int)((4.0f/3.0f)*(float)cluster_size);
	    min_nodes_consider = (int)((2.0f/3.0f)*(float)cluster_size);
	    geo_ratio = alpha * (1.0 - alpha) / 2.0;
	    
	    initialize_degrees();

	    // we automatically generate clusters on instantiation
	    clusterize();
	}
    ~WeightedPageRank()
	{
	  //   delete p;
// 	    delete r;
// 	    delete color_map;
// 	    delete weighted_degree;
	    delete [] cluster_map;
	    delete p_comp;
	    delete cluster_consider;
	}
    int* get_clusters() const
	{
	    return cluster_map;
	}

private:


    // graph
    graph_t& g;
    int num_verts;
    // profiler
    Profiler& profiler;
    // ordered nodes
    std::vector<vertex_descriptor>& ordered_nodes;

    // private methods
    void initialize_degrees();
    void clusterize();
    void page_rank(vertex_descriptor start_node);
    int score_clusters();
    void clean_up_pr();
  void eliminate_singletons(int cluster_num);
    bool check_connectedness(int num_clusters, bool note_singletons= false, bool note_empties=false);  // slow test function / not to be used for real

    // parameters
    int cluster_size;
    float alpha;
    float epsilon_divider;
    float epsilon;
    float max_conductance;
    int max_reclusters;
    int max_nodes_consider;
    int min_nodes_consider;
    float geo_ratio;

    // data structs
    std::vector<float> p;//(num_verts,0.0);  
    std::vector<float> r;//(num_verts,0.0); 
    std::vector<vertex_descriptor> processing_stack;
    std::vector<vertex_descriptor> touched;  // helps us clean everything up - really should be set
    std::vector<int> color_map;//(num_verts, WHITE);
    std::vector<float> cluster_score; //(num_vertices(g));
    std::vector<float> weighted_degree;
    std::vector<vertex_descriptor> singletons;
    distance_map_t degreemap;
    weight_map_t weightmap;
    int* cluster_map;//[num_verts];
    std::vector<vertex_descriptor> cluster_consider_array;
    std::vector<vertex_descriptor> neighboring_greys;

    identity_property_map ident;
    p_comp_class* p_comp;// = new p_comp_class(&p);
    PropertyMaxQueue* cluster_consider;// = new PropertyMaxQueue(num_verts, *p_comp, ident);
};
