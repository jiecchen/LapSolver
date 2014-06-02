// graphcluster.cpp
// Aaron Barnet 2007-07-09
// Clustering alg based on the ACL PageRank algorithm

// See char* DOC_STRING for usage

#include "graph.hpp"
#include "dijkstra.hpp"
//extern "C" {
//#include "../oldc/st_defs.h"
//#include "../oldc/st_basic.h"
//}
//#include <boost/config.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//#include <queue>


// #include <boost/graph/graph_traits.hpp>
// #include <boost/graph/adjacency_list.hpp>
// #include <boost/pending/mutable_queue.hpp>

using namespace boost;

const char DOC_STRING[] = "use: graphclusterpr CLUSTER_SIZE ALPHA EPSILON_DIVIDER IJV_FILE OUT_FILE\n";


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


// see DOC_STRING for expected argv parameters
int main(int argc, char** argv)
{
    if (argc < 6)
    {
	printf(DOC_STRING);
	exit(1);
    }
    
    int CLUSTER_SIZE(0);
    float ALPHA(0.0);
    float EPSILON_DIVIDER(0.0);

    if (1 != sscanf(argv[1], "%d", &CLUSTER_SIZE))
	exit(1);
    if (1 != sscanf(argv[2], "%f", &ALPHA))
	exit(1);
    if (1 != sscanf(argv[3], "%f", &EPSILON_DIVIDER))
	exit(1);
    char* IN_FILE = argv[4];
    char* OUT_FILE = argv[5];
    
    // get graph from file

    IJVclass< graph_t > * ijv = new IJVclass< graph_t >(IN_FILE);

    clock_t clock_ticks = clock();
    double seconds = (double)clock_ticks / (double)CLOCKS_PER_SEC;
    printf("%f seconds elapsed: done reading file into memory.\n", seconds);

    
    graph_t* g_ptr = ijv->getBoostGraph();
    graph_t& g = *g_ptr;
    
    int num_verts = num_vertices(g);

    delete ijv;

    distance_map_t distancemap = get(vertex_distance, g);
    weight_map_t weightmap = get(edge_weight, g);
    
    clock_ticks = clock();
    seconds = (double)clock_ticks / (double)CLOCKS_PER_SEC;
    printf("%f seconds elapsed: finished creating graph.\n", seconds);
    double time_started_alg = seconds;

    // set up Q
    dist_comp_class mycomp(distancemap);
    //PropertyNodeComp pNodeComp;
    identity_property_map ident;
    MutableQueue Q(num_verts, mycomp, ident);
    
    clock_ticks = clock();
    seconds = (double)clock_ticks / (double)CLOCKS_PER_SEC;
    printf("%f seconds elapsed: created Q.\n", seconds);
    

    // find random vertex; then find furthest node from that
    srand(time(NULL));
    int random = (int)floor((double)num_verts * (double)rand()/((double)(RAND_MAX) + 1.0));
    vertex_descriptor rand_v = random;
    
    //printf("randomly selected node %d \n", rand_v);


    // dijkstra it to find furthest node
    std::vector<vertex_descriptor> ordered_nodes(num_verts);
    dijkstra(rand_v, g_ptr, &Q, &ordered_nodes);
    
    graph_traits < graph_t >::vertex_iterator i, end;
    
    // I am under assumption we only deal with connected graphs!!
    
    vertex_descriptor start_node = ordered_nodes.back();
    

    clock_ticks = clock();
    seconds = (double)clock_ticks / (double)CLOCKS_PER_SEC;
    printf("%f seconds elapsed: ran dijkstra and found furthest node.\n", seconds);
    

    // Find out what order to cluster nodes in
    //printf("will start at node %d \n", start_node);
    dijkstra(start_node, g_ptr, &Q, &ordered_nodes);

    clock_ticks = clock();
    seconds = (double)clock_ticks / (double)CLOCKS_PER_SEC;
    printf("%f seconds elapsed: ran dijkstra again.\n", seconds);
    

    // Do actual clustering
    int max_nodes_consider = (int)((4.0/3.0)*(double)CLUSTER_SIZE);    
    float EPSILON = 1.0 / (EPSILON_DIVIDER * (double)CLUSTER_SIZE);
    int cluster_num = 0;
    int cur_node_idx = 0;     // num_vertices() returns unsigned int annoyingly
 
    // r and s values as in paper
    std::vector<float> p(num_verts,0.0);  
    std::vector<float> r(num_verts,0.0); 
    std::vector<vertex_descriptor> processing_stack;
    //std::vector<vertex_descriptor> cluster_consider(num_vertices(g));
    std::vector<vertex_descriptor> touched;  // helps us clean everything up - really should be set
    std::vector<int> color_map(num_verts, WHITE);
    std::vector<float> cluster_score; //(num_vertices(g));
    //std::vector<int> cluster_map(num_vertices(g), -1);
    int cluster_map[num_verts];
    for (int x = 0; x < num_verts; x++)
    {
	cluster_map[x] = -1;
    }
    std::vector<vertex_descriptor> cluster_consider_array(max_nodes_consider);
    // set up priority queue for consideration
    p_comp_class* p_comp = new p_comp_class(&p);
    //identity_property_map ident;

    PropertyMaxQueue* cluster_consider = new PropertyMaxQueue(num_verts, *p_comp, ident);
    // std::priority_queue<vertex_descriptor, std::vector<vertex_descriptor>, p_comp_class>* cluster_consider = new std::priority_queue<vertex_descriptor, std::vector<vertex_descriptor>, p_comp_class>(*p_comp);
    //std::priority_queue<vertex_descriptor>* cluster_consider = new std::priority_queue<vertex_descriptor>();


    // calculate degree of all nodes and store it in the distance property map
    distance_map_t degreemap = get(vertex_distance, g);
    for (tie(i, end) = vertices(g); i != end; ++i)
    {
	vertex_descriptor v = *i;
	degreemap[v] = out_degree(v, g);
    }

// ///////// TESTING
//    for (int x = num_vertices(g) - 1; x >= 0; x--)
//     {
// 	p[x] = (double)rand()/((double)(RAND_MAX) + 1.0);
// 	cluster_consider->push(x);
//     }
//     while (!cluster_consider->empty())
//     {
// 	cluster_consider->top();
// 	cluster_consider->pop();
//     }
//     for (int x = 0; x < num_vertices(g); x++)
//     {
// 	p[x] = 0.0;
//     }
// // ///////// 
    
    //int deleteme = 0;
    float geo_ratio = ALPHA * (1.0 - ALPHA) / 2.0;
    while (cur_node_idx < num_verts)
    {
//	num_nodes_considering = 0;
//	num_inside_edges = 0;

	vertex_descriptor cur_node = ordered_nodes[cur_node_idx];
	r[cur_node] = 1.0;
	//cluster_consider->push(cur_node);
	processing_stack.push_back(cur_node);
	touched.push_back(cur_node);

	while (processing_stack.size() > 0)
	{
	   
	    vertex_descriptor processing = processing_stack.back();
	    processing_stack.pop_back();

	    //num_nodes_considering++;
	    //cluster_consider[num_nodes_considering - 1] = processing;

// 	    int unclust_degree = 0;  // degree not including already clustered nodes
 	    graph_traits < graph_t >::adjacency_iterator neighbor_i, neighbor_end;
// 	    for (tie(neighbor_i, neighbor_end) = adjacent_vertices(processing, g); neighbor_i != neighbor_end; ++neighbor_i)
// 	    {
// 		if (cluster_map[*neighbor_i] < 0)
// 		    unclust_degree++;
// 	    }

	    float unclust_degree = degreemap[processing];
	    if (unclust_degree < .5)
	    {
		p[processing] = r[processing];
		r[processing] = 0.0;
		continue;
	    }

	    float t = ceil((log(EPSILON * unclust_degree)- log(r[processing]))/log((1.0 - ALPHA)/2.0)); // num iterations to avoid excessive pops
	    float new_r = r[processing] * pow((1.0 - ALPHA)/2.0, t); //EPSILON * unclust_degree; //  // should be EPSILON * d
	    
	    float add_to_p = (ALPHA * r[processing] * ( 1.0 - pow(geo_ratio, t)))/(1.0 - geo_ratio);
	    
	    float add_to_neighboring_r = (r[processing] - new_r - add_to_p)/unclust_degree;
	    //printf("node=%d t=%f   cur_r=%f add_to_p=%f,   add_to_neighbor_r=%f degree=%f\n",processing, t ,r[processing], add_to_p, add_to_neighboring_r, unclust_degree);
	    //deleteme++;
	    //if (deleteme > 20)
		//exit(1);
	    for (tie(neighbor_i, neighbor_end) = adjacent_vertices(processing, g); neighbor_i != neighbor_end; ++neighbor_i)
	    {
		vertex_descriptor v = *neighbor_i;
		if (cluster_map[v] < 0)
		{
		    if (r[v] < SMIDGEN && p[v] < SMIDGEN)
			touched.push_back(v);//cluster_consider->push(v);
		    //r[v] = r[v] + ((1.0 - ALPHA) * r[processing]) / (2.0 * degreemap[processing]);
		    r[v] += add_to_neighboring_r;
		    // int v_unclust_degree = 0;
// 		    graph_traits < graph_t >::adjacency_iterator v_neighbor_i, v_neighbor_end;
// 		    for (tie(v_neighbor_i, v_neighbor_end) = adjacent_vertices(v, g); v_neighbor_i != v_neighbor_end; ++v_neighbor_i)
// 		    {
// 			if (cluster_map[*v_neighbor_i] < 0)
// 			{
// 			    v_unclust_degree++;
// 			}
// 		    }
		    if (r[v] >= EPSILON * degreemap[v])
			processing_stack.push_back(v);
		}
	    }

//	    if (p[processing] < SMIDGEN)
//	    {
	    //p[processing] = p[processing] + ALPHA * r[processing];
	    p[processing] += add_to_p;
	    
	    //cluster_consider->update(processing);
//	    }
//	    else
//	    {
//		p[processing] = p[processing] + ALPHA * r[processing];
//		cluster_consider->update(processing);
//	    }

	    //r[processing] = (1.0 - ALPHA) * r[processing] / 2.0;
	    r[processing] = new_r;
	    // float unclust_degree = degreemap[processing];
// 	    if (unclust_degree == 0)
// 		unclust_degree++;
// 	    if (r[processing] >= EPSILON * unclust_degree)
// 		processing_stack.push_back(processing);
	    //else
	    //printf("did not put on stack");
	    //printf("just finished processing node: %d with p=%9f, r=%9f \n", processing, p[processing], r[processing]);
	}
    
	//std::vector<vertex_descriptor> cluster_consider_array;
	for (unsigned int x = 0; x < touched.size(); x++)
	{
	    vertex_descriptor v = touched[x];
	    p[v] += r[v];
	    cluster_consider->push(v);
	}
	int num_considering = 0;
	while (!cluster_consider->empty())
	{
	    vertex_descriptor top_node = cluster_consider->top();
	    if (num_considering < max_nodes_consider)
		cluster_consider_array[num_considering++] = top_node;
	    cluster_consider->pop();
	    r[top_node] = 0.0;
	    p[top_node] = 0.0;
	}

	// need to score cluster
	int num_inside_edges = 0;
	for (int x = 0; x < num_considering; x++)
	{
	    vertex_descriptor u = cluster_consider_array[x];
	    color_map[u] = BLACK;
	    graph_traits < graph_t >::adjacency_iterator neighbor_i, neighbor_end;
	    for (tie(neighbor_i, neighbor_end) = adjacent_vertices(u, g); neighbor_i != neighbor_end; ++neighbor_i)
	    {
		vertex_descriptor v = *neighbor_i;
		if (color_map[v] == BLACK)
		    num_inside_edges++;
	    }
	    cluster_score.push_back((float)num_inside_edges / ((float)x + 1));
	}
	float best_score = 0;
	int where_best = 0;

	//printf("Have %d elements in array to consider, and we want to look at %f many\n", cluster_consider_array.size(), (4.0/3.0)*(double)CLUSTER_SIZE);

	for (int x = (int)((2.0/3.0)*CLUSTER_SIZE); x < num_considering; x++)
	{
	    if (cluster_score[x] > best_score)
	    {
		best_score = cluster_score[x];
		where_best = x;
	    }
	}
	// in case we are forced to have a cluster smaller than the normal min
	if (where_best == 0)
	    where_best = num_considering - 1;

	// now assign the new cluster
	for (int x = 0; x <= where_best; x++)
	{
	    vertex_descriptor u = cluster_consider_array[x];
	    cluster_map[u] = cluster_num;
	    
	    // decrease degree of neighboring nodes
	    graph_traits < graph_t >::adjacency_iterator neighbor_i, neighbor_end;
	    for (tie(neighbor_i, neighbor_end) = adjacent_vertices(u, g); neighbor_i != neighbor_end; ++neighbor_i)
	    {
		vertex_descriptor v = *neighbor_i;
		//if (color_map[v] == BLACK)
		degreemap[v]--;
	    }

	}

	cluster_num++;
    
	// clean up
 	for (int x = 0; x < num_considering; x++)
 	{
	    color_map[cluster_consider_array[x]] = WHITE;
	    //[cluster_consider_array[x]] = 0.0;
	    //[cluster_consider_array[x]] = 0.0;
 	}
 	touched.clear();
	cluster_score.clear();

	// find the clost node to start_node that is not in a cluster
	while (cur_node_idx < num_verts && cluster_map[ordered_nodes[cur_node_idx]] >= 0)
	{
	    //printf("in curnode loop\n");
	    cur_node_idx++;
	}
    }

    clock_ticks = clock();
    seconds = (double)clock_ticks / (double)CLOCKS_PER_SEC;
    printf("%f seconds elapsed: done with alg.\n", seconds);
    double time_ended_alg = seconds;


    binWriteArray(cluster_map, num_verts, OUT_FILE);

// output file in simple format: one line for each node (in order of numbering). Each line has string giving its cluster number.
//     FILE* cluster_out;
//     if ((cluster_out = fopen(argv[5], "wb")) == NULL)
// 	exit(1);
   
    
//     for (tie(i, end) = vertices(g); i != end; ++i)
//     {
// 	vertex_descriptor v = *i;
// 	fprintf(cluster_out, "%d\n", cluster_map[v]);
// 	//printf("vertex, %d, of %d is in clust %d\n", v , *(end-1), cluster_map[v]);
//     }
//     int n = num_vertices(g);
//     fwrite(&n, sizeof(int), 1, cluster_out);
//     fwrite(cluster_map, sizeof(int), n, cluster_out);
//     fclose(cluster_out);


    // Time to free memory!
    //delete Q;
    //delete mycomp;
    delete g_ptr;
    delete p_comp;
    delete cluster_consider;
    //delete [] cluster_map;
    clock_t end_clustering = clock();
    seconds = (double)end_clustering / (double)CLOCKS_PER_SEC;
    printf("%f seconds elapsed: done outputting to file.\n", seconds);


    printf("***Whole algorithm took %f seconds.***\n", time_ended_alg - time_started_alg);

    return 0;
}
