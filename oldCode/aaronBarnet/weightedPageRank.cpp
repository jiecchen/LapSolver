#include "weightedPageRank.hpp"
#include "profiler.hpp"
#include <math.h>

void WeightedPageRank::page_rank(vertex_descriptor start_node)
{
    r[start_node] = 1.0;
    processing_stack.push_back(start_node);
    touched.push_back(start_node);

    while (processing_stack.size() > 0)
    {
	   
	vertex_descriptor processing = processing_stack.back();
	processing_stack.pop_back();

	graph_traits < graph_t >::adjacency_iterator neighbor_i, neighbor_end;

	float unclust_degree = degreemap[processing];
	if (unclust_degree < SMIDGEN)
	{
	    p[processing] = r[processing];
	    r[processing] = 0.0;
	    continue;
	}

	float t = ceil((log(epsilon * unclust_degree)- log(r[processing]))/log((1.0 - alpha)/2.0)); // num iterations to avoid excessive pops
	float new_r = r[processing] * pow((1.0 - alpha)/2.0, t); //epsilon * unclust_degree; //  // should be epsilon * d
	float add_to_p = (alpha * r[processing] * ( 1.0 - pow(geo_ratio, t)))/(1.0 - geo_ratio);	    
	float add_to_neighboring_rs = (r[processing] - new_r - add_to_p);//* PUSH_DAMPENING;

	for (tie(neighbor_i, neighbor_end) = adjacent_vertices(processing, g); neighbor_i != neighbor_end; ++neighbor_i)
	{
	    vertex_descriptor v = *neighbor_i;
	    if (cluster_map[v] < 0)
	    {
		if (r[v] < SMIDGEN && p[v] < SMIDGEN)
		    touched.push_back(v);//cluster_consider->push(v);
		    
		edge_descriptor e;
		bool is_edge;
		tie(e, is_edge) = edge(processing, v, g);
		r[v] += add_to_neighboring_rs * weightmap[e] / weighted_degree[processing];

		if (r[v] >= epsilon * degreemap[v])
		    processing_stack.push_back(v);
	    }
	}

	p[processing] += add_to_p;// + ((1.0 - SELF_DAMPENING) * new_r);
	    
	r[processing] = new_r;// * SELF_DAMPENING;
    }
    
    //printf(" ... leaving pagerank\n");
}

// want to be scroing conductance: some of weights leaving cluster / (sum of weighted degrees of nodes in cluster)
int WeightedPageRank::score_clusters()
{
	float weight_inside = 0;
	float weight_on_edge = 0;
	while (!cluster_consider->empty())
	{
	    vertex_descriptor u = cluster_consider->top();
	    cluster_consider->pop();
	    if (cluster_consider_array.size() <= (unsigned int)max_nodes_consider)
	    {
		//cluster_consider_array[num_considering++] = top_node;
		bool connected = false;
		float new_weight_inside = 0.0;
		float new_weight_on_edge = 0.0;
		graph_traits < graph_t >::adjacency_iterator neighbor_i, neighbor_end;
		for (tie(neighbor_i, neighbor_end) = adjacent_vertices(u, g); neighbor_i != neighbor_end; ++neighbor_i)
		{
		    vertex_descriptor v = *neighbor_i;
		    
		    edge_descriptor e;
		    bool is_edge;
		    tie(e, is_edge) = edge(u, v, g);
		    
		    if (color_map[v] == BLACK)
		    {
			connected = true;
			new_weight_inside += weightmap[e];
			new_weight_on_edge -= weightmap[e];
		    }
		    else
		    {
			new_weight_inside += weightmap[e];
			new_weight_on_edge += weightmap[e];
		    }

		    if (color_map[v] == GREY)
		    neighboring_greys.push_back(v);
		}
		
		
		// NEGATIVE cluster score -> not connected to component
		// ignore nodes not connected to first starting component
		if (!connected && cluster_consider_array.size() > 0)
		{
		    //cluster_score.push_back(-1.0);
		    color_map[u] = GREY;    // NOt connected now but maybe in the future
		    neighboring_greys.clear();
		    continue;
		}
		
		// reconsider neighboring greys which would now be attached to main cluster
		while (neighboring_greys.size() > 0)
		{
		    vertex_descriptor grey_node = neighboring_greys.back();		    
		    cluster_consider->push(grey_node);
		    color_map[grey_node] = WHITE;
		    neighboring_greys.pop_back();
		}

		weight_inside += new_weight_inside;
		weight_on_edge += new_weight_on_edge;
		color_map[u] = BLACK;
		
		if (weight_inside < SMIDGEN)
		{
		    cluster_score.push_back(weight_on_edge / SMIDGEN);
		   
		}
		else
		{
		    cluster_score.push_back(weight_on_edge / weight_inside);
		}
		cluster_consider_array.push_back(u);
	    }
	}

	float best_score = BIG_FLOAT; // we want a low score (conductance)
	int where_best = -1;
	int valid_consider = 0;
	for (unsigned int x = min_nodes_consider; x < cluster_consider_array.size(); x++)
	{
	    if (cluster_score[x] < best_score)
	    {
		valid_consider++;
		best_score = cluster_score[x];
		where_best = x;
	    }
	}
	// in case we are forced to have a cluster smaller than the normal min
	if (where_best == -1)
	{	   
	    where_best = cluster_consider_array.size() - 1;
	}
	return where_best;
}


void WeightedPageRank::initialize_degrees()
{
    // calculate degree of all nodes and store it in the distance property map
    // do the same with weighted degreemap
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

	weighted_degree[v] = total_weight;
	degreemap[v] = out_degree(v, g);
    }
}

// clean up page rank
void WeightedPageRank::clean_up_pr()
{
    for (unsigned int x = 0; x < touched.size(); x++)
    {
	color_map[touched[x]] = WHITE;
	p[touched[x]] = 0.0;
	r[touched[x]] = 0.0;
    }
    touched.clear();
    cluster_score.clear();
    cluster_consider_array.clear();
}

void WeightedPageRank::clusterize()
{
    

    // Do actual clustering
    int cluster_num = 0;
    int cur_node_idx = 0;     // num_vertices() returns unsigned int annoyingly
 

    while (cur_node_idx < num_verts)
    {
	//printf("going to start with new node...");
	vertex_descriptor cur_node = ordered_nodes[cur_node_idx];
	
	page_rank(cur_node);
	
	for (unsigned int x = 0; x < touched.size(); x++)
	{
	    vertex_descriptor v = touched[x];
	    
	    //printf("node: %d in cluster %d\n",v, cluster_map[v]);
	    
	    p[v] += r[v];
	    cluster_consider->push(v);
	}

	
	int where_best = score_clusters();
	profiler.new_cluster(where_best + 1, cluster_score[where_best]);

	// if we have bad conductance retry clustering from node with highest weighted degree
	int num_reclusters = 0;
	//printf("new recluster while loop\n");
	//printf("cluster_score where best= %f\n", cluster_score[where_best]);
	
	while (cluster_score[where_best] > max_conductance && num_reclusters < max_reclusters)
	{
	    //printf("reclustering, num_reclusters=%d\n", num_reclusters);
	    float high_weight = -1.0;
	    vertex_descriptor high_weight_node = 0;
	    for (unsigned int x = 0; x < touched.size(); x++)
	    {
		if (weighted_degree[touched[x]] > high_weight)
		{
		    high_weight = weighted_degree[touched[x]];
		    high_weight_node = touched[x];
		}
	    }

	    // clean up
	    clean_up_pr();
	   
	    page_rank(high_weight_node);
	
	    for (unsigned int x = 0; x < touched.size(); x++)
	    {
		vertex_descriptor v = touched[x];
		p[v] += r[v];
		cluster_consider->push(v);
	    }
	    
	    where_best = score_clusters();
	    profiler.new_recluster(where_best + 1, cluster_score[where_best]);
	    //printf("where best=%d\n", where_best);
	    num_reclusters++;
	}

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

		edge_descriptor e;
		bool is_edge;
		tie(e, is_edge) = edge(u, v, g);
		weighted_degree[v] -= weightmap[e];
		degreemap[v]--;
	    }
	}

	cluster_num++;

	// keep track of singleton nodes
	if (where_best == 0)
	    singletons.push_back(cluster_consider_array[0]);

	// clean up
	clean_up_pr();

	// find the clost node to start_node that is not in a cluster
	while (cur_node_idx < num_verts && cluster_map[ordered_nodes[cur_node_idx]] >= 0)
	{
	    //printf("in curnode loop\n");
	    cur_node_idx++;
	}
    }

    if (AVOID_SINGLETONS)
      {
	printf("eliminating singletons...\n");
	eliminate_singletons(cluster_num);

      }
}
 
void WeightedPageRank::eliminate_singletons(int cluster_num)
{   
// //// temp code to test for disconnections  // DOING CONNECTEDNESS TEST WRONG!!!!! 
//     printf("beginning connectedness test...\n");
//     check_connectedness(cluster_num, false, true);
//     printf("done testing connectedness bofroe singleton elimination code\n");
// ///// End of temp code

    
//     // slow check on singletons
//     for(int x = 0; x < singletons.size(); x++)
//     {
// 	vertex_descriptor v = singletons[x];
// 	for (int xx = 0; xx < num_verts; xx++)
// 	{
// 	    if (cluster_map[v] == cluster_map[xx] && xx != v)
// 	    {
// 		printf("Non singleton listed as singleton\n");
// 		exit(43);
// 	    }
// 	}
//     }
//     printf("finished check to make sure all singletons are really singletons\n");

    // add singletons to new other clusters
    std::vector<int> eliminated_clusters;
    std::vector<bool> eliminated_clusters_hash(cluster_num, false);
    std::vector<int> neighbor_clusters;
    std::vector<vertex_descriptor> unsingled;

    std::vector<vertex_descriptor> copy_singletons; //(deleteme)
    for (unsigned int x = 0; x < singletons.size(); x++)
    {
	copy_singletons.push_back(singletons[x]);
    }

    while (singletons.size() > 0)
    {
	while (singletons.size() > 0)
	{
	    //printf("alpha\n");
	    vertex_descriptor u = singletons.back();
	    if (u == 245)
	    {
		printf("working with vert 60457\n");
	    }
	    unsingled.push_back(u);
	    singletons.pop_back();
	    graph_traits < graph_t >::adjacency_iterator neighbor_i, neighbor_end;
	    for (tie(neighbor_i, neighbor_end) = adjacent_vertices(u, g); neighbor_i != neighbor_end; ++neighbor_i)
	    {
		vertex_descriptor v = *neighbor_i;
		if (!eliminated_clusters_hash[cluster_map[v]])	// what about case where all neighbors have invalid options?    
		    neighbor_clusters.push_back(cluster_map[v]);
	    }
	    std::sort(neighbor_clusters.begin(), neighbor_clusters.end());

	    int num_most_edges = -1;
	    int clust_most_edges = -1;
	    int cur_cluster = -1;
	    int cur_num_edges = 0;
	    if (neighbor_clusters.size() == 0)
	    {
		printf("could not find cluster to attach singleton to.\n");
		exit(77);
	    }
	    for (unsigned int x = 0; x < neighbor_clusters.size(); x++)
	    {
		if (neighbor_clusters[x] == cur_cluster)
		    cur_num_edges++;
		else
		{
		    // new best cluster
		    if (cur_num_edges > num_most_edges)
		    {
			num_most_edges = cur_num_edges;
			clust_most_edges = cur_cluster;
		    }
		    cur_cluster = neighbor_clusters[x];
		    cur_num_edges = 1;
		}
	    }
	    if (cur_num_edges > num_most_edges)
	    {
		num_most_edges = cur_num_edges;
		clust_most_edges = cur_cluster;
	    }

	    if (!eliminated_clusters_hash[cluster_map[u]])
	    {
		eliminated_clusters.push_back(cluster_map[u]);
		eliminated_clusters_hash[cluster_map[u]] = true;
	    }
	    cluster_map[u] = clust_most_edges;
	    neighbor_clusters.clear();
	}
		
	while (unsingled.size() > 0)
	{
	    vertex_descriptor u = unsingled.back();
	    unsingled.pop_back();
	    if (eliminated_clusters_hash[cluster_map[u]])
		singletons.push_back(u);
	}
    }

    // remap clusters to fill gaps
    if (eliminated_clusters.size() > 0)
    {
	std::vector<int> cluster_function(cluster_num);
	std::sort(eliminated_clusters.begin(), eliminated_clusters.end());
	int cur_gap_cluster_idx = 0;
	int cur_gap_cluster = eliminated_clusters[cur_gap_cluster_idx];

	//cluster_function.push_back(cluster_num);  // avoid buffer overrun
	eliminated_clusters.push_back(cluster_num);
	for(int x = 0; x < cluster_num; x++)
	{
	    if (x >= cur_gap_cluster)
	    {
		cur_gap_cluster_idx++;
// 		if (cur_gap_cluster_idx >= eliminated_clusters.size())
// 		  {
// 		    printf("***buffer overflow!****\n");
// 		  }
		cur_gap_cluster = eliminated_clusters[cur_gap_cluster_idx];
	    }
	    cluster_function[x] = x - cur_gap_cluster_idx;
	}

// 	printf("making sure all eliminated clusters are really eliminated\n");
// 	for (unsigned int x = 0; x < eliminated_clusters.size(); x++)
// 	  {
// 	    for (int v = 0; v < num_verts; v++)
// 	      {
// 		if (cluster_map[v] == eliminated_clusters[x])
// 		  printf("cluster %d not eliminatad!\n", x);

// 	      }

//	  }


	for(int x = 0; x < num_verts; x++)
	{
	  cluster_map[x] = cluster_function[cluster_map[x]];
	}

    }
}


bool WeightedPageRank::check_connectedness(int num_clusters, bool note_singletons, bool note_empties)
{
    int new_num_clusts = num_clusters;// - 1 - cur_gap_cluster_idx;
    for (int x = 0; x < new_num_clusts; x++)
    {
	std::vector<vertex_descriptor> clust_verts;
	for (int v = 0; v < num_verts; v++)
	{
	    if (cluster_map[v] == x)
	    {
		clust_verts.push_back(v);
	    }
	}

	if (clust_verts.size() == 0)
	{
	    if (note_empties)
		printf("cluster %d is empty\n", x);
	    continue;
	}
	else if (clust_verts.size() == 1)
	{
	    if (note_singletons)
		printf("cluster %d is a singleton\n", x);
	    continue;
	    
	}
	std::vector<vertex_descriptor> connected_comp;
	connected_comp.push_back(clust_verts.back());
	clust_verts.pop_back();
	while (clust_verts.size() > 0)
	{
	    bool connected = false;
	    vertex_descriptor v = 0;
	    std::vector<vertex_descriptor>::iterator iter;
	    for (iter = clust_verts.begin(); iter != clust_verts.end(); iter++)
	    {
		v = *iter;//clust_verts[ii];//clust_verts.back();
//			clust_verts.pop_back();
		    
		for (unsigned int i = 0; i < connected_comp.size(); i++)
		{
		    vertex_descriptor u = connected_comp[i];
		    graph_traits < graph_t >::adjacency_iterator neighbor_i, neighbor_end;
		    for (tie(neighbor_i, neighbor_end) = adjacent_vertices(u, g); neighbor_i != neighbor_end; ++neighbor_i)
		    {
			if (*neighbor_i == v)
			{
			    connected = true;
//				clust_verts_idx_connected == ii;
			    connected_comp.push_back(v);

			    clust_verts.erase(iter);
			    break;
			}
		    }
		    if (connected)
			break;
		    
		}
		if (connected)
		    break;
	    }
	    if (!connected)
	    {
		printf("disconnection!\n");
		printf("vertex: %d in cluster %d\n", int(v), x);
		bool in_singletons = false;
		for (unsigned int i = 0; i < singletons.size(); i++)
		{
		    if (singletons[i] == v)
		    {
			printf("%d was the %dth node in singletons\n", int(v), i);
			in_singletons = true;
		    }
		}
		if (!in_singletons)
		    printf("%d was not in singletons\n", int(v));
		exit(12);
	    }
	    else
	    {
		connected_comp.push_back(v);
	    }
		
	}
	    
	    
    }

    return true;
}
