//#define PETSCKSP_DLL

#include <math.h>
#include <queue>
//#include "private/pcimpl.h"   /*I "petscpc.h" I*/
//#include "graph.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/pending/mutable_queue.hpp>

#include <cstdio>
using namespace boost;

// from graph.hpp because had to change def of graph
// boost typedef shortcuts
typedef subgraph< adjacency_list <listS, vecS, undirectedS, property <vertex_distance_t, float >,
			property <edge_weight_t, float, property<edge_index_t, bool> > > > graph_t;

typedef property_map<graph_t, edge_index_t>::type
    f_edges_t; 
typedef property_map<graph_t, vertex_distance_t>::type
    vert_dist_t; 

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


/*
  Type definitions
*/

// functor: used by priority queue to organize heap by an external property map map
class p_comp_min_class
{
public:
    p_comp_min_class(vert_dist_t* prop_map) : p_map(prop_map) {}
    
    bool operator()(const vertex_descriptor& x,const vertex_descriptor& y) const
	{
	    return (*p_map)[x] < (*p_map)[y];
	}
private:
     vert_dist_t* p_map;
};
 
// the priority queue type
typedef mutable_queue<vertex_descriptor, std::vector<vertex_descriptor>, p_comp_min_class, identity_property_map> PropertyMinQueue;

// typedef int PetscInt;
// typedef float PetscScalar;
// typedef int PetscErrorCode;

// struct vertex_weight_t {
//     typedef vertex_property_tag kind;
// };
// typedef subgraph<adjacency_list<vecS, vecS, undirectedS, 
//  				property<vertex_weight_t, float, property<vertex_index_t, int> >, 
//  				property<edge_weight_t, float, property<edge_index_t, int> > > > Graph;

// typedef graph_traits< Graph >::edge_descriptor edge_descriptor;

// typedef property_map<Graph, vertex_weight_t>::type VertexWeight;
// typedef property_map<Graph, edge_weight_t>::type EdgeWeight;
// typedef std::pair<int, int> Edge;

// struct PQNode {
//     int vertex;
//     int pred;
//     float dist;

//     PQNode() {}

//     PQNode(const int v,const float d) {
// 	vertex = v;
// 	dist = d;
//     }

//     PQNode(const int v,const int p,const float d) {
// 	vertex = v;
// 	pred = p;
// 	dist = d;
//     }

//     bool operator<( const PQNode &a ) const {
// 	return dist > a.dist;
//     }
// };
// typedef std::priority_queue<PQNode> ShortestPathPriorityQueue;

/*
  Function headers
*/
graph_t* LowStretchSpanningTree(graph_t& g);//Mat mat,Mat *pre);
void LowStretchSpanningTreeHelper(graph_t& g,const int root,const float alpha,
				 graph_t& h);//,int perm[]);
int StarDecomp(const graph_t g, const int root,const float delta,const float epsilon,
	       int& k,std::vector<int>& size,std::vector<std::vector<int> >& idx,
	       std::vector<int>& x,std::vector<int>& y);



const char DOC_STRING[] = "Run as lowstretch FILE_IN FILE_OUT\n";

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
	printf("%s", DOC_STRING);
	exit(1);
    }

    //printf("about to load graph from file\n");
    IJVclass<graph_t> ijv(argv[1]);
    graph_t* g_ptr = ijv.getBoostGraph();
    graph_t& g = *g_ptr;

//     graph_t g(4);
//     add_edge(0,1,1.0,g);
//     add_edge(0,2,1.0,g);
//     add_edge(0,3,1.0,g);
//     add_edge(1,2,1.0,g);
//     add_edge(1,3,1.0,g);
//     add_edge(2,3,1.0,g);
    //printf("about to do lowstretch\n");
    graph_t* h_ptr = LowStretchSpanningTree(g);

//    printf("h has %d eges and %d verts\n", num_edges(*h_ptr), num_vertices(*h_ptr));

    graph_t& h = *h_ptr;
    
//    binWritegraph_t<graph_t, graph_traits<graph_t>::edge_iterator >(h, argv[2]);

    int n_verts = num_vertices(h);
    int n_edges = num_edges(h);
    int* i = new int[n_edges];
    int* j = new int[n_edges];
    double* v = new double[n_edges];
    int edge_ctr = 0;
    weight_map_t weightmap = get(edge_weight, h);
    graph_traits < graph_t >::edge_iterator e, end;
    for (tie(e, end) = edges(*h_ptr); e != end; e++)
    {
	graph_traits < graph_t >::vertex_descriptor unode, vnode;
	unode = source(*e, h);
	vnode = target(*e, h);
	float weight = weightmap[*e];
	i[edge_ctr] = unode + 1;
	j[edge_ctr] = vnode + 1;
	v[edge_ctr] = weight;
	edge_ctr++;
    }

    printf("edge_ctr: %d, n_verts: %d, n_edges: %d\n", edge_ctr, n_verts, n_edges);

    FILE* fp;
    //printf("%s\n", argv[2]);
    if ((fp = fopen(argv[2], "wb")) == NULL)
    {
	printf("failed to open output file\n");
	exit(1);
    }    
    fwrite(&n_verts, sizeof(int), 1, fp);
    fwrite(&n_edges, sizeof(int), 1, fp);
    fwrite(i, sizeof(int), n_edges, fp);
    fwrite(j, sizeof(int), n_edges, fp);
    fwrite(v, sizeof(double), n_edges, fp);
    fclose(fp);


    delete [] i;
    delete [] j;
    delete [] v;
    delete g_ptr;
    delete h_ptr;

    double elapsed = (double)clock()/(double)CLOCKS_PER_SEC;
    printf("finished in %lf seconds\n", elapsed);

    return 0;

}


// need to replace PetscFunctionReturn / matrix crap

/* -------------------------------------------------------------------------- */
/*
  LowStretchSpanningTree - Applies EEST algorithm to construct 
  low-stretch spanning tree preconditioner
  unpetsc'd version just takes BGL graph and returns BGL spanning tree
  NO cholesky factorization!


  Input Parameters:
  .  mat - input matrix

  Output Parameter:
  .  pre - preconditioner matrix with cholesky factorization precomputed in place
*/
//#undef __FUNCT__  
//#define __FUNCT__ "LowStretchSpanningTree"

graph_t* LowStretchSpanningTree(graph_t& g)// Mat mat,Mat *prefact)
{
//   PetscErrorCode    ierr;
//   PetscInt          *idx;
//   PetscInt          n,ncols,i,k;
//   MatFactorInfo     info;
//   IS                perm;
//   const PetscInt    *cols;
//   const PetscScalar *vals;
//   Mat               pre;

//  PetscFunctionBegin;

//  ierr = MatGetSize(mat,PETSC_NULL,&n);CHKERRQ(ierr);

    //graph_t g(n),h(n);
    int n = num_vertices(g);
    graph_t* h_ptr = new graph_t(n);
    graph_t& h = *h_ptr;
    //   VertexWeight vertex_weight_h = get(vertex_weight_t(), h);
//     weight_map_t edge_weight_h = get(edge_weight_t(), h);
//     weight_map_t edge_weight_g = get(edge_weight_t(), g);

//   for (i=0; i<n; i++) {
//     PetscScalar sum = 0;
//     ierr = MatGetRow(mat,i,&ncols,&cols,&vals);CHKERRQ(ierr);
//     for (k=0; k<ncols; k++) {
//       if (cols[k] == i || vals[k] < 0) { /****** fix this ******/
// 	sum += vals[k];
// 	if (cols[k] > i) {
// 	  add_edge(i,cols[k],-vals[k],g);
// 	}
//       }
//     }
//     graph_traits < graph_t >::vertex_iterator i, end;
//     for (tie(i, end) = vertices(g); i != end; ++i)
//     {
// 	graph_traits < graph_t >::vertex_descriptor v = *i;
// 	//graph_traits < graph_t >::adjacency_iterator neighbor_i, neighbor_end;
// 	float sum(0);
// 	graph_traits<graph_t>::out_edge_iterator out_i, out_end;
// 	for (tie(out_i, out_end) = out_edges(v, g); out_i != out_end; out_i++)
// 	    //for (tie(neighbor_i, neighbor_end) = adjacent_vertices(v, g); neighbor_i != neighbor_end; neighbor_i++)
// 	{
// 	    sum += edge_weight_g[*out_i];
// 	}
// //	vertex_weight_h[v] = sum;
// 	//put(vertex_weight_h,*i,sum);
// 	//ierr = MatRestoreRow(mat,i,&ncols,&cols,&vals);//CHKERRQ(ierr);
//     }


    //ierr = PetscMalloc(n*sizeof(PetscInt),&idx);CHKERRQ(ierr);
//    int* idx = new int[n];
    LowStretchSpanningTreeHelper(g,0,log(4.0/3)/(2.0*log(n)),h); //,idx);//CHKERRQ(ierr);

//   ierr = ISCreateGeneral(PETSC_COMM_WORLD,n,idx,&perm);CHKERRQ(ierr);
//   ierr = ISSetPermutation(perm);CHKERRQ(ierr);
//   ierr = PetscFree(idx);CHKERRQ(ierr);

//   /*
//   printf ("\nperm:\n");
//   ISView(perm,PETSC_VIEWER_STDOUT_SELF);
//   */

//   ierr = MatCreate(PETSC_COMM_WORLD,&pre);CHKERRQ(ierr);
//   ierr = MatSetSizes(pre,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
//   ierr = MatSetType(pre,MATAIJ);CHKERRQ(ierr);

//   for (i=0; i<n; i++) {
//     ierr = MatSetValue(pre,i,i,get(vertex_weight_h,i),INSERT_VALUES);CHKERRQ(ierr);
//   }
  
//   graph_traits<graph_t>::edge_iterator e, e_end;
//   for (tie(e, e_end) = edges(h); e != e_end; e++) {
//     ierr = MatSetValue(pre,source(*e,h),target(*e,h),-get(edge_weight_h,*e),INSERT_VALUES);CHKERRQ(ierr);
//     ierr = MatSetValue(pre,target(*e,h),source(*e,h),-get(edge_weight_h,*e),INSERT_VALUES);CHKERRQ(ierr);
//   }
//   ierr = MatAssemblyBegin(pre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//   ierr = MatAssemblyEnd(pre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  

//   /*
//   printf("\n----------\nOriginal matrix:\n"); 
//   ierr = MatView(mat,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
//   printf("\n----------\nPreconditioner:\n\n"); 
//   ierr = MatView(pre,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
//   printf("\n----------\n"); 
//   */

//   ierr = MatFactorInfoInitialize(&info);CHKERRQ(ierr);
//   ierr = MatCholeskyFactorSymbolic(pre,perm,&info,prefact);CHKERRQ(ierr);
//   ierr = MatCholeskyFactorNumeric(pre,&info,prefact);CHKERRQ(ierr);
//   ierr = MatDestroy(pre);CHKERRQ(ierr);

    return h_ptr;
}

/* -------------------------------------------------------------------------- */
/*
  LowStretchSpanningTreeHelper

  Input Parameters:
  .  g - input graph
  .  alpha - parameter
  .  h - empty graph with same number of vertices as g, and vertex weights equal row sums
  .  perm - preallocated array in which to store vertex ordering

  Output Parameter:
  .  h - low-stretch spanning tree of input graph
  .  perm - vertex ordering   ??? why do I care about this
*/
//#undef __FUNCT__  
//#define __FUNCT__ "LowStretchSpanningTreeHelper"
void LowStretchSpanningTreeHelper(graph_t& g,const int root,const float alpha,
				 graph_t& h)//,int perm[])
{
    int n,i,j,k;
    std::vector<int> size,x,y;
    std::vector<std::vector<int> > idx;
    int ierr;

//  PetscFunctionBegin;
    
    weight_map_t edge_weight_g = get(edge_weight_t(),g);
//    VertexWeight vertex_weight_h = get(vertex_weight_t(),h);
    n = num_vertices(g);

    if (n > 2) {
	ierr = StarDecomp(g,root,1.0/3.0,alpha,k,size,idx,x,y);
	j = 0;
	for (i=1;i<=k;i++) {
	    graph_t& g1 = g.create_subgraph(idx[i].begin(),idx[i].end());
	    graph_t& h1 = h.create_subgraph(idx[i].begin(),idx[i].end());
	    LowStretchSpanningTreeHelper(g1,g1.global_to_local(g.local_to_global(x[i-1])),alpha,h1);//,perm+j);
	    j += size[i];
	}
	graph_t& g1 = g.create_subgraph(idx[0].begin(),idx[0].end());
	graph_t& h1 = h.create_subgraph(idx[0].begin(),idx[0].end());
	LowStretchSpanningTreeHelper(g1,g1.global_to_local(g.local_to_global(root)),alpha,h1);//,perm+j);
	for (i=0;i<k;i++) {
	    float w = get(edge_weight_g,edge(x[i],y[i],g).first);
	    add_edge(x[i],y[i],w,h);
//	    put(vertex_weight_h,x[i],get(vertex_weight_h,x[i])+w);
//	    put(vertex_weight_h,y[i],get(vertex_weight_h,y[i])+w);
	}
    } else if (n == 2) {
	graph_traits<graph_t>::edge_descriptor e = *(out_edges(root,g).first);
	int t = target(e,g);

	float w = get(edge_weight_g,e);
	add_edge(root,t,w,h);
//	put(vertex_weight_h,root,get(vertex_weight_h,root)+w);
//	put(vertex_weight_h,t,get(vertex_weight_h,t)+w);
//	perm[0] = g.local_to_global(t);
//	perm[1] = g.local_to_global(root);
    } else /* n == 1 */ {
//	perm[0] = g.local_to_global(root);
    }

  

    //return 0;
}



/* -------------------------------------------------------------------------- */
/*
  StarDecomp - calculate a star decomposition of the graph

  Input Parameters:
  .  g - input graph
  .  root - center of star decomposition
  .  delta, epsilon - parameters of star decomposition

  Output Parameter:
  .  k - number of partitions, not-including central one.
  .  size[i] - number of vertices in partition i
  .  idx[i] - list of vertices in partition i; partition 0 contains root
  .  x[i-1] - "anchor" vertex of non-central partition i
  .  y[i-i] - vertex in partition 0 forming "bridge" with x[i-1]
*/
//#undef __FUNCT__  
//#define __FUNCT__ "LowStretchSpanningTreeHelper"
int StarDecomp(graph_t g,const int root,const float delta,const float epsilon,
	       int& k,std::vector<int>& size,std::vector<std::vector<int> >& idx,
	       std::vector<int>& x,std::vector<int>& y)
{
    int n,m,edgesLeft;
    //PetscErrorCode ierr;
    //ShortestPathPriorityQueue pq;
//    float radius;

    std::vector<int> centerIdx;
    


    //PQNode node;

//  PetscFunctionBegin;

    //printf("going to into ball cut segment\n");

    graph_t& root_graph = g.is_root() ? g : g.root();
    weight_map_t edge_weight_g = get(edge_weight,root_graph);  // actually all property maps are same! so don't worry!
    f_edges_t f_edges_g = get(edge_index, g);
    vert_dist_t root_dist = get(vertex_distance, g);
    n = num_vertices(g);
    m = num_edges(g);
    edgesLeft = m;

//    std::vector<int> 
//    std::vector<int> pred(n,-1);
    //std::vector<int> succ[n]; 
    std::vector<int>::iterator i;
//    float dist[n];
//    std::vector<bool> taken(n,false);

    // running dijkstra
    float radius;
    std::vector<vertex_descriptor> ordered_nodes(n);
    int cntr = 0;
    //std::vector<float> root_dist(n);
    std::vector<int> pred(n, -1);
    std::vector<int> color_map(n, WHITE);
    //int num_root_edges = g.is_root() ? num_edges(g) : num_edges(g.root());
    //printf("num_root_edges: %d\n", num_root_edges);
    //std::vector<bool> fedges(num_root_edges);// turns out edge_index not local fedges(m, false);
    p_comp_min_class comp(&root_dist);
    identity_property_map ident;
    PropertyMinQueue pq(n, comp, ident);
    root_dist[root] = 0.0;
    pred[root] = root;
    pq.push(root);
    //printf("about to go to dijkstra while loop\n");
    while(!pq.empty())
    {
	//printf("top of while\n");
	vertex_descriptor u = pq.top();
	pq.pop();
	color_map[u] = BLACK;
	// put the edge that got us to u into fedges
	vertex_descriptor pred_u = pred[u];
	//radius = root_dist[u];
	ordered_nodes[cntr++] = u;
	if (pred_u != u)
	{
	    edge_descriptor e = edge(pred_u, u, g).first;
	    put(f_edges_g, e, true);
	}
	//printf("after fedges\n");
	adjacency_it i, i_end;
	for (tie(i, i_end) = adjacent_vertices(u, g); i != i_end; i++)
	{
	    vertex_descriptor v = *i;
	    if (color_map[v] == WHITE)
	    {
		//printf("before pred[v]\n");
		pred[v] = u;
		//printf("before edge_weight_g; pred_u: %d, u: %d\n", pred_u, u);
		root_dist[v] = root_dist[u] + 1.0/get(edge_weight_g, g.local_to_global(edge(u,v, g).first));
		//printf("before color_map\n");
		color_map[v] = GREY;
		pq.push(v);
	    }
	    else if (color_map[v] == GREY)
	    {
		float new_dist = root_dist[u] + 1.0/get(edge_weight_g, g.local_to_global(edge(u, v, g).first));
		if (new_dist < root_dist[v])
		{
		    root_dist[v] = new_dist;
		    pred[v] = u;
		    pq.update(v);
		}
	    }
	}
    }

    //printf("did dijkstra\n");

    radius = root_dist[ordered_nodes.back()];
    float min_radius = delta*radius;
    int center_size = 0;
    float boundary = 0.0;
    int edge_count = 0;
    
    // repurpose color_map to tell whether nodes are in ball or not
    color_map.clear();             // does clear cause mem alloc slowness
//    root_dist.clear();

    //color_map.resize(n, WHITE); 
    std::vector<int> partition_map(n, -1); // -1: not in a partition yet , 0 -> n: the partion number (0= center ball)
    int cur_part = 0;
    while (root_dist[ordered_nodes[center_size]] < min_radius
	|| boundary > (edge_count+1)*log(m)/(log(2.0)*(1.0-2.0*delta)*radius))
    {
	vertex_descriptor u = ordered_nodes[center_size++];
	partition_map[u] = cur_part;
	centerIdx.push_back(g.local_to_global(u));
	adjacency_it i, i_end;
	for (tie(i, i_end) = adjacent_vertices(u, g); i != i_end; i++)
	{
	    vertex_descriptor v = *i;
	    if (partition_map[v] != cur_part)
	    {
		edge_count++;
		boundary += get(edge_weight_g, g.local_to_global(edge(u,v,g).first));
	    }
	    else
		boundary -= get(edge_weight_g, g.local_to_global(edge(u,v,g).first));
	}
    }
    size.push_back(centerIdx.size());
    idx.push_back(centerIdx);
    
    // NOTE we found center region but have not done anything about returning it yet!

    // center region should now just be the first center_size nodes in ordered_nodes
    // find "shell" around ball that make up potential anchor
    std::queue<vertex_descriptor> anchor_q;
    for (int i = center_size; i < n; i++)
    {
	vertex_descriptor u = ordered_nodes[i];
	adjacency_it neighbor_i, neighbor_end;
	for (tie(neighbor_i, neighbor_end) = adjacent_vertices(u, g); neighbor_i != neighbor_end; neighbor_i++)
	{
	    vertex_descriptor v = *neighbor_i;
	    if (partition_map[v] == cur_part)
	    {
		anchor_q.push(u);
		break;
	    }
	}
    }

    //printf("done with ball cut segment\n");
    
    // now we make cones
    while (!anchor_q.empty())
    {
	vertex_descriptor anchor = anchor_q.front();
	anchor_q.pop();
	if (partition_map[anchor] >= 0)   // have we already included node in new cone
	    continue;

	cur_part++; // cones numbered 1,2,3,...

	std::vector<int> cone_idx;
	std::vector<int> cone_color_map(n, WHITE);
	//std::vector<float> cone_dist(n);
	vert_dist_t cone_dist = get(vertex_distance, g);
	//std::vector<vertex_descriptor> cone_ordered_nodes(n);
	int cone_edge_count = 0;
	int cone_internal_edge_count = 0;
	float cone_boundary = 0.0;
	
	p_comp_min_class cone_comp(&cone_dist);
	//identity_property_map cone_ident;
	PropertyMinQueue cone_pq(n, cone_comp, ident);

	
	cone_dist[anchor] = 0.0;
	cone_color_map[anchor] = GREY;
	cone_pq.push(anchor);
	
	//int node_cntr = 0;
	bool initialized_cone = false; // i.e. have we included all nodes within dist 0.0
	float cone_r = 0.0;
	float init_r_limit = 0.0;
	float boundary_limit = 0.0;
	while (!cone_pq.empty())
	{
	    vertex_descriptor u = cone_pq.top();
	    float new_r = cone_dist[u];
	    if (!initialized_cone && new_r > init_r_limit + SMIDGEN)
	    {
		initialized_cone = true;
		if (cone_edge_count == 0)
		    boundary_limit = log(m+1)*2.0/(log(2.0) * epsilon);
		else
		    boundary_limit = cone_edge_count * log(m / cone_internal_edge_count) / (log(2.0) * epsilon); // (edgeCount)*log(edgesLeft*1.0/initialInternalConeEdges)*2.0/(log(2.0)*epsilon*radius))
	    }

	    if (initialized_cone
		&& new_r > cone_r + SMIDGEN
		&& cone_boundary <= boundary_limit)
		break;
	    cone_r = new_r;

	    cone_pq.pop();
	    cone_color_map[u] = BLACK;
	    partition_map[u] = cur_part;
	    cone_idx.push_back(g.local_to_global(u));

	    adjacency_it i, i_end;
	    for (tie(i, i_end) = adjacent_vertices(u, g); i != i_end; i++)
	    {
		// calculate boundary/volume props
		vertex_descriptor v = *i;
		if (partition_map[v] != cur_part)
		{
		    cone_edge_count++;
		    cone_boundary += get(edge_weight_g, g.local_to_global(edge(u,v,g).first));
		}
		else
		{
		    cone_boundary -= get(edge_weight_g, g.local_to_global(edge(u,v,g).first));
		    cone_internal_edge_count++;
		}

		// consider more nodes
		if (partition_map[v] < 0)
		{
		    //int e_local_index = get(edge_index_g, g.global_to_local(edge(u,v,g).first));
		    //int e_index = get(edge_index_g, edge(u,v,g).first);
		    
		    //printf("e_index:%d  num_edges:%d  num_edges_root:%d \n", num_edges(g), num_edges(g.root()));
		    if (cone_color_map[v] == WHITE
			&& get(f_edges_g, edge(u,v,g).first))
		    {
			cone_color_map[v] = GREY;
			cone_dist[v] = cone_dist[u];
			cone_pq.push(v);
		    }
		    else if (cone_color_map[v] == WHITE
			     && !get(f_edges_g, edge(u,v,g).first))
		    {
			cone_color_map[v] = GREY;
			cone_dist[v] = cone_dist[u] + get(edge_weight_g, g.local_to_global(edge(u,v,g).first));
			cone_pq.push(v);
		    }
		    else if (cone_color_map[v] == GREY
			     && get(f_edges_g, edge(u,v,g).first))
		    {
			float new_dist = cone_dist[u];  // note this should always be at least as small as old dist
			if (new_dist < cone_dist[v])
			{
			    cone_dist[v] = new_dist;
			    cone_pq.update(v);
			}
		    }
		    else if (cone_color_map[v] == GREY
			     && !get(f_edges_g, edge(u,v,g).first))
		    {
			float new_dist = cone_dist[u] + get(edge_weight_g, g.local_to_global(edge(u,v,g).first));
			if (new_dist < cone_dist[v])
			{
			    cone_dist[v] = new_dist;
			    cone_pq.update(v);
			}
		    }
		}		    
	    }
	}
	x.push_back(anchor);
	y.push_back(pred[anchor]);
	size.push_back(cone_idx.size());
	idx.push_back(cone_idx);
    }
    k = cur_part;
// still have to pass back the new partionings etc.




//     /** form tree of shortest paths to root **/
//     graph_traits<graph_t>::out_edge_iterator e, e_end;  
//     for (tie(e,e_end)=out_edges(root,g); e!=e_end; e++) {
// 	int t = target(*e,g);
// 	pq.push(PQNode(t,root,1.0/get(edge_weight_g,*e)));
//     }
//     pred[root] = root;
//     while (!pq.empty()) {
// 	node = pq.top();pq.pop();
// 	if (pred[node.vertex] == -1) {
// 	    succ[node.pred].push_back(node.vertex);
// 	    pred[node.vertex] = node.pred;
// 	    dist[node.vertex] = node.dist;
// 	    for (tie(e,e_end)=out_edges(node.vertex,g); e!=e_end; e++) {
// 		int t = target(*e,g);
// 		if (pred[t] == -1) {
// 		    pq.push(PQNode(t,node.vertex,node.dist+1.0/get(edge_weight_g,*e)));
// 		}
// 	    }
// 	    radius = node.dist;
// 	}
//     }

//     /** BALL CUT **/
//     for (i=succ[root].begin();i!=succ[root].end();i++) {
// 	pq.push(PQNode(*i,dist[*i]));
//     }
//     float boundary = 0;
//     int edgeCount = 0; // all edges not just edges on boundary?
//     centerIdx.push_back(g.local_to_global(root));
//     taken[root] = true;//PETSC_TRUE;
//     centerSize = 1;
//     for (tie(e,e_end)=out_edges(root,g); e!=e_end; e++) {
// 	boundary += get(edge_weight_g,*e);
// 	edgeCount++;
//     }
//     const float minRadius = delta*radius;
//     while (dist[pq.top().vertex] < minRadius) {
// 	assert(!pq.empty());
// 	node = pq.top();pq.pop();
// 	centerIdx.push_back(g.local_to_global(node.vertex));
// 	taken[node.vertex] = true;//PETSC_TRUE;
// 	centerSize++;
// 	for (tie(e,e_end)=out_edges(node.vertex,g); e!=e_end; e++) {
// 	    if (taken[target(*e,g)]) {
// 		boundary -= get(edge_weight_g,*e);
// 	    } else {
// 		boundary += get(edge_weight_g,*e);
// 		edgeCount++;
// 	    }
// 	}
// 	for (i=succ[node.vertex].begin();i!=succ[node.vertex].end();i++) {
// 	    pq.push(PQNode(*i,dist[*i]));
// 	}
//     }
//     while (boundary > (edgeCount+1)*log(m)/(log(2)*(1-2*delta)*radius)) {
// 	assert(!pq.empty());
// 	node = pq.top();pq.pop();
// 	centerIdx.push_back(g.local_to_global(node.vertex));
// 	taken[node.vertex] = true;//PETSC_TRUE;
// 	centerSize++;
// 	for (tie(e,e_end)=out_edges(node.vertex,g); e!=e_end; e++) {
// 	    if (taken[target(*e,g)]) {
// 		boundary -= get(edge_weight_g,*e);
// 	    } else {
// 		boundary += get(edge_weight_g,*e);
// 		edgeCount++;
// 	    }
// 	}
// 	for (i=succ[node.vertex].begin();i!=succ[node.vertex].end();i++) {
// 	    pq.push(PQNode(*i,dist[*i]));
// 	}
//     }
//     size.push_back(centerSize);
//     idx.push_back(centerIdx);
//     edgesLeft -= edgeCount;

//     k = 0;
//     assert(!pq.empty());
//     std::queue<int> anchor_q;
//     ShortestPathPriorityQueue cone_pq;
//     std::vector<int> cone_succ[n]; 
//     std::vector<bool> cone_found(n,false);//PETSC_FALSE);

//     /** form tree of shortest paths to an anchor **/
//     while (!pq.empty()) {
// 	node = pq.top();pq.pop();
// 	cone_found[node.vertex] = true;//PETSC_TRUE;
// 	anchor_q.push(node.vertex);
// 	for (tie(e,e_end)=out_edges(node.vertex,g); e!=e_end; e++) {
// 	    int t = target(*e,g);
// 	    if (!taken[t]) {
// 		cone_pq.push(PQNode(t,node.vertex,1.0/get(edge_weight_g,*e)));
// 	    }
// 	}
//     }
//     while (!cone_pq.empty()) {
// 	node = cone_pq.top();cone_pq.pop();
// 	if (!cone_found[node.vertex]) {
// 	    cone_succ[node.pred].push_back(node.vertex);
// 	    cone_found[node.vertex] = true;//PETSC_TRUE;
// 	    for (tie(e,e_end)=out_edges(node.vertex,g); e!=e_end; e++) {
// 		int t = target(*e,g);
// 		if (!taken[t] && !cone_found[t]) {
// 		    cone_pq.push(PQNode(t,node.vertex,node.dist+1.0/get(edge_weight_g,*e)));
// 		}
// 	    }
// 	}
//     }

//     while (!anchor_q.empty()) {
// 	/** CONE CUT **/
// 	int anchor = anchor_q.front();anchor_q.pop();
// 	if (!taken[anchor]) {
// 	    int v;
// 	    int thisSize = 0;
// 	    std::vector<int> thisIdx;
// 	    std::queue<int> q;
// 	    ShortestPathPriorityQueue mycone_pq;
// 	    std::vector<bool> mycone_taken(n,false);//);PETSC_FALSE);
// 	    int initialInternalConeEdges = 0;

// 	    boundary = 0;
// 	    edgeCount = 0;
// 	    q.push(anchor);
// 	    while (!q.empty()) {
// 		v = q.front();q.pop();
// 		taken[v] = true;//PETSC_TRUE;
// 		mycone_taken[v] = true;//PETSC_TRUE;
// 		thisIdx.push_back(g.local_to_global(v));
// 		thisSize++;
// 		for (i=cone_succ[v].begin();i!=cone_succ[v].end();i++) {
// 		    q.push(*i);
// 		}
// 		for (tie(e,e_end)=out_edges(v,g); e!=e_end; e++) {
// 		    int t = target(*e,g);
// 		    if (!taken[t]) {
// 			mycone_pq.push(PQNode(t,v,1.0/get(edge_weight_g,*e)));
// 			boundary += get(edge_weight_g,*e);
// 			edgeCount++;
// 		    } else if (mycone_taken[t]) {
// 			boundary -= get(edge_weight_g,*e);
// 			initialInternalConeEdges++;
// 		    }
// 		}
// 	    }
// 	    if (initialInternalConeEdges < edgesLeft) {
// 		while (initialInternalConeEdges == 0 ?
// 		       boundary > (edgeCount+1)*log(edgesLeft+1)*2.0/(log(2.0)*epsilon*radius) : 
// 		       boundary > (edgeCount)*log(edgesLeft*1.0/initialInternalConeEdges)*2.0/(log(2.0)*epsilon*radius))
// 		{
// 		    assert(!mycone_pq.empty());
// 		    node = mycone_pq.top();mycone_pq.pop();
// 		    if (!mycone_taken[node.vertex]) {
// 			q.push(node.vertex);
// 			while (!q.empty()) {
// 			    v = q.front();q.pop();
// 			    taken[v] = true;//PETSC_TRUE;
// 			    mycone_taken[v] = true;//PETSC_TRUE;
// 			    thisIdx.push_back(g.local_to_global(v));
// 			    thisSize++;
// 			    for (i=cone_succ[v].begin();i!=cone_succ[v].end();i++) {
// 				q.push(*i);
// 			    }
// 			    for (tie(e,e_end)=out_edges(v,g); e!=e_end; e++) {
// 				int t = target(*e,g);
// 				if (!taken[t]) {
// 				    mycone_pq.push(PQNode(t,v,node.dist+1.0/get(edge_weight_g,*e)));
// 				    boundary += get(edge_weight_g,*e);
// 				    edgeCount++;
// 				} else if (mycone_taken[t]) {
// 				    boundary -= get(edge_weight_g,*e);
// 				}
// 			    }
// 			}
// 		    }
// 		}
// 	    }
// 	    edgesLeft -= edgeCount;
// 	    size.push_back(thisSize);
// 	    idx.push_back(thisIdx);
// 	    x.push_back(anchor);
// 	    y.push_back(pred[anchor]);
// 	    k++;
// 	}
//     }
    
  

//     /*
//     // pseudo cone cut
//     while (!pq.empty()) {
//     node = pq.top();pq.pop();

//     PetscInt thisSize = 1;
//     std::vector<PetscInt> thisIdx;
//     std::queue<PetscInt> q;

//     thisIdx.push_back(g.local_to_global(node.vertex));
//     for (i=succ[node.vertex].begin();i!=succ[node.vertex].end();i++) {
//     q.push(*i);
//     }

//     PetscInt v;
//     while (!q.empty()) {
//     v = q.front();q.pop();
//     thisSize++;
//     thisIdx.push_back(g.local_to_global(v));
//     for (i=succ[v].begin();i!=succ[v].end();i++) {
//     q.push(*i);
//     }
//     }
//     size.push_back(thisSize);
//     idx.push_back(thisIdx);
//     x.push_back(node.vertex);
//     y.push_back(pred[node.vertex]);
//     k++;
//     }
//     */

  


    return 0;
}
