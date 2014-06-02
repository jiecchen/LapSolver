//#define PETSCKSP_DLL

#include <math.h>
#include <queue>
//#include "private/pcimpl.h"   /*I "petscpc.h" I*/
#include "boost/graph/adjacency_list.hpp"
#include "graph.hpp"
#include "boost/graph/subgraph.hpp"
#include <time.h>
using namespace boost;


/*
  Type definitions
*/

// typedef int PetscInt;
// typedef float PetscScalar;
// typedef int PetscErrorCode;

struct vertex_weight_t {
    typedef vertex_property_tag kind;
};
typedef subgraph<adjacency_list<vecS, vecS, undirectedS, 
 				property<vertex_weight_t, float, property<vertex_index_t, int> >, 
 				property<edge_weight_t, float, property<edge_index_t, int> > > > Graph;

typedef graph_traits< Graph >::edge_descriptor edge_descriptor;

typedef property_map<Graph, vertex_weight_t>::type VertexWeight;
typedef property_map<Graph, edge_weight_t>::type EdgeWeight;
typedef std::pair<int, int> Edge;

struct PQNode {
    int vertex;
    int pred;
    float dist;

    PQNode() {}

    PQNode(const int v,const float d) {
	vertex = v;
	dist = d;
    }

    PQNode(const int v,const int p,const float d) {
	vertex = v;
	pred = p;
	dist = d;
    }

    bool operator<( const PQNode &a ) const {
	return dist > a.dist;
    }
};
typedef std::priority_queue<PQNode> ShortestPathPriorityQueue;

/*
  Function headers
*/
Graph* LowStretchSpanningTree(Graph& g);//Mat mat,Mat *pre);
int LowStretchSpanningTreeHelper(Graph& g,const int root,const float alpha,
				 Graph& h);//,int perm[]);
int StarDecomp(const Graph g,const int root,const float delta,const float epsilon,
	       int& k,std::vector<int>& size,std::vector<std::vector<int> >& idx,
	       std::vector<int>& x,std::vector<int>& y);



const char DOC_STRING[] = "Run as lowstretch FILE_IN FILE_OUT\n";
void binWriteGraph2IJV(Graph& h, const char* fileName)
{
    int n_verts = num_vertices(h);
    int n_edges = num_edges(h);
    int* i = new int[n_edges];
    int* j = new int[n_edges];
    double* v = new double[n_edges];
    int edge_ctr = 0;
    EdgeWeight weightmap = get(edge_weight, h);
    graph_traits < Graph >::edge_iterator e, end;
    for (tie(e, end) = edges(h); e != end; e++)
    {
	graph_traits < Graph >::vertex_descriptor unode, vnode;
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
    if ((fp = fopen(fileName, "wt")) == NULL)
	exit(1);

    fwrite(&n_verts, sizeof(int), 1, fp);
    fwrite(&n_edges, sizeof(int), 1, fp);
    fwrite(i, sizeof(int), n_edges, fp);
    fwrite(j, sizeof(int), n_edges, fp);
    fwrite(v, sizeof(double), n_edges, fp);
    fclose(fp);

    delete [] i;
    delete [] j;
    delete [] v;
}

// we assume h is a tree!
void binWriteGraph2pArray(Graph& h, const char* fileName)
{
    int n_verts = num_vertices(h);

    if (n_verts < 1)
    {
	fprintf(stderr, "h as no vertices to output!\n");
	exit(99);
    }
    // c style queue for bfs of tree
    int* nodes_Q = new int[n_verts];
    int nodes_front = 0;
    int nodes_back = 0;

    bool* touched = new bool[n_verts];// whether node has been processed
    for (int i = 0; i < n_verts; i++)
    {
	touched[i] = false;
    }
    int* parray = new int[n_verts];       // what we will outpub

    nodes_Q[nodes_back++] = 0;
    parray[0] = 0;
    touched[0] = true;
    while (nodes_front < n_verts)
    {
	int cur_vert = nodes_Q[nodes_front++];
	graph_traits<Graph>::adjacency_iterator i, i_end;
	for (tie(i, i_end) = adjacent_vertices(cur_vert, h); i != i_end; i++)
	{
	    vertex_descriptor v = *i;
	    if (touched[v])
		continue;
	    
	    nodes_Q[nodes_back++] = v;
	    parray[v] = cur_vert;
	    touched[v] = true;
	}
    }

    FILE* fp;
    if ((fp = fopen(fileName, "wt")) == NULL)
	exit(1);
    
    
    fwrite(&n_verts, sizeof(int), 1, fp);
    fwrite(parray, sizeof(int), n_verts, fp);
    fclose(fp);


    delete [] parray;
    delete [] touched;
    delete [] nodes_Q;

}

int main(int argc, char* argv[])
{

    if (argc < 3)
    {
	printf("%s", DOC_STRING);
	exit(1);
    }

    IJVclass<Graph> ijv(argv[1]);
    Graph* g_ptr = ijv.getBoostGraph();
    Graph& g = *g_ptr;

    Graph* h_ptr = LowStretchSpanningTree(g);

//    printf("h has %d eges and %d verts\n", num_edges(*h_ptr), num_vertices(*h_ptr));

    Graph& h = *h_ptr;
    

    // for binIJV out
    // binWriteGraph2IJV(h, argv[2]);

    // for pArray
    binWriteGraph2pArray(h, argv[2]);
    

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

Graph* LowStretchSpanningTree(Graph& g)// Mat mat,Mat *prefact)
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

    //Graph g(n),h(n);
    int n = num_vertices(g);
    Graph* h_ptr = new Graph(n);
    Graph& h = *h_ptr;
    //   VertexWeight vertex_weight_h = get(vertex_weight_t(), h);
    EdgeWeight edge_weight_h = get(edge_weight_t(), h);
    EdgeWeight edge_weight_g = get(edge_weight_t(), g);

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
    graph_traits < Graph >::vertex_iterator i, end;
    for (tie(i, end) = vertices(g); i != end; ++i)
    {
	graph_traits < Graph >::vertex_descriptor v = *i;
	//graph_traits < Graph >::adjacency_iterator neighbor_i, neighbor_end;
	float sum(0);
	graph_traits<Graph>::out_edge_iterator out_i, out_end;
	for (tie(out_i, out_end) = out_edges(v, g); out_i != out_end; out_i++)
	    //for (tie(neighbor_i, neighbor_end) = adjacent_vertices(v, g); neighbor_i != neighbor_end; neighbor_i++)
	{
	    sum += edge_weight_g[*out_i];
	}
//	vertex_weight_h[v] = sum;
	//put(vertex_weight_h,*i,sum);
	//ierr = MatRestoreRow(mat,i,&ncols,&cols,&vals);//CHKERRQ(ierr);
    }


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
  
//   graph_traits<Graph>::edge_iterator e, e_end;
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
int LowStretchSpanningTreeHelper(Graph& g,const int root,const float alpha,
				 Graph& h)//,int perm[])
{
    int n,i,j,k;
    std::vector<int> size,x,y;
    std::vector<std::vector<int> > idx;
    int ierr;

//  PetscFunctionBegin;

    EdgeWeight edge_weight_g = get(edge_weight_t(),g);
//    VertexWeight vertex_weight_h = get(vertex_weight_t(),h);
    n = num_vertices(g);

    if (n > 2) {
	ierr = StarDecomp(g,root,1.0/3,alpha,k,size,idx,x,y);
	j = 0;
	for (i=1;i<=k;i++) {
	    Graph& g1 = g.create_subgraph(idx[i].begin(),idx[i].end());
	    Graph& h1 = h.create_subgraph(idx[i].begin(),idx[i].end());
	    LowStretchSpanningTreeHelper(g1,g1.global_to_local(g.local_to_global(x[i-1])),alpha,h1);//,perm+j);
	    j += size[i];
	}
	Graph& g1 = g.create_subgraph(idx[0].begin(),idx[0].end());
	Graph& h1 = h.create_subgraph(idx[0].begin(),idx[0].end());
	LowStretchSpanningTreeHelper(g1,g1.global_to_local(g.local_to_global(root)),alpha,h1);//,perm+j);
	for (i=0;i<k;i++) {
	    float w = get(edge_weight_g,edge(x[i],y[i],g).first);
	    add_edge(x[i],y[i],w,h);
//	    put(vertex_weight_h,x[i],get(vertex_weight_h,x[i])+w);
//	    put(vertex_weight_h,y[i],get(vertex_weight_h,y[i])+w);
	}
    } else if (n == 2) {
	graph_traits<Graph>::edge_descriptor e = *(out_edges(root,g).first);
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

  

    return 0;
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
int StarDecomp(Graph g,const int root,const float delta,const float epsilon,
	       int& k,std::vector<int>& size,std::vector<std::vector<int> >& idx,
	       std::vector<int>& x,std::vector<int>& y)
{
    int n,m,edgesLeft;
    //PetscErrorCode ierr;
    ShortestPathPriorityQueue pq;
    float radius = 0.0f;
    int centerSize = 0;
    std::vector<int> centerIdx;
    PQNode node;

//  PetscFunctionBegin;

    EdgeWeight edge_weight_g = get(edge_weight_t(),g);
    n = num_vertices(g);
    m = num_edges(g);
    edgesLeft = m;

    std::vector<int> pred(n,-1);
    std::vector<int> succ[n]; 
    std::vector<int>::iterator i;
    float dist[n];
    std::vector<bool> taken(n,false);

    /** form tree of shortest paths to root **/
    graph_traits<Graph>::out_edge_iterator e, e_end;  
    for (tie(e,e_end)=out_edges(root,g); e!=e_end; e++) {
	int t = target(*e,g);
	pq.push(PQNode(t,root,1.0/get(edge_weight_g,*e)));
    }
    pred[root] = root;
    while (!pq.empty()) {
	node = pq.top();pq.pop();
	if (pred[node.vertex] == -1) {
	    succ[node.pred].push_back(node.vertex);
	    pred[node.vertex] = node.pred;
	    dist[node.vertex] = node.dist;
	    for (tie(e,e_end)=out_edges(node.vertex,g); e!=e_end; e++) {
		int t = target(*e,g);
		if (pred[t] == -1) {
		    pq.push(PQNode(t,node.vertex,node.dist+1.0/get(edge_weight_g,*e)));
		}
	    }
	    radius = node.dist;
	}
    }

    /** BALL CUT **/
    for (i=succ[root].begin();i!=succ[root].end();i++) {
	pq.push(PQNode(*i,dist[*i]));
    }
    float boundary = 0;
    int edgeCount = 0;
    centerIdx.push_back(g.local_to_global(root));
    taken[root] = true;//PETSC_TRUE;
    centerSize = 1;
    for (tie(e,e_end)=out_edges(root,g); e!=e_end; e++) {
	boundary += get(edge_weight_g,*e);
	edgeCount++;
    }
    const float minRadius = delta*radius;
    while (dist[pq.top().vertex] < minRadius) {
	assert(!pq.empty());
	node = pq.top();pq.pop();
	centerIdx.push_back(g.local_to_global(node.vertex));
	taken[node.vertex] = true;//PETSC_TRUE;
	centerSize++;
	for (tie(e,e_end)=out_edges(node.vertex,g); e!=e_end; e++) {
	    if (taken[target(*e,g)]) {
		boundary -= get(edge_weight_g,*e);
	    } else {
		boundary += get(edge_weight_g,*e);
		edgeCount++;
	    }
	}
	for (i=succ[node.vertex].begin();i!=succ[node.vertex].end();i++) {
	    pq.push(PQNode(*i,dist[*i]));
	}
    }
    while (boundary > (edgeCount+1)*log(m)/(log(2)*(1-2*delta)*radius)) {
	assert(!pq.empty());
	node = pq.top();pq.pop();
	centerIdx.push_back(g.local_to_global(node.vertex));
	taken[node.vertex] = true;//PETSC_TRUE;
	centerSize++;
	for (tie(e,e_end)=out_edges(node.vertex,g); e!=e_end; e++) {
	    if (taken[target(*e,g)]) {
		boundary -= get(edge_weight_g,*e);
	    } else {
		boundary += get(edge_weight_g,*e);
		edgeCount++;
	    }
	}
	for (i=succ[node.vertex].begin();i!=succ[node.vertex].end();i++) {
	    pq.push(PQNode(*i,dist[*i]));
	}
    }
    size.push_back(centerSize);
    idx.push_back(centerIdx);
    edgesLeft -= edgeCount;

    k = 0;
    assert(!pq.empty());
    std::queue<int> anchor_q;
    ShortestPathPriorityQueue cone_pq;
    std::vector<int> cone_succ[n]; 
    std::vector<bool> cone_found(n,false);//PETSC_FALSE);

    /** form tree of shortest paths to an anchor **/
    while (!pq.empty()) {
	node = pq.top();pq.pop();
	cone_found[node.vertex] = true;//PETSC_TRUE;
	anchor_q.push(node.vertex);
	for (tie(e,e_end)=out_edges(node.vertex,g); e!=e_end; e++) {
	    int t = target(*e,g);
	    if (!taken[t]) {
		cone_pq.push(PQNode(t,node.vertex,1.0/get(edge_weight_g,*e)));
	    }
	}
    }
    while (!cone_pq.empty()) {
	node = cone_pq.top();cone_pq.pop();
	if (!cone_found[node.vertex]) {
	    cone_succ[node.pred].push_back(node.vertex);
	    cone_found[node.vertex] = true;//PETSC_TRUE;
	    for (tie(e,e_end)=out_edges(node.vertex,g); e!=e_end; e++) {
		int t = target(*e,g);
		if (!taken[t] && !cone_found[t]) {
		    cone_pq.push(PQNode(t,node.vertex,node.dist+1.0/get(edge_weight_g,*e)));
		}
	    }
	}
    }

    while (!anchor_q.empty()) {
	/** CONE CUT **/
	int anchor = anchor_q.front();anchor_q.pop();
	if (!taken[anchor]) {
	    int v;
	    int thisSize = 0;
	    std::vector<int> thisIdx;
	    std::queue<int> q;
	    ShortestPathPriorityQueue mycone_pq;
	    std::vector<bool> mycone_taken(n,false);//);PETSC_FALSE);
	    int initialInternalConeEdges = 0;

	    boundary = 0;
	    edgeCount = 0;
	    q.push(anchor);
	    while (!q.empty()) {
		v = q.front();q.pop();
		taken[v] = true;//PETSC_TRUE;
		mycone_taken[v] = true;//PETSC_TRUE;
		thisIdx.push_back(g.local_to_global(v));
		thisSize++;
		for (i=cone_succ[v].begin();i!=cone_succ[v].end();i++) {
		    q.push(*i);
		}
		for (tie(e,e_end)=out_edges(v,g); e!=e_end; e++) {
		    int t = target(*e,g);
		    if (!taken[t]) {
			mycone_pq.push(PQNode(t,v,1.0/get(edge_weight_g,*e)));
			boundary += get(edge_weight_g,*e);
			edgeCount++;
		    } else if (mycone_taken[t]) {
			boundary -= get(edge_weight_g,*e);
			initialInternalConeEdges++;
		    }
		}
	    }
	    if (initialInternalConeEdges < edgesLeft) {
		while (initialInternalConeEdges == 0 ?
		       boundary > (edgeCount+1)*log(edgesLeft+1)*2.0/(log(2.0)*epsilon*radius) : 
		       boundary > (edgeCount)*log(edgesLeft*1.0/initialInternalConeEdges)*2.0/(log(2.0)*epsilon*radius))
		{
		    assert(!mycone_pq.empty());
		    node = mycone_pq.top();mycone_pq.pop();
		    if (!mycone_taken[node.vertex]) {
			q.push(node.vertex);
			while (!q.empty()) {
			    v = q.front();q.pop();
			    taken[v] = true;//PETSC_TRUE;
			    mycone_taken[v] = true;//PETSC_TRUE;
			    thisIdx.push_back(g.local_to_global(v));
			    thisSize++;
			    for (i=cone_succ[v].begin();i!=cone_succ[v].end();i++) {
				q.push(*i);
			    }
			    for (tie(e,e_end)=out_edges(v,g); e!=e_end; e++) {
				int t = target(*e,g);
				if (!taken[t]) {
				    mycone_pq.push(PQNode(t,v,node.dist+1.0/get(edge_weight_g,*e)));
				    boundary += get(edge_weight_g,*e);
				    edgeCount++;
				} else if (mycone_taken[t]) {
				    boundary -= get(edge_weight_g,*e);
				}
			    }
			}
		    }
		}
	    }
	    edgesLeft -= edgeCount;
	    size.push_back(thisSize);
	    idx.push_back(thisIdx);
	    x.push_back(anchor);
	    y.push_back(pred[anchor]);
	    k++;
	}
    }
    
  

    /*
    // pseudo cone cut
    while (!pq.empty()) {
    node = pq.top();pq.pop();

    PetscInt thisSize = 1;
    std::vector<PetscInt> thisIdx;
    std::queue<PetscInt> q;

    thisIdx.push_back(g.local_to_global(node.vertex));
    for (i=succ[node.vertex].begin();i!=succ[node.vertex].end();i++) {
    q.push(*i);
    }

    PetscInt v;
    while (!q.empty()) {
    v = q.front();q.pop();
    thisSize++;
    thisIdx.push_back(g.local_to_global(v));
    for (i=succ[v].begin();i!=succ[v].end();i++) {
    q.push(*i);
    }
    }
    size.push_back(thisSize);
    idx.push_back(thisIdx);
    x.push_back(node.vertex);
    y.push_back(pred[node.vertex]);
    k++;
    }
    */

  


    return 0;
}
