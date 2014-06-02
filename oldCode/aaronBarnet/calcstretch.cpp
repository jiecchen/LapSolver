// calculate stretch of a graph

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <queue>

const char DOC_STRING[] = "run like ./calcstretch GRAPH TREE\n";


// template <class T>
// SpielGraph* ijv2graph(int n, int nnz, int* i, int* j, T* v);


class SpielGraph
{
public:
    SpielGraph() : n(0), nnz(0), deg(NULL), nbrs(NULL), nbrsBlock(NULL), wts(NULL), wtsBlock(NULL), backBlock(NULL), back(NULL)
	{
	    //
	}
    ~SpielGraph()
	{
	    if (deg != NULL)
		delete [] deg;
	    if (nbrs != NULL)
	    	delete [] nbrs;
	    if (nbrsBlock != NULL)
	    	delete [] nbrsBlock;
	    if (wts != NULL)
	    	delete [] wts;
	    if (wtsBlock != NULL)
		delete [] wtsBlock;
	    if (backBlock != NULL)
		delete [] backBlock;
	    if (back != NULL)
		delete [] back;
	}
    
    int n;
    int nnz;
    
    int *deg; /* up to n */
    int **nbrs; /* up to n */
    int *nbrsBlock; /* holds the nbrs, nnz */
    float **wts; /* up to n */
    float *wtsBlock; /* up to nnz (has duplicates) */
    int *backBlock; /* up to nnz */
    int **back; /* up to n */  
};

class IJVMatClass
{
public:
    IJVMatClass(const char* fileName) : n(0), nnz(0), i(NULL), j(NULL), v(NULL)
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
    ~IJVMatClass()
	{
	    if (i != NULL)
		delete [] i;
	    if (j != NULL)
		delete [] j;
	    if (v != NULL)
		delete [] v;
	}
    int n;
    int nnz;
    int *i;
    int *j;
    double *v;
};

template <class T>
SpielGraph* ijv2graph(int n, int nnz, int* i, int* j, T* v)
{
    int x;
    
    int ind;
    
    SpielGraph *G;
    
    G = new SpielGraph();
    
    G->n = n;
    G->nnz = 2 * nnz;
    
    G->deg = new int[n]; //(int *) cCalloc(G->n, sizeof(int),"G->deg in ijv2sparse");
    G->nbrs = new int*[n]; //(int **) cCalloc(G->n, sizeof(int *),"G->nbrs in ijv2sparse");
    G->wts = new float*[n]; //(double **) cCalloc(G->n, sizeof(double *), "G->wts in ijv2sparse");
    G->back = new int*[n]; //(int **) cCalloc(G->n, sizeof(int *),"G->back in ijv2sparse");
    
    if (nnz == 0)
    {
	for (x = 0; x < G->n; x++)
	{
	    G->deg[x] = 0;
	    G->nbrs[x] = NULL;
	    G->wts[x] = NULL;
	    G->back[x] = NULL;
	    
	    G->nbrsBlock = NULL;
	    G->wtsBlock = NULL;
	    G->backBlock = NULL;
	}
	return (G);
    }
    
    /* allocate the data structures we will need to 
       store the graph */
    
    //printf("about to allocate mem for nbrsBlock\n");
    G->nbrsBlock = new int[G->nnz]; //(int *) cCalloc(G->nnz, sizeof(int),"nbrsBlock in ijv2sparse");
    G->wtsBlock = new float[G->nnz]; //(double *) cCalloc(G->nnz, sizeof(double),"wtsBlock in ijv2sparse");
    G->backBlock = new int[G->nnz]; //(int *) cCalloc(G->nnz, sizeof(int),"backBlock in ijv2sparse");
    
    for (x = 0; x < G->n; x++)
	G->deg[x] = 0;
    
    /* first, set all the degrees */
    for (x = 0; x < nnz; x++) {
	G->deg[i[x]]++;
	G->deg[j[x]]++;
    }
    
    /* then, allocate space for the data for each vertex */
    ind = 0;
    for (x = 0; x < G->n; x++)
    {
	G->nbrs[x] = &(G->nbrsBlock[ind]);
	G->wts[x] = &(G->wtsBlock[ind]);
	G->back[x] = &(G->backBlock[ind]);
	ind += G->deg[x];
    }
    
    /* and now, populate the graph with edges,
       temporarily resetting the degrees
       (don't worry: they will return to their prev values */

    //printf("preparing to populate graph with edges\n");

    for (x = 0; x < G->n; x++)
	G->deg[x] = 0;
    
    int row; /* i */
    int col; /* j */
    
    for (x = 0; x < nnz; x++)
    {
	row = i[x];
	col = j[x];


	
	G->nbrs[col][G->deg[col]] = row;
	G->nbrs[row][G->deg[row]] = col;

//  	if (row == 14 || col == 14)
//  	{
//  	    printf("forming graph: row: %d col: %d  deg[14]=%d, nbrs[14][G->deg[14]]=%d\n", row, col, G->deg[14], G->nbrs[col][G->deg[col]]);
	    
//  	}
	G->wts[col][G->deg[col]] = v[x];
	G->wts[row][G->deg[row]] = v[x];
	
	G->back[col][G->deg[col]] = G->deg[row];
	G->back[row][G->deg[row]] = G->deg[col];
	
	G->deg[row]++;
	G->deg[col]++;
    }

  //   printf("neighbors of node 14\n");
//     for (int ctr = 0; ctr < G->deg[14]; ctr++)
//     {
// 	printf("neighor %d is %d\n", ctr, G->nbrs[14][ctr]);
//     }
	    
    return G;
}


struct TreeNode
{
    int v;  // node
    TreeNode* p;    // parent node
    TreeNode** children;
    int num_children;
    float d;  // distance to root
    int c;    // which component it's in
};

class Tree
{
public:
    Tree(const char* filename)
	{
	    FILE* fpr;
	    if ((fpr = fopen(filename, "rb")) == NULL)
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
    
	    if (!fread(v, sizeof(float), nnz, fpr))
		exit(1); // cError ("error reading v from file\n");
  	    fclose(fpr);
	    
	    
	    for (int ctr = 0; ctr < nnz; ctr++)
	    {
		i[ctr]--;
		j[ctr]--;
	    }


	    SpielGraph* g = ijv2graph(n, nnz, i, j, v);
	    
	    // let node 0 be root


	    components.resize(n);
	    nodes.resize(n);
	    root = new TreeNode();
	    root->v = 0;
	    root->p = root;
	    root->d = 0.0;
	    root->c = 0;
	    
	    
	    generate_tree(root, g);

	    delete [] i;
	    i = NULL;
	    delete [] j;
	    j = NULL;
	    delete [] v;
	    v = NULL;
	    delete g;
	}
    ~Tree()
	{
	    std::vector<TreeNode*>::iterator it;
	    for (it = nodes.begin(); it != nodes.end(); it++)
	    {
		delete *it;
	    }
	}

    float get_stretch(const SpielGraph* h)
	{
	   return get_stretch_h(root, h);
	}

    std::vector<std::vector<TreeNode*> > components;
    std::vector<TreeNode*> nodes;
    TreeNode* root;
    int n;
    int nnz;
    int* i;
    int* j;
    double* v;
    

private:
    void generate_tree(TreeNode* anchor, const SpielGraph* g)
	{
	    std::queue<TreeNode*> nodeQ;
 	    nodeQ.push(anchor);

 	    while(!nodeQ.empty())
 	    {
		TreeNode* curNode = nodeQ.front();
		nodeQ.pop();
		int vert = curNode->v;

		components[vert].push_back(curNode);
		nodes[vert] = curNode;
		
		if (g->deg[vert] == 1 && curNode->p != curNode)
		    curNode->children = NULL;
		if (curNode->p != curNode)
		    curNode->children = new TreeNode*[g->deg[vert] - 1];
		else
		    curNode->children = new TreeNode*[g->deg[vert]];
		
		int child_ctr = 0;
		for (int ctr = 0; ctr < g->deg[vert]; ctr++)
		{
		    //printf("top of for loop\n");
		    int neighbor_vert = g->nbrs[vert][ctr];
		    //printf("neighbor %d\n", neighbor_vert);
		    if ((curNode->p)->v == neighbor_vert)
			continue;
		    
		    //printf("new tree node...\n");
		    TreeNode* node = new TreeNode();
		    node->v = neighbor_vert;
		    node->p = curNode;
		    node->d = curNode->d + g->wts[vert][ctr];
		    node->c = neighbor_vert;
		    
		    curNode->children[child_ctr++] = node;
		    
		    //printf("gen tree...\n");
		    //generate_tree(node, g);
		    nodeQ.push(node);
		    //printf("comp\n");
		    
		}
		curNode->num_children = child_ctr;
 	    }
	    
	    //printf( "root->v:%d\n", root->v);
	    
	}
    float get_stretch_h(TreeNode* anchor, const SpielGraph* h)
	{
	    // make a vector of all nodes root first/ leaves last (i.e. every child must come after it's parent)
	    std::vector<TreeNode*> nodeList;
	    nodeList.push_back(anchor);
	    int nodeListPos = 0;

	    while (nodeListPos < h->n)
	    {
		TreeNode* curNode = nodeList[nodeListPos++]; // note if not correct tree for graph could have seg fault!
		for (int ctr = 0; ctr < curNode->num_children; ctr++)
		{
		    nodeList.push_back(curNode->children[ctr]);
		}
	    }

	    float total_stretch = 0.0;
	    while (!nodeList.empty())
	    {
		TreeNode* top = nodeList.back();
		nodeList.pop_back();

		if (top->num_children == 0)
		    continue;
		
// 		for (int ctr = 0; ctr < top->num_children; ctr++)
// 		{
// 		    total_stretch += get_stretch_h(top->children[ctr], h);
// 		}
	    
		// calc stretech between first child and top
		int top_vert = top->v;
		int child_comp = top->children[0]->c;
		for (int ctr = 0; ctr < h->deg[top_vert]; ctr++)
		{
		    TreeNode* neighbor = nodes[h->nbrs[top_vert][ctr]];
		    if (neighbor->c == child_comp)
		    {
			total_stretch += (neighbor->d - top->d)/h->wts[top_vert][ctr];
		    }
		}
		top->c = child_comp;
		components[top_vert].pop_back(); // should be alone in comp
		components[top->children[0]->c].push_back(top);
	    

		// now we deal with rest of children
		int old_component = child_comp;
		for (int ctr = 1; ctr < top->num_children; ctr++)
		{
		    int new_comp = (top->children[ctr])->c;
		    int big_comp = 0;
		    int small_comp = 0;
		    if (components[new_comp].size() > components[old_component].size())
		    {
			big_comp = new_comp;
			small_comp = old_component;
		    }
		    else
		    {
			big_comp = old_component;
			small_comp = new_comp;
		    }
		
		    std::vector<TreeNode*>::iterator it;
		    for (it = components[small_comp].begin(); it != components[small_comp].end(); it++)
		    {
			TreeNode* node = *it;
			for (int nbr_ctr = 0; nbr_ctr < h->deg[node->v]; nbr_ctr++)
			{
			    if (nodes[h->nbrs[node->v][nbr_ctr]]->c == big_comp)
				total_stretch += (node->d + nodes[h->nbrs[node->v][nbr_ctr]]->d - 2.0*top->d)/h->wts[node->v][nbr_ctr];
			}
		    }

		    // relabel small_comp as part of big_comp
		    for (it = components[small_comp].begin(); it != components[small_comp].end(); it++)
		    {
			(*it)->c = big_comp;
		    }
		    while (!components[small_comp].empty())
		    {
			components[big_comp].push_back(components[small_comp].back());
			components[small_comp].pop_back();
		    }
		    old_component = big_comp;
		}
	    }
	    return total_stretch;
	}
	
};


int main(int argc, char* argv[])
{
    if (argc < 3)
    {
	printf("%s", DOC_STRING);
	exit(1);
    }

    IJVMatClass* ijv = new IJVMatClass(argv[1]);

    //printf("created ijv\n");

    SpielGraph* graph = ijv2graph(ijv->n, ijv->nnz, ijv->i, ijv->j, ijv->v);

    //printf("finished creating SpielGraph\n");

    Tree* tree = new Tree(argv[2]);
    
    printf("finished getting tree and graph at %lf seconds\n", (double)clock()/(double)CLOCKS_PER_SEC);

    float stretch = tree->get_stretch(graph); 

    printf("the stretch is: %f, the average stretch is: %f\n", stretch, stretch/graph->nnz);
    printf("finished calc stretch at %lf seconds\n", (double)clock()/(double)CLOCKS_PER_SEC);
    delete tree;

    delete ijv;

    delete graph;

    printf("exiting at %lf seconds\n", (double)clock()/(double)CLOCKS_PER_SEC);

    return 0;
    

}
