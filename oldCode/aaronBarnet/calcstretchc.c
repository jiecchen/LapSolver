// calculate stretch of a graph

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "graph.h"
const char DOC_STRING[] = "run like ./calcstretchc GRAPH TREE\n";


typedef struct TreeNode_s
{
    int v;  // node
    struct TreeNode_s* p;    // parent node
    struct TreeNode_s** children; 
    int num_children;
    double d;  // distance to root
    int c;    // which component it's in
} TreeNode;


// linked list for specifying components
typedef struct Cell_s
{
    TreeNode* treeNode;
    struct Cell_s* next;
} Cell;

typedef struct Tree_s
{
    int n;                   // num nodes
    TreeNode* root;          
    Cell** components;       // each component is a linked list of Cell's that point to treeNodes
    Cell* raw_cells;         // raw cells (so can malloc/free them in one block) (not for actually accessing stuff
    int* components_sizes;   // an array of how big all the components are
    TreeNode* nodes;         // where all the nodes are actually stored (parents are always stored before children
    TreeNode** verts;        // table of where nodes corresponding to verts are
} Tree;


Tree* graph2tree(myGraph* g)
{
    Tree* tree = NULL;
    TreeNode* curNode = NULL;
    int vert = 0;
    int nodes_front = 0;  // use nodes as a queue
    int nodes_back = 0;
    int ctr = 0;

    //printf("top of graph2tree\n");


    tree = (Tree*)malloc(sizeof(Tree));
    tree->n = g->n;
    tree->nodes = (TreeNode*)malloc(sizeof(TreeNode)*tree->n);
    tree->verts = (TreeNode**)malloc(sizeof(TreeNode*)*tree->n);
    
    tree->components = (Cell**)malloc(sizeof(Cell*)*tree->n);
    tree->raw_cells = (Cell*)malloc(sizeof(Cell)*tree->n);
    tree->components_sizes = (int*)malloc(sizeof(int)*tree->n);

    //printf("done iwth mallocs\n");

    for (ctr = 0; ctr < tree->n; ctr++)
    {
	//tree->components[ctr] = (TreeNode**)malloc(sizeof(TreeNode*));
	tree->components_sizes[ctr] = 1;
    }
    //printf("done with for loop\n");
    tree->root = &tree->nodes[0];
    
    tree->root->v = 0;
    tree->root->p = tree->root;
    tree->root->d = 0.0;
    tree->root->c = 0;
    nodes_back++;

    //printf("before make tree's while loop\n");
    //printf("don't forget to remove check on node_back\n");
    while (nodes_front < tree->n)
    {
	int child_ctr = 0;
	int neighbor_vert = 0;
	
	curNode = &tree->nodes[nodes_front++];
	vert = curNode->v;
	tree->verts[vert] = curNode;
	tree->components[vert] = &tree->raw_cells[vert];
	tree->components[vert]->treeNode = curNode;
	tree->components[vert]->next = NULL;

	//printf("about to setup children for the %d'th node: %d\n", nodes_front, vert);

	//tree->components[vert][0] = curNode;
       	
	if (g->deg[vert] == 1 && curNode->p != curNode)
	    curNode->children = NULL;
	if (curNode->p != curNode)
	    curNode->children = (TreeNode**)malloc(sizeof(TreeNode*) * g->deg[vert] - 1);
	else
	    curNode->children = (TreeNode**)malloc(sizeof(TreeNode*) * g->deg[vert]);
		
	for (ctr = 0; ctr < g->deg[vert]; ctr++)
	{   
	    TreeNode* node = NULL;

	    neighbor_vert = g->nbrs[vert][ctr];
	    if ((curNode->p)->v == neighbor_vert)
		continue;
	    //if (nodes_back > tree->n)
	    //printf("nodes_back > tree->n\n");

	    node = &tree->nodes[nodes_back++];
	    node->v = neighbor_vert;
	    node->p = curNode;
	    node->d = curNode->d + g->wts[vert][ctr];
	    node->c = neighbor_vert;
	    
	    curNode->children[child_ctr++] = node;
		    
	}

	curNode->num_children = child_ctr;
    }
    return tree;
}

void freeTreeNodeInternals(TreeNode* node)
{
    if (node->num_children > 0)
    {
	free(node->children);
    }
}

void freeTree(Tree* tree)
{
    int nodeItr = 0;
    free(tree->components);
    free(tree->raw_cells);
    free(tree->components_sizes);
    for (nodeItr = 0; nodeItr < tree->n; nodeItr++)
    {
	freeTreeNodeInternals(&tree->nodes[nodeItr]);
    }
    free(tree->nodes);
    free(tree->verts);
}

double calcStretch(Tree* tree, const myGraph* h)
{
    //printf("entered calcStretch\n");
    int nodes_pos = tree->n - 1;
    double total_stretch = 0.0;

    while (nodes_pos >= 0)
    {
	int top_vert = 0;
	int child_comp = 0;
	int ctr = 0;
	int old_component = 0;

	TreeNode* top = &tree->nodes[nodes_pos--];

	if (top->num_children == 0)
	    continue;
    
	// calc stretech between first child and top
	top_vert = top->v;
	child_comp = top->children[0]->c;
	for (ctr = 0; ctr < h->deg[top_vert]; ctr++)
	{
	    TreeNode* neighbor = tree->verts[h->nbrs[top_vert][ctr]];
	    if (neighbor->c == child_comp)
	    {
		total_stretch += (neighbor->d - top->d)/h->wts[top_vert][ctr];
	    }
	}
	top->c = child_comp;
	//free(tree->components[top_vert]);
	tree->components[top_vert]->next = tree->components[child_comp];
	tree->components[child_comp] = tree->components[top_vert];
 	tree->components[top_vert] = NULL;
	tree->components_sizes[top_vert] = 0;
	tree->components_sizes[child_comp]++;
	

	// now we deal with rest of children
	old_component = child_comp;
	for (ctr = 1; ctr < top->num_children; ctr++)
	{
	    int new_comp = (top->children[ctr])->c;
	    int big_comp = 0;
	    int small_comp = 0;
	    Cell* itr = NULL;
	    
	    if (tree->components_sizes[new_comp] > tree->components_sizes[old_component])
	    {
		big_comp = new_comp;
		small_comp = old_component;
	    }
	    else
	    {
		big_comp = old_component;
		small_comp = new_comp;
	    }
	    
	    itr = tree->components[small_comp];
	    while (1)
	    {
		TreeNode* node = itr->treeNode;
		int nbr_ctr = 0;
		for (nbr_ctr = 0; nbr_ctr < h->deg[node->v]; nbr_ctr++)
		{
		    if (tree->verts[h->nbrs[node->v][nbr_ctr]]->c == big_comp)
			total_stretch += (node->d + tree->verts[h->nbrs[node->v][nbr_ctr]]->d - 2.0*top->d)/h->wts[node->v][nbr_ctr];
		}
		if (itr->next != NULL)
		    itr = itr->next;
		else
		    break;
	    }

   	    // relabel small_comp as part of big_comp	    
	    itr = tree->components[small_comp];
	    while (1)
	    {
		TreeNode* node = itr->treeNode;
		node->c = big_comp;
		if (itr->next != NULL)
		    itr = itr->next;
		else
		    break;
	    }
	    
	    // combine the linked lists of components
	    itr->next = tree->components[big_comp];
	    tree->components[big_comp] = tree->components[small_comp];
	    tree->components[small_comp] = NULL;
	    tree->components_sizes[big_comp]+= tree->components_sizes[small_comp];
	    tree->components_sizes[small_comp] = 0;

	    old_component = big_comp;
	}
    }
    return total_stretch;
}


int main(int argc, char* argv[])
{

    FILE* graph_fp = NULL;
    FILE* tree_fp = NULL;
    ijvType* ijv = NULL;
    myGraph* graph = NULL;
    myGraph* treeGraph = NULL;
    Tree* tree = NULL;
    
    double stretch = 0.0;


    if (argc < 3)
    {
	printf("%s", DOC_STRING);
	exit(1);
    }
    // deal with graph
    if ((graph_fp = fopen(argv[1], "rt")) == NULL)
	exit(1);
    ijv = binReadIJV(graph_fp);
    fclose(graph_fp);
    graph = ijv2graph(ijv);
    freeIJV(ijv);
    
    // deal with tree
    if ((tree_fp = fopen(argv[2], "rt")) == NULL)
	exit(1);
    ijv = binReadIJV(tree_fp);
    fclose(tree_fp);
    treeGraph = ijv2graph(ijv);
    freeIJV(ijv);
    tree = graph2tree(treeGraph);
    freeGraph(treeGraph);

    printf("finished getting graph and tree at %lf seconds\n", (double)clock()/(double)CLOCKS_PER_SEC);
    
    stretch = calcStretch(tree, graph);
    printf("the stretch is: %f, the average stretch is: %f\n", stretch, stretch/graph->nnz);
    
    freeGraph(graph);
    freeTree(tree);

    printf("exiting at %lf seconds\n", (double)clock()/(double)CLOCKS_PER_SEC);

    return 0;
    

}
