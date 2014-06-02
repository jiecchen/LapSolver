// common tree structs and functions for treeimprove and treeimprove2
#include "treeimprovetree.h"

// parray determines structure of tree
// g determines weights of edges in tree
Tree* pArray2tree(pArray* parray, myGraph* g)
{
    COMP = 0;  // will be a problem if there is more than one tree!

    Tree* tree = NULL;
    TreeNode* curNode = NULL;
    int nodes_front = 0;  // use nodes as a queue
    int nodes_back = 0;
    int ctr = 0;

    tree = (Tree*)malloc(sizeof(Tree));
    tree->n = parray->n;
    tree->raw_nodes = (TreeNode*)malloc(sizeof(TreeNode)*tree->n);
    tree->nodes = (TreeNode**)malloc(sizeof(TreeNode*)*tree->n);
    tree->verts = (TreeNode**)malloc(sizeof(TreeNode*)*tree->n);
    
    // have to know how many children each nodes has for malloc'ing
    int* num_children = (int*)malloc(sizeof(int) * tree->n);
    for (ctr = 0; ctr < tree->n; ctr++)
    {
	num_children[ctr] = 0;
    }
    for (ctr = 0; ctr < tree->n; ctr++)
    {
	if (parray->array[ctr] == ctr)
	    tree->root = &(tree->raw_nodes[ctr]);

	num_children[parray->array[ctr]]++;
    }

    // set up treeNodes (except for assigning children)
    for (ctr = 0; ctr < tree->n; ctr++)
    {
	curNode = &(tree->raw_nodes[ctr]);
	tree->verts[ctr] = curNode;
	curNode->v = ctr;
	curNode->p = &(tree->raw_nodes[parray->array[ctr]]);
	curNode->c = COMP;

	curNode->num_children = 0; // we will increment this in the next loop
	if (num_children[ctr] > 0)
	    curNode->children = (TreeNode**)malloc(sizeof(TreeNode*) * num_children[ctr]);
	else
	    curNode->children = NULL;
	
	// don't need to initialize these values, but we will anyways
	curNode->topDist = 0.0;
	curNode->topVol = 0;
	curNode->totChildDist = 0.0;
	curNode->totChildVol = 0;
    }

    // now give each node its children
    // and assign weight to edges
    for (ctr = 0; ctr < tree->n; ctr++)
    {
	curNode = tree->verts[ctr];
	TreeNode* parentNode = curNode->p;
	if (tree->root == curNode)
	    continue;
	parentNode->children[parentNode->num_children++] = curNode;

	// must get weight of edge
	int nbrItr = 0;
	for (nbrItr = 0; nbrItr < g->deg[curNode->v]; nbrItr++)
	{
	    TreeNode* nbrNode = tree->verts[g->nbrs[curNode->v][nbrItr]];
	    if (nbrNode == curNode->p)
	    {
		curNode->weight = g->wts[curNode->v][nbrItr];
	    }
	}

    }

    // setup nodes (As a BFS Q struct)
    tree->nodes[nodes_back++] = tree->root;
    while (nodes_front < tree->n)
    {
	curNode = tree->nodes[nodes_front++];
	
	for (ctr = 0; ctr < curNode->num_children; ctr++)
	{
	    TreeNode* childNode = curNode->children[ctr];
	    tree->nodes[nodes_back++] = childNode;
	}
    }
    
    free(num_children);

    return tree;
}

// should be called on every TreeNode in tree (raw_nodes)
void freeTreeNodeInternals(TreeNode* node)
{
    if (node->children != NULL)
    {
	free(node->children);
    }
}

// assumes g is a valid tree to begin with (use when tree is passed as a binijv)
Tree* graph2tree(myGraph* g)
{
    Tree* tree = NULL;
    TreeNode* curNode = NULL;
    int vert = 0;
    int nodes_front = 0;  // use nodes as a queue
    int nodes_back = 0;
    int ctr = 0;

    tree = (Tree*)malloc(sizeof(Tree));
    tree->n = g->n;
    tree->nodes = (TreeNode**)malloc(sizeof(TreeNode*)*tree->n);
    tree->raw_nodes = (TreeNode*)malloc(sizeof(TreeNode) * tree->n);
    tree->verts = (TreeNode**)malloc(sizeof(TreeNode*)*tree->n);
   

    tree->root = &tree->raw_nodes[0];
    
    tree->root->v = 0;
    tree->root->p = tree->root;
    tree->root->c = COMP;
    nodes_back++;

    //printf("before make tree's while loop\n");
    while (nodes_front < tree->n)
    {
	int child_ctr = 0;
	int neighbor_vert = 0;
	int num_children = 0;
	curNode = &tree->raw_nodes[nodes_front];
	tree->nodes[nodes_front++] = curNode;
	vert = curNode->v;
	tree->verts[vert] = curNode;
       	
	if (g->deg[vert] == 1 && curNode->p != curNode)
	    curNode->children = NULL;
	if (curNode->p != curNode)
	    num_children = g->deg[vert] - 1;
	else
	    num_children = g->deg[vert];

	curNode->children = (TreeNode**)malloc(sizeof(TreeNode*) * num_children);
	curNode->num_children = num_children;
	for (ctr = 0; ctr < g->deg[vert]; ctr++)
	{   
	    TreeNode* node = NULL;

	    neighbor_vert = g->nbrs[vert][ctr];
	    if ((curNode->p)->v == neighbor_vert)
		continue;
	    
	    node = &tree->raw_nodes[nodes_back++];
	    node->v = neighbor_vert;
	    node->p = curNode;
	    node->weight = g->wts[vert][ctr];
	    node->c = COMP;
	    curNode->children[child_ctr++] = node;
	}
	
    }
    return tree;
}

// can be called on the masterTree or any subtree (subtrees don't keep a verts or raw_nodes array but
// are otherwise identical
void freeTree(Tree* tree)
{
    int i = 0;
    if (tree->raw_nodes != NULL)
    {
	for (i = 0; i < tree->n; i++)
	{
	    freeTreeNodeInternals(&tree->raw_nodes[i]);
	}
	free(tree->raw_nodes);  // wait need to free nodes individually (malloc'd mem b/c variable num children!)
    }
    free(tree->nodes);
    if (tree->verts != NULL)
	free(tree->verts);
    free(tree);
    tree=NULL;
}


// same output format as binWriteIJV but we have different input
void binWriteTree2IJV(Tree* tree, const char* fileName)
{
    FILE* fp = NULL;
    int num_edges = tree->n - 1;
    int nodeItr = 0;
    if ((fp = fopen(fileName, "wt")) == NULL)
	exit(1);
    
    fwrite(&tree->n, sizeof(int), 1, fp);
    fwrite(&num_edges, sizeof(int), 1, fp);

    int* i = malloc(sizeof(int)*num_edges);
    int* j = malloc(sizeof(int)*num_edges);
    double* v = malloc(sizeof(double)*num_edges);

    int edgeCtr = 0;
    for (nodeItr = 0; nodeItr < tree->n; nodeItr++)
    {
	if (tree->nodes[nodeItr]->p == tree->nodes[nodeItr])
	    continue;
	i[edgeCtr] = tree->nodes[nodeItr]->v + 1;    // matlab start indices at 1 not 0!
	j[edgeCtr] = tree->nodes[nodeItr]->p->v + 1;
	v[edgeCtr] = 1.0/tree->nodes[nodeItr]->weight;   // we recipricals in graph.h - need to undo this 
	edgeCtr++;
    }

    fwrite(i, sizeof(int), num_edges, fp);
    fwrite(j, sizeof(int), num_edges, fp);
    fwrite(v, sizeof(double), num_edges, fp);

    fclose(fp);

    free(i);
    free(j);
    free(v);
}

void binWriteTree2pArray(Tree* tree, const char* fileName)
{
    FILE* fp;
    if ((fp = fopen(fileName, "wt")) == NULL)
	exit(1);

    fwrite(&tree->n, sizeof(int), 1, fp);
    
    int vertCtr = 0;
    for (vertCtr = 0; vertCtr < tree->n; vertCtr++)
    {
	int parentVert = tree->verts[vertCtr]->p->v;
	fwrite(&parentVert, sizeof(int), 1, fp);
    }
    fclose(fp);
}
