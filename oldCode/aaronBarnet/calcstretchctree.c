#include "calcstretchctree.h"
#include <stdlib.h>
// parray determines structure of tree
// g determines distances on tree
Tree* pArray2tree(pArray* parray, myGraph* g)
{
    Tree* tree = NULL;
    TreeNode* curNode = NULL;
    int nodes_front = 0;  // use nodes as a queue
    int nodes_back = 0;
    int ctr = 0;

    //printf("top of graph2tree\n");


    tree = (Tree*)malloc(sizeof(Tree));
    tree->n = parray->n;
    tree->raw_nodes = (TreeNode*)malloc(sizeof(TreeNode)*tree->n);
    tree->nodes = (TreeNode**)malloc(sizeof(TreeNode*)*tree->n);
    tree->verts = (TreeNode**)malloc(sizeof(TreeNode*)*tree->n);
    
    tree->components = (Cell**)malloc(sizeof(Cell*)*tree->n);
    tree->raw_cells = (Cell*)malloc(sizeof(Cell)*tree->n);
    tree->components_sizes = (int*)malloc(sizeof(int)*tree->n);

    //printf("done iwth mallocs\n");

    for (ctr = 0; ctr < tree->n; ctr++)
    {
	tree->components_sizes[ctr] = 1;
    }

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
	curNode->c = ctr;

	curNode->num_children = 0; // we will increment this in the next loop
	if (num_children[ctr] > 0)
	    curNode->children = (TreeNode**)malloc(sizeof(TreeNode*) * num_children[ctr]);
	else
	    curNode->children = NULL;

	tree->components[ctr] = &tree->raw_cells[ctr];
	tree->components[ctr]->treeNode = curNode;
	tree->components[ctr]->next = NULL;
    }

    // now give its node its children
    for (ctr = 0; ctr < tree->n; ctr++)
    {
	curNode = tree->verts[ctr];
	TreeNode* parentNode = curNode->p;
	if (tree->root ==  curNode)
	    continue;
	parentNode->children[parentNode->num_children++] = curNode;
    }

    // setup nodes (As a BFS Q struct and assign dists)
    tree->nodes[nodes_back++] = tree->root;
    tree->root->d = 0.0;
    while (nodes_front < tree->n)
    {
	curNode = tree->nodes[nodes_front++];
	
	for (ctr = 0; ctr < curNode->num_children; ctr++)
	{
	    TreeNode* childNode = curNode->children[ctr];
	    tree->nodes[nodes_back++] = childNode;
	    // must find weight of edge in graph
	    int edgeItr = 0;
	    for (edgeItr = 0; edgeItr < g->deg[curNode->v]; edgeItr++)
	    {
		if (g->nbrs[curNode->v][edgeItr] == childNode->v)
		{
		    childNode->d = curNode->d + g->wts[curNode->v][edgeItr];
		}
	    }
	}
    }
    
    free(num_children);

    return tree;
}

/* Tree* graph2tree(myGraph* g) */
/* { */
/*     Tree* tree = NULL; */
/*     TreeNode* curNode = NULL; */
/*     int vert = 0; */
/*     int nodes_front = 0;  // use nodes as a queue */
/*     int nodes_back = 0; */
/*     int ctr = 0; */

/*     //printf("top of graph2tree\n"); */


/*     tree = (Tree*)malloc(sizeof(Tree)); */
/*     tree->n = g->n; */
/*     tree->nodes = (TreeNode*)malloc(sizeof(TreeNode)*tree->n); */
/*     tree->verts = (TreeNode**)malloc(sizeof(TreeNode*)*tree->n); */
    
/*     tree->components = (Cell**)malloc(sizeof(Cell*)*tree->n); */
/*     tree->raw_cells = (Cell*)malloc(sizeof(Cell)*tree->n); */
/*     tree->components_sizes = (int*)malloc(sizeof(int)*tree->n); */

/*     //printf("done iwth mallocs\n"); */

/*     for (ctr = 0; ctr < tree->n; ctr++) */
/*     { */
/* 	//tree->components[ctr] = (TreeNode**)malloc(sizeof(TreeNode*)); */
/* 	tree->components_sizes[ctr] = 1; */
/*     } */
/*     //printf("done with for loop\n"); */
/*     tree->root = &tree->nodes[0]; */
    
/*     tree->root->v = 0; */
/*     tree->root->p = tree->root; */
/*     tree->root->d = 0.0; */
/*     tree->root->c = 0; */
/*     nodes_back++; */

/*     //printf("before make tree's while loop\n"); */
/*     while (nodes_front < tree->n) */
/*     { */
/* 	int child_ctr = 0; */
/* 	int neighbor_vert = 0; */
	
/* 	curNode = &tree->nodes[nodes_front++]; */
/* 	vert = curNode->v; */
/* 	tree->verts[vert] = curNode; */
/* 	tree->components[vert] = &tree->raw_cells[vert]; */
/* 	tree->components[vert]->treeNode = curNode; */
/* 	tree->components[vert]->next = NULL; */

/* 	//tree->components[vert][0] = curNode; */
       	
/* 	if (g->deg[vert] == 1 && curNode->p != curNode) */
/* 	    curNode->children = NULL; */
/* 	if (curNode->p != curNode) */
/* 	    curNode->children = (TreeNode**)malloc(sizeof(TreeNode*) * g->deg[vert] - 1); */
/* 	else */
/* 	    curNode->children = (TreeNode**)malloc(sizeof(TreeNode*) * g->deg[vert]); */
		
/* 	for (ctr = 0; ctr < g->deg[vert]; ctr++) */
/* 	{    */
/* 	    TreeNode* node = NULL; */

/* 	    neighbor_vert = g->nbrs[vert][ctr]; */
/* 	    if ((curNode->p)->v == neighbor_vert) */
/* 		continue; */
	    
/* 	    node = &tree->nodes[nodes_back++]; */
/* 	    node->v = neighbor_vert; */
/* 	    node->p = curNode; */
/* 	    node->d = curNode->d + g->wts[vert][ctr]; */
/* 	    node->c = neighbor_vert; */
	    
/* 	    curNode->children[child_ctr++] = node; */
		    
/* 	} */
/* 	curNode->num_children = child_ctr; */
/*     } */
/*     return tree; */
/* } */

void freeTreeNodeInternals(TreeNode* node)
{
    if (node->children != NULL)
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
	freeTreeNodeInternals(tree->nodes[nodeItr]);
    }
    free(tree->raw_nodes);
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

	TreeNode* top = tree->nodes[nodes_pos--];

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

pArray* tree2pArray(Tree* tree)
{
    int* q = malloc(sizeof(int) * tree->n);
    int q_front = 0;
    int q_back = 0;
    pArray* parray = newPArray(tree->n);
    int ctr = 0;
    for (ctr = 0; ctr < tree->n; ctr++)
	parray->array[ctr] = -1; // negative value means unvisited
    // will let vert 0 be root
    parray->array[0] = 0;
    q[q_back++] = 0;
    while (q_front < tree->n)
    {
	int curVert = q[q_front++];
	int childItr = 0;
	for (childItr = 0;childItr < tree->verts[curVert]->num_children; childItr++)
	{
	    int childVert = tree->verts[curVert]->children[childItr]->v;
	    parray->array[childVert] = curVert;
	    q[q_back++] = childVert;
	}
    }
    free(q);
    return parray;
}

int checkTree(Tree* tree)
{
    TreeNode** q = malloc(sizeof(TreeNode*) * tree->n);
    int q_front = 0;
    int q_back = 0;
    
    double maxDist = 0.0;

    q[q_back++] = tree->root;
    while (q_front < tree->n)
    {
	TreeNode* curNode = q[q_front++];
	int childItr = 0;
	for (childItr = 0; childItr < curNode->num_children; childItr++)
	{
	    if (curNode->children[childItr]->d > maxDist)
		maxDist = curNode->children[childItr]->d;
	    q[q_back++] = curNode->children[childItr];
	}
	    
    }
    if (q_front != tree->n || q_back != tree->n)
    {
	printf("error in tree\n");
	exit(1);
    }
    printf("check OK! maxDist: %f\n", maxDist);
    return 1;
}
