// take in a low stretch spanning tree and try to improve it (this one focuses on lowstretch)
// this version creates virtual edges that allow us to guarantee that stretch is always improved
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "graph.h"
#include "treeimprovetree.h"

// Globals and consts
const char* DOC_STRING = "run as ./treeimprove4 GRAPH.BINIJV TREE.PARRAY IMPROVED_TREE.PARRAY\n";
//int COMP;                                    // remakeTree first "COMP++" and then labels all nodes as COMP
enum SPLIT_MODES { RANDOM, VOLUME, DIST};   // RANDOM = choose rand edge | VOLUME = edge evenly splitting nodes | DIST = sum of all paths to node is smallest | ALL_EDGES is where we take every single edge in tree breaking and reforming it
enum TREE_FORMAT_TYPE {PARRAY, BINIJV};
const enum SPLIT_MODES SPLIT_MODE = DIST;    // how we will split up tree
const enum TREE_FORMAT_TYPE TREE_FORMAT = PARRAY;
const int NUM_ITERATIONS = 10;

typedef struct Edge_s
{
    int v;
    double w;
    struct Edge_s* next;
} Edge;


typedef struct
{
    int n; // num vertices in graph
    Edge** verts; // each vert has a list of attached edges
} VirtGraph;

VirtGraph* newVirtGraph(int n)
{
    VirtGraph* out = malloc(sizeof(VirtGraph));
    out->n = n;
    out->verts = (Edge**)malloc(sizeof(Edge*) * n);
    int ctr = 0;
    for (ctr = 0; ctr < n; ctr++)
	out->verts[ctr] = NULL;
    
    return out;
}

void freeEdgeList(Edge* in)
{
    Edge* cur = in;
    while (cur)
    {
	Edge* next = cur->next;
	free(cur);
	cur = next;
    }

    in = NULL;
}

void freeVirtGraph(VirtGraph* in)
{
    int ctr = 0;
    for (ctr = 0; ctr < in->n; ctr++)
	freeEdgeList(in->verts[ctr]);
    free(in->verts);
    free(in);
    in = NULL;
}

// will insert edge twice (for each direction)
// we take care of duplicate edges (combine weights like resistors in parallel)
void vgInsertEdge(VirtGraph* vg, int i, int j, double weight)
{
    if (i == j)   // don't add self-loops
	return;

    // insert forward edge
    Edge* e = malloc(sizeof(Edge));
    if (e == NULL)
	fprintf(stderr, "malloc failed\n");
    e->v = j;
    e->w = weight;
    Edge* itr = vg->verts[i];
    Edge* last = NULL;
    while (itr && (itr->v < e->v))
    {
	last = itr;
	itr = itr->next;
    }
    // deal with dup
    if (itr && itr->v == e->v)
    {
	itr->w = 1.0/((1.0/itr->w) + (1.0/e->w));
	free(e);
    }
    else
    {
	if (last)
	    last->next = e;
	else
	    vg->verts[i] = e;
	e->next = itr;
    }
    // insert backwards edge
    e = malloc(sizeof(Edge));
    if (e == NULL)
	fprintf(stderr, "malloc failed\n");
    e->v = i;
    e->w = weight;
    itr = vg->verts[i];
    last = NULL;
    while (itr && (itr->v < e->v))
    {
	last = itr;
	itr = itr->next;
    }
    if (itr && itr->v == e->v)
    {
	itr->w = 1.0/((1.0/itr->w) + (1.0/e->w));
	free(e);
    }
    else
    {
	if (last)
	    last->next = e;
	else
	    vg->verts[i] = e;
	e->next = itr;
    }

/*     e->next = vg->verts[j]; */
/*     vg->verts[j] = e; */

}


// take the anme of a binijv file and return it as a tree struct (must be tree!)
Tree* getBinIJVTree(const char* fileName)
{
    myGraph* graph = NULL;
    Tree* tree = NULL;

    graph = getGraph(fileName);
    tree = graph2tree(graph);

    freeGraph(graph);
    
    return tree;
}

void stretchDistDownTree(Tree* tree, Tree* masterTree, myGraph* graph, VirtGraph* vg, int c)
{
    int nodeItr = 0;
    for (nodeItr = 0; nodeItr < tree->n; nodeItr++)
    {
	TreeNode* curNode = tree->nodes[nodeItr];

	// root is special case
	if (curNode->p == curNode)
	{
	    curNode->topDist = 0.0;
	    curNode->topVol = 0;
	}
	
	int nbrItr = 0;
	double curNodeVolContrib = 0.0;
	// iterate through main graph
	for (nbrItr = 0; nbrItr < graph->deg[curNode->v]; nbrItr++)
	{
	    TreeNode* nbrNode = masterTree->verts[graph->nbrs[curNode->v][nbrItr]];
	    if (nbrNode->c == c)
		curNodeVolContrib += 1.0/graph->wts[curNode->v][nbrItr];
	}
	// iterate through virt graph
	Edge* edge = vg->verts[curNode->v];
	while (edge)
	{
	    TreeNode* nbrNode = masterTree->verts[edge->v];
	    if (nbrNode->c == c)
		curNodeVolContrib += 1.0/edge->w;
	    edge = edge->next;
	}
	curNodeVolContrib += curNode->topVol;

	int childItr = 0;
	for (childItr = 0; childItr < curNode->num_children; childItr++)
	{
	    TreeNode* childNode = curNode->children[childItr];
	    double childVolContrib = childNode->totChildVol;
	    // only increase volume for edges from child to other tree
	    // iterate through main graph
	    int childNbrItr = 0;
	    for (childNbrItr = 0; childNbrItr < graph->deg[childNode->v]; childNbrItr++)
	    {
		TreeNode* childNbrNode = masterTree->verts[graph->nbrs[childNode->v][childNbrItr]];
		if (childNbrNode->c == c)
		    childVolContrib += 1.0/graph->wts[childNode->v][childNbrItr];
	    }
	    // iterate through virt graph
	    Edge* edge = vg->verts[childNode->v];
	    while (edge)
	    {
		TreeNode* childNbrNode = masterTree->verts[edge->v];
		if (childNbrNode->c == c)
		    childVolContrib += 1.0/edge->w;
		edge = edge->next;
	    }
	    
	    double childDistContrib = childNode->totChildDist + (childNode->weight * childVolContrib);
	    childNode->topVol = curNodeVolContrib + curNode->totChildVol - childVolContrib;
	    childNode->topDist = curNode->topDist + curNode->totChildDist - childDistContrib + (curNode->topVol * childNode->weight);
	}
    }
}

void stretchDistUpTree(Tree* tree, Tree* masterTree, myGraph* graph, VirtGraph* vg, int c)
{
    int nodeItr = 0;
    for (nodeItr = tree->n - 1; nodeItr >= 0; nodeItr--)
    {
	TreeNode* curNode = tree->nodes[nodeItr];
	double totalVol = 0;
	double totalDist = 0.0;
	int childItr = 0;
	for (childItr = 0; childItr < curNode->num_children; childItr++)
	{
	    TreeNode* childNode = curNode->children[childItr];
	    
	    double childVol = childNode->totChildVol;
	    // only increase volume for edges from child to other tree
	    // iterate through main graph
	    int childNbrItr = 0;
	    for (childNbrItr = 0; childNbrItr < graph->deg[childNode->v]; childNbrItr++)
	    {
		TreeNode* childNbrNode = masterTree->verts[graph->nbrs[childNode->v][childNbrItr]];
		if (childNbrNode->c == c)
		    childVol += 1.0/graph->wts[childNode->v][childNbrItr];
	    }
	    // iterate through virtual graph
	    Edge* edge = vg->verts[childNode->v];
	    while (edge)
	    {
		TreeNode* childNbrNode = masterTree->verts[edge->v];
		if (childNbrNode->c == c)
		    childVol += 1.0/edge->w;
		edge = edge->next;
	    }
	    totalDist += childNode->totChildDist + childNode->weight * childVol;
	    totalVol += childVol;
	}
	curNode->totChildDist = totalDist;
	curNode->totChildVol = totalVol;
    }
}

void stretchDistTree(Tree* tree, Tree* masterTree, myGraph* graph, VirtGraph* vg, int c)
{
    //printf("start of stretchDistTree...");
    stretchDistUpTree(tree, masterTree, graph, vg, c);
    stretchDistDownTree(tree, masterTree, graph, vg, c);
    //printf("end of stretchDistTree...");
}

void distDownTree(Tree* tree)
{
    int nodeItr = 0;
    for (nodeItr = 0; nodeItr < tree->n; nodeItr++)
    {
	TreeNode* curNode = tree->nodes[nodeItr];

	// root is special case
	if (curNode->p == curNode)
	{
	    curNode->topDist = 0.0;
	    curNode->topVol = 0;
	}
	
	int childItr = 0;
	for (childItr = 0; childItr < curNode->num_children; childItr++)
	{
	    TreeNode* childNode = curNode->children[childItr];
	    double childVolContrib = childNode->totChildVol + 1;
	    double childDistContrib = childNode->totChildDist + (childNode->weight * childVolContrib);
	    childNode->topVol = curNode->topVol + curNode->totChildVol - childVolContrib + 1;
	    childNode->topDist = curNode->topDist + curNode->totChildDist - childDistContrib + (childNode->topVol * childNode->weight);
	}
    }
}

void distUpTree(Tree* tree)
{
    int nodeItr = 0;
    for (nodeItr = tree->n - 1; nodeItr >= 0; nodeItr--)
    {
	TreeNode* curNode = tree->nodes[nodeItr];
	int totalVol = 0;
	double totalDist = 0.0;
	int childItr = 0;
	for (childItr = 0; childItr < curNode->num_children; childItr++)
	{
	    TreeNode* childNode = curNode->children[childItr];
	    totalVol += childNode->totChildVol + 1;
	    totalDist += childNode->totChildDist + childNode->weight * (childNode->totChildVol + 1);
	}
	curNode->totChildDist = totalDist;
	curNode->totChildVol = totalVol;
    }
}
void distTree(Tree* tree)
{
    distUpTree(tree);
    distDownTree(tree);
}

// for debugging only
void checkTree(Tree* tree, myGraph* graph, const char* str)
{
    TreeNode* curNode;
    TreeNode** nodes = (TreeNode**)malloc(sizeof(TreeNode*) * tree->n);
    int nodes_front = 0;
    int nodes_back = 0;
    nodes[nodes_back++] = tree->root;
    
    // check that we have right number of nodes in tree
    while (nodes_front < tree->n)
    {
	curNode = nodes[nodes_front++];
	int childItr = 0;
	for (childItr = 0; childItr < curNode->num_children; childItr++)
	{
	    TreeNode* childNode = curNode->children[childItr];
	    if (childNode->c != tree->root->c)
	    {
		fprintf(stderr,"child comp does not match root comp\n");
	    }
	    if (nodes_back >= tree->n)
 	    { 
 		fprintf(stderr, "CHECKTREE FAILED %s : nodes_back:%d,   tree->n: %d\n",str, nodes_back, tree->n);
		exit(1);
 	    }
	    nodes[nodes_back++] = childNode;
	}
    }

    // check that edge weights are correct
    double small_double = .00001;
    int nodeItr = 0;
    for (nodeItr = 0; nodeItr < tree->n; nodeItr++)
    {
	curNode = nodes[nodeItr];
	int curVert = curNode->v;
	int parentVert = curNode->p->v;
	if (curNode == tree->root)
	    continue;
	
	int graphItr = 0;
	int foundNbr = 0;
	for (graphItr = 0; graphItr < graph->deg[curVert]; graphItr++)
	{
	    if (graph->nbrs[curVert][graphItr] == parentVert)
	    {
		if (graph->wts[curVert][graphItr] < curNode->weight - small_double
		    || graph->wts[curVert][graphItr] > curNode->weight + small_double)
		{
		    fprintf(stderr, "CHECKTREE FAILED %s : weights on tree did not match graph\n", str);
		    exit(1);
		}
		foundNbr = 1;
		break;
	    }
	}
	if (!foundNbr)
	{
	    fprintf(stderr, "CHECKTREE FAILED %s : edge in tree not found in graph\n", str);
	    exit(1);
	}
    }

    free(nodes);
    printf("CHECKTREE %s:passed\n", str);
}

// redesign it to take new root (so don't have to traverse up tree and n precalced)
// return a new tree given node (i.e. create new tree/ set up root and nodes)
// leave raw_nodes and verts blank
Tree* remakeTree(TreeNode* root, int n)
{
    COMP++;
    //printf("in remake tree\n");
    Tree* tree = (Tree*)malloc(sizeof(Tree));
    TreeNode* curNode = NULL;
    tree->raw_nodes = NULL;
    tree->verts = NULL;
    tree->root = root;
    tree->n = n;
    int nodes_front = 0;
    int nodes_back = 0;

    //printf("before tree->nodes\n");
    tree->nodes = (TreeNode**)malloc(sizeof(TreeNode*) * tree->n);
    tree->nodes[nodes_back++] = tree->root;
    while (nodes_front < tree->n)
    {
	curNode = tree->nodes[nodes_front++];
	curNode->c = COMP;
	int childItr = 0;
	for (childItr = 0; childItr < curNode->num_children; childItr++)
	{
	    TreeNode* childNode = curNode->children[childItr];
	    if (nodes_back >= tree->n)
 	    { 
 		printf("nodes_back:%d,   tree->n: %d\n", nodes_back, tree->n);
 	    }
	    tree->nodes[nodes_back++] = childNode;
	}
    }
    //printf("end of remake tree\n");
    return tree;
}

// redo tree so node is the root of tree   // don't forget to deal with weights!
void rehangTree(Tree* tree, TreeNode* node)
{
    TreeNode* curNode;
    TreeNode* lastNode;
    double lastWeight;
    if (tree->root == node)
    {
	//printf("rehanging tree on root!\n");
	return;
    }
    

    lastNode = node;
    curNode = node->p;
    node->num_children++;
    if (node->children != NULL)
	node->children = (TreeNode**)realloc(node->children, sizeof(TreeNode*) * node->num_children);
    else
	node->children = (TreeNode**)malloc(sizeof(TreeNode*) * node->num_children);
    node->children[node->num_children - 1] = node->p;
    //node->p->p = node;
    lastWeight = node->weight;
    node->p = node;
    
    // must keep moving up parent path until reach old root
    while (curNode->p != curNode)
    {
	// must make old parent a child
	
	//curNode->num_children++;
	//curNode->children = (TreeNode**)realloc(curNode->children, sizeof(TreeNode*) * curNode->num_children);
	//curNode->children[curNode->num_children - 1] = curNode->p;
	int childItr = 0;
	for (childItr = 0; childItr < curNode->num_children; childItr++)
	{
	    if (curNode->children[childItr] == lastNode)
	    {
		double tempWeight;
		curNode->children[childItr] = curNode->p;
		curNode->p = lastNode;
		//curNode->weight = lastNode->weight;
		tempWeight = curNode->weight; // after reassigning curNode/lastNode?
		curNode->weight = lastWeight;
		lastWeight = tempWeight;

		lastNode = curNode;
		curNode = curNode->children[childItr];
		break;
	    }
	}

	//curNode->p = lastNode;

	//lastNode = curNode;
	//curNode = curNode->children[curNode->num_children - 1];
    }
    
/*     if (curNode != tree->root) */
/*     { */
/* 	printf("node not in tree\n"); */
/* 	exit(1); */
/*     } */
/*     else */
/* 	printf("node is in tree\n"); */

    // curNode is now root
    curNode->p = lastNode;
    // want to eliminate lastNode as a child
    int childItr = 0;
    int childCtr = 0;
    for (childItr = 0; childItr < curNode->num_children; childItr++)
    {
	if (curNode->children[childItr] == lastNode)
	    continue;
	curNode->children[childCtr++] = curNode->children[childItr];
    }
    curNode->num_children--;
    curNode->weight = lastWeight;  //??
    tree->root = node;
}

// choose rand edge
void divideRandom(Tree* tree, TreeNode** nodeTop, TreeNode** nodeBottom)
{
    int randNodeItr = 0;
    while (1)
    {
     randNodeItr = floor((double)rand()/(double)RAND_MAX * (double)tree->n);
     if (tree->nodes[randNodeItr] == tree->root)
         continue;
     //int score = tree->nodes[randNodeItr]->totChildVol;
     //if (score > tree->n / 4 && score < 3 * tree->n / 4)
     break;
     // pass
    }
    *nodeBottom = tree->nodes[randNodeItr];
    *nodeTop = (*nodeBottom)->p;
}

// choose edge evenly splitting nodes
void divideByVolume(Tree* tree, TreeNode** nodeTop, TreeNode** nodeBottom)
{
    int best_split = 999999999;
    int half_size = tree->n / 2;
    int nodeItr = 0;
    for (nodeItr = 0; nodeItr < tree->n; nodeItr++)
    {
	TreeNode* curNode = tree->nodes[nodeItr];
	
	if (curNode->p == curNode)   // this line should not be necessary
	    continue;
	
	int score = half_size - curNode->totChildVol;
	score = (score < 0 ? -1 * score : score);
	if (score < best_split)
	{
	    best_split = score;
	    *nodeBottom = curNode;
	}
    }
    *nodeTop = (*nodeBottom)->p;
}

// choose nodeBottom s.t. sum of all paths to node is smallest
void divideByDist(Tree* tree, TreeNode** nodeTop, TreeNode** nodeBottom)
{
    double best_dist = 99999999999999999999.9;
    int nodeItr = 0;
    for (nodeItr = 0; nodeItr < tree->n; nodeItr++)
    {
        TreeNode* curNode = tree->nodes[nodeItr];

        if (curNode->p == curNode)   // this line should not be necessary                                                        
            continue;

        double score = curNode->totChildDist + curNode->topDist;
        if (score < best_dist)
        {
            best_dist = score;
            *nodeBottom = curNode;
        }
    }
    *nodeTop = (*nodeBottom)->p;
}


// recursively break up graph and join back together to improve stretch
// tree: tree/subtree that is to improved
// graph:  we need to know what alternative edges we can use to rejoin trees
// virtgraph: every time we break up tree we add virtual edges so that we can properaly calc stretch
// masterTree: this is the tree that has a verts array (maybe make verts a seperate struct
Tree* improveTree(Tree* treeIn, myGraph* graph, VirtGraph* vg, Tree* masterTree)
{
    Tree* tree = remakeTree(treeIn->root, treeIn->n);  // why do I need this (get seg fault w/out)?
    //checkTree(tree, graph, "tree");

    // base case
    if (tree->n < 2) //masterTree->n)
	return tree;

    //printf("top of improve tree\n");

    distTree(tree);

    TreeNode* nodeTop = NULL;     // the edge nodeBottom->nodeTop is where we will split tree
    TreeNode* nodeBottom = NULL;
    Tree* treeTop = NULL;         // tree minuse treeBottom
    Tree* treeBottom = NULL;      // this is the half the tree with nodeBottom as it's root
    Tree* treeNew = NULL;         // what we will return

    
    // decide how to split up tree
    if (SPLIT_MODE == RANDOM)
	divideRandom(tree, &nodeTop, &nodeBottom);
    else if (SPLIT_MODE == VOLUME)
	divideByVolume(tree, &nodeTop, &nodeBottom);
    else if (SPLIT_MODE == DIST)
	divideByDist(tree, &nodeTop, &nodeBottom);
    else
    {
	fprintf(stderr, "SPLIT_MODE invalid\n");
	exit(2);
    }
    
    //printf("made decisiong about splitting tree\n");

    // actual splitting of tree between nodeTop and nodeBottom
    int childCtr = 0;
    int nodeItr = 0;    
    for (nodeItr = 0; nodeItr < nodeTop->num_children; nodeItr++)
    {
	TreeNode* childNode = nodeTop->children[nodeItr];
	if (childNode == nodeBottom)
	    continue;

	nodeTop->children[childCtr++] = childNode;
    }
    
    //printf("before remake tree\n");
    nodeTop->num_children--;
    treeTop = remakeTree(tree->root, nodeTop->topVol + nodeTop->totChildVol - nodeBottom->totChildVol);
    //printf("after remake tree\n");

    //printf("before remake tree\n");
    nodeBottom->p = nodeBottom;
    treeBottom = remakeTree(nodeBottom, nodeBottom->totChildVol + 1);
    //printf("after remake tree\n");

    // find all edges that go from treeBottom to treeTop and add them as edges to vg
    //printf("before adding virt edge...\n");
    int top_comp = nodeTop->c;
    for (nodeItr = 0; nodeItr < treeBottom->n; nodeItr++)
    {
	TreeNode* curNode = treeBottom->nodes[nodeItr];
	
	int vert = curNode->v;

	// look for edges in virtGraph
	Edge* edge = vg->verts[vert];
	while (edge)
	{
	    int nbrVert = edge->v;
	    TreeNode* nbrNode = masterTree->verts[nbrVert];
	    if (nbrNode->c == top_comp);
	    {
		vgInsertEdge(vg, vert, nodeBottom->v, edge->w);
		vgInsertEdge(vg, nodeTop->v, nbrVert, edge->w);
	    }
	    edge = edge->next;
	}

	// look for edges in original graph
	int nbrItr = 0;
	for (nbrItr = 0; nbrItr < graph->deg[vert]; nbrItr++)
	{
	    int nbrVert = graph->nbrs[vert][nbrItr];
	    TreeNode* nbrNode = masterTree->verts[nbrVert];
	    if (nbrNode->c == top_comp)
	    {
		// be careful of self loops 
		vgInsertEdge(vg, vert, nodeBottom->v, graph->wts[vert][nbrItr]);
		vgInsertEdge(vg, nodeTop->v, nbrVert, graph->wts[vert][nbrItr]);
	    }
	}

    }
    //printf("after adding virt edges\n");

    // recursion
    Tree* treeTopTemp = improveTree(treeTop, graph, vg, masterTree);
    freeTree(treeTop);
    treeTop = treeTopTemp;

    // recursion
    Tree* treeBottomTemp = improveTree(treeBottom, graph, vg, masterTree);
    freeTree(treeBottom);
    treeBottom = treeBottomTemp;

    // recalc dists on our new trees (using stretch version
    stretchDistTree(treeBottom, masterTree, graph, vg, treeTop->root->c);
    stretchDistTree(treeTop, masterTree, graph, vg, treeBottom->root->c);
    /* distUpTree(treeBottom); */
/*     distDownTree(treeBottom); */
/*     distUpTree(treeTop); */
/*     distDownTree(treeTop); */

    // look for best edge from treeTop to treeBottom to rejoin trees
    int target_comp = treeBottom->root->c;
    double best_score = 99999999999999999999.9;
    TreeNode* n1 = NULL;
    TreeNode* n2 = NULL;
    double weight = 0.0;
    for (nodeItr = 0; nodeItr < treeTop->n; nodeItr++)
    {
	// should nodes also keep track of edges not in tree?
	TreeNode* curNode = treeTop->nodes[nodeItr];
	int vert = curNode->v;
	int nbrItr = 0;
	for (nbrItr = 0; nbrItr < graph->deg[vert]; nbrItr++)
	{
	    int nbrVert = graph->nbrs[vert][nbrItr];
	    TreeNode* nbrNode = masterTree->verts[nbrVert];
	    if (nbrNode->c == target_comp)
	    {
		double score = curNode->totChildDist + curNode->topDist + ((curNode->totChildVol + curNode->topVol) * graph->wts[vert][nbrItr])
		    + nbrNode->totChildDist + nbrNode->topDist;

/* 		if (curNode->totChildVol + curNode->topVol != nbrNode->totChildVol + nbrNode->topVol) */
/* 		{ */
/* 		    printf("nbr volumes don't match\n"); */
/* 		} */

		if (score < best_score)
		{
		    best_score = score;
		    n1 = curNode;
		    n2 = nbrNode;
		    weight = graph->wts[vert][nbrItr];
		}
	    }
	}
    }
    
    // join the trees
    rehangTree(treeBottom, n2);
    // now make n2 a new child of n1
    n1->num_children++;
    if (n1->children != NULL)
	n1->children = (TreeNode**)realloc(n1->children, sizeof(TreeNode*) * n1->num_children);
    else
	n1->children = (TreeNode**)malloc(sizeof(TreeNode*) * n1->num_children);
    n1->children[n1->num_children - 1] = n2;
    n2->p = n1;
    n2->weight = weight;

    // deal with output and free mem
    treeNew = remakeTree(treeTop->root, tree->n);

    freeTree(tree);
    freeTree(treeTop);
    freeTree(treeBottom);

    //checkTree(treeNew, graph, "treeNdew");
    return treeNew;
}

int main(int argc, char* argv[])
{
    srand(5); // keep results repeatable for testing

    myGraph* graph = NULL;
    Tree* tree = NULL;
    Tree* improvedTree = NULL;
    COMP = 0; // always label nodes with next comp number when making new tree
    if (argc != 4)
    {
	printf("%s", DOC_STRING);
	exit(1);
    }

    // get are input data structs
    graph = getGraph(argv[1]);
    if (TREE_FORMAT == PARRAY)
    {
	pArray* parray = binReadPArray(argv[2]);
	tree = pArray2tree(parray, graph);
 	freePArray(parray);
    }
    else
	tree = getBinIJVTree(argv[2]);
    
    VirtGraph* vg = newVirtGraph(graph->n);
    int ctr = 0;
    improvedTree = improveTree(tree, graph, vg, tree);
    ctr++;
    while (ctr++ < NUM_ITERATIONS)
    {
	freeVirtGraph(vg);
	vg = newVirtGraph(graph->n);
	Tree* tmp = improveTree(improvedTree, graph, vg, tree);
	freeTree(improvedTree);
	improvedTree = tmp;
    }

    // output
    if (TREE_FORMAT == PARRAY)
    {
	// we need verts, which imrovedTree doesn't have (old tree now has internals of new tree)
	binWriteTree2pArray(tree, argv[3]);
    }
    else
	binWriteTree2IJV(improvedTree, argv[3]);

    //free mem
    freeGraph(graph);
    freeTree(tree);
    freeTree(improvedTree);
    freeVirtGraph(vg);
    // tic; ...; toc
    printf("exiting at %lf seconds\n", (double)clock()/(double)CLOCKS_PER_SEC);

    return 0;
}
