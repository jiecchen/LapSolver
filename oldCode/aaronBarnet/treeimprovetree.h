// common tree structs and functions for treeimprove and treeimprove2
#ifndef TREEIMPROVETREE_H
#define TREEIMPROVETREE_H
#include "graph.h"

static int COMP;                  // remakeTree first "COMP++" and then labels all nodes as COMP
typedef struct TreeNode_s
{
    int v;                 // node
    struct TreeNode_s* p;  // parent node
    double weight;         // weight of edge between this node and parent
    double topDist;        // sum of path lengths to every node above this node
    double topVol;            // how many nodes are in tree connected by parents edge?
    int num_children;
    struct TreeNode_s** children;
    double totChildDist;   // sum of path lengths to every node below this node
    double totChildVol;       // number of nodes below this node
    int c;                 // which component it's in
} TreeNode;


typedef struct Tree_s
{
    int n;                    // num nodes
    TreeNode* root;          
    TreeNode** nodes;         // all the nodes ordered so that parents are always before children
    TreeNode* raw_nodes;      // NULL for subtree -- where all the nodes are actually stored
    TreeNode** verts;         // NULL for subtree -- table of where nodes corresponding to verts are
} Tree;

void freeTreeNodeInternals(TreeNode* node);
Tree* graph2tree(myGraph* g);
void freeTree(Tree* tree);
void binWriteTree2IJV(Tree* tree, const char* fileName);
void binWriteTree2pArray(Tree* tree, const char* fileName);
Tree* pArray2tree(pArray* parray, myGraph* g);
#endif // TREEIMPROVETREE_H
