#ifndef CALCSTRETCHCTREE_H
#define CALCSTRETCHCTREE_H
// module that does all the actual work for calculating stretch
// calcstretch also has its own tree datatypes that are incompatible with treeimproves
#include "graph.h"
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
    TreeNode** nodes;        // list of all nodes with parents always before children
    TreeNode* raw_nodes;     // where nodes are actually stored
    TreeNode** verts;        // table of where nodes corresponding to verts are
} Tree;

// parray determines structure of tree
// g determines distances on tree
Tree* pArray2tree(pArray* parray, myGraph* g);
pArray* tree2pArray(Tree* tree);
void freeTreeNodeInternals(TreeNode* node);
void freeTree(Tree* tree);
double calcStretch(Tree* tree, const myGraph* h);
#endif // CALCSTRETCHCTREE_H
