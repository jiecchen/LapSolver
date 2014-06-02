/*******************
   treedist.c

   compile with mex treedist.c

*******************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h> 
#include <sys/stat.h>
#include <sys/times.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <dirent.h>

#include "st_defs.h"

/* #include "datatype.h" */

#define DEBUG_FLAG 0

/********************************************************
 *   function headers
 *******************************************************/

void freeGraph(myGraph *G);
void printGraph(myGraph *G);




/********************************************************
  *   treeDist part
  *******************************************************/

static struct TreeNode {
    int name;
    double load;
    double cost;
    int numberOfChildren;
    struct TreeNode * parent;
    struct TreeNode * children;
    double *w;
};

static struct TreeWorkSpace{
  struct TreeNode * first_treeNode;
  double * first_w;
  struct TreeNode *current_treeNode;
  double *current_w;
};



static struct Edge {
   int e1;
   int e2;
   double w;
};

static struct GraphByTwoArrays {
   int n;
   int m;
   int *E1;
   int *E2;
   double *w;
};


static struct GraphByOneArray {
   int n;
   int m;
   struct Edge *E;
};




/* The following function is Needed for calling
  *   qsort on edges by the first endpoint
  */

static int compareEdgeByFirst(a,b)
char *a, *b;
{
   struct Edge * e1 = ( struct Edge *) a;
   struct Edge * e2 = ( struct Edge *) b;
   return (e1->e1 - e2->e1);
}


static two2OneArray (struct GraphByTwoArrays *G1, struct GraphByOneArray *G2){
   int k;
   G2->n = G1->n;
   G2->m = G1->m;
   G2->E  = (struct Edge *) cCalloc(G2->m+1,sizeof(struct Edge), "G2->E in two2OneArray");

   for (k = 0; k < G2->m; k++){
     G2->E[k].e1 = G1->E1[k];
     G2->E[k].e2 = G1->E2[k];
     G2->E[k].w = G1->w[k];
   }
}

static int RootTreeFromOneArray(struct GraphByOneArray *G, struct TreeNode *root, 
			    int chosenRoot, struct TreeWorkSpace *workspace,
			    struct TreeNode ** BFS_Ordering){
   int i,j,k, n, m, counts;
   int *vertexIndex, *nbrcounts, *covered,  qhead, qtail, current,qlocation, *currentnbrs;
   struct Edge *symE;
   struct TreeNode *node, **queue, nodeTemp;
   int bfs_count = 0;

   struct tms t1 ;
   struct tms t2 ;

   n = G->n;
   m = G->m;
   if (m/2 != n - 1) {
     fprintf(stderr,"It is not a tree \n");
     return 1;
   }
   else {
     symE = G->E;
   }

   vertexIndex = (int *) cCalloc(n,sizeof(int),"vertexIndex in RootTreeFromOneArray");
   nbrcounts = (int *) cCalloc(n,sizeof(int),"nbrcounts in RootTreeFromOneArray");
   covered = (int *) cCalloc(n,sizeof(int),"covered in RootTreeFromOneArray");
   currentnbrs = (int *) cCalloc(n,sizeof(int),"currentnbrs in RootTreeFromOneArray");
   queue = (struct TreeNode **) cCalloc(n,sizeof(struct TreeNode *),"queue in RootTreeFromOneArray");

   for (i = 0; i< n; i++){ covered[i] = 0; }
   k = 0;
   counts = 0;
   vertexIndex[k] = 0;
   for (i = 0; i< m; i++){
     if (k == symE[i].e1) counts++;
     else {
	nbrcounts[k] = counts;
        counts = 0;
        k = k+1;
	vertexIndex[k] = i;
        i = i-1;
     }
     nbrcounts[k] = counts;
   }

/*   root the tree with chosenRoot as root   */
     current = chosenRoot;
     qhead = 0;
     qtail =0;
     root->name = current;
     root->parent = NULL;
     queue[qtail++] = root;
     covered[current] = 1;
     while (qtail - qhead > 0){
        node = queue[qhead++];
/* Building the BFS Ordering */
        BFS_Ordering[bfs_count++] = node;
        current = node->name;
        node->load = 0;
        node->cost = 0;
        qlocation = qtail;
        node->numberOfChildren = 0;

        k = 0;
        for (i = 0; i< nbrcounts[current];i++){
          if (covered[symE[vertexIndex[current]+i].e2] == 0){
            currentnbrs[k++] = vertexIndex[current]+i;
            node->numberOfChildren++;
            covered[symE[vertexIndex[current]+i].e2] = 1;
	 }
        }

/* Allocate space for children and edge weights */

        node->children = workspace->current_treeNode;  
        workspace->current_treeNode+=node->numberOfChildren;

        node->w = workspace->current_w;  
        workspace->current_w+=node->numberOfChildren;



/* Assign children and edge weights */

        for (i = 0; i< node->numberOfChildren; i++){
          node->children[i].name = symE[currentnbrs[i]].e2;
          node->children[i].parent = node;
          node->w[i] = symE[currentnbrs[i]].w;
          queue[qtail++] = &node->children[i];
        }

     }

/* Free the temp space allocated */
     cFree(vertexIndex);
     cFree(nbrcounts);
     cFree(covered);
     cFree(currentnbrs);
     cFree(queue);


     return 0;
}

static printTreeBFS(int n, struct TreeNode **BFS_Ordering){
   int numberChildren;
   int i, j;
   struct TreeNode *node, *current;
   for (i=0; i<n;i++){
     node = BFS_Ordering[i];
     if (node->parent != NULL){
       fprintf (stderr,"name,parent,load,cost,numberOfChildren= %d, %d,%f, %f,%d \n",
              node->name, node->parent->name, node->load,
	     node->cost,node->numberOfChildren);
     }
     else{
       fprintf(stderr,"name, parent,load,cost,numberOfChildren= %d, NULL, %f, %f,%d \n",
              node->name, node->load, node->cost,node->numberOfChildren);
     }
     numberChildren = node->numberOfChildren;
     current = node;
     node = node->children;
     for (j=0; j<numberChildren;j++){
       fprintf(stderr,"weight of (%d,%d) is %f \n",current->name,node->name,current->w[j]);
       node++;
     }
   }
}

static ExtractCostBFS(int n, struct TreeNode **BFS_Ordering, double *cost){
   int j;
   struct TreeNode *node;
   for (j=0; j< n; j++){
     node = BFS_Ordering[j];
     if (node != NULL){
       cost[node->name] = node->cost;
     }
   }
}

static AssignLoadBFS(int n,struct TreeNode **BFS_Ordering, double *load){
   int i;
   struct TreeNode *node;
   for (i=0;i < n; i++){
     node= BFS_Ordering[i];
     node->load = load[node->name];
   }
}

static bottomUpBFS(int n, struct TreeNode **BFS_Ordering)
{
   int numberChildren;
   struct TreeNode *node, *child;
   int i,j;

   for (j=n-1; j> -1; j--){
     node = BFS_Ordering[j];
     if (node->numberOfChildren != 0){
       numberChildren = node->numberOfChildren;
       child = node->children;
       for (i=0; i<numberChildren;i++){
         node->load += child->load;
         node->cost += child->cost + child->load*node->w[i];
         child++;
       }
     }
   }
}

static topDownBFS(int n,struct TreeNode **BFS_Ordering)
{
   int numberChildren;
   struct TreeNode *child, *node;
   int i,j;
   double downwardLoad, downwardCost;
   for (j=0; j<n; j++){
     node = BFS_Ordering[j];
     if (node->numberOfChildren != 0){
       numberChildren = node->numberOfChildren;
       child = node->children;
       for (i=0; i<numberChildren;i++){
         downwardLoad = node->load - child->load;
         downwardCost = node->cost - node->w[i]*child->load - child->cost;
         child->load += downwardLoad;
         child->cost += downwardCost + downwardLoad*node->w[i];
         child++;
       }
     }
   }
}



/**
  *
  *In treedist, it uses the assumption that (E1,E2) are expressing
  *   each edge twice, one for each direction. In addition, E1 is
  *   sorted from small to large. Also in the input, m is equal to
  *   2 times of the number of the edge in the tree.
  */
int st_treedist(int n, int m, int *E1, int *E2,
		double *w, double *load, double *cost)
{
   struct TreeNode *root, *node;
   struct GraphByTwoArrays G2;
   struct GraphByOneArray G1;
   struct TreeWorkSpace  workspace;
   struct TreeNode ** BFS_Ordering;

   
   int i;

   G2.n = n;
   G2.m = m;
   G2.E1 = E1;
   G2.E2 = E2;
   G2.w =  w;

/*   Allocate the workspace so that it can be freed easily */

   workspace.first_treeNode = (struct TreeNode *) cCalloc(n,sizeof(struct TreeNode),"workspace,first_treeNode in st_treeDist");
   workspace.first_w = (double *) cCalloc(n,sizeof(double),"workspace.first_w in st_treeDist");
   workspace.current_treeNode = workspace.first_treeNode;
   workspace.current_w = workspace.first_w;


/* Convert from two array representation to one array representation */
   two2OneArray ( &G2, &G1);

/* Allocate space for the root */
   root = (struct TreeNode *) cCalloc(1,sizeof(struct TreeNode),"root in st_treeDist");

/* Build the rooted tree with a specfied root */
   int err;
   
   BFS_Ordering = (struct TreeNode **) cCalloc(n, sizeof(struct TreeNode *),"BFS_Ordering in st_treedist");
   err = RootTreeFromOneArray(&G1, root,0,&workspace, BFS_Ordering);

/*    printTreeBFS(n, BFS_Ordering); */
   if (err == 1) {
     cFree(G1.E);
     cFree(workspace.first_treeNode);
     cFree(workspace.first_w);
     cFree(root);
     cFree(BFS_Ordering);
     return err;
   }
   

   /* Assign load to tree vertices */
 
   AssignLoadBFS(n,BFS_Ordering,load);
   bottomUpBFS(n,BFS_Ordering); 
   topDownBFS(n,BFS_Ordering); 
   ExtractCostBFS(n, BFS_Ordering,cost); 

   
   /* Free memory allocated  */
  
   cFree(G1.E);
   cFree(workspace.first_treeNode);
   cFree(workspace.first_w);
   cFree(root);
   cFree(BFS_Ordering);
   return 0;

};


/* two to one array
   children
   w
 */

