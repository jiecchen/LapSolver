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
   G2->E  = (struct Edge *) cMalloc(G2->m*sizeof(struct Edge));
   for (k = 0; k < G2->m; k++){
     G2->E[k].e1 = G1->E1[k];
     G2->E[k].e2 = G1->E2[k];
     G2->E[k].w = G1->w[k];
   }
}


static RootTreeFromOneArray(struct GraphByOneArray *G, struct TreeNode *root, 
			    int chosenRoot, struct TreeWorkSpace *workspace){
   int i,j,k, n, m, counts;
   int *vertexIndex, *nbrcounts, *covered,  qhead, qtail, current,qlocation, *currentnbrs;
   struct Edge *symE;
   struct TreeNode *node, **queue, nodeTemp;

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

if (DEBUG_FLAG == 1) {
   for (i =0; i< m; i++){
      printf("symE [%d] is [%d,%d] \n", i, symE[i].e1,symE[i].e2);
   }
}


     vertexIndex = (int *) cMalloc(n*sizeof(int));
     nbrcounts = (int *) cMalloc(n*sizeof(int));
     covered = (int *) cMalloc(n*sizeof(int));
     currentnbrs = (int *) cMalloc(n*sizeof(int));
     queue = (struct TreeNode **) cMalloc(n*sizeof(struct TreeNode *));
     for (i = 0; i< n; i++){
       covered[i] = 0;
     }
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

if (DEBUG_FLAG == 1) {
   for (i =0; i< n; i++){
     printf("nbrcount and index of %d is %d, %d \n", i, 
nbrcounts[i],vertexIndex[i]);
   }

   for (i =0; i< n; i++){
     printf("number of nbrs of %d is %d \n", i, nbrcounts[i]);
   }
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

if (DEBUG_FLAG == 1) {
    printf("the queue is \n");
    for (j=qhead;j< qtail; j++){
      printf("%d \n", queue[j]->name);
    }
}
        node = queue[qhead++];
        current = node->name;
if (DEBUG_FLAG == 1) {
    printf("current = %d \n", current);
}
        node->load = 0;
        node->cost = 0;
        qlocation = qtail;
        node->numberOfChildren = 0;
if (DEBUG_FLAG == 1) {
     printf("number of nbrs = %d \n", nbrcounts[current]);
     for (j =0; j< n; j++){
       printf("cover of %d is %d \n", j, covered[j]);
     }
}

        k = 0;
        for (i = 0; i< nbrcounts[current];i++){
if (DEBUG_FLAG == 1) {
    printf("process vertex %d \n", symE[vertexIndex[current]+i].e2);
}
          if (covered[symE[vertexIndex[current]+i].e2] == 0){
/*            currentnbrs[k++] = symE[vertexIndex[current]+i].e2; */
            currentnbrs[k++] = vertexIndex[current]+i;
            node->numberOfChildren++;
            covered[symE[vertexIndex[current]+i].e2] = 1;
	 }
        }
if (DEBUG_FLAG == 1) {
    printf("number of Children = %d \n", node->numberOfChildren);
}


/* Allocate space for children and edge weights */

/*         node->children = (struct TreeNode *)  */
/*                cMalloc(node->numberOfChildren*sizeof(struct TreeNode)); */
/*         node->w = (double *) cMalloc(node->numberOfChildren*sizeof(double)); */

        node->children = workspace->current_treeNode;  
        workspace->current_treeNode+=node->numberOfChildren;

        node->w = workspace->current_w;  
        workspace->current_w+=node->numberOfChildren;



/* Assign children and edge weights */

        for (i = 0; i< node->numberOfChildren; i++){
/*          node->children[i].name = currentnbrs[i]; */
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

static printTree(struct TreeNode *node){
   int numberChildren;
   int i;
   struct TreeNode *current;
   if (node != NULL){
     if (node->parent != NULL){
       printf("name,parent,load,cost,numberOfChildren= %d, %d,%f, %f,%d \n",
              node->name, node->parent->name, node->load,
	     node->cost,node->numberOfChildren);
     }
     else{
       printf("name, parent,load,cost,numberOfChildren= %d, NULL, %f, %f,%d \n",
              node->name, node->load, node->cost,node->numberOfChildren);
     }
     numberChildren = node->numberOfChildren;
     current = node;
     node = node->children;
     for (i=0; i<numberChildren;i++){
       printf("weight of (%d,%d) is %f \n",current->name,node->name,current->w[i]);
       printTree(node);
       node++;
     }
   }

}

static ExtractCost(struct TreeNode *node, double *cost){
   int numberChildren;
   int i;
   struct TreeNode *current;
   if (node != NULL){
     cost[node->name] = node->cost;
     /*    mexPrintf("%d -> %g\n",node->name, node->cost); */
     numberChildren = node->numberOfChildren;
     node = node->children;
     for (i=0; i<numberChildren;i++){
       ExtractCost(node,cost);
       node++;
     }
   }
}


static AssignLoad(struct TreeNode *node, double *load){
   int numberChildren;
   int i;
   struct TreeNode *current;
   if (node != NULL){
     node->load = load[node->name];
     numberChildren = node->numberOfChildren;
     current = node;
     node = node->children;
     for (i=0; i<numberChildren;i++){
       AssignLoad(node,load);
       node++;
     }
   }
}


static bottomUp(struct TreeNode *node)
{
   int numberChildren;
   struct TreeNode *child;
   int i;

   if (node->numberOfChildren != 0){
     numberChildren = node->numberOfChildren;
     child = node->children;
     for (i=0; i<numberChildren;i++){
       bottomUp(child);
       node->load += child->load;
       node->cost += child->cost + child->load*node->w[i];
       child++;
     }
   }
}

static topDown(struct TreeNode *node)
{
   int numberChildren;
   struct TreeNode *child;
   int i;
   double downwardLoad, downwardCost;

   if (node->numberOfChildren != 0){
     numberChildren = node->numberOfChildren;
     child = node->children;
     for (i=0; i<numberChildren;i++){
       downwardLoad = node->load - child->load;
       downwardCost = node->cost - node->w[i]*child->load - child->cost;
       child->load += downwardLoad;
       child->cost += downwardCost + downwardLoad*node->w[i];
       topDown(child);
       child++;
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

   int i;

   G2.n = n;
   G2.m = m;
   G2.E1 = E1;
   G2.E2 = E2;
   G2.w =  w;

/*   Allocate the workspace so that it can be freed easily */
     workspace.first_treeNode = (struct TreeNode *) cMalloc(n*sizeof(struct TreeNode));
     workspace.first_w = (double *) cMalloc(n*sizeof(double));
     workspace.current_treeNode = workspace.first_treeNode;
     workspace.current_w = workspace.first_w;



/* Convert from two array representation to one array representation */
   two2OneArray ( &G2, &G1);

/* Allocate space for the root */
   root = (struct TreeNode *) cMalloc(sizeof(struct TreeNode));

/* Build the rooted tree with a specfied root */
   int err;
   
   err = RootTreeFromOneArray(&G1, root,0,&workspace);

   if (err == 1) {
     cFree(G1.E);
     cFree(workspace.first_treeNode);
     cFree(workspace.first_w);
     cFree(root);
     return err;
   }
   
   /*  printTree(root); */
   /*  printf("\n"); */
   
   /* Assign load to tree vertices */
   AssignLoad(root,load);
   
   bottomUp(root);
   topDown(root);
   ExtractCost(root,cost);
   /*   printTree(root);  */
   
   /* Free memory allocated  */
   

   cFree(G1.E);
   cFree(workspace.first_treeNode);
   cFree(workspace.first_w);
   cFree(root);

};


/* two to one array
   children
   w
 */

