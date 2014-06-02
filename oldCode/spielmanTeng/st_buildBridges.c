/**
 * \file st_buildBridges.c
 *
 * \brief choose which edges will go between
 * blocks in the parition
 *
 *
 *******************/


#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>

#include "st_defs.h"


#define min(A,B) ((A) < (B) ? (A) : (B))
#define max(A,B) ((A) > (B) ? (A) : (B))

/********************************************************
 *   function headers
 *******************************************************/

dubIjvType *doubleIJV(ijvType *ijvIn);
void *cCalloc(size_t n, size_t size, char *errMsg);
int st_treedist(int n, int m, int *E1, int *E2,
		double *w, double *load, double *cost);
ijvType *makeIJV(int nnz);




/**
 *  compLoad
 *  \brief computes the loads
 *
 *  the inputs are the ijv block, and the number
 *   of items on each side 
 *  (which really should be in the block)
 */
void static compLoad(const ijvType *ijv, const int size1, const int size2,
		double *load1, double *load2)
{

  int x;

  for (x = 0; x < ijv->nnz; x++) {
    load1[ijv->i[x]] += ijv->v[x];
    load2[ijv->j[x]] += ijv->v[x];
  }
  
}


/**
 *  st_buildBridge
 *  \brief builds one bridge
 *
 * the output is the index of the bridge edge
 *    in the input ijv.
 *
 * in case of error, returns -1
 *
 * \todo go more directly to treeDist input form
 */
static int buildBridge(ijvType *tree1,
		       ijvType *tree2,
		       ijvType *ijv)
{

  int i;

  int bridgeIndex; /* the output */
  
  double *load1, *load2;
  double *cost1, *cost2;

  dubIjvType *dubTree1, *dubTree2;

  dubTree1 = doubleIJV(tree1);
  /* claim is that st_treeDist only needs first sorted */
  /*  ijvSortj((ijvType *) dubTree1); */
  ijvSorti((ijvType *) dubTree1);


  dubTree2 = doubleIJV(tree2);
  /*  ijvSortj(dubTree2);   */
  ijvSorti(dubTree2);

  load1 = (double *) cCalloc(tree1->n, sizeof(double),
			     "load1 in buildBridge");
  load2 = (double *) cCalloc(tree2->n, sizeof(double),
			     "load2 in buildBridge");

  cost1 = (double *) cCalloc(tree1->n, sizeof(double),
			     "cost1 in buildBridge");
  cost2 = (double *) cCalloc(tree2->n, sizeof(double),
			     "cost2 in buildBridge");

  compLoad(ijv, tree1->n, tree2->n,
	   load1, load2);

  int err;
  
  /* treedist can make errors on nnz = 0,
     so protect for now */
  if (dubTree1->nnz == 0) {
    for (i = 0; i < tree1->n; i++)
      cost1[i] = load1[i];

  } else {
    
    err =  st_treedist(dubTree1->n, dubTree1->nnz, dubTree1->i, dubTree1->j,
		       dubTree1->v, load1, cost1);
    if (err == 1) {
      fprintf(stderr,"error in buildBridge\n");
      printIJV(dubTree1);
      return -1;
    }
  }
  freeIJV(dubTree1);

  if (dubTree2->nnz == 0) {
    for (i = 0; i < tree2->n; i++)
      cost2[i] = load2[i];

  } else {
    err =  st_treedist(dubTree2->n, dubTree2->nnz, dubTree2->i, dubTree2->j,
		       dubTree2->v, load2, cost2);
    if (err == 1) {
      fprintf(stderr,"error in buildBridge\n");
      printIJV(dubTree2);
      return -1;
    }
  }
  freeIJV(dubTree2);
    
  /*-------------------------------

     now, decide which edge should
     be the bridge

     --------------------------------*/

  bridgeIndex = -1;


  double bridgeCost = DBL_MAX;

  int x;
  for (x = 0; x < ijv->nnz; x++) {
    double cost = (cost1[ijv->i[x]] + cost2[ijv->j[x]])/(ijv->v[x]);
    if (cost < bridgeCost) {
      bridgeIndex = x;
      bridgeCost = cost;
    }
  }

  cFree(cost1);
  cFree(cost2);

  cFree(load1);
  cFree(load2);

  return bridgeIndex;
}


/**
 *  st_buildBridges
 *  \brief builds all the bridges
 *
 * note that the xth edge in the ijv list meta
 * corresponds to the xth ijv in ijvArray,
 *   which contains the edges between the componenents
 *   indexed in meta.
 *
 * returns null on error
 */
ijvType *buildBridges(int n, /**< the number of verts in graph */
		      partType *p, /**< partition, needed to figure
				      edges */
		 ijvType **ijvArray, /**< the off-diag blocks */
		 ijvType *meta, /**< the graph on components */
		 ijvType **trees)
{

  int k;

  int x;
  

  ijvType *bridgeEdges; /* the output */
  bridgeEdges = makeIJV(meta->nnz);
  bridgeEdges->n = n;

  /* enum over the meta edges */
  for (x = 0; x < meta->nnz; x++) {
    int ind;
    int comp1;
    int comp2;
    comp1 = meta->i[x];
    comp2 = meta->j[x];

    ind = buildBridge(trees[comp1], trees[comp2],
		      ijvArray[x]);
    if (ind == -1) {
      fprintf(stderr,"error in buildBridges\n");
      return NULL;
    }

    bridgeEdges->i[x] = max(p->partMap[comp1][ijvArray[x]->i[ind]],p->partMap[comp2][ijvArray[x]->j[ind]]);
    bridgeEdges->j[x] = min(p->partMap[comp1][ijvArray[x]->i[ind]],p->partMap[comp2][ijvArray[x]->j[ind]]);
    bridgeEdges->v[x] = ijvArray[x]->v[ind];

  }
  
  return bridgeEdges;
}
