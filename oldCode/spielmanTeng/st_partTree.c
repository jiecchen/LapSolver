/*******************
    st_partTree.c

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


/* #define ST_DEBUG 1 */

/********************************************************
  *   function headers
  *******************************************************
  */


ijvType *mergeIJVs(int n, /** the size of the output graph */
		   ijvType **ijvs, /** the array of ijvs */
		   int numIJVs);

ijvType *buildBridges(int n, /**< the number of verts in graph */
		      partType *p, /**< partition, needed to figure
				      edges */
		 ijvType **ijvArray, /**< the off-diag blocks */
		 ijvType *meta, /**< the graph on components */
		      ijvType **trees);

ijvType **st_offDiags(ijvType *ijv, partType *p, ijvType *meta);
ijvType *st_metaGraph(ijvType *ijv, int *compVec);

ijvBlockDiag *st_rePart(ijvType *ijv, int *compVec);

ijvType *copyIJV(ijvType *ijv_in);
void freeGraph(myGraph *G);
void printGraph(myGraph *G);
int * CountingSort(int * A, int n, int k);
ijvType *makeIJV(int nnz);


/* It maps each vertex of the original graph to
 *  its component vertex. It also update the n entry.
 */



ijvType *st_partTree(ijvType *ijv, int *opts)
{ 

  ijvType *ijvOut;

  int *comp;

  ijvBlockDiag *bd;
  
  /* if it is a tree, return original */
  if (ijv->nnz == ((ijv->n)-1)) {
    return copyIJV(ijv);
  }

  comp = (int *) cCalloc(ijv->n, sizeof(int),"comp in mex_metis");
  st_metis(ijv, 2, comp);

  /* modifies comp */
  bd = st_rePart(ijv, comp);

  int numParts;
  numParts = bd->part->numParts;


  /* for each part, recursively compute a tree */
  ijvType **trees;
  trees = (ijvType **) cCalloc(numParts, sizeof(ijvType *), "trees in st_partTree");

  int i;
  for (i = 0; i < numParts; i++) {
    trees[i] = st_partTree(bd->blocks[i]);
  }
    
  ijvType *meta, *metaTree;
  meta = st_metaGraph(ijv, bd->part->compVec);

  int numTrees;
  numTrees = numParts;

  metaTree = st_partTree(meta);
  freeIJV(meta);



  ijvType **ijvArray;
  ijvArray = st_offDiags(ijv, bd->part, metaTree);

  
  ijvType *ijvBridges;
  ijvBridges = 
    buildBridges(ijv->n, bd->part,
		 ijvArray, 
		 metaTree, 
		 trees);

  if (!ijvBridges) {
    fprintf(stderr,"error in st_partTree\n");
    printIJV(ijv);
    fprintf(stderr,"\n COMP: \n");
    printIntVec(bd->part->compVec,ijv->n,"compnew");
    printIntVec(comp,ijv->n,"compOrig");
    for (i = 0; i < numParts; i++) {
      fprintf(stderr,"block %d\n",i);
      printIJV(bd->blocks[i]);
      fprintf(stderr,"\ntree %d\n",i);
      printIJV(trees[i]);
    }
    return NULL;
  }


  cFree(comp);



  /* now, merge the trees with the bridges,
     first, allocating space for the new array of
     pointers, and then performing the merge

     before we merge, we need to re-map the vertex
     labels in the trees
  */
  ijvType **treesAndBridges;
  treesAndBridges = (ijvType **)
    cCalloc(numTrees+1,
	    sizeof(ijvType *),
	    "treesAndBridges in mex_partTreeB");

  for (i = 0; i < numTrees; i++) {
    remapIJV(trees[i], bd->part, i);
    treesAndBridges[i] = trees[i];
  }

  treesAndBridges[numTrees] = ijvBridges;

  ijvOut = mergeIJVs(ijv->n, treesAndBridges, numTrees+1);
  /*  ijvSorti(ijvOut); */
  ijvSortij(ijvOut);

  
  freeIjvBlockDiagAndPart(bd);

  for (i = 0; i < numTrees; i++) {
    freeIJV(trees[i]);
  }
  cFree(trees);

  for (i = 0; i < metaTree->nnz; i++)
    freeIJV(ijvArray[i]);
  cFree(ijvArray);

  freeIJV(metaTree);
  
  freeIJV(ijvBridges);

  cFree(treesAndBridges);

  return ijvOut;
  
}

