/*******************
    st_part.c
 
 @file handles trees, partPre,
       and whatever
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
#include <sys/resource.h>

#include "st_defs.h"


/* #define ST_DEBUG 1 */

/********************************************************
  *   function headers
  *******************************************************
  */

ijvType * vertexWise_Sampling(myGraph *G, int min_deg, float minfrac_deg, float minfrac_wt);
ijvType *compressIJV(ijvType *ijvIn);
int *copyOpts(int *opts);
myGraph *ijv2graph(ijvType *ijv);

ijvType *multi_Prim (myGraph *G, int k);

int *partTreeOpts();

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



#define XTIMER taucs_sttime


/*DAS*/
double taucs_sttime()
{
  struct rusage a;
  
  getrusage(RUSAGE_SELF,&a);

  return (double) 
    (double) ((a.ru_utime).tv_sec )+
    (double) ((a.ru_utime).tv_usec )* 1.0e-6;
  /*#endif*/
}


double wtime_total;
double wtime_bridges; 
double wtime_offdiags; 
double wtime_merge; 
double wtime_metis; 



ijvType *st_part(ijvType *ijv, int *opts)
{ 

  ijvType *ijvOut;

  int *comp;

  ijvBlockDiag *bd;

  struct tms t1 ;
  struct tms t2 ;


  if (!(opts[stOptNumTree] == stOptValTreeYes) && (opts[stOptNumTime2])) {
    wtime_bridges = XTIMER();
  }

  
  /* debug stuff */
  if (opts[stOptNumDebug]) {
    fprintf(stderr,"debug, n: %d, nnz: %d\n",ijv->n,ijv->nnz);
  }
  
  /* if it is a tree, return original */
  if (ijv->nnz == ((ijv->n)-1)) {
    return copyIJV(ijv);
  }


  /* DAS: 10/30/04: alter the balance) */
  int saven;
  saven = ijv->n;

  ijv->n = (int) (ijv->n * metisBalance);
  
  if (opts[stOptNumTime]) {
    printf("metis start\n");

    wtime_metis =  XTIMER();
  }

  comp = (int *) cCalloc(ijv->n, sizeof(int),"comp in mex_metis");
  if (opts[stOptNumTree] == stOptValTreeYes)
    st_metis(ijv, 2, comp, opts[stOptNumMetis]);
  else
    st_metis(ijv, opts[stOptNumParts], comp, opts[stOptNumMetis]);

  (void) times (&t2) ;
  if (opts[stOptNumTime]) {

    wtime_metis = XTIMER() - wtime_metis;
    
    printf("metis time: %10.3f seconds\n",wtime_metis);
    opts[stOptNumTime] = 0;
  }


  ijv->n = saven;

  
  /* modifies comp */
  bd = st_rePart(ijv, comp);

  int numParts;
  numParts = bd->part->numParts;


  /* for each part, recursively compute a tree */
  ijvType **trees;
  trees = (ijvType **) cCalloc(numParts, sizeof(ijvType *), "trees in st_part");

  int i;
  int *treeOpts;
  if (opts[stOptNumTree] == stOptValTreeYes)
    treeOpts = opts;
  else {
    treeOpts = copyOpts(opts);
    treeOpts[stOptNumTree] = stOptValTreeYes;
  }

  
  for (i = 0; i < numParts; i++) {
    if (bd->part->partSizes[i] > 0) {
      trees[i] = st_part(bd->blocks[i],treeOpts);
    } else {
      trees[i] = NULL;
    }
  }

  if (opts[stOptNumTree] != stOptValTreeYes)
    cFree(treeOpts);

  
  ijvType *meta, *metaSparse;
  meta = st_metaGraph(ijv, bd->part->compVec);

  int numTrees;
  numTrees = numParts;


  /* handling the meta graph */
  ijvType *metaSparse1, *metaSparse2;
  metaSparse1 = NULL;
  metaSparse2 = NULL;



  if (opts[stOptNumTree] == stOptValTreeYes) {
    metaSparse = st_part(meta,opts);
    freeIJV(meta);
  } else {

    if (opts[stOptNumSpars] == stOptValSparsNo)
      metaSparse = meta;
    else {
      if (opts[stOptNumSpars] == stOptValSparsPrim) {

	if (opts[stOptNumDebug] || opts[stOptNumData]) {
	  fprintf(stderr,"multi-Prim on meta, n: %d, nnz: %d\n",meta->n,meta->nnz);
	}

	
	myGraph *G;
	G = ijv2graph(meta);
	metaSparse1 =  multi_Prim (G, opts[stOptNumSparsParam]);
	freeGraph(G);

	if (opts[stOptNumDebug] || opts[stOptNumData]) {
	  fprintf(stderr,"Multi-Prim output: metaSparse1, n: %d, nnz: %d\n",metaSparse1->n,metaSparse1->nnz);
	}

      }

      if (opts[stOptNumSparsMinDeg] ||
	  opts[stOptNumSparsMinFracDeg] ||
	  opts[stOptNumSparsMinFracWt]) {

	if (opts[stOptNumDebug] || opts[stOptNumData]) {
	  fprintf(stderr,"Min Degrees on meta, n: %d, nnz: %d\n",meta->n,meta->nnz);
	}

	
	myGraph *G;
	G = ijv2graph(meta);

	metaSparse2 =
	  vertexWise_Sampling(G,opts[stOptNumSparsMinDeg], sparsMinFracDeg, sparsMinFracWt);
	freeGraph(G);

	if (opts[stOptNumDebug] || opts[stOptNumData]) {
	  fprintf(stderr,"Min Degrees output: metaSparse2, n: %d, nnz: %d\n",metaSparse2->n,metaSparse2->nnz);
	}	

	
      }

      if (!metaSparse2)
	metaSparse = metaSparse1;
      else
	if (!metaSparse1) {
	  
	  metaSparse = compressIJV(metaSparse2);
	  freeIJV(metaSparse2);
	  
	  if (opts[stOptNumDebug] || opts[stOptNumData]) {
	    fprintf(stderr,"Min Degrees output: metaSparse, n: %d, nnz: %d\n",metaSparse->n,metaSparse->nnz);
	}	
  
	}      

        else {
	  /* both sparse1 and sparse2 exist */
	  ijvType *both[2];
	  both[0] = metaSparse1;
	  both[1] = metaSparse2;

	  ijvType *metaSparse12;
	  metaSparse12 = mergeIJVs(meta->n, both, 2);
	  freeIJV(metaSparse1);
	  freeIJV(metaSparse2);

	  metaSparse = compressIJV(metaSparse12);
	  freeIJV(metaSparse12);


	  if (opts[stOptNumDebug] || opts[stOptNumData]) {
	    fprintf(stderr,"Combined sparsifier: metaSparse, n: %d, nnz: %d\n",metaSparse->n,metaSparse->nnz);
	  }

	}

    }
  }
  
  ijvType **ijvArray;


  if (!(opts[stOptNumTree] == stOptValTreeYes) && (opts[stOptNumTime2])) {
    wtime_offdiags = XTIMER();
  }
  

  ijvArray = st_offDiags(ijv, bd->part, metaSparse);


  if (!(opts[stOptNumTree] == stOptValTreeYes) && (opts[stOptNumTime2])) {
    wtime_offdiags = XTIMER() - wtime_offdiags;

    printf("offdiags time: %10.3f seconds\n",wtime_offdiags);

  }

    
  
  if (!(opts[stOptNumTree] == stOptValTreeYes) && (opts[stOptNumTime2])) {
    wtime_bridges = XTIMER();
  }

  ijvType *ijvBridges;
  ijvBridges = 
    buildBridges(ijv->n, bd->part,
		 ijvArray, 
		 metaSparse, 
		 trees);

  if (!(opts[stOptNumTree] == stOptValTreeYes) && (opts[stOptNumTime2])) {
    wtime_bridges = XTIMER() - wtime_bridges;

    printf("bridges time: %10.3f seconds\n",wtime_bridges);
  }
  
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



  if (!(opts[stOptNumTree] == stOptValTreeYes) && (opts[stOptNumTime2])) {
    wtime_merge = XTIMER();
  }

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

  /* issue is that some trees may be null,
     so we don't put those on this list */
  int ind;
  int numGoodTrees;
  ind = 0;
  for (i = 0; i < numTrees; i++) {
    if (bd->part->partSizes[i] > 0) {
      remapIJV(trees[i], bd->part, i);
      treesAndBridges[ind++] = trees[i];
    }
  }
  numGoodTrees = ind;
  
  treesAndBridges[numGoodTrees] = ijvBridges;

  ijvOut = mergeIJVs(ijv->n, treesAndBridges, numGoodTrees+1);
  /*  ijvSorti(ijvOut); */
  ijvSortij(ijvOut);


  if (!(opts[stOptNumTree] == stOptValTreeYes) && (opts[stOptNumTime2])) {
    wtime_merge = XTIMER() - wtime_merge;

    printf("merge time: %10.3f seconds\n",wtime_merge);
  }
  
  freeIjvBlockDiagAndPart(bd);

  for (i = 0; i < numTrees; i++) {
    if (trees[i]) freeIJV(trees[i]);
  }
  cFree(trees);

  for (i = 0; i < metaSparse->nnz; i++)
    freeIJV(ijvArray[i]);
  cFree(ijvArray);

  freeIJV(metaSparse);
  
  freeIJV(ijvBridges);

  cFree(treesAndBridges);

  if (!(opts[stOptNumTree] == stOptValTreeYes) && (opts[stOptNumTime2])) {
    wtime_total = XTIMER() - wtime_total;

    printf("total time: %10.3f seconds\n",wtime_total);
  }

  
  return ijvOut;
  
}

