/*******************
   st_offDiags.c

   make graphs of off-diagonal
   blocks according to edges in
   a meta-graph

*******************/


#include <stdio.h>
#include <stdlib.h>

#include "st_defs.h"


/********************************************************
 *   function headers
 *******************************************************/

void *cCalloc(size_t n, size_t size, char *errMsg);
void ijvSorti(ijvType *ijv);
void ijvSortj(ijvType *ijv);
int * CountingSort(int * A, int n, int k);
ijvType *makeIJV(int nnz);



/********************************************************
 *   main routine
 *******************************************************/


/** st_offDiags has inputs:
  ijv: the original graph,
  p  : a partition of that graph
  meta : a meta graph on the blocks of that partition,

  it's outputs are 
  meta : modified so that meta.ex points to a graph in 
  the returned value, ijvArray : the array of the needed off-Diag graphs

  the mx-th item in the ijvArray corresponds to the mx-th edge
    in the ijv list of meta.
  
  Requires that meta
    be sorted with i most significant
  
  FOR NOW: IT SORTS JUST TO BE SURE
*/
ijvType **st_offDiags(ijvType *ijv, partType *p, ijvType *meta)
{
  ijvType *ijvElem;


  int *perm;
  int *perm2;
  int *sortBy;
  

  int n; /* num verts */
  int nnz; /* num nonzeros in ijv */
  int k; /* num parts in partition */

  int x; /* for for loops on ijv */
  int mx; /* for for loops on meta */

  ijvType *tmp_ijv;

  ijvType **ijvArray; /* the output */

  int *counts;
  
  n = ijv->n;
  nnz = ijv->nnz;
  k = p->numParts;


  /* ELIMINATE EVENTUALLY! */
  ijvSortj(meta); 
  ijvSorti(meta); 

  /*------------------------------------------
   *
   * now, sort ijv according to p->compVec,
   *  placing result into tmp_ijv 
   *
   */
  sortBy = (int *) cCalloc(nnz, sizeof(int),
			   "sortBy in st_offDiags");
  tmp_ijv = makeIJV(nnz);
  
  /* first, sort according to j */
  for (x = 0; x < nnz; x++)
    sortBy[x] = min(p->compVec[ijv->j[x]],p->compVec[ijv->i[x]]);

  
  perm = CountingSort(sortBy, nnz, k);

  
  for (x = 0; x < nnz; x++)
    sortBy[perm[x]] = max(p->compVec[ijv->j[x]],p->compVec[ijv->i[x]]);

  perm2 = CountingSort(sortBy, nnz, k);

  tmp_ijv->n = ijv->n;
  for (x = 0; x < nnz; x++) {
    tmp_ijv->i[perm2[perm[x]]] = ijv->i[x];
    tmp_ijv->j[perm2[perm[x]]] = ijv->j[x];
    tmp_ijv->v[perm2[perm[x]]] = ijv->v[x];
  }


  /**********
     DEBUG
   check if all sorted right
  ijvType *tmp;
  tmp = makeIJV(tmp_ijv->nnz);
  for (x = 0; x < tmp->nnz; x++) {
    tmp->i[x] = max((p->compVec[tmp_ijv->i[x]]),(p->compVec[tmp_ijv->j[x]]));
    tmp->j[x] = min((p->compVec[tmp_ijv->i[x]]),(p->compVec[tmp_ijv->j[x]]));
  }
  if (!ijvCheckSortji(tmp)) {
    fprintf(stderr,"tmp not sorted\n");
  }
  tmp->n = 2;

  if (!ijvCheckSortji(meta)) {
    fprintf(stderr,"meta not sorted\n");
  }
  **********/
 


  /*----------------------------

  pass through the ijv, 
  computing the number edges in each off-diag block,
  and putting this number in counts

  */

  
  counts = (int *) cCalloc(meta->nnz, sizeof(int), "counts in offDiags");
  for (mx = 0; mx < meta->nnz; mx++)
    counts[mx] = 0;


  x = 0;
  mx = 0;
  while ((x < nnz) && (mx < meta->nnz)) {
    
    while ((x < nnz) &&
	   ((meta->j[mx] != min((p->compVec[tmp_ijv->i[x]]),(p->compVec[tmp_ijv->j[x]])))  || 
	    (meta->i[mx] != max((p->compVec[tmp_ijv->i[x]]),(p->compVec[tmp_ijv->j[x]]))))) {
      x++;
    }
    if (x == nnz) {
      fprintf(stderr,"error in offDiags, %d, %d, %d\n",mx,meta->nnz, nnz);
      printIJV(meta);
      fprintf(stderr,"\n\n");
      /*  printIJV(tmp); */
      break;
    }

    while ((x < nnz) &&
	   ((meta->j[mx] == min((p->compVec[tmp_ijv->i[x]]),(p->compVec[tmp_ijv->j[x]])))  &&
	    (meta->i[mx] == max((p->compVec[tmp_ijv->i[x]]),(p->compVec[tmp_ijv->j[x]]))))) {
      counts[mx]++;
      x++;
    }
    
    mx++;
  }

  if (mx < meta->nnz) {
    fprintf(stderr,"error2 in offDiags\n");
    return;
  }

  /* part of DEBUG stuff 
  freeIJV(tmp);
  */
  
  /*-----------------------------------------

  allocate the array of off-diagonal matrices

  */

  ijvArray = (ijvType **) cCalloc(meta->nnz,sizeof(ijvType *),
				  "ijvArray in st_offDiags");
  for (mx = 0; mx < meta->nnz; mx++) {
    ijvArray[mx] = makeIJV(counts[mx]);
  }

  /*-----------------------------------------

  pass through the same loop as before,
  this time populating the arrays,
  and keeping indices in counts

  */
  for (mx = 0; mx < meta->nnz; mx++)
    counts[mx] = 0;

  x = 0;
  mx = 0;
  while ((x < nnz) && (mx < meta->nnz)) {

    
    while ((x < nnz) &&
	   ((meta->j[mx] != min((p->compVec[tmp_ijv->i[x]]),(p->compVec[tmp_ijv->j[x]])))  || 
	    (meta->i[mx] != max((p->compVec[tmp_ijv->i[x]]),(p->compVec[tmp_ijv->j[x]]))))) {
      x++;
    }


    /* set the size of ijvArray[mx] here */
    ijvArray[mx]->n = max(p->partSizes[meta->j[mx]],
			  p->partSizes[meta->i[mx]]);

    while ((x < nnz) &&
	   ((meta->j[mx] == min((p->compVec[tmp_ijv->i[x]]),(p->compVec[tmp_ijv->j[x]])))  &&
	    (meta->i[mx] == max((p->compVec[tmp_ijv->i[x]]),(p->compVec[tmp_ijv->j[x]]))))) {


      /* make first element from component meta->i[mx] */
      if (meta->i[mx] == p->compVec[tmp_ijv->i[x]]) {
	ijvArray[mx]->i[counts[mx]] = p->compMap[tmp_ijv->i[x]];
	ijvArray[mx]->j[counts[mx]] = p->compMap[tmp_ijv->j[x]];
      } else {
	ijvArray[mx]->i[counts[mx]] = p->compMap[tmp_ijv->j[x]];
	ijvArray[mx]->j[counts[mx]] = p->compMap[tmp_ijv->i[x]];
      }
      ijvArray[mx]->v[counts[mx]] = tmp_ijv->v[x];

      counts[mx]++;
      x++;
    }

    mx++;
  }

  cFree(sortBy);
  cFree(counts);
  freeIJV(tmp_ijv);
  cFree(perm);
  cFree(perm2);
  
  
  return ijvArray;


  
}
