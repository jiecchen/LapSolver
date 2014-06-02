
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
#include "mex.h"

#include "st_defs.h"

/**
 * \defgroup mex
 *
 */

/**
 *
 * input must be symmetric
 * \ingroup mex
*/
myGraph *sparse2graph(const mxArray *array_ptr)
{
  double  *pr;
  int     *ir, *jc;
  int      col, row;
  int      starting_row_index, stopping_row_index, current_row_index;
  int      n, nnz;

  int deg;
  int cnt;

  int *counts;

  int i;
  
  myGraph *G;


  /* do some basic type checking */
  if (mxGetNumberOfDimensions(array_ptr) > 2)
    mexErrMsgTxt("Must be two-dimensional array.");

  if (mxGetClassID(array_ptr) != mxDOUBLE_CLASS) 
    mexErrMsgTxt("Must be array of doubles.");

  if (!mxIsSparse(array_ptr)) 
    mexErrMsgTxt("Must be sparse.");


  
  G = (myGraph *) cCalloc(1,sizeof(myGraph),"");

  /* Get the starting positions of all four data arrays. */ 
  pr = mxGetPr(array_ptr);
  ir = mxGetIr(array_ptr);
  jc = mxGetJc(array_ptr);

  n = mxGetN(array_ptr);
  nnz = jc[n];

  
  G->n = n;
  G->nnz = nnz;

  G->deg = (int *) cCalloc(n, sizeof(int),"");
  G->nbrs = (int **) cCalloc(n, sizeof(int *),"");
  G->wts = (double **) cCalloc(n, sizeof(double *),"");
  G->back = (int **) cCalloc(n, sizeof(int *),"");

  if (nnz == 0) {

    for (i = 0; i < n; i++) {
      G->deg[i] = 0;
      G->nbrs[i] = NULL;
      G->wts[i] = NULL;
      G->back[i] = NULL;

      G->nbrsBlock = NULL;
      G->wtsBlock = NULL;
      G->backBlock = NULL;
      
    }

    return (G);
  }

  /* allocate the data structures we will need to
     store the graph */

  G->nbrsBlock = (int *) cCalloc(nnz, sizeof(int ),"");
  G->wtsBlock = (double *) cCalloc(nnz, sizeof(double ),"");
  G->backBlock = (int *) cCalloc(nnz, sizeof(int),"");

  counts = (int *) cCalloc(n, sizeof(int),"counts");

  /* first, set all the degrees */
  for (col = 0; col < n; col++) {
    starting_row_index = jc[col]; 
    stopping_row_index = jc[col+1];
    deg = (stopping_row_index - starting_row_index);
    G->deg[col] = deg;
    G->nbrs[col] = &(G->nbrsBlock[starting_row_index]);
    G->wts[col] = &(G->wtsBlock[starting_row_index]);
    G->back[col] = &(G->backBlock[starting_row_index]);

    counts[col] = 0;
  }    


  /* Get the nonzero elements of the sparse array.
     Only looking at the sub-diagonal elements */
  for (col=0; col<n; col++)  { 
    starting_row_index = jc[col]; 
    stopping_row_index = jc[col+1];

    deg = stopping_row_index - starting_row_index;

    if (deg == 0)
      continue;
    else {
      for (current_row_index = starting_row_index; 
	   current_row_index < stopping_row_index; 
	   current_row_index++)  {

	row = ir[current_row_index];

	if (col < row) {

	  if ((counts[col] >= G->deg[col]) ||
	      (counts[row] >= G->deg[row]))
	    mexErrMsgTxt("Must be symmetric.");


	  
	  G->nbrs[col][counts[col]] = row;
	  G->nbrs[row][counts[row]] = col;

	  G->wts[col][counts[col]] = pr[current_row_index];
	  G->wts[row][counts[row]] = pr[current_row_index];

	  G->back[col][counts[col]] = counts[row];
	  G->back[row][counts[row]] = counts[col];

	  counts[col]++;
	  counts[row]++;

	}
      }
    }
  }

  /* check symmetry once again */
  for (col = 0; col < n; col++)
    if (counts[col] != G->deg[col])
      mexErrMsgTxt("Must be symmetric with zero diagonal.");


  return G;
}


/**
 *
 * \ingroup mex
*/
mxArray *graph2sparse(myGraph *G)
{
  int n, nnz;
  mxArray *mxA;

  double  *pr;
  int     *ir, *jc;
  
  int col, row, i;
  int ind;
  
  n = G->n;
  nnz = G->nnz;

  mxA = mxCreateSparse(n, n, nnz, mxREAL); 

  pr = mxGetPr(mxA);
  ir = mxGetIr(mxA);
  jc = mxGetJc(mxA);

  ind = 0;
  
  for (col = 0; col < n; col++) {

    jc[col] = ind;
    
    for (i = 0; i < G->deg[col]; i++) {
      ir[ind] = G->nbrs[col][i];
      pr[ind] = G->wts[col][i];

      ind++;
    }
  }

  jc[col] = ind;

  mxSetPr(mxA, pr);
  mxSetIr(mxA, ir);
  mxSetJc(mxA, jc);

  
  return mxA;
  
}


/**
 * needs a symmetric matrix!!!
 * \ingroup mex
*/
ijvType *sparse2ijv(const mxArray *array_ptr)
{

  double  *pr;
  int     *ir, *jc;
  int      col, row;
  int      starting_row_index, stopping_row_index, current_row_index;
  int      n, nnz;

  int ind;

  ijvType *ijv;


  /* do some basic type checking */
  if (mxGetNumberOfDimensions(array_ptr) > 2)
    mexErrMsgTxt("Must be two-dimensional array.");

  if (mxGetClassID(array_ptr) != mxDOUBLE_CLASS) 
    mexErrMsgTxt("Must be array of doubles.");

  if (!mxIsSparse(array_ptr)) 
    mexErrMsgTxt("Must be sparse.");

  
  ijv = (ijvType *) cCalloc(1,sizeof(ijvType),"");
  
  /* Get the starting positions of all four data arrays. */ 
  pr = mxGetPr(array_ptr);
  ir = mxGetIr(array_ptr);
  jc = mxGetJc(array_ptr);

  n = mxGetN(array_ptr);
  nnz = jc[n]/2;

  ijv->n = n;
  ijv->nnz = nnz;

  ijv->i = (int *) cCalloc(nnz, sizeof(int),"");
  ijv->j = (int *) cCalloc(nnz, sizeof(int),"");
  ijv->v = (double *) cCalloc(nnz, sizeof(double),"");
  
  ind = 0;
  for (col=0; col<n; col++)  { 
    starting_row_index = jc[col]; 
    stopping_row_index = jc[col+1];

    if (stopping_row_index > starting_row_index)
      for (current_row_index = starting_row_index; 
	   current_row_index < stopping_row_index; 
	   current_row_index++)  {

	row = ir[current_row_index];

	if (col < row) {
	  ijv->i[ind] = row;
	  ijv->j[ind] = col;
	  ijv->v[ind] = pr[current_row_index];
	  ind++;
	}
      }
    }

  return ijv;
}

/**
 * Note: only outputs lower-triangular part!
 * needs data to be sorted by column, the j coord
 * \ingroup mex
*/
mxArray *ijv2sparse(ijvType *ijv)
{
  int n, nnz;
  mxArray *mxA;

  double  *pr;
  int     *ir, *jc;
  
  int col, row, i;
  int ind;
  int prevCol;
  
  n = ijv->n;
  nnz = ijv->nnz;

  if (ijvCheckSortj(ijv) == 0) {
    ijvSortj(ijv);
  }

  mxA = mxCreateSparse(n, n, nnz, mxREAL); 

  pr = mxGetPr(mxA);
  ir = mxGetIr(mxA);
  jc = mxGetJc(mxA);

  ind = 0;
  col = 0;

  for (col = 0; col < n; col++) {
    jc[col] = ind;
    
    while ((ind < nnz) && (ijv->j[ind] == col))
      ind++;
  }
  jc[n] = nnz;


  for (i = 0; i < nnz; i++) {
    ir[i] = ijv->i[i];
    pr[i] = ijv->v[i];
  }


  mxSetPr(mxA, pr);
  mxSetIr(mxA, ir);
  mxSetJc(mxA, jc);

  return mxA;
}


/**
 * mainly for diagnostics
 * \ingroup mex
*/
void printSparse(const mxArray *array_ptr)
{
  double  *pr;
  int     *ir, *jc;
  int      n, nnz;

  int i;
  
  /* do some basic type checking */
  if (mxGetNumberOfDimensions(array_ptr) > 2)
    mexErrMsgTxt("Must be two-dimensional array.");

  if (mxGetClassID(array_ptr) != mxDOUBLE_CLASS) 
    mexErrMsgTxt("Must be array of doubles.");

  if (!mxIsSparse(array_ptr)) 
    mexErrMsgTxt("Must be sparse.");


  /* Get the starting positions of all four data arrays. */ 
  pr = mxGetPr(array_ptr);
  ir = mxGetIr(array_ptr);
  jc = mxGetJc(array_ptr);

  n = mxGetN(array_ptr);
  nnz = jc[n];

  fprintf(stderr,"n: %d, nnz: %d\n",n,nnz);

  fprintf(stderr,"jc: ");
  for (i = 0; i <= n; i++)
    fprintf(stderr,"%d ",jc[i]);
  fprintf(stderr,"\n");


  fprintf(stderr,"ir: ");
  for (i = 0; i < nnz; i++)
    fprintf(stderr,"%d ",ir[i]);
  fprintf(stderr,"\n");


  fprintf(stderr,"pr: ");
  for (i = 0; i < nnz; i++)
    fprintf(stderr,"%g ",pr[i]);
  fprintf(stderr,"\n");

}
    
/** 
 * 
 * 
 * @param p the partition to convert
 * 
 * @return a cell array with one array for each part
 */
mxArray *part2cell (partType *p)
{
  mxArray *cell;
  mxArray *cellElem;
  double *cellElemData;

  int i,j;
  
  int dims[2];

  dims[0] = 1;
  dims[1] = p->numParts;
  
  cell = mxCreateCellArray(2, dims);

  for (j = 0; j < p->numParts; j++) {
    cellElem = mxCreateDoubleMatrix(1, p->partSizes[j] , mxREAL);
    cellElemData = mxGetPr(cellElem);
    for (i = 0; i < p->partSizes[j]; i++)
      cellElemData[i] = (double) p->partMap[j][i] + 1;
    mxSetCell(cell, j, cellElem);
  }

  return cell;
}

/**
 * \ingroup mex
 */
mxArray *ijvBlockDiag2cell (ijvBlockDiag *bd)
{

  mxArray *cell;
  mxArray *cellElem;
  
  partType *p;

  int i,j;
  
  int dims[2];

  p = bd->part;
  
  dims[0] = 1;
  dims[1] = p->numParts;
  
  cell = mxCreateCellArray(2, dims);

  for (j = 0; j < p->numParts; j++) {
    cellElem = ijv2sparse(bd->blocks[j]);
    mxSetCell(cell, j, cellElem);
  }

  return cell;
}

/**
 * \ingroup mex
 * takes as input a cell array of sparse graphs,
 * and outputs an array of ijvs
 *
 * the input output is the number of graphs found
 */
int cell2array(const mxArray *cellIn,  ijvType*** arrayOut) 
{
  int       numCells;
  int       i;
  const mxArray *cellElem;

  ijvType **aOut;

  mxClassID  category;
  category = mxGetClassID(cellIn);
  if (category != mxCELL_CLASS)
    mexErrMsgTxt("Input must be a cell array.");

  numCells = mxGetNumberOfElements(cellIn); 
  
  aOut = (ijvType **) cCalloc(numCells, sizeof(ijvType *),"aOut in cell2array");

  for (i = 0; i < numCells; i++) {
    aOut[i] = sparse2ijv(mxGetCell(cellIn,i));
    if (!aOut[i])
      mexErrMsgTxt("Error in cell2array.");
  }

  *arrayOut = aOut;

  return numCells;

}
