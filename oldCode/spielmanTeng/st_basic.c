/*******************
   st_basic.c

   basic operations on our data structures

*******************/


#include <stdio.h>
#include <stdlib.h>

#include "st_defs.h"


#define min(A,B) ((A) < (B) ? (A) : (B))
#define max(A,B) ((A) > (B) ? (A) : (B))

void *cCalloc(size_t n, size_t size, char *errMsg)
{

  void *obj;

  if (n == 0)
    return NULL;

  
#ifdef mex_h
  fprintf(stderr,"xxx\n");
  obj = (void *) mxCalloc(n,size);
#else  
  obj = calloc(n,size);
#endif
  
  if (!(obj)) {
    fprintf(stderr,"calloc failed: %s\n, mem %d\n", errMsg,n*size);
    exit(1);
  }

  return obj; 
}

void *cMalloc(size_t n, char *errMsg)
{

  void *obj;
  

  if (n == 0)
    return NULL;


#ifdef mex_h
  fprintf(stderr,"xxx\n");
  obj = (void *) mxMalloc(n);
#else  
  obj = (void *) malloc(n);
#endif
  
  if (!(obj)) {
    fprintf(stderr,"malloc failed: %s, mem %d\n", errMsg, n);
    exit(1);
  }

  return obj; 
}


void cFree(void *obj)
{

  if (!obj)
    fprintf(stderr,"cFree error!\n");
  
#ifdef mex_h
  if (obj) mxFree(obj);
#else  
  if (obj) free(obj);
#endif
}



void printGraph(myGraph *G)
{
  int col;
  int i;

  fprintf(stderr,"\n");

  fprintf(stderr,"G->n %d\n",G->n);
  
  for (col = 0; col < G->n; col++) {

    fprintf(stderr,"%d: %d \n",col, G->deg[col]);
    
    for (i = 0; i < G->deg[col]; i++) {
      fprintf(stderr,"%d (%g), ", G->nbrs[col][i], G->wts[col][i]); 
    }
    fprintf(stderr,"\n");

  }

}

/**
 * \relatesalso myGraph
 *
 */
void freeGraph(myGraph *G)
{
  int i;

  for (i = 0; i < G->n; i++) {
    G->nbrs[i] = NULL;
    G->wts[i] = NULL;
  }
    
  if (G->nbrs) cFree(G->nbrs);
  if (G->wts) cFree(G->wts);
  if (G->back) cFree(G->back);
  if (G->deg) cFree(G->deg);
  
  if (G->nnz > 0) {
    cFree(G->backBlock);
    cFree(G->wtsBlock);
    cFree(G->nbrsBlock);
  }
  
  cFree(G);
  
}

int * CountingSortOld(int * A, int n, int k){
  int *perm, *C;
  int i;
  C = cCalloc(k, sizeof(int),"");  
  perm = cCalloc(n, sizeof(int),"");

  /* not needed as calloc makes 0
  for (i = 0; i< k; i++){
    C[i] = 0;
  }
  */

  /* replace with below
  for (i = 0; i< n; i++){
     C[A[i]] = C[A[i]] + 1;
  }
  */


  for (i = 0; i< n; i++){
    C[A[i]]++;
  }

  for (i = 1; i< k; i++){
    C[i] = C[i] + C[i-1];
  }

  for (i = n-1; i>-1; i--){
    perm[i] = C[A[i]] - 1;
    C[A[i]] = C[A[i]] -1;
  }

  if (C) cFree(C);

  return perm;
}
 

int * CountingSort(int * A, int n, int k){
  int *perm, *C;
  int i;
  C = cCalloc(k, sizeof(int),"");  
  perm = cCalloc(n, sizeof(int),"");

  /*
  for (i = 0; i< k; i++){
    C[i] = 0;
  }
  */

  /*
  for (i = 0; i< n; i++){
     C[A[i]] = C[A[i]] + 1;
  }
  */

  int limit;
  limit = 3*(n/3);
  register int a1, a2, a3;

  for (i = 0; i< limit; i+= 3){
    a1 = A[i];
    a2 = A[i+1];
    a3 = A[i+2];
    C[a1]++;
    C[a2]++;
    C[a3]++;
  }
  if( i < n ) 
    { 
      switch( n - i ) 
	{ 
        case 2 : C[A[i++]]++;
        case 1 : C[A[i++]]++;
	}
    } 

  
  for (i = 1; i< k; i++){
    C[i] = C[i] + C[i-1];
  }

  register int b1;
  for (i = n-1; i>-1; i--){
    b1 = A[i];
    perm[i] = C[b1] - 1;
    C[b1]--;
  }

  if (C) cFree(C);

  return perm;
}
 


/*
 * \ingroup ijv
 *
 */
void printIJV(ijvType *ijv)
{
  int col;
  int i;

  fprintf(stderr,"IJV: \n");

  fprintf(stderr,"ijv->n: %d, ijv->nnz: %d \n",ijv->n, ijv->nnz);
  
  for (i = 0; i < ijv->nnz; i++) {
    fprintf(stderr,"(%d,%d) : %g \n",ijv->i[i], ijv->j[i], ijv->v[i]);
  }
    
  fprintf(stderr,"\n");
}


/**
 * \ingroup ijv
 *
 */
  
ijvType *makeIJV(int nnz)
{

  ijvType *ijv;

  ijv = (ijvType *) cCalloc(1,sizeof(ijvType),
			    "makeIJV");

  ijv->nnz = nnz;
  
  ijv->i = (int *) cCalloc(nnz,sizeof(int),"makeIJV");
  ijv->j = (int *) cCalloc(nnz,sizeof(int),"makeIJV");
  ijv->v = (double *) cCalloc(nnz,sizeof(double),"makeIJV");

  return ijv;
}



#define int32 int

ijvType *binReadIJV(FILE* fpr) 
{
  int n, nnz;
  
  if (!fread(&(n), sizeof(int32), 1, fpr))
    cError ("error reading n from file\n");

  if (!fread(&(nnz), sizeof(int32), 1, fpr))
    cError ("error reading nnz from file\n");

  ijvType *ijv;
  ijv = makeIJV(nnz);
  ijv->n = n;

  int x;
  int i, j;
  float v;

  if (!fread(ijv->i, sizeof(int32), nnz, fpr))
    cError ("error reading i from file\n");

  if (!fread(ijv->j, sizeof(int32), nnz, fpr))
    cError ("error reading j from file\n");

  if (!fread(ijv->v, sizeof(double), nnz, fpr))
    cError ("error reading v from file\n");

  for (x = 0; x < nnz; x++) {
    ijv->i[x] = ijv->i[x]-1;
    ijv->j[x] = ijv->j[x]-1;
  }

  return ijv;
}
		    

/**
 * \ingroup ijv
 *
 */
  
ijvType *copyIJV(ijvType *ijv_in)
{

  ijvType *ijv;

  ijv = (ijvType *) cCalloc(1,sizeof(ijvType),
			    "copyIJV");

  ijv->nnz = ijv_in->nnz;
  ijv->n = ijv_in->n;

  
  ijv->i = (int *) cCalloc(ijv->nnz,sizeof(int),"copyIJV");
  ijv->j = (int *) cCalloc(ijv->nnz,sizeof(int),"copyIJV");
  ijv->v = (double *) cCalloc(ijv->nnz,sizeof(double),"copyIJV");

  int x;
  for (x = 0; x < ijv->nnz; x++) {
    ijv->i[x] = ijv_in->i[x];
    ijv->j[x] = ijv_in->j[x];
    ijv->v[x] = ijv_in->v[x];
  }		  
  
  return ijv;
}


/**
 * \ingroup ijv
 *
 */
void freeIJV(ijvType *ijv)
{
  int i;

  if (ijv->i) cFree(ijv->i);
  if (ijv->j) cFree(ijv->j);
  if (ijv->v) cFree(ijv->v);

  cFree(ijv);
  
}




/**

  \brief
   prints s to the standard error and halts execution
*/
cError(char *s)
{
    fprintf (stderr, s);
    exit (1);
}


/**
 * \ingroup ijv
 * \brief sort an ijv on the i entries *
 *
 */
void ijvSorti(ijvType *ijv)
{

  int *perm;
  ijvType *ijv_temp;
  int e;
  
  ijv_temp = makeIJV(ijv->nnz);
  ijv_temp->n = ijv->n;

  perm = CountingSort(ijv->i, ijv->nnz, ijv->n);

  for (e=0; e < ijv->nnz; e++){
    ijv_temp->i[perm[e]] = ijv->i[e];
    ijv_temp->j[perm[e]] = ijv->j[e];
    ijv_temp->v[perm[e]] = ijv->v[e];
  }
  if (perm) cFree(perm);

  for (e=0; e < ijv->nnz; e++){
    ijv->i[e] = ijv_temp->i[e];
    ijv->j[e] = ijv_temp->j[e];
    ijv->v[e] = ijv_temp->v[e];
  }
  
  freeIJV(ijv_temp);
}

/**
 * \ingroup ijv
 * \brief sort an ijv on the j entries *
 *
 */
void ijvSortj(ijvType *ijv)
{

  int *perm;
  ijvType *ijv_temp;
  int e;

  ijv_temp = makeIJV(ijv->nnz);

  ijv_temp->n = ijv->n;
 

  perm = CountingSort(ijv->j, ijv->nnz, ijv->n);

  for (e=0; e < ijv->nnz; e++){
    ijv_temp->i[perm[e]] = ijv->i[e];
    ijv_temp->j[perm[e]] = ijv->j[e];
    ijv_temp->v[perm[e]] = ijv->v[e];
  }
  if (perm) cFree(perm);


  for (e=0; e < ijv->nnz; e++){
    ijv->i[e] = ijv_temp->i[e];
    ijv->j[e] = ijv_temp->j[e];
    ijv->v[e] = ijv_temp->v[e];
  }

  freeIJV(ijv_temp);

}


/**
 * \ingroup ijv
 * \brief sort an ijv on the i entries ,
 *         and then on the j entries
 *
 * not yet working!!!
 */
void ijvSortij(ijvType *ijv_orig)
{

  int *perm;
  int e;

  ijvType *ijv_temp;
  
  perm = CountingSort(ijv_orig->i, ijv_orig->nnz, ijv_orig->n);
  ijv_temp = makeIJV(ijv_orig->nnz);
  ijv_temp->n = ijv_orig->n;

  for (e=0; e < ijv_orig->nnz; e++){
    ijv_temp->i[perm[e]] = ijv_orig->i[e];
    ijv_temp->j[perm[e]] = ijv_orig->j[e];
    ijv_temp->v[perm[e]] = ijv_orig->v[e];
  }
  cFree(perm);
  
  perm = CountingSort(ijv_temp->j, ijv_temp->nnz, ijv_temp->n);

  for (e=0; e < ijv_temp->nnz; e++){
    ijv_orig->i[perm[e]] = ijv_temp->i[e];
    ijv_orig->j[perm[e]] = ijv_temp->j[e];
    ijv_orig->v[perm[e]] = ijv_temp->v[e];
  }
  cFree(perm);
}

/**
 * \ingroup ijv
 * \brief return 1 if sorted on j, 0 otherwise
 *
 */
int ijvCheckSortj(ijvType *ijv)
{
  int i;
  int v;

  /* if empty, is sorted */
  if (ijv->nnz == 0)
    return 1;
  
  v = ijv->j[0];
  for (i = 0; i < ijv->nnz; i++) {
    if (ijv->j[i] < v)
      return 0;
    else
      v = ijv->j[i];
  }
    
  return 1;
}

/**
 * \ingroup ijv
 * \brief return 1 if sorted on (i,j), with
 *  i primary, 0 otherwise
 *
 */
int ijvCheckSortji(ijvType *ijv)
{
  int x;

  int previ, prevj;

  previ = 0;
  prevj = 0;

  /* if empty, is sorted */
  if (ijv->nnz == 0)
    return 1;
  
  for (x = 0; x < ijv->nnz; x++) {
    if (ijv->i[x] < previ)
      return 0;
    else
      if (ijv->i[x] == previ && (ijv->j[x] < prevj))
	return 0;
      else {
	previ = ijv->i[x];
	prevj = ijv->j[x];
      }
    
  }
    
  return 1;
}


/** 
 * 
 * 
 * @param ijv the input ijv
 * 
 * @return a myGraph type
 *
 * @todo  make a light version of that that just
 *        uses adjacency structure,
 *        as the only place we use it now only
 *        requires that.
 */
myGraph *ijv2graph(ijvType *ijv)
{
  int x;
  
  int ind;
  
  myGraph *G;

  G = (myGraph *) cCalloc(1,sizeof(myGraph),"G in ijv2sparse");

  G->n = ijv->n;
  G->nnz = 2 * ijv->nnz;


  G->deg = (int *) cCalloc(G->n, sizeof(int),"G->deg in ijv2sparse");
  G->nbrs = (int **) cCalloc(G->n, sizeof(int *),"G->nbrs in ijv2sparse");
  G->wts = (double **) cCalloc(G->n, sizeof(double *), "G->wts in ijv2sparse");
  G->back = (int **) cCalloc(G->n, sizeof(int *),"G->back in ijv2sparse");

  if (ijv->nnz == 0) {
    for (x = 0; x < G->n; x++) {
      G->deg[x] = 0;
      G->nbrs[x] = NULL;
      G->wts[x] = NULL;
      G->back[x] = NULL;

      G->nbrsBlock = NULL;
      G->wtsBlock = NULL;
      G->backBlock = NULL;
    }
    return (G);
  }

  /* allocate the data structures we will need to 
     store the graph */

  G->nbrsBlock = (int *) cCalloc(G->nnz, sizeof(int),"nbrsBlock in ijv2sparse");
  G->wtsBlock = (double *) cCalloc(G->nnz, sizeof(double),"wtsBlock in ijv2sparse");
  G->backBlock = (int *) cCalloc(G->nnz, sizeof(int),"backBlock in ijv2sparse");

  for (x = 0; x < G->n; x++)
    G->deg[x] = 0;

  /* first, set all the degrees */
  for (x = 0; x < ijv->nnz; x++) {
     G->deg[ijv->i[x]]++;
     G->deg[ijv->j[x]]++;
  }

  /* then, allocate space for the data for each vertex */
  ind = 0;
  for (x = 0; x < G->n; x++) {
    G->nbrs[x] = &(G->nbrsBlock[ind]);
    G->wts[x] = &(G->wtsBlock[ind]);
    G->back[x] = &(G->backBlock[ind]);
    ind += G->deg[x];
  }

  /* and now, populate the graph with edges,
     temporarily resetting the degrees
     (don't worry: they will return to their prev values */
  for (x = 0; x < G->n; x++)
    G->deg[x] = 0;

  int row; /* i */
  int col; /* j */

  for (x = 0; x < ijv->nnz; x++) {

    row = ijv->i[x];
    col = ijv->j[x];

    G->nbrs[col][G->deg[col]] = row;
    G->nbrs[row][G->deg[row]] = col;

    G->wts[col][G->deg[col]] = ijv->v[x];
    G->wts[row][G->deg[row]] = ijv->v[x];
    
    G->back[col][G->deg[col]] = G->deg[row];
    G->back[row][G->deg[row]] = G->deg[col];
	  
    G->deg[row]++;
    G->deg[col]++;
  }

  return G;
}

/* n is the length of compVec */
partType *compVec2part(int *compVec, int n)
{
  int i;

  int ind;
  
  partType *p;


  if (n < 1) {
    fprintf(stderr,"compVec2part passed length 0 compVec\n");
    return NULL;
  }

  p = (partType *) cCalloc(1,sizeof(partType),"p in compVec2part");

  p->compVec = (int *) cCalloc(n,sizeof(int),"compVec in  compVec2part");
  p->compMap = (int *) cCalloc(n,sizeof(int),"compMap in  compVec2part");

  /* copy in compVec */
  for (i = 0; i < n; i++)
    p->compVec[i] = compVec[i];

  
  /* compute the number of parts */
  /* (note compVec is 0 based ) */
  p->numParts = 0;
  for (i = 0; i < n; i++)
    if (compVec[i] >= p->numParts)
      p->numParts = compVec[i] + 1;

  /* allocate the partSizes */
  p->partSizes = (int *) cCalloc(p->numParts, sizeof(int),
				    "partSizes in compVec2part");
  for (i = 0; i < p->numParts; i++)
    p->partSizes[i] = 0;


  /* compute the partSizes */
  for (i = 0; i < n; i++) 
    p->partSizes[compVec[i]]++;


  /* allocate and divvy up partMapBlock */
  p->partMapBlock = (int *) cCalloc(n, sizeof(int),
				    "partMapBlock in compVec2part");
  p->partMap = (int **) cCalloc(p->numParts, sizeof(int),
				    "partMap in compVec2part");
  
  ind = 0;
  for (i = 0; i < p->numParts; i++) {
    p->partMap[i] = &(p->partMapBlock[ind]);
    ind += p->partSizes[i];
  }

  /* in final pass, create partMap and compMap */
  for (i = 0; i < p->numParts; i++) 
    p->partSizes[i] = 0;

  for (i = 0; i < n; i++) {
    p->compMap[i] = p->partSizes[compVec[i]];
    p->partMap[compVec[i]][p->partSizes[compVec[i]]] = i;
    p->partSizes[compVec[i]]++;
  }

  
  return p;
}

void freePart(partType *p) 
{
  if (p != NULL) {

    if (p->compVec) cFree(p->compVec);
    if (p->compMap) cFree(p->compMap);
    if (p->partSizes) cFree(p->partSizes);
    if (p->partMap) cFree(p->partMap);
    if (p->partMapBlock) cFree(p->partMapBlock);

    cFree(p);
  }
}

ijvBlockDiag *formDiagBlocks(ijvType *ijv, partType *p)
{
  ijvBlockDiag *bd;

  int i,x;

  int *nnzs;
  
  /* make sure part and ijv exist */
  if (!p || !ijv)
    return NULL;

  /* allocate bd and the list of ijvs */
  bd = (ijvBlockDiag *) cCalloc(1,sizeof(ijvBlockDiag),"bd in formDiagBlocks");
  bd->n = ijv->n;
  bd->blocks = (ijvType **) cCalloc(p->numParts,sizeof(ijvType *),
				    "bd->blocks in formDiagBlocks");
  bd->part = p;

  /* find the number of non-zeros in each graph */
  nnzs = (int *) cCalloc(p->numParts, sizeof(int),
			 "nnzs in formDiagBlocks");
  for (i = 0; i < p->numParts; i++)
    nnzs[i] = 0;

  for (x = 0; x < ijv->nnz; x++) {
    if (p->compVec[ijv->i[x]] == p->compVec[ijv->j[x]])
      nnzs[p->compVec[ijv->i[x]]]++;
  }

  /* now, we can alloc the blocks */
  for (i = 0; i < p->numParts; i++) {
    bd->blocks[i] =  makeIJV(nnzs[i]);
    bd->blocks[i]->n = p->partSizes[i];
  }

  /* in a final pass, populate the blocks */
  for (i = 0; i < p->numParts; i++)
    nnzs[i] = 0;

  int bl;
  for (x = 0; x < ijv->nnz; x++) 
    if (p->compVec[ijv->i[x]] == p->compVec[ijv->j[x]]) {
      bl = p->compVec[ijv->i[x]];
      
      bd->blocks[bl]->i [ nnzs[bl] ] = 
	p->compMap[ijv->i[x]];
      
      bd->blocks[bl]->j [ nnzs[bl] ] = 
	p->compMap[ijv->j[x]];

      bd->blocks[bl]->v [ nnzs[bl] ] = 
	ijv->v[x];

      nnzs[bl]++;
    }

  return bd;

}

void freeIjvBlockDiagAndPart(ijvBlockDiag *bd) 
{
  int i;
  
  if (bd) {
    if (bd->blocks)
      for (i = 0; i < bd->part->numParts; i++)
	freeIJV(bd->blocks[i]);

    if (bd->part) freePart(bd->part);

    cFree(bd);
  }
}

void freeIjvBlockDiagNotPart(ijvBlockDiag *bd) 
{
  int i;
  
  if (bd) {
    if (bd->blocks)
      for (i = 0; i < bd->part->numParts; i++)
	freeIJV(bd->blocks[i]);

    cFree(bd);
  }
}


void printIntVec(int *vec, int n, char* msg)
{
  int i;

  fprintf(stderr,"%s\n",msg);
  for (i = 0; i < n; i++)
    fprintf(stderr,"%d ",vec[i]);
  fprintf(stderr,"\n");
  
}

/**
 *  \ingroup ijv
 *  \brief combine the ijvs tail to tail, but does not
 *         sort or remove duplicate edges.
 */
ijvType *mergeIJVs(int n, /** the size of the output graph */
		   ijvType **ijvs, /** the array of ijvs */
		   int numIJVs) /** the length of the list */
{

  int x;
  int nnz;

  ijvType *ijv;

  /* scan down the list, collecting nnz,
     and finding the maximum n */
  nnz = 0;
  for (x = 0; x < numIJVs; x++) {
    n = max(n, ijvs[x]->n);
    nnz += ijvs[x]->nnz;
  }

  /* allocate the new ijv */
  ijv = makeIJV(nnz);
  ijv->n = n;

  /* re-scan, populating the ijv */
  int ind;
  ind = 0;
  for (x = 0; x < numIJVs; x++) {
    int y;

    for (y = 0; y < ijvs[x]->nnz; y++) {
      ijv->i[ind] = ijvs[x]->i[y];
      ijv->j[ind] = ijvs[x]->j[y];
      ijv->v[ind] = ijvs[x]->v[y];
      ind++;
    }
  }

  return ijv;

}


/**
 *  \ingroup ijv
 *  \brief removed duplicate edges in an IJV
 *         sorts its input, and returns input sorted
 *         with j major
 */
ijvType *compressIJV(ijvType *ijvIn)
{
  char *newEntry;
  
  int x;

  ijvType *ijvOut;

  int prevI, prevJ;
  int count;
  
  ijvSorti(ijvIn);
  ijvSortj(ijvIn);

  newEntry = (char *) cCalloc(ijvIn->nnz, sizeof(char),"compressIJV");
  /* note all set to 0 */


  prevI = -1;
  prevJ = -1;
  count = 0;

  for (x = 0; x < ijvIn->nnz; x++) {
    
    /* skip until get a new one */
    while ((x < ijvIn->nnz) &&
	   ((ijvIn->i[x] == prevI) && (ijvIn->j[x] == prevJ)))
      x++;

    if (x < ijvIn->nnz) {
      count++;
      newEntry[x] = 1;
      prevI = ijvIn->i[x];
      prevJ = ijvIn->j[x];
    }
  }

  ijvOut = makeIJV(count);
  ijvOut->n = ijvIn->n;

  int ind;
  ind = 0;
  
  for (x = 0; x < ijvIn->nnz; x++) 
    if (newEntry[x]) {
      ijvOut->i[ind] = ijvIn->i[x];
      ijvOut->j[ind] = ijvIn->j[x];
      ijvOut->v[ind] = ijvIn->v[x];
      ind++;
    }

  cFree(newEntry);
  
  return ijvOut;
}



/**
 * \ingroup ijv
 *
 * \brief repeats each edge both ways,
 *    unsorted
 *
 * converts an ordinary ijv
 * to a dubIjvtype: one in which
 * each edge appears twice
 */
dubIjvType *doubleIJV(const ijvType *ijvIn)
{

  int x;

  ijvType *ijv;

  ijv = (ijvType *) cCalloc(1,sizeof(ijvType),
			    "doubleIJV");
  ijv->nnz = 2* (ijvIn->nnz);
  ijv->n = ijvIn->n;
  ijv->i = (int *) cCalloc(ijv->nnz,sizeof(int),"doubleIJV");
  ijv->j = (int *) cCalloc(ijv->nnz,sizeof(int),"doubleIJV");
  ijv->v = (double *) cCalloc(ijv->nnz,sizeof(double),"doubleIJV");


  for (x = 0; x < ijvIn->nnz; x++) {
    ijv->i[x] = ijvIn->i[x];
    ijv->j[x] = ijvIn->j[x];
    ijv->v[x] = ijvIn->v[x];

    ijv->i[x + ijvIn->nnz] = ijvIn->j[x];
    ijv->j[x + ijvIn->nnz] = ijvIn->i[x];
    ijv->v[x + ijvIn->nnz] = ijvIn->v[x];
  }

  /*
  memcpy(ijv->i, ijvIn->i, ijv->nnz * sizeof(int));
  memcpy(ijv->i + ijvIn->nnz, ijvIn->j, ijv->nnz * sizeof(int));

  memcpy(ijv->j, ijvIn->j, ijv->nnz * sizeof(int));
  memcpy(ijv->j + ijvIn->nnz, ijvIn->i, ijv->nnz * sizeof(int));

  memcpy(ijv->v, ijvIn->v, ijv->nnz * sizeof(double));
  memcpy(ijv->v + ijvIn->nnz, ijvIn->v, ijv->nnz * sizeof(double));
  */
  
  return ijv;

}

/**
 * \ingroup ijv
 *
 * \brief remaps the vertex labels in an ijv
 *   according to a partition
 *
 */
void remapIJV(ijvType *ijv, /*< the ijv */
	      partType *p,  /*< the partition */
	      int comp)    /*< the index into the partion */
{
  int x;

  for (x = 0; x < ijv->nnz; x++) {
    ijv->i[x] = p->partMap[comp][ijv->i[x]];
    ijv->j[x] = p->partMap[comp][ijv->j[x]];
  }
}
