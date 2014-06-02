/*******************
    st_metaGraph.c

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
  *******************************************************/

void freeGraph(myGraph *G);
void printGraph(myGraph *G);
int * CountingSort(int * A, int n, int k);
ijvType *makeIJV(int nnz);


/* It maps each vertex of the original graph to
 *  its component vertex. It also update the n entry.
 */

CompMap(ijvType *ijv_orig, int *compVec){
  int k, i, e;
  k = 0;
  
/*  compute the number of components */
  for (i=0; i < ijv_orig->n; i++){
    if (k < compVec[i]){
      k = compVec[i];
    }
  }
  k++;

  int tmp;
  
 /* map the vertices */
  for (e=0; e < ijv_orig->nnz; e++){
    if (compVec[ijv_orig->i[e]] > compVec[ijv_orig->j[e]]){
      ijv_orig->i[e] = compVec[ijv_orig->i[e]];
      ijv_orig->j[e] = compVec[ijv_orig->j[e]];
    } else{
      tmp = ijv_orig->i[e];
      ijv_orig->i[e] = compVec[ijv_orig->j[e]];
      ijv_orig->j[e] = compVec[tmp];
    }
  }
  ijv_orig->n = k;
}


ijvType *st_metaGraph(ijvType *ijv, int *compVec)
{ 
  int * perm;
  ijvType *ijv_temp, *ijv_comp;
  int k, e;
  int *index;
  double * value;
  int current_i, current_j;

  ijvType *ijv_orig;

  ijv_orig = makeIJV(ijv->nnz);
  ijv_orig->n = ijv->n;
  /* copy ijv into ijv_orig, since we destroy it */
  for (e = 0; e < ijv->nnz; e++) {
    ijv_orig->i[e] = ijv->i[e];
    ijv_orig->j[e] = ijv->j[e];
    ijv_orig->v[e] = ijv->v[e];
  }
  
  ijv_temp = (ijvType *) cCalloc(1,sizeof(ijvType),"ijv_temp in st_metaGraph");
  ijv_comp = (ijvType *) cCalloc(1,sizeof(ijvType),"ijv_comp in st_metaGraph");

  #ifdef ST_DEBUG
  int i;
  for (i = 0; i < ijv_orig->n; i++)
    fprintf(stderr,"%d:%d ",i,compVec[i]);
  fprintf(stderr,"\n");
  #endif
  
  #ifdef ST_DEBUG
  fprintf(stderr,"AAA\n");
  printIJV(ijv_orig);
  fprintf(stderr,"\n");
  #endif


  CompMap(ijv_orig, compVec);

  #ifdef ST_DEBUG
  fprintf(stderr,"000\n");
  printIJV(ijv_orig);
  fprintf(stderr,"\n");
  #endif
  
/* sort the j entry and permute the structure and save it to ijv_temp*/

  perm = CountingSort(ijv_orig->j, ijv_orig->nnz, ijv_orig->n);
  ijv_temp->n = ijv_orig->n;
  ijv_temp->nnz = ijv_orig->nnz;

  ijv_temp->i = (int *) cCalloc(ijv_orig->nnz,sizeof(int),"ijv_temp->i in st_metaGraph");
  ijv_temp->j = (int *) cCalloc(ijv_orig->nnz,sizeof(int),"ijv_temp->j in st_metaGraph");
  ijv_temp->v = (double *) cCalloc(ijv_orig->nnz,sizeof(double),"ijv_temp->v in st_metaGraph");

  for (e=0; e < ijv_orig->nnz; e++){
    ijv_temp->i[perm[e]] = ijv_orig->i[e];
    ijv_temp->j[perm[e]] = ijv_orig->j[e];
    ijv_temp->v[perm[e]] = ijv_orig->v[e];
  }
  cFree(perm);

  #ifdef ST_DEBUG
  fprintf(stderr,"111\n");
  printIJV(ijv_temp);
  fprintf(stderr,"\n");
  #endif


  
/* sort the i entry and permute the structure and save it back to ijv_orig*/

  perm = CountingSort(ijv_temp->i, ijv_temp->nnz, ijv_temp->n);

  for (e=0; e < ijv_temp->nnz; e++){
    ijv_orig->i[perm[e]] = ijv_temp->i[e];
    ijv_orig->j[perm[e]] = ijv_temp->j[e];
    ijv_orig->v[perm[e]] = ijv_temp->v[e];
  }
  cFree(perm);

  

  #ifdef ST_DEBUG
  fprintf(stderr,"222\n");
  printIJV(ijv_orig);
  fprintf(stderr,"\n");
  #endif


/* compress and sum */
  index = (int *) cCalloc(ijv_orig->nnz,sizeof(int),"index in st_metaGraph");
  value = (double *) cCalloc(ijv_orig->nnz,sizeof(double),"value in st_metaGraph");

  freeIJV(ijv_temp);
  
  k = 0;
  current_i =-1;
  for (e=0; e < ijv_orig->nnz; e++){
    if (ijv_orig->i[e] != ijv_orig->j[e]){
      if (current_i == -1){
	index[k] = e;
        value[k++] = ijv_orig->v[e];
        current_i = ijv_orig->i[e];
        current_j = ijv_orig->j[e];
      }else{
	if ((current_i==ijv_orig->i[e]) && (current_j==ijv_orig->j[e])){
          value[k-1] += ijv_orig->v[e];
	} else{
          index[k] = e;
          value[k++] = ijv_orig->v[e];
          current_i = ijv_orig->i[e];
          current_j = ijv_orig->j[e];
	}
      }
    }
  }
  ijv_comp->n = ijv_orig->n;
  ijv_comp->nnz = k;
  ijv_comp->i = (int *) cCalloc(k,sizeof(int),"ijv_comp->i in st_metaGraph");
  ijv_comp->j = (int *) cCalloc(k,sizeof(int),"ijv_comp->j in st_metaGraph");
  ijv_comp->v = (double *) cCalloc(k,sizeof(double),"ijv_comp->v in st_metaGraph");
  for(e = 0; e< k; e++){
    ijv_comp->i[e] = ijv_orig->i[index[e]];
    ijv_comp->j[e] = ijv_orig->j[index[e]];
    ijv_comp->v[e] = value[e];
  }
  
  cFree(value);
  cFree(index);

  
  
  #ifdef ST_DEBUG
  fprintf(stderr,"333\n");
  printIJV(ijv_comp);
  fprintf(stderr,"\n");
  #endif


  return ijv_comp;
}

