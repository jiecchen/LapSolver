/*******************
   Prim.c

   compile with mex treedist.c

*******************/

#include <stdio.h>
#include <tgmath.h>
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

#define DEBUG_FLAG 0
#define RAND_MAX    2147483647      /*  2^31 - 1 */
#define SEED 95012928

/********************************************************
 *   function headers
 *******************************************************/

void freeGraph(myGraph *G);
void printGraph(myGraph *G);
ijvType  *makeIJV(int nnz);
int rand();
void *cCalloc(size_t n, size_t size, char *errMsg);

/******************************************************
 *  Heap operations and Prim's algorithm for maximum spanning
 *  tree
 * 
 ******************************************************/

  struct heapElement{
    double value;
    int vertex;
    int time;
  };

  struct Heap {
    struct heapElement * heap;
    int *index;
    int size;
  };


/* Return a random integer from 0 -- n */
int st_random(int n){  
 int r;
 int rnum;
 r = rand();
 // rnum = floor( (double) r / RAND_MAX*n );
 rnum = ( (double) r / RAND_MAX*n );
 if (rnum == n)
   rnum--;
 return(rnum);
}

int heap_ElementComapre(struct heapElement *h1, struct heapElement *h2 ){
  if ( (h1->value > h2->value) || (h1->value == h2->value) && (h1->time < h2->time) ){
    return 1;
  }else{
    return -1;
  }
}


heap_Update(struct Heap *H, int vertex, double value){
  int p, q;
  struct heapElement temp;

  p = H->index[vertex];
  H->heap[p].value = value;
  //q = floor((p-1)/2);
  q = ((p-1)/2);
  while ( (p > 0) && (heap_ElementComapre(&H->heap[p], &H->heap[q]) > 0)) {
    H->index[H->heap[p].vertex] = q;
    H->index[H->heap[q].vertex] = p;
    temp.value = H->heap[p].value;
    temp.vertex = H->heap[p].vertex;
    temp.time = H->heap[p].time;
    H->heap[p].value = H->heap[q].value;
    H->heap[p].vertex = H->heap[q].vertex;
    H->heap[p].time = H->heap[q].time;
    H->heap[q].value = temp.value;
    H->heap[q].vertex = temp.vertex;
    H->heap[q].time = temp.time;
    p = q;
    //q = floor((p-1)/2);
    q = ((p-1)/2);
  } 

}


int heap_ExtractMax(struct Heap *H){
  int i;
  int u, p, c, completed, tail;
  struct heapElement temp;
  int DONE = -11;
  if (H->size >0){
    u = H->heap[0].vertex;
    H->index[u] = DONE;
/*     fprintf(stderr,"%d is done\n", u); */
    H->size -= 1;
    tail = H->size;
    if (tail > 0){
      H->index[H->heap[tail].vertex] = 0;
      H->heap[0].value = H->heap[tail].value;
      H->heap[0].vertex = H->heap[tail].vertex;
      H->heap[0].time = H->heap[tail].time;
    }

    p = 0;
    completed = 0;
    while ((2*p + 1 < tail) && (completed == 0)){
      c =2*p+1;
      if (c+1< tail){
        if ((heap_ElementComapre(&H->heap[c+1], &H->heap[c]) > 0)){
          c++;
	}
      }
      if ( heap_ElementComapre(&H->heap[c], &H->heap[p]) > 0){
         H->index[H->heap[p].vertex] = c;
         H->index[H->heap[c].vertex] = p;

         temp.value = H->heap[p].value;
         temp.vertex = H->heap[p].vertex;
         temp.time = H->heap[p].time;
         H->heap[p].value = H->heap[c].value;
         H->heap[p].vertex = H->heap[c].vertex;
         H->heap[p].time = H->heap[c].time;
 
         H->heap[c].value = temp.value;
         H->heap[c].vertex = temp.vertex;
         H->heap[c].time = temp.time;
         p = c;
      }else{
	completed = 1;
      }
    }
    return u;
  }else {
    return -1;
  }
}

heap_Insert(struct Heap *H, double value, int vertex, int time){
  int size, i;
  size = H->size++;
  H->index[vertex] = size;
  H->heap[size].vertex = vertex;

  H->heap[size].time = time;

/*
fprintf(stderr,"after insert\n"); 
for (i=0; i< H->size; i++){ 
  fprintf(stderr,"heap[%d]  = %d, %f, %d\n", i, H->heap[i].vertex, H->heap[i].value, H->heap[i].time); 
} 

*/
  heap_Update(H,vertex,value);

/* fprintf(stderr,"after update\n"); 
for (i=0; i< H->size; i++){ 
  fprintf(stderr,"heap[%d]  = %d, %f, %d\n", i, H->heap[i].vertex, H->heap[i].value, H->heap[i].time); 
} 
 */
}


/* Construct the Prim Maximum Spanning Tree from source s */

int * Prim (myGraph *G, int s){
  struct Heap * H;
/*  Note: H->index[u] = -2 (done), (-1) (new),  >-1 in the heap */

  int NEW = -2;
  int DONE = -11;
  int MASKEDOUT = -17;

  int * parent;
  int n,k,i,t;
  int *nbrs, degree, u, v;
  int time;

  n = G->n;
  H = (struct Heap *) cCalloc(1, sizeof(struct Heap),"H in Prim");
  H->heap = (struct heapElement *) cCalloc(n, sizeof(struct heapElement),"H->heap in Prim");
  H->index = (int *) cCalloc(n, sizeof(int),"H->hindex in Prim");
  parent = (int *) cCalloc(2*n, sizeof(int),"parent in Prim");
  
/* Initialize   */
  H->size = 0;
  time = 1;
  for (i=0;i < n; i++){
    H->index[i] = NEW;
    parent[i] = -1;
  }


  for (i=-1;i<n;i++){
    if (i == -1) {
      heap_Insert(H, 0, s,time++);
    } else{
      if (H->index[i]== NEW){
        heap_Insert(H, 0, i,time++);
      }
    }

/* Tree construction, if the graph is not connected, it will produce
 *  a tree for each component
 */

    while (H->size > 0){

/* fprintf(stderr,"to start\n");
for (t=0; t< n; t++){ 
  fprintf(stderr,"index[%d]  = %d\n", t, H->index[t]);
}
 */

      u = heap_ExtractMax(H);

/* fprintf(stderr,"after extract %d \n",u);
for (t=0; t< n; t++){ 
  fprintf(stderr,"index[%d]  = %d\n", t, H->index[t]);
}
*/
/*       fprintf(stderr,"processing %d\n", u); */


      for (k = 0; k < G->deg[u]; k++){
        v = G->nbrs[u][k];
        if ( (H->index[v]!= DONE)){
/* 	     fprintf(stderr,"%d  is not done yet\n", v); */
          if (G->back[u][k] != MASKEDOUT){
            if (H->index[v] == NEW){
       	      parent[v] = u;
	      parent[n+v] = k;
/*  	      fprintf(stderr,"make %d the parent of  %d\n", u,v); */
              heap_Insert(H,G->wts[u][k],v,time++);
/* fprintf(stderr,"after insert %d \n",v);
for (t=0; t< n; t++){ 
  fprintf(stderr,"index[%d]  = %d\n", t, H->index[t]);
}
 */
	    } else{ 
	      if (G->wts[u][k] > H->heap[H->index[v]].value){
        	 parent[v] = u;
   	         parent[n+v] = k;
/* 	         fprintf(stderr,"make %d the parent of  %d\n", u,v); */
                 heap_Update(H,v,G->wts[u][k]);
/* fprintf(stderr,"after update %d \n",v);
for (t=0; t< n; t++){ 
  fprintf(stderr,"index[%d]  = %d\n", t, H->index[t]);
}
 */

	      }
	    }
	  }else{
/* 	    fprintf(stderr,"edge (%d,%d) is used in the previous tree\n",u,v); */
	  }
	}
      }
    }
  }

/* Free dynamic variables */
  free(H->heap);
  free(H->index);
  free(H);
  return parent;
}


ijvType * multi_Prim (myGraph *G, int k){ 
  int MASKEDOUT = -17;
  int i, n, itrs, s, nnz;
  int *backBlockSave;
  int * parent;
  ijvType * multi_PrimTrees;
  int j;
 
  backBlockSave = (int *) cCalloc(G->nnz, sizeof(int),"backBlockSave in multi_Prim");
  memcpy(backBlockSave,G->backBlock, G->nnz*sizeof(int));


  multi_PrimTrees = makeIJV( k*(G->n - 1) );
  multi_PrimTrees->n = G->n;
  n = G->n;
  nnz = 0;
  srand(SEED);
/*   for (i=0; i< 10; i++){ 
     fprintf(stderr,"random [%d] = %d\n", i, st_random(G->n)); 
   } */

  for (itrs=0; itrs< k; itrs++){
    s = st_random(n);
/*     fprintf(stderr,"starting at %d\n", s); */
    parent = Prim(G,s);

/*     for (j=0; j< n; j++){ */
/*       fprintf(stderr,"parent[%d] = %d\n", j, parent[j]); */
/*     } */


    for (i = 0; i< G->n; i++){
      if (parent[i] != -1){ 
        G->back[i][G->back[parent[i]][parent[n+i]]] = MASKEDOUT;
        G->back[parent[i]][parent[n+i]] = MASKEDOUT;
	if (parent[i] > i){
           multi_PrimTrees->i[nnz] = parent[i];
           multi_PrimTrees->j[nnz] = i;
   	   multi_PrimTrees->v[nnz] = G->wts[parent[i]][parent[n+i]];
	}else{
           multi_PrimTrees->i[nnz] = i;
           multi_PrimTrees->j[nnz] = parent[i];
	   multi_PrimTrees->v[nnz] = G->wts[parent[i]][parent[n+i]];
	}
        nnz++;
      }
    }
    free(parent);
  }
  multi_PrimTrees->nnz = nnz;
  
  memcpy(G->backBlock, backBlockSave, G->nnz*sizeof(int));

  free(backBlockSave);
  return multi_PrimTrees;
};



struct WT_Index {
  double wt;
  int index;
};

static int compare(a,b)
char *a, *b;
{
  struct WT_Index * x = ( struct WT_Index *) a;
  struct WT_Index * y = ( struct WT_Index *) b;
  if (y->wt> x->wt){
    return(1);
  }else{
    return (-1);
  }
}


void insert_ijv( ijvType ** ijv, int i,  int j, double  v, int * capacity, int nnz ){
  ijvType *current;
  int k;

  current = *ijv;
  if (nnz < *capacity){
    current->i[nnz] = i;
    current->j[nnz] = j;
    current->v[nnz] = v;
  }else{
    *capacity = 2* (*capacity);

    
    ijvType *new;
    new = makeIJV( *capacity );
    for (k = 0; k < nnz; k++){
      new->i[k] = current->i[k];
      new->j[k] = current->j[k];
      new->v[k] = current->v[k];
    }

    
    new->i[nnz] = i;
    new->j[nnz] = j;
    new->v[nnz] = v;
    new->n = current->n;
    *ijv = new;
    freeIJV(current);
  }
}


ijvType * vertexWise_Sampling(myGraph *G, int min_deg, float minfrac_deg, float minfrac_wt){ 

  ijvType * sample;
  struct WT_Index * weights;
  int row, col,n, nnz;
  double swts, total_wt;
  int sdegs;
  int capacity;
  
  if ( (double)  (min_deg*G->n) > ((double) G->nnz)* ((double ) minfrac_deg)){
    capacity = min_deg*G->n; 
  } else{
    capacity = (int) ((double) G->nnz)* ((double ) minfrac_deg);
  }
  
  if (capacity < 2*G->n){
    capacity = 2*G->n;
  }
  sample = makeIJV( capacity );
  sample->n = G->n;

  nnz = 0;
  for (row = 0; row< G->n; row++){
    weights = (struct WT_Index *) cCalloc(G->deg[row], sizeof(struct WT_Index),"weights in vertexWise_Sampling");
    total_wt = 0.0;
    for (col = 0; col< G->deg[row]; col++){
      weights[col].wt = G->wts[row][col];
      weights[col].index = G->nbrs[row][col];
      total_wt+= G->wts[row][col];
    }
    qsort( (char *) weights, G->deg[row], sizeof(struct WT_Index), compare);

/*
int i;
 fprintf(stderr, "row = %d \n",row);
for (i = 0; i< G->deg[row]; i++){
  fprintf(stderr, "wts[%d] = %f \n", i, weights[i].wt);
}
fprintf(stderr, "\n");
*/
    sdegs= 0;
    swts = 0.0;
    for (col = 0; col<G->deg[row]; col++){
      if ( (sdegs < min_deg) || ( (double) sdegs < (double) minfrac_deg * (double) G->deg[row]) 
	   ||  (swts < (double) minfrac_wt*total_wt) ){ 
        if (row > weights[col].index){
  	  insert_ijv( &sample, row, weights[col].index, weights[col].wt, &capacity, nnz );
	}else{
  	  insert_ijv( &sample, weights[col].index, row, weights[col].wt, &capacity, nnz );
	}
	/*
	  fprintf(stderr,"added: %d - %d, %d : %g\n",nnz,row, weights[col].index, weights[col].wt);
	*/

        nnz++;
        sdegs++;
        swts+= weights[col].wt;
      }
    }
    cFree(weights);
    
  }

  sample->nnz = nnz;

  return sample;
};
