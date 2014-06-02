#ifndef GRAPH_H
#define GRAPH_H
// c functions/structs for dealing with graphs (from oldc folder)
#include <stdlib.h>
#include <stdio.h>

typedef struct myGrapht {
  int n;
  int nnz;

  int *deg; /* up to n */
  int **nbrs; /* up to n */
  int *nbrsBlock; /* holds the nbrs, nnz */
  double **wts; /* up to n */
  double *wtsBlock; /* up to nnz (has duplicates) */
  int *backBlock; /* up to nnz */
  int **back; /* up to n */

} myGraph;

typedef struct ijvTypet {
  int n;
  int nnz;
  int *i;
  int *j;
  double *v;
} ijvType;

typedef struct pArray_s {
    int n;
    int* array;    
} pArray;


void cError(const char* msg);
void *cCalloc(size_t n, size_t size, char *errMsg);
void cFree(void *obj);
void freePArray(pArray* array);
pArray* binReadPArray(const char* fileName);
pArray* graph2pArray(myGraph* graph);
pArray* newPArray(int n);
ijvType *makeIJV(int nnz);
ijvType *binReadIJV(FILE* fpr);
ijvType *binReadIJVcosts(FILE* fpr);
void binWriteIJV(ijvType* ijv, const char* fileName); // will modify ijv!
void textWriteIJV(ijvType* ijv, FILE* fp); // will modify ijv!
void freeIJV(ijvType *ijv);
myGraph *ijv2graph(ijvType *ijv);
myGraph* newGraph(int n, int nnz);
void freeGraph(myGraph *G);
myGraph* getGraph(const char* fileName);  // file is binIJV
void binWritePArray(pArray* parray, const char* fileName); // vert ints must be zero based (this function will increment them)
void textWritePArray(pArray* parray, FILE* fp); // vert ints must be zero based (this function will increment them)
// pArray is modified!
#endif  // GRAPH_H
