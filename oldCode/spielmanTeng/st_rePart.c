/*******************
    st_rePart.c

    given a partition of a graph,
    it breaks each component into
    its connected components

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

partType *compVec2part(int *compVec, int n);
ijvBlockDiag *formDiagBlocks(ijvType *ijv, partType *p);
myGraph *ijv2graph(ijvType *ijv);




ijvBlockDiag *st_rePart(ijvType *ijv, int *compVec)
{
  int x;
  int i;
  int j;
  
  int n;
  int *comp_new;
  int *comp_this;

  int curr, nextcurr;
  
  ijvBlockDiag *bd_orig, *bd_new;
  partType *p_orig, *p_new;

  myGraph *G;
  
  n = ijv->n;

  /* is freed with bd_orig */
  p_orig = compVec2part(compVec, n);

  bd_orig = formDiagBlocks(ijv,p_orig);

  /* do not free: goes into bd_new */
  p_new = (partType *) cCalloc(1,sizeof(partType),"p_new in st_rePart");
  
  comp_new = (int *) cCalloc(n, sizeof(int),"comp_new in st_rePart");
  for (x = 0; x < n; x++)
    comp_new[x] = -1;
  

  /* the current number of components */
  curr = p_orig->numParts;

  /* enum over the original blocks,
     to set up comp_new */
  for (j = 0; j < p_orig->numParts; j++) {

    /* only do non-empty blocks */
    if (bd_orig->part->partSizes[j] > 0) {
    
      /* convert to a graph format */
      G = ijv2graph(bd_orig->blocks[j]);
      
      /* gets the output of components */
      comp_this = (int *) cCalloc(G->n, sizeof(int),"comp_this in st_rePart");
      st_componentsIter(G, comp_this);
      
      nextcurr = curr;
      /* re-map by these */
      for (i = 0; i < G->n ; i++) 
	if (comp_this[i] == 0)
	  comp_new[p_orig->partMap[j][i]] = j;
	else {
	  comp_new[p_orig->partMap[j][i]] = curr-1 + comp_this[i];
	  if (curr-1 + comp_this[i] >= nextcurr)
	    nextcurr = curr-1 + comp_this[i] + 1;
	}

      curr = nextcurr;
      
      freeGraph(G);

      cFree(comp_this);
    }
  }


  freeIjvBlockDiagAndPart(bd_orig); /*p_orig goes with it */

  
  /* now, create the new block graph */
  /* CAN SAVE A FACTOR OF 2 IN TIME HERE!!! 
     BY GOING THROUGH BD2IJV FIRST */

  
  /* FOR DEBUGGING FOR NOW */
  /* make sure all entries assigned */
  for (x = 0; x < n; x++)
    if (comp_new[x] == -1) {
      fprintf(stderr,"error in comp_new\n");
      return NULL;
    }

  
  p_new = compVec2part(comp_new, n);
  cFree(comp_new);

  bd_new = formDiagBlocks(ijv,p_new);


  /* check if this is correct */

  /*
  int y;
  for (y = 0; y < bd_new->part->numParts; y++) {
    G = ijv2graph(bd_new->blocks[y]);
    if (bd_new->part->partSizes[y] != G->n)
      fprintf(stderr,"disc: %d\n",G->n);
    comp_this = (int *) cCalloc(G->n, sizeof(int),"comp_this in st_rePart");

    st_components(G, comp_this);
    freeGraph(G);

    for (x = 0; x < bd_new->part->partSizes[y]; x++)
      if (comp_this[x] > 0) {
	fprintf(stderr,"error in st_rePart:\n");
	printIJV(ijv);
	printIntVec(compVec,ijv->n, "compVec");
	return NULL;
      }
    cFree (comp_this);
  }
  */


  return bd_new;
  
}


