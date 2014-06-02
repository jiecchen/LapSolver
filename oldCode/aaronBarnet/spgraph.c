// create a simple graph to output to file "spgraph.binijv"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "graph.h"

int main(int argc, char* argv[])
{
  myGraph* g = newGraph(4, 5);
  int edgeCtr = 0;
  g->deg[0] = 3;
  g->nbrs[0] = g->nbrsBlock + edgeCtr;
  g->wts[0] = g->wtsBlock + edgeCtr;
  g->nbrs[0][0] = 1; edgeCtr++;
  g->nbrs[0][1] = 2; edgeCtr++;
  g->nbrs[0][2] = 3; edgeCtr++;

  g->deg[1] = 2;
  g->nbrs[1] = g->nbrsBlock + edgeCtr;
  g->wts[1] = g->wtsBlock + edgeCtr;
  g->nbrs[1][0] = 0; edgeCtr++;
  g->nbrs[1][1] = 2; edgeCtr++;

  g->deg[2] = 3;
  g->nbrs[2] = g->nbrsBlock + edgeCtr;
  g->wts[2] = g->wtsBlock + edgeCtr;
  g->nbrs[2][0] = 1; edgeCtr++;
  g->nbrs[2][1] = 0; edgeCtr++;
  g->nbrs[2][2] = 3; edgeCtr++;

  g->deg[3] = 2;
  g->nbrs[3] = g->nbrsBlock + edgeCtr;
  g->wts[3] = g->wtsBlock + edgeCtr;
  g->nbrs[3][0] = 0; edgeCtr++;
  g->nbrs[3][1] = 2; edgeCtr++;

  // make all wts 1.0
  int nodeItr = 0;
  for (nodeItr = 0; nodeItr < g->n; nodeItr++)
    {
      int nbrItr = 0;
      for (nbrItr = 0; nbrItr < g->deg[nodeItr]; nbrItr++)
	{
	  g->wts[nodeItr][nbrItr] = 1.0;
	}
    }

  ijvType* ijv = makeIJV(g->nnz);
  ijv->n = g->n;
  int ijvItr = 0;
  for (nodeItr = 0; nodeItr < g->n; nodeItr++)
    {
      int nbrItr = 0;
      for (nbrItr = 0; nbrItr < g->deg[nodeItr]; nbrItr++)
	{
	  int nbrNode = g->nbrs[nodeItr][nbrItr];	  
	  if (nbrNode < nodeItr)
	    continue;    // avoids writing down edge twice!
	  ijv->i[ijvItr] = nodeItr;
	  ijv->j[ijvItr] = nbrNode;
	  ijv->v[ijvItr] = 1.0; //g->wts[nodeItr][nbrItr];
	  ijvItr++;
	}
    }

  
  // write it all out
  binWriteIJV(ijv, "spgraph.binijv");
  
  freeGraph(g);
  freeIJV(ijv);
  return 0;
}
