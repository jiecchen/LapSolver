// use pseudo random trees to measure distance between nodes in graphs
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

#include "components.h"
#include "weightedselection.h"
#include "graph.h"
//#include "priorityq.h"

// Changed to use txtijv files for output: important to know the order in which edges were added

const char* DOC_STRING = "run like ./graphdist GRAPH_IN.BINIJV TREE_OUT.TXTIJV NUM_TREES\n";

typedef struct List_s
{
    int u;     // keep u smaller than v
    int v;
    struct List_s* next;
} List;



pArray* getPRandTree(ijvType* ijv)
{
  Components* comps = initGraphComponents(ijv->n);

  WeightedElements* we = newWeightedElements(ijv->nnz);   // the elements are the indexes to the edges in ijv
  int ctr = 0;
  for (ctr = 0; ctr < ijv->nnz; ctr++)
    {
      weInsert(we, ctr, 1.0/ijv->v[ctr]);
    }

  ijvType* treeOut = makeIJV(ijv->n - 1);
  treeOut->n = ijv->n;
  int treeEdgeCtr = 0;
    
  // randomly select edges until we have a spanning tree (only one comp)
  int numComponents = ijv->n;
  while (numComponents > 1)
    {
      int u = -1;
      int v = -1;
      double weight = 0.0;

      int randEdgeIdx = weSelect(we);
      u = ijv->i[randEdgeIdx];
      v = ijv->j[randEdgeIdx];
      weight = ijv->v[randEdgeIdx];
	
      int uComp = comps->vertComps[u];
      int vComp = comps->vertComps[v];
      if (uComp == vComp)
	{
	  continue;
	}
      // let u be the smaller component 
      if (comps->compSizes[vComp] < comps->compSizes[uComp])
	{
	  int tmp = u;
	  u = v;
	  v = tmp;
	  tmp = uComp;
	  uComp = vComp;
	  vComp = tmp;
	}

      // make uComp part of vComp
      CompCell* curUcell = comps->comps[uComp];
      while (1)
	{
	  comps->vertComps[curUcell->vert] = vComp;
	  if (curUcell->next)
	    curUcell = curUcell->next;
	  else
	    break;
	}
      curUcell->next = comps->comps[vComp];
      comps->comps[vComp] = comps->comps[uComp];
      comps->comps[uComp] = NULL;
      comps->compSizes[vComp] += comps->compSizes[uComp];
      comps->compSizes[uComp] = 0;

      // add edge to tree
      treeOut->i[treeEdgeCtr] = u;
      treeOut->j[treeEdgeCtr] = v;
      treeOut->v[treeEdgeCtr++] = weight;
      numComponents--;
    }
  //myGraph* treeOutGraph = ijv2graph(treeOut);
  //pArray* parray = graph2pArray(treeOutGraph);
  
  //freeGraph(treeOutGraph);
  //freeIJV(treeOut);

  freeWeightedElements(we);
  freeComponents(comps);
  return treeOut;
}


int main(int argc, char* argv[])
{
    srand(time(NULL));
    if (argc != 4)
	{
	    fprintf(stderr, "%s", DOC_STRING);
	    exit(1);
	}

    int numTrees = atoi(argv[3]);

    FILE* fpIn = NULL;
    ijvType* ijv = NULL;
    //myGraph* graph = NULL;
    if ((fpIn = fopen(argv[1], "rt")) == NULL)
        {
	    fprintf(stderr, "Error: could not read %s\n", argv[1]);
	    exit(1);
	}
    FILE* fpOut = NULL;
    if ((fpOut = fopen(argv[2], "wt")) == NULL)
      {
	fprintf(stderr, "Error: could not open %s\n", argv[2]);
	exit(1);
      }

    fprintf(fpOut, "%d\n", numTrees);

    ijv = binReadIJV(fpIn);
    fclose(fpIn);

    int i = 0;
    for (i = 0; i < numTrees; i++)
      {
	//fprintf(fpOut, "\n");
	ijvType* treeOut = getPRandTree(ijv);    
	textWriteIJV(treeOut, fpOut);
	freeIJV(treeOut);
      }
    freeIJV(ijv);

    fclose(fpOut);

    printf("exiting at %lf seconds\n", (double)clock()/(double)CLOCKS_PER_SEC);
    
    return 0;
}
