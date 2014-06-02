// Calculate the effective resitance for given edge on graph
// assumes we have series parallel graph
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "graph.h"

const char* DOC_STRING = "run like: ./calcer GRAPH_IN.BIN_IJV v1 v2\n";

void die(int e)
{
  fprintf(stderr, DOC_STRING);
  exit(e);
}

bool isEdge(myGraph* graph, int v1, int v2)
{
  int nbrItr = 0;
  for (nbrItr = 0; nbrItr < graph->deg[v1]; nbrItr++)
    {
      if (graph->nbrs[v1][nbrItr] == v2)
	return true;
    }
  return false;
}


void printGraph(myGraph* graph)
{
  int tmp2 = 0;
  printf("\t");
  for (tmp2 = 0; tmp2 < graph->n; tmp2++)
    {
      int tmp3 = 0;
      printf("\t%d: ", tmp2);
      for (tmp3 = 0; tmp3 < graph->deg[tmp2]; tmp3++)
	{
	  // check back
	  int nbr = graph->nbrs[tmp2][tmp3];
	  int back = graph->back[tmp2][tmp3];
	  if (graph->nbrs[nbr][back] != tmp2)
	    printf("ERROR BACK BROKEN!!!!!\n");

	  printf("%d ", graph->nbrs[tmp2][tmp3]);
	}     
    }
  printf("\n");
}

void checkWts(myGraph* graph)
{
  int nodeItr = 0;
  for (nodeItr = 0; nodeItr < graph->n; nodeItr++)
    {
      int nbrItr = 0;
      for (nbrItr = 0; nbrItr < graph->deg[nodeItr]; nbrItr++)
	{
	  printf("weight: %lf ", graph->wts[nodeItr][nbrItr]);
	}     
    }
}



// take double pointer to allow us to replace graph with new one
// look for two edge that can be combined (a vert of deg 2 that can be removed or a double edge between two 
// verts that can be turned into 1)
void reduceGraph(myGraph** graph_p, int v1, int v2)
{

  // remember not to ever eliminate v1 or v2
  myGraph* graph = *graph_p;
  // iterate over every vertex (then adjacent edge) to eliminate edges
  int nodeItr = 0;
  for (nodeItr = 0; nodeItr < graph->n; nodeItr++)
    {
      printf("at node: %d with deg: %d\n", nodeItr, graph->deg[nodeItr]);
      // check for vert of deg 2
      if (graph->deg[nodeItr] == 2 && nodeItr != v1 && nodeItr != v2)
	{
	  printf("found vert of deg 2 (vert: %d between %d and %d) \n", nodeItr, graph->nbrs[nodeItr][0], graph->nbrs[nodeItr][1]);
	  printGraph(graph);

	  // elimanate the vertex
	  double newWt = graph->wts[nodeItr][0] + graph->wts[nodeItr][1];
	  int firstNode = graph->nbrs[nodeItr][0];
	  int firstBack = graph->back[nodeItr][0];
	  int secondNode = graph->nbrs[nodeItr][1];
	  int secondBack = graph->back[nodeItr][1];

	  // check that backs line up
	  if (graph->nbrs[firstNode][firstBack] != nodeItr
	      || graph->nbrs[secondNode][secondBack] != nodeItr)
	    printf("error with backs!\n");


	  graph->nbrs[firstNode][firstBack] = secondNode;
	  graph->wts[firstNode][firstBack] = newWt;
	  graph->back[firstNode][firstBack] = secondBack;
	  graph->nbrs[secondNode][secondBack] = firstNode;
	  graph->wts[secondNode][secondBack] = newWt;
	  graph->back[secondNode][secondBack] = firstBack;
	  graph->deg[nodeItr] = 0;
	  graph->nnz -= 2;
	  
	  printGraph(graph);

	  printf("done with vert of deg 2\n");
	}

      // look for a double edge
      int nbrItr = 0;
      for (nbrItr = 0; nbrItr < graph->deg[nodeItr]; nbrItr++)
	{
	  int nbrItr2 = 0;
	  for (nbrItr2 = nbrItr + 1; nbrItr2 < graph->deg[nodeItr]; nbrItr2++)
	  {	    
	    if (graph->nbrs[nodeItr][nbrItr] == graph->nbrs[nodeItr][nbrItr2])
	      {
		printf("found double edge verts (%d & %d) \n", nodeItr, graph->nbrs[nodeItr][nbrItr]);
		printGraph(graph);
		// move nbrItr2 intoItr
		double newWt = 1.0 / ((1.0 / graph->wts[nodeItr][nbrItr]) + (1.0 / graph->wts[nodeItr][nbrItr2]));
		int nbrNode = graph->nbrs[nodeItr][nbrItr];
		int firstBack = graph->back[nodeItr][nbrItr];
		int secondBack = graph->back[nodeItr][nbrItr2];
		graph->wts[nodeItr][nbrItr] = newWt;
		graph->wts[nbrNode][firstBack] = newWt;
		
		// fill in gaps
		int ctr = 0;
		for (ctr = nbrItr2 + 1; ctr < graph->deg[nodeItr]; ctr++)
		  {
		    graph->nbrs[nodeItr][ctr - 1] = graph->nbrs[nodeItr][ctr];
		    graph->wts[nodeItr][ctr - 1] = graph->wts[nodeItr][ctr];
		    graph->back[nodeItr][ctr - 1] = graph->back[nodeItr][ctr];
		    graph->back[graph->nbrs[nodeItr][ctr]][graph->back[nodeItr][ctr]]--;
									       
		  }
		graph->deg[nodeItr]--;
		for (ctr = secondBack + 1;  ctr < graph->deg[nbrNode]; ctr++)
		  {
		    graph->nbrs[nbrNode][ctr - 1] = graph->nbrs[nbrNode][ctr];
		    graph->wts[nbrNode][ctr - 1] = graph->wts[nbrNode][ctr];
		    graph->back[nbrNode][ctr - 1] = graph->back[nbrNode][ctr];		   
		    graph->back[graph->nbrs[nbrNode][ctr]][graph->back[nbrNode][ctr]]--;
		  }
		graph->deg[nbrNode]--;
		graph->nnz -= 2;
		printGraph(graph);
		printf("done with double edcge\n");
		break;
	      }
	  }
	}
    }
}

int main(int argc, char* argv[])
{
  if (argc != 4)
    die(1);

  printf("going to get graph\n");
  myGraph* graph = getGraph(argv[1]);
  printf("got graph\n");

  int v1, v2 = -1;
  sscanf(argv[2], "%d", &v1);
  sscanf(argv[3], "%d", &v2);
  
  // check that we were given an actual edge
  if (!isEdge(graph, v1, v2))
    die(1);
  
  checkWts(graph);

  int lastNNZ = -1;
  while(graph->nnz > 2)   // note each edge is counted twice (forwards/backwards)
    {
      printf("graph->nnz: %d\n", graph->nnz);     
      reduceGraph(&graph, v1, v2);
      if (lastNNZ == graph->nnz)
	{
	  fprintf(stderr, "ERROR: reduceGraph failed to do anything!\n");
	  exit(2);
	}
      lastNNZ = graph->nnz;
    }
  
  printf("The effective resistence is: %lf\n", graph->wts[v1][0]); // assume all other edges gone, so it has to be first listing!

  freeGraph(graph);
  
  return 0;
}
