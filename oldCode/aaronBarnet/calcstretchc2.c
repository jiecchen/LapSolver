// calculate stretch of a graph

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "graph.h"
#include "calcstretchctree.h"

const char DOC_STRING[] = "run like ./calcstretchc2 GRAPH.BINIJV TREE.PARRAY\n";

int main(int argc, char* argv[])
{

    FILE* graph_fp = NULL;
    FILE* tree_fp = NULL;
    ijvType* ijv = NULL;
    myGraph* graph = NULL;    
    //myGraph* treeGraph = NULL;
    pArray* parray = NULL;
    Tree* tree = NULL;
    
    double stretch = 0.0;


    if (argc < 3)
    {
	printf("%s", DOC_STRING);
	exit(1);
    }
    // deal with graph
    if ((graph_fp = fopen(argv[1], "rt")) == NULL)
	exit(1);
    ijv = binReadIJV(graph_fp);
    fclose(graph_fp);
    graph = ijv2graph(ijv);
    freeIJV(ijv);
    
    // deal with tree
    parray = binReadPArray(argv[2]);
    tree = pArray2tree(parray, graph);
    freePArray(parray);

    printf("finished getting graph and tree at %lf seconds\n", (double)clock()/(double)CLOCKS_PER_SEC);
    
    //checkTree(tree);
    stretch = calcStretch(tree, graph);
    printf("the stretch is: %f, the average stretch is: %f\n", stretch, stretch/graph->nnz);
    
    freeGraph(graph);
    freeTree(tree);

    printf("exiting at %lf seconds\n", (double)clock()/(double)CLOCKS_PER_SEC);

    return 0;
    

}
