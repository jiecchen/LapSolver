/*******************
    st_components.c

    compile with mex st_components.c or make

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


/********************************************************
  *   st_components part
  *******************************************************/


void DFS_Visit(myGraph *G, int col, int * color, int comNumber, int *comIndex){
   int i;
   int row;
   color[col] = 1;
   for (i = 0; i < G->deg[col]; i++) {
       row = G->nbrs[col][i];
       if (color[row] == 0){
         comIndex[row] = comNumber;
         color[row] = 1;
         DFS_Visit(G,row,color,comNumber,comIndex);
       }
   }
   color[col] = 2;
}


void st_components(myGraph *G, int *comIndex)
{
   int col, row;
   int ind, comNumber;
   int *color;

   color = (int *) cCalloc(G->n, sizeof(int), "color in st_components");
   for (col = 0; col < G->n ; col++) {
     color[col] = 0;
   }

/*  color: 0 while, 1 gray, 2 black
  * I used the DFS language for CLRS just in case
  *   in the future we need a full version of DFS;
  */

   comNumber = 0;
   for (col = 0; col < G->n ; col++) {
     if (color[col] == 0) {
        comIndex[col] = comNumber;
        DFS_Visit(G,col,color,comNumber,comIndex);
        comNumber++;
     }
   }

   cFree(color);

}


