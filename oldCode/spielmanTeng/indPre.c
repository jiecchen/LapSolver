/* generates an error: ./indPre ../jan10.gr  16352 jan10.pre
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>

#include "st_defs.h"

ijvType *makeIJV(int nnz);
ijvType *st_part(ijvType *ijv, int *opts);


/* creates a new string by allocating space and concatentating s1 and s2 */
char *mystrcat (char *s1, char *s2)
{
    char *s;
    char *p;

    p = s = (char *) malloc ((strlen(s1)+strlen(s2)+1)*sizeof(char));

    while (*s1 != '\0') {
	*p = *s1;
	p++; s1++;
    }
    while (*s2 != '\0') {
	*p = *s2;
	p++; s2++;
    }
    *p = '\0';

    return s;
}


void usage()
{
  fprintf(stderr, "Usage: indPre graphin k graphout\n");
  fprintf(stderr, "        [-tree]: just produce a tree\n");
  fprintf(stderr, "        [-debug]: print all sorts of stuff\n");
  fprintf(stderr, "        [-metisGreedy]: a way for metis to refine boundaries \n");
  fprintf(stderr, "        [-metisRand]: a way for metis to refine boundaries \n");
  fprintf(stderr, "        [-metisMin]: the default way for metis to refine boundaries \n");
  fprintf(stderr, "        [-metisBal f]: cut with balance f >= 1 \n");
  fprintf(stderr, "        [-prim k]: precondition metagraph with k max spanning trees \n");
  fprintf(stderr, "        [-minDeg j]: metagraph sparsifier must have min deg j\n");
  fprintf(stderr, "        [-minFracDeg x]: metagraph sparsifier must have min frac deg x\n");
  fprintf(stderr, "        [-minFracWt x]: metagraph sparsifier must have min frac wt x\n");
  fprintf(stderr, "        [-data]: output data on sparsification success\n");
  fprintf(stderr, "        [-time]: print timing data\n");
  fprintf(stderr, "        [-time2]: print other timing data\n");
}


main(int argc, char **argv)
{
  
 char inName[256];
 char outName[256];


 int k;
 
 if (argc < 4) {
   usage();
   exit(1);
 }

 sprintf(inName,"%s",argv[1]);
 k = atoi(argv[2]);

 
 sprintf(outName,"%s",argv[3]);
 
 int *opts;
 opts = (int *) partPreOpts(k);

 
 /* default values */
 int prim_trees = 3;

 
 int i;


 /* set this default */
 metisBalance = 1;

 for (i=4; i<argc; i++) {


   if (!strcmp(argv[i],"-tree")) 
     opts[stOptNumTree] = stOptValTreeYes;
   if (!strcmp(argv[i],"-debug"))
     opts[stOptNumDebug] = 1;
   if (!strcmp(argv[i],"-metisGreedy"))
     opts[stOptNumMetis] = metisStrategyGreedy;
   if (!strcmp(argv[i],"-metisRand"))
     opts[stOptNumMetis] = metisStrategyRandom;
   if (!strcmp(argv[i],"-metisMin"))
     opts[stOptNumMetis] = metisStrategyMinConn;

   if (!strcmp(argv[i],"-time"))
     opts[stOptNumTime] = 1;
   
   if (!strcmp(argv[i],"-time2"))
     opts[stOptNumTime2] = 1;
   
   
   if (!strcmp(argv[i],"-data"))
     opts[stOptNumData] = 1;
   
   if (!strcmp(argv[i],"-prim") && i <= argc-1) {
      i++;
      if (sscanf(argv[i],"%d",&prim_trees) != 1) {
	fprintf(stderr, "numTrees follows -prim argument\n");
	exit(1);
      }
      opts[stOptNumSpars] = stOptValSparsPrim;
      opts[stOptNumSparsParam] = prim_trees;
    }

   if (!strcmp(argv[i],"-minDeg") && i <= argc-1) {
     i++;
     if (sscanf(argv[i],"%d",&(opts[stOptNumSparsMinDeg])) != 1) {
       fprintf(stderr, "integer follows -minDeg argument\n");
       exit(1);
     }
     if (opts[stOptNumSpars] == stOptValSparsNo)
       opts[stOptNumSpars] = stOptValSparsMin;
   }

   if (!strcmp(argv[i],"-minFracDeg") && i <= argc-1) {
     i++;
     if (sscanf(argv[i],"%g",&(sparsMinFracDeg)) != 1) {
       fprintf(stderr, "float follows -minFracDeg argument\n");
       exit(1);
     }
     opts[stOptNumSparsMinFracDeg] = 1;
     if (opts[stOptNumSpars] == stOptValSparsNo)
       opts[stOptNumSpars] = stOptValSparsMin;
   }

   if (!strcmp(argv[i],"-minFracWt") && i <= argc-1) {
     i++;
     if (sscanf(argv[i],"%g",&(sparsMinFracWt)) != 1) {
       fprintf(stderr, "float follows -minFracWt argument\n");
       exit(1);
     }
     opts[stOptNumSparsMinFracWt] = 1;
     if (opts[stOptNumSpars] == stOptValSparsNo)
       opts[stOptNumSpars] = stOptValSparsMin;
   }

   if (!strcmp(argv[i],"-metisBal") && i <= argc-1) {
     i++;
     if (sscanf(argv[i],"%g",&(metisBalance)) != 1) {
       fprintf(stderr, "float follows -metisBal argument\n");
       exit(1);
     }
   }

   
 }

 if (opts[stOptNumSparsMinFracDeg] == 0)
   sparsMinFracDeg = 0;
 if (opts[stOptNumSparsMinFracWt] == 0)
   sparsMinFracWt = 0;
 
 printf("opts: ");
 for (i = 0; i < numOpts; i++)
   printf("%d ",opts[i]);
 printf("\n");

 FILE *fpr;
 
 if ((fpr = fopen(inName,"r")) == NULL)
   cError (mystrcat(inName," <- can't create this file\n"));

 int n, nnz;
 fscanf(fpr,"%d %d\n",&n, &nnz);

 fprintf(stderr,"n: %d, nnz: %d\n",n,nnz);

 ijvType *ijv;
 ijv = makeIJV(nnz);
 ijv->n = n;

 int x;
 int j;
 float v;
 for (x = 0; x < nnz; x++) {
   fscanf(fpr,"%d %d %g\n",&i,&j,&v);
   ijv->i[x] = i-1;
   ijv->j[x] = j-1;
   ijv->v[x] = (double) v;
   /*   fprintf(stderr,"%d, %d : %g\n",i,j,v); */
 }
 fclose(fpr);

 ijvType *ijv_out;

 ijv_out = st_part(ijv,opts);

 FILE *fpw;
 
 if ((fpw = fopen(outName,"w")) == NULL)
   cError (mystrcat(outName," <- can't create this file\n"));
 for (x = 0; x < ijv_out->nnz; x++) {
   fprintf(fpw,"%d %d %g\n",ijv_out->i[x],ijv_out->j[x],ijv_out->v[x]);
 }
 fclose(fpw);

 
}


