/* readAndDumpIJV.c
 
   just exists to test our ability to read in an ijv

   compile it with:

   gcc readAndDumpIJV.c st_basic.c -o readAndDumpIJV
   
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "st_defs.h"

ijvType *makeIJV(int nnz);
ijvType *binReadIJV(FILE* fpr) ;


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
  fprintf(stderr, "Usage: readAndDumpIJV graphin\n");
}


main(int argc, char **argv)
{
  
 char inName[256];


 if (argc < 2) {
   usage();
   exit(1);
 }

 sprintf(inName,"%s",argv[1]);

 

 FILE *fpr;
 
 if ((fpr = fopen(inName,"r")) == NULL)
   cError (mystrcat(inName," <- can't create this file\n"));

 ijvType *ijv;
 ijv = binReadIJV(fpr);

 fclose(fpr);

 printIJV(ijv);

 
}


