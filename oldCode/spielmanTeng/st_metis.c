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



dubIjvType *doubleIJV(const ijvType *ijvIn);
int ijvCheckSortj(ijvType *ijv);
void ijvSortj(ijvType *ijv);
void ijvSorti(ijvType *ijv);


/**
 * \defgroup metis
 * 
 */

/**
 * \ingroup metis
 * \todo take out the sort, as it above
*/
void ijv2metis(const ijvType *ijv_in,
	       int **xadj, int **adjncy, int **adjwgt)
{
  int n, nnz;
  
  int col, i;
  int ind;

  dubIjvType *ijv;

  ijv = doubleIJV(ijv_in);
  
  n = ijv->n;
  nnz = ijv->nnz;

  ijvSorti(ijv);
  ijvSortj(ijv);

  /*
  if (ijvCheckSortj(ijv) == 0) {
    ijvSortj(ijv);
  }
  */

  
  *xadj = (int *) cCalloc(n+1, sizeof(int),"xadj in ijv2metis");
  *adjncy = (int *) cCalloc(nnz, sizeof(int),"adjncy in ijv2metis");
  *adjwgt = (int *) cCalloc(nnz, sizeof(int),"adjwgt in ijv2metis");

  ind = 0;
  col = 0;

  for (col = 0; col < n; col++) {
    (*xadj)[col] = ind;
    
    while ((ind < nnz) && (ijv->j[ind] == col))
      ind++;
  }
  (*xadj)[n] = nnz;


  for (i = 0; i < nnz; i++) {
    (*adjncy)[i] = ijv->i[i];
    (*adjwgt)[i] = (int) ijv->v[i];
  }
}

/**
 * is for call to metisV.
 * \ingroup metis
 * \todo take out the sort, as it above
*/
void ijv2metisV(const ijvType *ijv_in,
		int **xadj, int **adjncy, int **vwgt, int **vsize)
{
  int n, nnz;
  
  int col, i;
  int ind;

  dubIjvType *ijv;

  ijv = doubleIJV(ijv_in);
  
  n = ijv->n;
  nnz = ijv->nnz;

  ijvSorti(ijv);
  ijvSortj(ijv);

  /*
  if (ijvCheckSortj(ijv) == 0) {
    ijvSortj(ijv);
  }
  */

  
  *xadj = (int *) cCalloc(n+1, sizeof(int),"xadj in ijv2metis");
  *adjncy = (int *) cCalloc(nnz, sizeof(int),"adjncy in ijv2metis");
  *vsize = (int *) cCalloc(n+1, sizeof(int),"vsize in ijv2metis");
  *vwgt = (int *) cCalloc(n+1, sizeof(int),"vsize in ijv2metis");

  ind = 0;
  col = 0;

  for (col = 0; col < n; col++) {
    (*xadj)[col] = ind;
    
    while ((ind < nnz) && (ijv->j[ind] == col))
      ind++;
  }
  (*xadj)[n] = nnz;


  for (i = 0; i < nnz; i++) {
    (*adjncy)[i] = ijv->i[i];
  }

  for (i = 0; i <= n; i++) {
    (*vsize)[i] = 1;
    (*vwgt)[i] = 1;
  }

}


/**
 *
 * \brief k is number of parts
 *  
 * the output goes into comp
 * \ingroup metis
 *
 */
void st_metis(const ijvType *ijv, const int k, int *comp, int metisStrat)
{
  
  int wgtflag;

  wgtflag = 1;

  int options[8] = {1, 3, 1, 1, 0, 0, 0, 0};
  int numflag = 0;
  int *vwgt = NULL;

  int edgecut;
  
  int *xadj;
  int *adjncy;
  int *adjwgt;
  int *vsize;

  struct tms t1 ;
  struct tms t2 ;

  ijv2metis(ijv, &(xadj), &(adjncy), &(adjwgt));

  /*
  ijv2metisV(ijv, &(xadj), &(adjncy), &(vwgt), &(vsize));
  */

  (void) times (&t1) ;

  
  if (k == 2) {
    /* partgraph recursive */
    
    /* Do the call */
    METIS_PartGraphRecursive (&(ijv->n), xadj, adjncy, vwgt, adjwgt, &wgtflag,
			      &numflag, &k, options, &edgecut, comp);
  } else {
    /* partgraph kway */

    options[3] = metisStrat;

    METIS_PartGraphKway  (&(ijv->n), xadj, adjncy, vwgt, adjwgt, &wgtflag,
			  &numflag, &k, options, &edgecut, comp);


    /* never worked!
    METIS_PartGraphVKway  (&(ijv->n), xadj, adjncy, vwgt, vsize, &wgtflag,
			  &numflag, &k, options, &edgecut, comp);
    */

    /*
    METIS_PartGraphRecursive (&(ijv->n), xadj, adjncy, vwgt, adjwgt, &wgtflag,
			      &numflag, &k, options, &edgecut, comp);

    */
    
  }


  (void) times (&t2) ;
  /*
  fprintf(stderr,"metis time: %g\n",(double) (t2.tms_utime +
				      t2.tms_stime - t1.tms_utime - t1.tms_stime));
  */


  

  cFree(xadj);
  cFree(adjwgt);
  cFree(adjncy);
  
}
