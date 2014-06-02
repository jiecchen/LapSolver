#ifdef mex_h
#include "mex.h" 
#endif

#define min(A,B) ((A) < (B) ? (A) : (B))
#define max(A,B) ((A) > (B) ? (A) : (B))


/**
 * \ingroup mex
 */
struct myGrapht {
  int n;
  int nnz;

  int *deg; /* up to n */
  int **nbrs; /* up to n */
  int *nbrsBlock; /* holds the nbrs, nnz */
  double **wts; /* up to n */
  double *wtsBlock; /* up to nnz (has duplicates) */
  int *backBlock; /* up to nnz */
  int **back; /* up to n */

};

#define myGraph struct myGrapht

/**
 * should maintain j[x] < i[x]
 *
 * \defgroup ijv
 */
struct ijvTypet {
  int n;
  int nnz;
  int *i;
  int *j;
  double *v;
};

#define ijvType struct ijvTypet

/**
 *
 * \brief dubIjvType is just like ijvType,
 * except that each edge appears in
 * each direction
 *
 */
typedef ijvType dubIjvType;



/**
 *
 *      partType
 *
 * A partition contains:
 *
 * 1. compVec: length n, 
 *
 *      compVec[i] is index of component
 *      of node i
 *
 *  2. compMap: length n,
 *       compMap[i] is the index of the node
 *       in component compVec[i] to which
 *       node i is mapped
 *
 *  3. numParts: the number of components
 *
 *  4. partSizes: the number of nodes in each component
 *
 *  5. partMap[][], partMap[p][i] is the
 *       index of the node in the original
 *       graph that is mapped to node i in part p.
 *
 *  6. partMapBlock[], the block of memory containing
 *       all of these.
*/

struct partTypet {
  int *compVec;
  int *compMap;
  int numParts;
  int *partSizes;
  int **partMap;
  int *partMapBlock;
};

#define partType struct partTypet

/*---------------------------------

    ijvBlockDiag

   is a type contain the diagonal blocks
   of a graph, specified by a partition.
   each block is given by an ijv.

   ----------------------------------*/

struct ijvBlockDiagt {
  int n; /* the number of nodes in the graph */
  partType *part; /* a partition giving the mapping to blocks */
  ijvType **blocks;
};

#define ijvBlockDiag struct ijvBlockDiagt

  
/* format for opts:
  stOptNumX is the index of the option, and
  stOptValXY is value Y for option X
*/



#define numOpts 12

#define stOptNumTree 0
#define stOptValTreeYes 1
#define stOptValTreeNo 0

/**
 * \def stOptNumParts
 * \brief should be an integer
 */
#define stOptNumParts 1

#define stOptNumSpars 2
#define stOptValSparsNo 0
#define stOptValSparsPrim 1
#define stOptValSparsMin 2

/** should be an integer */
#define stOptNumSparsParam 3

/** \brief 1 for yes and 0 for no */
#define stOptNumDebug 4

/** \brief one of the metisStrategyX values */
#define stOptNumMetis 5


/** \brief the one that uses memory to minimize*/
#define metisStrategyMinConn 3

/** \brief reducing metis memory*/
#define metisStrategyRandom 1

/** \brief reducing metis memory*/
#define metisStrategyGreedy 2

/**
 * \def stOptNumSparsMinDeg
 * \brief an integer, 0 if not being used
 */
#define stOptNumSparsMinDeg 6

/**
 * \def stOptNumSparsMinFracDeg
 * \brief 0 if unused, 1 if used, in which case
 *    the global sparseMinFracDeg contains the value
 */
#define stOptNumSparsMinFracDeg 7
/**
 *
 * \brief a sparsifier parameter, global
 */
float sparsMinFracDeg;

/**
 * \def stOptNumSparsMinFracWt
 * \brief 0 if unused, 1 if used, in which case
 *    the global sparseMinFracWt contains the value
 */
#define stOptNumSparsMinFracWt 8
/**
 *
 * \brief a sparsifier parameter, global
 */
float sparsMinFracWt;

/**
 * \def stOptNumData
 * \brief if 1, outputs a little data on how doing
 */
#define stOptNumData 9

#define stOptNumTime 10

#define stOptNumTime2 11

/**
 *
 * \brief a param for Metis, global
 */
float metisBalance;
