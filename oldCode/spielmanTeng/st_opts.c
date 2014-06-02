/***************
   st_opts.c

  \newgroup opts
*****************/

#include "st_defs.h"


/** 
 * \brief default opts for partTree
 * \ingroup opts
 * 
 * @return a pointer to the defaults opts for partTree
 */
int *partTreeOpts() {
  int *opts;

  opts = (int *) cCalloc(numOpts, sizeof(int),"opts");
  opts[stOptNumTree] = stOptValTreeYes;
  opts[stOptNumParts] = 2;
  
  opts[stOptNumMetis] = metisStrategyMinConn;

  return opts;
}

/** 
 * \brief default opts for partPre
 * \ingroup opts
 * 
 * @return a pointer to the defaults opts for partPre
 * @param k the number of parts requested
 */
int *partPreOpts(int k) {
  int *opts;

  opts = (int *) cCalloc(numOpts, sizeof(int),"opts");
  opts[stOptNumTree] = stOptValTreeNo;
  opts[stOptNumParts] = k;
  opts[stOptNumSpars] = stOptValSparsNo;

  opts[stOptNumMetis] = metisStrategyMinConn;

  opts[stOptNumTime] = 0;

  opts[stOptNumTime2] = 0;
  
  return opts;
}


/** 
 * \ingroup opts
 * 
 */
int *copyOpts(int *opts)
{
 int *newOpts;
 int i;

 newOpts = (int *) cCalloc(numOpts, sizeof(int),"copyOpts");
 for (i = 0; i < numOpts; i++) 
   newOpts[i] = opts[i];

 return newOpts;

}
