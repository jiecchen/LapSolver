/* save an ijv graph to an ascii file */
#include "math.h"
#include "stdio.h"
#include "mex.h"   /*--This one is required*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* declare crap b/c this is c */
    const mxArray *iData;
    double *iValues;
    int iRowLen, iColLen;

    const mxArray *jData;
    double *jValues;
    int jRowLen, jColLen;

    const mxArray *vData;
    double *vValues;
    int vRowLen, vColLen;

    const mxArray *nData;
    double *nValues;
    int nRowLen, nColLen;

    const mxArray *fileData;
    int fileLength;
    char* fileString;
    FILE* fOut;

    int rowCntr = 0;

    if (nrhs != 5)
    {
	printf("call with parameters i, j, v, NUM_NODES, 'FILE_OUT_NAME'!\n");
	return;
    }
    /* get inputs */
    iData = prhs[0];
    iValues = mxGetPr(iData);
    iRowLen = mxGetN(iData);
    iColLen = mxGetM(iData);

    jData = prhs[1];
    jValues = mxGetPr(jData);
    jRowLen = mxGetN(jData);
    jColLen = mxGetM(jData);

    vData = prhs[2];
    vValues = mxGetPr(vData);
    vRowLen = mxGetN(vData);
    vColLen = mxGetM(vData);

    nData = prhs[3];
    nValues = mxGetPr(nData);
    nRowLen = mxGetN(nData);
    nColLen = mxGetM(nData);
    
    if (nRowLen != 1 || nColLen != 1)
    {
	printf("fourth arg must be a scalar giving number of vertices are in graph");
	return;
    }

    
    fileData = prhs[4];
    fileLength = mxGetN(fileData) + 1;
    fileString = mxCalloc(fileLength, sizeof(char));
    mxGetString(fileData, fileString,fileLength);

    if ((iColLen != jColLen || jColLen != vColLen)
	|| (iRowLen != 1 || jRowLen != 1 || vRowLen != 1))
    {
	printf("matrices are different sizes or don't all have 1 col\n");
	return;
    }
    
    fOut = fopen(fileString, "wt");

    fprintf(fOut, "%d %d\n", (int)nValues[0], iColLen);

    for(rowCntr = 0; rowCntr < iColLen; rowCntr++)
    {
	fprintf(fOut, "%d %d %f\n", (int)iValues[rowCntr], (int)jValues[rowCntr], vValues[rowCntr]);
    }
    
    fclose(fOut);

    return;
}
        
