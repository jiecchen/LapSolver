// c functions/structs for dealing with graphs (from oldc folder)
#include "graph.h"

void cError(const char* msg)
{
    fprintf(stderr, "%s", msg);
}

void *cCalloc(size_t n, size_t size, char *errMsg)
{

    void *obj;

    if (n == 0)
	return NULL;

 
    obj = calloc(n,size);
  
    if (!(obj)) {
	fprintf(stderr,"calloc failed: %s\n, mem %d\n", errMsg,n*size);
	exit(1);
    }

    return obj; 
}

void cFree(void *obj)
{

    if (!obj)
	fprintf(stderr,"cFree error!\n");

    if (obj)
	free(obj);
}

void freePArray(pArray* array)
{
    free(array->array);
    free(array);
    array = NULL;
}

pArray* binReadPArray(const char* fileName)
{
    pArray* out = malloc(sizeof(pArray));
    FILE* fp;
    if ((fp = fopen(fileName, "rb")) == NULL)
    {
	fprintf (stderr, "could not open the file: %s\n", fileName);
	exit(1);
    }
    int n;
    if (!fread(&n, sizeof(int), 1, fp))
    {
	fprintf (stderr, "error reading n from pArray file\n");
	exit(1);
    }
    out->n = n;
    out->array = malloc(sizeof(int) * n);
    if (!fread(out->array, sizeof(int), n, fp))
    {
	fprintf (stderr, "error reading pArray from file\n");
	exit(1);
    }
    return out;
}

// graph must be a spanning tree
pArray* graph2pArray(myGraph* graph)
{
    pArray* parray = newPArray(graph->n);
    int* q = malloc(sizeof(int) * graph->n);
    int q_front = 0;
    int q_back = 0;
    int ctr = 0;
    for (ctr = 0; ctr < graph->n; ctr++)
	parray->array[ctr] = -1;         // negative value means unvisited
    
    // make vert 0 root and then do bfs
    parray->array[0] = 0;
    q[q_back++] = 0;
    while (q_front < graph->n)
    {
	int curVert = q[q_front++];
	int nbrItr =  0;
	for (nbrItr = 0; nbrItr < graph->deg[curVert]; nbrItr++)
	{
	    int nbrVert = graph->nbrs[curVert][nbrItr];
	    if (parray->array[nbrVert] < 0)
	    {
		parray->array[nbrVert] = curVert;
		q[q_back++] = nbrVert;
	    }
	}
    }

    free(q);
    return parray;
}

ijvType *makeIJV(int nnz)
{

    ijvType *ijv;

    ijv = (ijvType *) cCalloc(1,sizeof(ijvType),
			      "makeIJV");

    ijv->nnz = nnz;
  
    ijv->i = (int *) cCalloc(nnz,sizeof(int),"makeIJV");
    ijv->j = (int *) cCalloc(nnz,sizeof(int),"makeIJV");
    ijv->v = (double *) cCalloc(nnz,sizeof(double),"makeIJV");

    return ijv;
}

#define int32 int

ijvType *binReadIJV(FILE* fpr) 
{
    int n, nnz;
  
    if (!fread(&(n), sizeof(int32), 1, fpr))
	cError ("error reading n from file\n");

    if (!fread(&(nnz), sizeof(int32), 1, fpr))
	cError ("error reading nnz from file\n");

    ijvType *ijv;
    ijv = makeIJV(nnz);
    ijv->n = n;

    int x;
    //int i, j;
    //float v;

    if (!fread(ijv->i, sizeof(int32), nnz, fpr))
	cError ("error reading i from file\n");

    if (!fread(ijv->j, sizeof(int32), nnz, fpr))
	cError ("error reading j from file\n");

    if (!fread(ijv->v, sizeof(double), nnz, fpr))
	cError ("error reading v from file\n");

    for (x = 0; x < nnz; x++) {
	ijv->i[x] = ijv->i[x]-1;
	ijv->j[x] = ijv->j[x]-1;
	ijv->v[x] = 1.0/ijv->v[x];
    }

    return ijv;
}

// this version just reads in the costs instead of taking reciprical for length
ijvType *binReadIJVcosts(FILE* fpr) 
{
    int n, nnz;
  
    if (!fread(&(n), sizeof(int32), 1, fpr))
	cError ("error reading n from file\n");

    if (!fread(&(nnz), sizeof(int32), 1, fpr))
	cError ("error reading nnz from file\n");

    ijvType *ijv;
    ijv = makeIJV(nnz);
    ijv->n = n;

    int x;
    //int i, j;
    //float v;

    if (!fread(ijv->i, sizeof(int32), nnz, fpr))
	cError ("error reading i from file\n");

    if (!fread(ijv->j, sizeof(int32), nnz, fpr))
	cError ("error reading j from file\n");

    if (!fread(ijv->v, sizeof(double), nnz, fpr))
	cError ("error reading v from file\n");

    for (x = 0; x < nnz; x++) {
	ijv->i[x] = ijv->i[x]-1;
	ijv->j[x] = ijv->j[x]-1;
	ijv->v[x] = ijv->v[x];
    }

    return ijv;
}

/**
 * \ingroup ijv
 *
 */
void freeIJV(ijvType *ijv)
{
    if (ijv->i) cFree(ijv->i);
    if (ijv->j) cFree(ijv->j);
    if (ijv->v) cFree(ijv->v);
    
    cFree(ijv);
  
}


/** 
 * 
 * 
 * @param ijv the input ijv
 * 
 * @return a myGraph type
 *
 * @todo  make a light version of that that just
 *        uses adjacency structure,
 *        as the only place we use it now only
 *        requires that.
 */
myGraph *ijv2graph(ijvType *ijv)
{
    int x;
  
    int ind;
  
    myGraph *G;

    G = (myGraph *) cCalloc(1,sizeof(myGraph),"G in ijv2sparse");

    G->n = ijv->n;
    G->nnz = 2 * ijv->nnz;


    G->deg = (int *) cCalloc(G->n, sizeof(int),"G->deg in ijv2sparse");
    G->nbrs = (int **) cCalloc(G->n, sizeof(int *),"G->nbrs in ijv2sparse");
    G->wts = (double **) cCalloc(G->n, sizeof(double *), "G->wts in ijv2sparse");
    G->back = (int **) cCalloc(G->n, sizeof(int *),"G->back in ijv2sparse");

    if (ijv->nnz == 0) {
	for (x = 0; x < G->n; x++) {
	    G->deg[x] = 0;
	    G->nbrs[x] = NULL;
	    G->wts[x] = NULL;
	    G->back[x] = NULL;

	    G->nbrsBlock = NULL;
	    G->wtsBlock = NULL;
	    G->backBlock = NULL;
	}
	return (G);
    }

    /* allocate the data structures we will need to 
       store the graph */

    G->nbrsBlock = (int *) cCalloc(G->nnz, sizeof(int),"nbrsBlock in ijv2sparse");
    G->wtsBlock = (double *) cCalloc(G->nnz, sizeof(double),"wtsBlock in ijv2sparse");
    G->backBlock = (int *) cCalloc(G->nnz, sizeof(int),"backBlock in ijv2sparse");

    for (x = 0; x < G->n; x++)
	G->deg[x] = 0;

    /* first, set all the degrees */
    for (x = 0; x < ijv->nnz; x++) {
	G->deg[ijv->i[x]]++;
	G->deg[ijv->j[x]]++;
    }

    /* then, allocate space for the data for each vertex */
    ind = 0;
    for (x = 0; x < G->n; x++) {
	G->nbrs[x] = &(G->nbrsBlock[ind]);
	G->wts[x] = &(G->wtsBlock[ind]);
	G->back[x] = &(G->backBlock[ind]);
	ind += G->deg[x];
    }

    /* and now, populate the graph with edges,
       temporarily resetting the degrees
       (don't worry: they will return to their prev values */
    for (x = 0; x < G->n; x++)
	G->deg[x] = 0;

    int row; /* i */
    int col; /* j */

    for (x = 0; x < ijv->nnz; x++) {

	row = ijv->i[x];
	col = ijv->j[x];

	G->nbrs[col][G->deg[col]] = row;
	G->nbrs[row][G->deg[row]] = col;

	G->wts[col][G->deg[col]] = ijv->v[x];
	G->wts[row][G->deg[row]] = ijv->v[x];
    
	G->back[col][G->deg[col]] = G->deg[row];
	G->back[row][G->deg[row]] = G->deg[col];
	  
	G->deg[row]++;
	G->deg[col]++;
    }

    return G;
}

// just do all the mem allocation, nothin else
// don't deal with back block
myGraph* newGraph(int n, int nnz)
{
    myGraph* h;

    h = malloc(sizeof(myGraph));

    h->n = n;
    h->nnz = nnz;

    h->deg = (int*)malloc(sizeof(int)* h->n);
    h->nbrs = (int**)malloc(sizeof(int)* h->n);
    h->wts = (double**)malloc(sizeof(double*) * h->n);
    h->back = NULL; //(int**)malloc(sizeof(int*) * h->n);   // I don't need this!
    h->nbrsBlock = (int*) malloc(sizeof(int) * h->nnz);
    h->wtsBlock = (double*) malloc(sizeof(double) * h->nnz);
    h->backBlock = NULL;
    
    return h;
}


void freeGraph(myGraph *G)
{
    int i;

    for (i = 0; i < G->n; i++) {
	G->nbrs[i] = NULL;
	G->wts[i] = NULL;
    }
    
    if (G->nbrs) cFree(G->nbrs);
    if (G->wts) cFree(G->wts);
    if (G->back) cFree(G->back);
    if (G->deg) cFree(G->deg);
  
    if (G->nnz > 0) {
	if (G->backBlock)       // (I don't always allocate backBlock)
	    cFree(G->backBlock);
	cFree(G->wtsBlock);
	cFree(G->nbrsBlock);
    }
  
    cFree(G);
  
}

// wrapper function to combine reading ijv and converting to myGraph
// take the name of a binijv file and return its graph
myGraph* getGraph(const char* fileName)
{
    FILE* fp = NULL;
    ijvType* ijv = NULL;
    myGraph* graph = NULL;
    
    if ((fp = fopen(fileName, "rt")) == NULL)
	exit(1);
    
    ijv = binReadIJV(fp);
    fclose(fp);
    graph = ijv2graph(ijv);
    freeIJV(ijv);
    return graph;
}

pArray* newPArray(int n)
{
    pArray* out = malloc(sizeof(pArray));
    out->n = n;
    out->array = malloc(sizeof(int) * n);
    return out;
}

// assume that we need to add 1 to all vertex numbers for matlab
// ie loading a pArray and then directly saving it again would not work!
// pArray is modified!
void binWritePArray(pArray* parray, const char* fileName)
{
    FILE* fp;
    if ((fp = fopen(fileName, "wt")) == NULL)
	exit(1);
    
    int ctr = 0;
    for (ctr = 0; ctr < parray->n; ctr++)
    {
	//parray->array[ctr]++;                                // wait did I forget matlab compatibility?
    }

    fwrite(&(parray->n), sizeof(int), 1, fp);
    fwrite(parray->array, sizeof(int), parray->n, fp);

    fclose(fp);

}

void textWritePArray(pArray* parray, FILE* fp)
{
    fprintf(fp, "%d", parray->n);

    int ctr = 0;
    for (ctr = 0; ctr < parray->n; ctr++)
    {
	fprintf(fp, "\n%d", parray->array[ctr]);
    }
}

// will modify ijv!
void binWriteIJV(ijvType* ijv, const char* fileName)
{
    // matlab format it and turn lengths to costs
    int ctr = 0;
    for (ctr = 0; ctr < ijv->nnz; ctr++)
    {
	ijv->i[ctr]++;
	ijv->j[ctr]++;
	ijv->v[ctr] = 1.0/ijv->v[ctr];
	
    }

    FILE* fp = NULL;
    if ((fp = fopen(fileName, "wt")) == NULL)
    {
	fprintf(stderr, "ERROR: failed to open %s\n", fileName);
	exit(2);
    }
    fwrite(&(ijv->n), sizeof(int), 1, fp);
    fwrite(&(ijv->nnz), sizeof(int), 1, fp);
    fwrite(ijv->i, sizeof(int), ijv->nnz, fp);
    fwrite(ijv->j, sizeof(int), ijv->nnz, fp);
    fwrite(ijv->v, sizeof(double), ijv->nnz, fp);
    fclose(fp);
}

// will modify ijv!
void textWriteIJV(ijvType* ijv, FILE* fp)
{
    // matlab format it and turn lengths to costs
    int ctr = 0;
    fprintf(fp, "%d %d\n", ijv->n, ijv->nnz);
    for (ctr = 0; ctr < ijv->nnz; ctr++)
    {
	ijv->i[ctr]++;
	ijv->j[ctr]++;
	ijv->v[ctr] = 1.0/ijv->v[ctr];
	fprintf(fp, "%d %d %lf\n", ijv->i[ctr], ijv->j[ctr], ijv->v[ctr]);
    }
}
