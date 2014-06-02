// Had to templatize IJVclass so that it could handle the lowstretch graph type. Templates have to be entirely included. Thus moved whole class to header

// #include "graph.hpp"
// #include <stdio.h>

// IJVclass::IJVclass(const char* fileName) : n(0), nnz(0), i(NULL), j(NULL), v(NULL)
// {


//     FILE* fpr;
//     if ((fpr = fopen(fileName, "rb")) == NULL)
// 	exit(1);

//     if (!fread(&(n), sizeof(int), 1, fpr))
// 	exit(1);//cError ("error reading n from file\n");
    
//     if (!fread(&(nnz), sizeof(int), 1, fpr))
// 	exit(1);//cError ("error reading nnz from file\n");
    
//     i = new int[nnz];
//     j = new int[nnz];
//     v = new double[nnz];

//     if (!fread(i, sizeof(int), nnz, fpr))
// 	exit(1); //cError ("error reading i from file\n");
    
//     if (!fread(j, sizeof(int), nnz, fpr))
// 	exit(1); //cError ("error reading j from file\n");
    
//     if (!fread(v, sizeof(double), nnz, fpr))
// 	exit(1); // cError ("error reading v from file\n");
    
//     for (int x = 0; x < nnz; x++)
//     {
// 	i[x] -= 1;
// 	j[x] -= 1;
//     }

//     fclose(fpr);
// }

// graph_t* IJVclass::getBoostGraph()
// { 
//     graph_t* g_ptr = new graph_t(n);
//     graph_t& g = *g_ptr;
//     weight_map_t weightmap = get(edge_weight, g);

//     for (int x = 0; x < nnz; x++)
//     {
// 	edge_descriptor e;
// 	bool inserted;
	
// 	tie(e, inserted) = add_edge(i[x], j[x], g);
// 	weightmap[e] = v[x];
//     }
    
//     return g_ptr;
// }
