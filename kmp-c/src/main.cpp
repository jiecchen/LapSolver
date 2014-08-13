/* Copyright 2014 Yale Institute for Network Science
 * Author: Alex Reinking
 *
 * Example benchmark for vector arithmetic
 */
#include <sys/time.h>
#include <mkl.h>
#include <vector>
#include <numeric>
#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <algorithm>
#include "util/aligned.h"
#include "structures/graph.h"

int main(int argc, char const *argv[])
{
    vmlSetMode(VML_EP);
    srand(time(NULL));

    aligned_vector<int> u = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1 };
    aligned_vector<int> v = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 5, 8 };
    aligned_vector<double> w = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

    EdgeList edges(u,v,w);
    Graph g(edges);

    printGraph(g);

    return 0;
}
