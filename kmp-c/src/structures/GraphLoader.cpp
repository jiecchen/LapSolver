#include "GraphLoader.h"
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

Graph GraphLoader::fromArrays(aligned_vector<int> u, aligned_vector<int> v, aligned_vector<double> w)
{
    return Graph(EdgeList(u, v, w));
}

Graph GraphLoader::fromFile(const string &filename)
{
    ifstream ijvFile(filename);
    if (!ijvFile)
        throw -1;

    int nv, ne;
    ijvFile >> nv >> ne;

    aligned_vector<int> u = aligned_vector<int>(ne);
    aligned_vector<int> v = aligned_vector<int>(ne);
    aligned_vector<double> w = aligned_vector<double>(ne);

    for (int i = 0; i < ne; i++) {
        ijvFile >> u[i] >> v[i] >> w[i];
    }

    ijvFile.close();

    return fromArrays(u, v, w);
}