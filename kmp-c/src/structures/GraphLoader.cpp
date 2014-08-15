#include "GraphLoader.h"
#include <iostream>
#include <fstream>
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

    return fromStream(ijvFile);
}

Graph GraphLoader::fromStdin()
{
    return fromStream(cin);
}

Graph GraphLoader::fromStream(istream &inStream)
{
    int nv, ne;
    inStream >> nv >> ne;

    aligned_vector<int> u = aligned_vector<int>(ne);
    aligned_vector<int> v = aligned_vector<int>(ne);
    aligned_vector<double> w = aligned_vector<double>(ne);

    for (int i = 0; i < ne; i++) {
        inStream >> u[i] >> v[i] >> w[i];
    }

    return fromArrays(u, v, w);
}