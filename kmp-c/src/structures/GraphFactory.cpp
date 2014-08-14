#include "GraphFactory.h"
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

Graph GraphFactory::fromArrays(aligned_vector<int> u, aligned_vector<int> v, aligned_vector<double> w)
{
    return Graph(EdgeList(u, v, w));
}

Graph GraphFactory::fromFile(const string &filename)
{
    aligned_vector<int> u = aligned_vector<int>();
    aligned_vector<int> v = aligned_vector<int>();
    aligned_vector<double> w = aligned_vector<double>();

    string line;
    ifstream ijvFile(filename);
    if (!ijvFile)
        throw - 1;

    while (getline(ijvFile, line))
    {
        istringstream iss(line);
        int a, b;
        double c;
        if (!(iss >> a >> b >> c))
            throw u.size();
        u.push_back(a);
        v.push_back(b);
        w.push_back(c);
    }

    return fromArrays(u, v, w);
}