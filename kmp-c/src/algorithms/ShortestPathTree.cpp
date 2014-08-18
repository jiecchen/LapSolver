#include <vector>
#include <limits>
#include <set>
#include "ShortestPathTree.h"
using namespace std;

class ArrayComparator
{
public:
    ArrayComparator(double * dist):dist(dist) { };
    
    bool operator() (const int& x, const int& y) const {
        if(dist[x] == dist[y])
            return x < y;
        return dist[x] < dist[y];
    }

private:
    double *dist;
};

template <class Key, class Comparator>
static set<Key, Comparator> make_set(Comparator c) {
    return set<Key, Comparator>(c);
}

ShortestPathTree::ShortestPathTree(const Graph &g, int source)
{
    int n = g.nv();
    dist = new double[n];
    parent = new int[n];
    parentIndex = new int[n];
    weight = new double[n];

    parent[source] = source;

    dist[0:n] = numeric_limits<double>::infinity();
    dist[source] = 0;
    
    auto nextNodes = make_set<int, ArrayComparator>(ArrayComparator(dist));

    for (int i = 0; i < n; ++i)
        nextNodes.insert(i);

    vector<bool> settled(n);

    while(!nextNodes.empty()) {
        int u = *nextNodes.begin();
        nextNodes.erase(u);
        auto nbrs_u = g.getNeighbors(u);
        auto wght_u = g.getWeights(u);

        for (int i = 0, deg = g.getDegree(u); i < deg; ++i) {
            int v = nbrs_u[i];
            if(!settled[v]) {
                double alt = dist[u] + 1/wght_u[i];
                if(alt < dist[v]) {
                    nextNodes.erase(v);
                    dist[v] = alt;
                    weight[v] = 1/wght_u[i];
                    parent[v] = u;
                    parentIndex[v] = i;
                    nextNodes.insert(v);
                }
            }
        }

        settled[u] = true;
    }
}

ShortestPathTree::~ShortestPathTree()
{
    delete [] dist;
    delete [] parent;
    delete [] parentIndex;
    delete [] weight;
}
