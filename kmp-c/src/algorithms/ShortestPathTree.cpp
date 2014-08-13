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
    dist = new double[g.nv];
    parent = new int[g.nv];
    weight = new double[g.nv];

    parent[source] = source;

    dist[0:g.nv] = numeric_limits<double>::infinity();
    dist[source] = 0;
    
    auto nextNodes = make_set<int, ArrayComparator>(ArrayComparator(dist));

    for (int i = 0; i < g.nv; ++i)
        nextNodes.insert(i);

    vector<bool> settled(g.nv);

    while(!nextNodes.empty()) {
        int u = *nextNodes.begin();
        nextNodes.erase(u);
        auto nbrs_u = g.getNeighbors(u);
        auto wght_u = g.getWeights(u);

        for (int i = 0, deg = g.getDegree(u); i < deg; ++i) {
            int v = nbrs_u[i];
            if(!settled[v]) {
                double alt = dist[u] + wght_u[i];
                if(alt < dist[v]) {
                    nextNodes.erase(v);
                    dist[v] = alt;
                    weight[v] = wght_u[i];
                    parent[v] = u;
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
    delete [] weight;
}
