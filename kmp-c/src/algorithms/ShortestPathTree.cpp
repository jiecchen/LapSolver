#include <vector>
#include <limits>
#include "structures/heap/Aligned4aryHeap.h"
#include "ShortestPathTree.h"
using namespace std;

ShortestPathTree::ShortestPathTree(const Graph &g, int source)
{
    size_t n = (size_t) g.nv();
    dist = new double[n];
    parent = new int[n];
    parentIndex = new int[n];
    weight = new double[n];

    parent[source] = source;

    dist[0:n] = numeric_limits<double>::infinity();
    dist[source] = 0;

    Aligned4aryHeap<double, int> pq(n << 1);
    pq.push(0, source);

    vector<int> settled(n, false);

    while (!pq.isEmpty())
    {
        double k; int u;
        pq.pop(&k, &u);
        for (int i = 0, deg = g.degree(u); i < deg; ++i)
        {
            int v = g.neighbor(u, i);
            if (!settled[v])
            {
                double edgeW = 1 / g.weight(u, i);
                double alt = dist[u] + edgeW;
                if (alt < dist[v])
                {
                    pq.push(alt, v);

                    dist[v] = alt;
                    weight[v] = edgeW;
                    parent[v] = u;
                    parentIndex[v] = i;
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
