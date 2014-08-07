/**
 * @file VoronoiTree.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Thu Aug 7 2014
 *
 * A sketchy LSST algorithm, using multiple Dijkstra balls.
 *
 */

package lapsolver.lsst;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.MappedEdgeList;
import lapsolver.algorithms.ShortestPathTree;
import lapsolver.algorithms.UnionFind;

import java.util.Arrays;
import java.util.Random;

public class VoronoiTree extends OrderedSpanningTreeStrategy {
    public int nSamples = 5;

    @Override
    public EdgeList getTreeEdges (Graph graph) {
        MappedEdgeList edgeMap = new MappedEdgeList(graph);

        double[] times = new double[graph.ne];
        double[] ecc = new double[graph.nv];

        Arrays.fill(times, Double.POSITIVE_INFINITY);
        Arrays.fill(ecc, Double.POSITIVE_INFINITY);

        Random rand = new Random();
        for (int i = 0; i < graph.nv; i++) {
            ecc[i] = 1e10 + rand.nextDouble();
        }

        for (int sample = 0; sample < nSamples; sample++) {
            int source = 0;
            for (int i = 0; i < graph.nv; i++) {
                if (ecc[i] > ecc[source]) {
                    source = i;
                }
            }

            ShortestPathTree spt = new ShortestPathTree(graph, source);
            int[] parentIndex = spt.getParentIndex();
            double[] dist = spt.getDist();

            for (int i = 0; i < graph.nv; i++) {
                if (dist[i] < ecc[i]) ecc[i] = dist[i];
            }

            for (int u = 0; u < graph.nv; u++) {
                if (u == source) continue;
                int treeEdge = edgeMap.index[u][parentIndex[u]];
                times[treeEdge] = Math.min(times[treeEdge], dist[u]);
            }
        }

        return UnionFindTreeBuilder.buildTree(edgeMap.edges, graph.nv, times);
    }
}
