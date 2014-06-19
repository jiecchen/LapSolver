/**
 * @file KruskalTree.java
 * @author Alex Reinking <alexander.reinking@yale.edu>
 * @date Tue Jun 3 2014
 *
 * A not-so-low-stretch spanning tree strategy, which just uses Kruskal's algorithm.
 * This will not produce good low stretch spanning trees.
 *
 */
package lapsolver.lsst;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.UnionFind;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class KruskalTree implements SpanningTreeStrategy {
    @Override
    public Tree getTree(Graph graph) {
        UnionFind disjointSet = new UnionFind(graph.nv);
        EdgeList treeEdges = new EdgeList(graph.nv - 1);

        // Compute auxiliarySize
        EdgeList inputEdges = new EdgeList(graph);

        final double[] weights = inputEdges.weight;
        final List<Integer> indices = new ArrayList<>(weights.length);
        for (int i = 0; i < weights.length; i++) {
            indices.add(i, i);
        }

        Collections.sort(indices, new Comparator<Integer>() {
            @Override
            public int compare(Integer o1, Integer o2) {
                return Double.compare(weights[o1], weights[o2]);
            }
        });

        int currentEdge = 0;
        for (Integer edge : indices) {
            int u = inputEdges.u[edge];
            int v = inputEdges.v[edge];

            if (disjointSet.find(u) != disjointSet.find(v)) {
                disjointSet.union(u, v);

                treeEdges.u[currentEdge] = u;
                treeEdges.v[currentEdge] = v;
                treeEdges.weight[currentEdge++] = weights[edge];
            }
        }
        return new Tree(treeEdges);
    }
}
