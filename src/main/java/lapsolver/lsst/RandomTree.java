/**
 * @file RandomTree.java
 * @author Alex Reinking <alexander.reinking@yale.edu>
 * @date Fri Jun 13 2014
 *
 * A random spanning tree. Very high stretch, probably ;)
 *
 */
package lapsolver.lsst;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.UnionFind;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class RandomTree implements SpanningTreeStrategy {
    @Override
    public Tree getTree(Graph graph) {
        UnionFind disjointSet = new UnionFind(graph.nv);
        EdgeList treeEdges = new EdgeList(graph.nv - 1);

        // Compute index
        EdgeList inputEdges = new EdgeList(graph);

        final double[] weights = inputEdges.weight;
        List<Integer> indices = new ArrayList<>(weights.length);
        for (int i = 0; i < weights.length; i++) {
            indices.add(i, i);
        }

        Collections.shuffle(indices);

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
