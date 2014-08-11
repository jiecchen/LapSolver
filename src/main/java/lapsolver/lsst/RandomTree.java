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

import java.util.Random;

public class RandomTree extends OrderedSpanningTreeStrategy {
    @Override
    public EdgeList getTreeEdges(Graph graph) {
        EdgeList edges = new EdgeList(graph);
        double[] randTimes = new double[graph.ne];

        Random rand = new Random();
        for (int i = 0; i < graph.ne; i++) {
            randTimes[i] = rand.nextDouble();
        }

        return UnionFindTreeBuilder.buildTree(edges, graph.nv, randTimes);
    }
}
