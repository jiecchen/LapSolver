/**
 * @file StarDecompositionTree.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Fri Jun 13 2014
 *
 * An implementation of the star decomposition tree algorithm in EEST05.
 * Just a skeleton for now.
 */

package lapsolver.lsst;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.ShortestPathTree;

import java.util.Arrays;

public class StarDecompositionTree implements SpanningTreeStrategy {
    public Graph graph;

    @Override
    public Tree getTree(Graph graph) {
       // TODO(Cyril): set beta/n contraction threshold her
       // TODO(Cyril): randomized centre picking strategy?
       return new Tree( getLowStretchTree(graph, 0) );
    }

    // generate the spanning tree
    // (LowStretchTree in EEST05)
    public EdgeList getLowStretchTree(Graph graph, int x0) {
        // below some threshold, just return the shortest path tree
        // same threshold as in Yu's code for now
        if (graph.nv < 32) {
            // TODO(Cyril): play around with this
            ShortestPathTree sptInstance = new ShortestPathTree(graph, x0);
            return new EdgeList( sptInstance.getTree() );
        }

        // contract small edges
        // TODO(Cyril): implement this! change weights to 0 (don't forget both directions)

        // obtain star coloring
        int[] colors = new int[graph.nv];
        EdgeList bridges = getStarColoring(graph, x0, colors);
        int nColors = bridges.ne + 1;

        // expand contracted edges
        // TODO(Cyril): implement this! change weights back to original

        // generate induced subgraphs
        int[][] relabelUp = new int[graph.nv][];
        int[] relabelDown = new int[graph.nv];
        Graph[] subgraphs = splitGraph(graph, colors, nColors, relabelUp, relabelDown);
        EdgeList[] childTreeEdges = new EdgeList[nColors];

        // get trees recursively
        // TODO(Cyril): parallelize here
        for (int color = 0; color < nColors; color++) {
            int xi = relabelDown[color == 0 ? x0 : bridges.v[color-1]];
            childTreeEdges[color] = getLowStretchTree(subgraphs[color], xi);
        }

        // merge results from child calls
        EdgeList parentTreeEdges = new EdgeList(graph.nv - 1);
        int edgePos = 0;
        for (int color = 0; color < nColors; color++) {
            // merge child with relabeled vertices
            for (int i = 0; i < childTreeEdges[color].ne; i++) {
                parentTreeEdges.u[edgePos] = relabelUp[color][childTreeEdges[color].u[i]];
                parentTreeEdges.v[edgePos] = relabelUp[color][childTreeEdges[color].v[i]];
                parentTreeEdges.weight[edgePos] = childTreeEdges[color].weight[i];
                edgePos++;
            }
        }

        // merge bridges into result
        for (int i = 0; i < bridges.ne; i++) {
            parentTreeEdges.u[edgePos] = bridges.u[i];
            parentTreeEdges.v[edgePos] = bridges.v[i];
            parentTreeEdges.weight[edgePos] = bridges.weight[i];
            edgePos++;
        }

        return parentTreeEdges;
    }

    // if you pass in an array of size nv, I write back:
    //     colors, where colors[v] = component label of v (0 for ball)
    // and I'll return an EdgeList with bridge edges (u[i] = x_(i+1), v[i] = y_(i+1))
    // (StarDecomp in EEST05)
    public EdgeList getStarColoring(Graph graph, int x0, int[] colors) {
        // initially nobody's part of the decomposition
        colors = new int[graph.nv];
        Arrays.fill(colors, -1);

        ShortestPathTree sptInstance = new ShortestPathTree(graph, x0);
        Tree shortestPathTree = sptInstance.getTree();
        double[] dist = sptInstance.getDist();

        // grow low-cut ball, build bridges from shell
        int[] roots = growBall(graph, shortestPathTree, dist, x0, colors, 0);
        EdgeList bridges = new EdgeList(roots.length);
        for (int i = 0; i < roots.length; i++) {
            bridges.u[i] = shortestPathTree.parent[roots[i]];
            bridges.v[i] = roots[i];
            bridges.weight[i] = shortestPathTree.weight[roots[i]];
        }

        // grow low-cut cones from shell
        for (int i = 0; i < roots.length; i++) {
            growCone(graph, shortestPathTree, roots[i], 0.0, colors, i+1);
        }

        return bridges;
    }

    // update colors (write back) with color for a low-cut ball
    // return vertices in shell
    // should be called with color=0 in this algorithm
    public int[] growBall(Graph graph, Tree shortestPathTree, double[] dist,
                          int source, int[] colors, int color) {
        // TODO(Cyril): implement this. might need other parameters

        colors[source] = color;
        return new int[0];
    }

    // update colors (write back) with color for low-cut cone
    // should be called with
    public void growCone(Graph graph, Tree shortestPathTree,
                         int source, double radius, int[] colors, int color) {
        // TODO(Cyril): implement this. might need other parameters

        colors[source] = color;
    }

    // split graph into induced subgraphs with same color
    // if you pass in arrays of length graph.nv, I write back:
    //     relabelUp, such that relabelUp[k][v] is the name in
    //          parent graph of vertex v in subgraph k
    //     relabelDown, such that relabelDown[i] is the name of
    //          vertex v in its subgraph
    public Graph[] splitGraph(Graph graph, int[] colors, int nColors,
                              int[][] relabelUp, int[] relabelDown) {
        Graph[] subgraphs = new Graph[nColors];

        // count frequency of each color
        int[] subgraphSizes = new int[nColors];
        for (int color : colors) {
            subgraphSizes[color]++;
        }

        for (int color = 0; color < nColors; color++) {
            subgraphs[color] = null;
            relabelUp[color] = new int[subgraphSizes[color]];
            // TODO(Cyril): compute induced subgraphs
            // TODO(Cyril): consider making GraphUtil for this? too specialized?
        }

        return subgraphs;
    }
}
