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

import java.util.*;

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
            // TODO(Cyril): play around with this constant (32)
            return new EdgeList(new ShortestPathTree(graph, x0).getTree());
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
        // TODO(Alex): parallelize here
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
        int[] ballShell = growBall(graph, dist, colors, 0);

        // grow low-cut cones from shell
        int nColors = 1;
        ArrayList<Integer> bridgeSources = new ArrayList<>();
        for (int i = 0; i < ballShell.length; i++) {
            if (colors[ballShell[i]] != -1) {
                // oops, another cone took this already
                continue;
            }

            growCone(graph, shortestPathTree, ballShell[i], 0.0, colors, nColors);
            bridgeSources.add(ballShell[i]);
            nColors++;
        }

        // build bridge vertices
        EdgeList bridges = new EdgeList(nColors);
        for (int i = 0; i < nColors; i++) {
            bridges.v[i] = bridgeSources.get(i);
            bridges.u[i] = shortestPathTree.parent[bridges.v[i]];
            bridges.weight[i] = shortestPathTree.weight[bridges.v[i]];
        }

        return bridges;
    }

    /**
     * @param graph  the containing graph
     * @param dist   the distance array from a shortest path tree
     * @param colors the colors of each vertex in graph
     * @param color  the color of the ball (should be 0)
     * @return the set of vertices directly outside the ball (in the shortestPathTree's vertex labeling)
     */
    public int[] growBall(Graph graph, final double[] dist, int[] colors, int color) {
        ArrayList<Integer> order = new ArrayList<>(graph.nv);
        for (int i = 0; i < graph.nv; i++) order.add(i, i);

        Collections.sort(order, new Comparator<Integer>() {
            @Override
            public int compare(Integer o1, Integer o2) {
                return Double.compare(dist[o1], dist[o2]);
            }
        });

        boolean[] inCut = new boolean[graph.nv];
        int bestCut = 0;
        double bestValue = Double.POSITIVE_INFINITY;

        double cutValue = 0;
        for (int i = 0; i < 2 * order.size() / 3; i++) {
            int u = order.get(i);
            inCut[u] = true;

            for (int j = 0; j < graph.deg[u]; j++) {
                int v = graph.nbrs[u][j];
                cutValue += ((inCut[v]) ? -1 : 1) * graph.weights[u][j];
            }

            if (i >= order.size() / 3 && cutValue < bestValue) {
                bestCut = i;
                bestValue = cutValue;
            }
        }

        Arrays.fill(inCut, false);
        boolean[] inShell = new boolean[graph.nv];
        for (int i = 0; i <= bestCut; i++) {
            int u = order.get(i);
            inCut[u] = true;
            inShell[u] = false;

            colors[u] = color; // mark the ball

            for (int v : graph.nbrs[u])
                if (!inCut[v])
                    inShell[v] = true; // mark the shell
        }

        // construct the shell
        ArrayList<Integer> shellList = new ArrayList<>();
        for (int i = 0; i < inShell.length; i++)
            if (inShell[i]) shellList.add(i);

        // convert it to a primitive array because lol Java lol
        int[] shell = new int[shellList.size()];
        for (int i = 0; i < shell.length; i++)
            shell[i] = shellList.get(i);

        return shell;
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
        // compute relabeling
        int[] subgraphSizes = new int[nColors];
        for (int i = 0; i < graph.nv; i++) {
            relabelUp[colors[i]][subgraphSizes[colors[i]]] = i;
            relabelDown[i] = subgraphSizes[colors[i]];
            subgraphSizes[colors[i]]++;
        }

        // traverse edges once to get subgraph edge counts
        EdgeList parentEdges = new EdgeList(graph);
        int[] subgraphEdgeCounts = new int[nColors];
        for (int i = 0; i < parentEdges.ne; i++) {
            int u = parentEdges.u[i], v = parentEdges.v[i];
            if (colors[u] == colors[v]) {
                // count edge in induced subgraph
                subgraphEdgeCounts[colors[u]]++;
            }
        }

        // build subgraph edge lists
        EdgeList[] subgraphEdges = new EdgeList[nColors];
        int[] edgePos = new int[nColors];
        for (int color = 0; color < nColors; color++) {
            subgraphEdges[color] = new EdgeList(color);
        }
        for (int i = 0; i < parentEdges.ne; i++) {
            int u = parentEdges.u[i], v = parentEdges.v[i];
            if (colors[u] == colors[v]) {
                // add edge to induced subgraph
                int color = colors[u];
                double weight = parentEdges.weight[i];
                subgraphEdges[color].u[edgePos[color]] = relabelDown[u];
                subgraphEdges[color].v[edgePos[color]] = relabelDown[v];
                subgraphEdges[color].weight[edgePos[color]] = weight;
                edgePos[colors[u]]++;
            }
            double w = parentEdges.weight[i];
        }

        // convert edge lists to graphs
        Graph[] subgraphs = new Graph[nColors];
        for (int color = 0; color < nColors; color++) {
            subgraphs[color] = new Graph(subgraphEdges[color]);
        }

        return subgraphs;
    }
}
