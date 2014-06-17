package lapsolver.lsst;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.ShortestPathTree;

import java.util.*;

public class StarDecompositionWorker {// scratch space for cut colorings

    public int[] getColors() {
        return colors;
    }

    private int[] colors;

    public StarDecompositionWorker(Graph graph) {
        colors = new int[graph.nv];
    }

    public static class Decomposition {
        public final Graph[] subgraphs;
        public final int[][] originalLabels;
        public final int[] newLabels;

        public Decomposition(Graph[] subgraphs, int[][] originalLabels, int[] newLabels) {
            this.subgraphs = subgraphs;
            this.originalLabels = originalLabels;
            this.newLabels = newLabels;
        }
    }

    // if you pass in an array of size nv, I write back:
    //     colors, where colors[v] = component label of v (0 for ball)
    // and I'll return an EdgeList with bridge edges (u[i] = x_(i+1), v[i] = y_(i+1))
    // (StarDecomp in EEST05)
    public EdgeList makeStarCut(Graph graph, int x0, ShortestPathTree sptInstance) {
        // initially nobody's part of the decomposition
        Arrays.fill(colors, -1);

        Tree shortestPathTree = sptInstance.getTree();
        double[] dist = sptInstance.getDist();

        // grow low-cut ball, build bridges from shell
        int[] ballShell = growBall(graph, dist, 0);

        // grow low-cut cones from shell
        int nColors = 1;
        ArrayList<Integer> bridgeSources = new ArrayList<Integer>();
        for (int aBallShell : ballShell) {
            if (colors[aBallShell] != -1)
                continue; // oops, another cone took this already

            growCone(graph, shortestPathTree, aBallShell, 0.0, nColors);
            bridgeSources.add(aBallShell);
            nColors++;
        }

        // build bridge vertices
        EdgeList bridges = new EdgeList(nColors - 1);
        for (int i = 0; i < nColors - 1; i++) {
            bridges.v[i] = bridgeSources.get(i);
            bridges.u[i] = shortestPathTree.parent[bridges.v[i]];
            bridges.weight[i] = shortestPathTree.weight[bridges.v[i]];
        }

        return bridges;
    }

    /**
     * Find a low-cut ball between 1/3 and 2/3 of the total vertices
     *
     * @param graph the containing graph
     * @param dist  the distance array from a shortest path tree
     * @param color the color of the ball (should be 0)
     * @return the set of vertices directly outside the ball (in the shortestPathTree's vertex labeling)
     */
    private int[] growBall(Graph graph, final double[] dist, double low, double high, int color) {
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
        ArrayList<Integer> shellList = new ArrayList<Integer>();
        for (int i = 0; i < inShell.length; i++)
            if (inShell[i]) shellList.add(i);

        // convert it to a primitive array because lol Java lol
        int[] shell = new int[shellList.size()];
        for (int i = 0; i < shell.length; i++)
            shell[i] = shellList.get(i);

        return shell;
    }

    // generate the spanning tree

    // (LowStretchTree in EEST05)
    // update colors (write back) with color for low-cut cone
    // should be called with
    public void growCone(Graph graph, Tree shortestPathTree,
                         int source, double radius, int color) {
        Queue<Integer> bfs = new LinkedList<>();
        bfs.add(source);
        while (!bfs.isEmpty()) {
            int u = bfs.poll();
            if (colors[u] != -1) continue;
            colors[u] = color;
            for (int v : shortestPathTree.children[u])
                bfs.add(v);
        }
    }

    // split graph into induced subgraphs with same color

    // if you pass in arrays of length graph.nv, I write back:
    //     originalLabels, such that originalLabels[k][v] is the name in
    //          parent graph of vertex v in subgraph k
    //     newLabels, such that newLabels[i] is the name of
    //          vertex v in its subgraph
    public Decomposition splitGraph(Graph graph, int nColors) {
        // compute relabeling
        int[][] relabelUp = new int[nColors][];
        int[] relabelDown = new int[graph.nv];

        int[] subgraphSizes = new int[nColors];
        for (int i = 0; i < graph.nv; i++) {
            relabelDown[i] = subgraphSizes[colors[i]];
            subgraphSizes[colors[i]]++;
        }
        for (int i = 0; i < nColors; i++) {
            relabelUp[i] = new int[subgraphSizes[i]];
        }
        for (int i = 0; i < graph.nv; i++) {
            relabelUp[colors[i]][relabelDown[i]] = i;
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
            subgraphEdges[color] = new EdgeList(subgraphEdgeCounts[color]);
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
                edgePos[color]++;
            }
        }

        // convert edge lists to graphs
        Graph[] subgraphs = new Graph[nColors];
        for (int color = 0; color < nColors; color++)
            subgraphs[color] = new Graph(subgraphEdges[color]);

        return new Decomposition(subgraphs, relabelUp, relabelDown);
    }
}