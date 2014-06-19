package lapsolver.lsst;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.ShortestPathTree;

import java.util.*;

public class StarDecompositionWorker {

    public int[] getColors() {
        return colors;
    }

    private int[] colors; // scratch space for cut colorings
    private double[] coneCost; // scratch space for growCone (we need total O(n))

    public StarDecompositionWorker(Graph graph) {
        colors = new int[graph.nv];
        coneCost = new double[graph.nv];
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

    /**
     *
     * @param graph The input graph.
     * @param x0 Central vertex of the cut (source of ball).
     * @param sptInstance ShortestPathTree object, made by parent with call to
     *                    constructor on graph and x0.
     * @return The set of bridge edges going from the ball to the shell vertices.
     *         Returns null if no ball was found.
     */
    public EdgeList makeStarCut(Graph graph, int x0, ShortestPathTree sptInstance) {
        // initially nobody's part of the decomposition
        Arrays.fill(colors, -1);
        Arrays.fill(coneCost, Double.POSITIVE_INFINITY);

        Tree shortestPathTree = sptInstance.getTree();
        double[] dist = sptInstance.getDist();
        double radius = sptInstance.getRadius();

        // grow low-cut ball, build bridges from shell
        int[] ballShell = growBall(graph, shortestPathTree, dist, radius/3, 2*radius/3, 0);
        if (ballShell == null) return null;

        // grow low-cut cones from shell
        int nColors = 1;
        ArrayList<Integer> bridgeSources = new ArrayList<>();
        for (int coneCenter : ballShell) {
            if (colors[coneCenter] != -1)
                continue; // oops, another cone took this already

            // find radius of x_i (coneCenter) in the ideal
            int[] ideal = getIdeal(shortestPathTree, coneCenter);
            double eccentricity = 0;
            for (int u : ideal) {
                eccentricity = Math.max(eccentricity, dist[u] - dist[coneCenter]);
            }

            // grow the cone
            growCone(graph, shortestPathTree, dist, coneCenter, 2*radius/3, nColors);
            bridgeSources.add(coneCenter);
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
     * Find a low-cut ball of radius in the open interval (low, high)
     *
     * @param graph the containing graph
     * @param shortestPathTree the shortest path tree of the graph
     * @param dist  a distance array from a shortest path tree
     *              the position marked 0 is considered the root
     * @param minRadius the lower bound of the target radius
     * @param maxRadius the upper bound of the target radius
     * @param color the color of the ball (should be 0)
     * @return the set of vertices directly outside the ball (in the shortestPathTree's vertex labeling)
     */
    private int[] growBall(Graph graph, Tree shortestPathTree, final double[] dist,
                           double minRadius, double maxRadius, int color) {
        ArrayList<Integer> order = new ArrayList<>(graph.nv);
        for (int i = 0; i < graph.nv; i++) order.add(i, i);

        Collections.sort(order, new Comparator<Integer>() {
            @Override
            public int compare(Integer o1, Integer o2) {
                return Double.compare(dist[o1], dist[o2]);
            }
        });

        boolean[] inCut = new boolean[graph.nv];
        int bestCut = -1;
        double bestRatio = Double.POSITIVE_INFINITY;

        double cutCost = 0, cutVolume = 0;
        for (int i = 0; i < order.size() && dist[order.get(i)] < maxRadius; i++) {
            int u = order.get(i);
            inCut[u] = true;

            for (int j = 0; j < graph.deg[u]; j++) {
                int v = graph.nbrs[u][j];
                if (inCut[v]) { // this edge is leaving the cut, entering the volume
                    cutCost -= graph.weights[u][j];
                    cutVolume += graph.weights[u][j];
                }
                else { // this edge is entering the cut
                    cutCost += graph.weights[u][j];
                }
            }

            if (dist[u] > minRadius && cutCost/cutVolume < bestRatio) {
                bestCut = i;
                bestRatio = cutCost/cutVolume;
            }
        }

        // If we couldn't find a cut, return null -- replace with exception?
        if (bestCut == -1) return null;

        // reset the array -- at least we don't have to reallocate memory
        Arrays.fill(inCut, false);

        boolean[] inShell = new boolean[graph.nv];
        for (int i = 0; i <= bestCut; i++) {
            int u = order.get(i);
            inCut[u] = true;
            inShell[u] = false;

            colors[u] = color; // mark the ball

            for (int v : shortestPathTree.children[u])
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

    /**
     * Grows a low-cut cone in the graph.
     *
     * @param graph The input graph.
     * @param shortestPathTree The shortest path tree on the graph's x0.
     * @param dist The array of distances from x0.
     * @param source The center of this cone (i.e. x_color, color > 0)
     * @param maxRadius The maximum allowed eccentricity `source` in this cone.
     * @param color The color (identifier) of this cone.
     */
    private void growCone(Graph graph, Tree shortestPathTree, double[] dist,
                          int source, double maxRadius, int color) {

        // expand node that would cause cone to have smallest radius
        // similar to Dijkstra's algorithm
        PriorityQueue<Integer> toGrow = new PriorityQueue<>(graph.nv, new Comparator<Integer>() {
            public int compare(Integer x, Integer y) {
                return Double.compare(coneCost[x], coneCost[y]);
            }
        });

        // start with ideal as cone of radius 0
        coneCost[source] = 0;
        toGrow.add(source);

        // cone will be a union of the first cutIndex ideals
        int cutIndex = 0, bestCut = -1;
        double cutCost = 0, cutVolume = 0, bestRatio = 0;
        ArrayList<int[]> ideals = new ArrayList<>();

        while (!toGrow.isEmpty()) {
            int next = toGrow.poll();
            if (colors[next] != -1) continue;

            int[] ideal = getIdeal(shortestPathTree, next);
            ideals.add(ideal);

            if (coneCost[next] > 0) {
                // see if we've exceeded max allowed radius
                double newRadius = 0;
                for (int u : ideal) {
                    double ecc = coneCost[next] + dist[u] - dist[source];
                    newRadius = Math.max(newRadius, ecc);
                }
                if (newRadius > maxRadius) break;
            }

            // use colors as scratch space to mark this ideal
            for (int u : ideal) {
                colors[u] = -2;
            }

            // adjust cost and volume
            for (int u : ideal) {
                for (int i = 0; i < graph.deg[u]; i++) {
                    int v = graph.nbrs[u][i];
                    double w = graph.weights[u][i];

                    if (colors[v] == color) {
                        // edge is leaving cut, entering volume
                        cutCost -= w;
                        cutVolume += w;
                    }
                    else if (colors[v] == -2) {
                        // edge will always be in volume (avoid double-counting)
                        if (u < v) cutVolume += w/2;
                    }
                    else {
                        // edge is entering cut
                        cutCost += w;
                    }
                }
            }

            // now mark the ideal as part of the cone
            for (int u : ideal) {
                colors[u] = color;
            }

            // expand frontier
            for (int u : ideal) {
                coneCost[u] = coneCost[next];
                for (int i = 0; i < graph.deg[u]; i++) {
                    int v = graph.nbrs[u][i];
                    double w = graph.weights[u][i];
                    if (colors[v] != -1) continue;

                    // add new vertex with priority equal to radius of cone
                    // if we were to add it
                    if (coneCost[v] == Double.POSITIVE_INFINITY) {
                        toGrow.add(v);
                    }
                    coneCost[v] = Math.min(coneCost[v], coneCost[u] + w);
                }
            }

            // identify best cut
            if (bestCut == -1 || cutCost/cutVolume < bestRatio) {
                bestCut = cutIndex;
                bestRatio = cutCost/cutVolume;
            }

            cutIndex++;
        }

        // uncolor vertices not in best conex
        for (int i = ideals.size()-1; i > bestCut; i--) {
            for (int u : ideals.get(i)) {
                colors[u] = -1;
            }
        }

        // reset scratch space
        for (int[] ideal : ideals) {
            for (int u : ideal) {
                // reset radius priorities
                coneCost[u] = Double.POSITIVE_INFINITY;
                for (int v : graph.nbrs[u]) {
                    coneCost[v] = Double.POSITIVE_INFINITY;
                }
            }
        }
    }

    /**
     * Gets a subtree of the rooted shortest path tree.
     *
     * @param shortestPathTree The shortest path tree of the graph.
     * @param source The root of the ideal.
     * @return The set of vertices in the ideal, in BFS order.
     */
    public int[] getIdeal(Tree shortestPathTree, int source) {
        Queue<Integer> bfsQueue = new LinkedList<>();
        ArrayList<Integer> ideal = new ArrayList<>();

        // do BFS
        bfsQueue.add(source);
        while (!bfsQueue.isEmpty()) {
            int u = bfsQueue.poll();
            if (colors[u] != -1) continue;
            ideal.add(u);
            for (int v : shortestPathTree.children[u])
                bfsQueue.add(v);
        }

        // convert to primitive array
        int[] idealArray = new int[ideal.size()];
        for (int i = 0; i < idealArray.length; i++) {
            idealArray[i] = ideal.get(i);
        }
        return idealArray;
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