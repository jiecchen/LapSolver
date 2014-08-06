/**
 * @file SimulPathTree.java
 * @author Daniel Spielman <spielman@cs.yale.edu>
 * @date Fri May 30 2014
 *
 * An attempt at a low-stretch spanning tree, formed by
 * growing shortest path trees at different rates.
 *
 *  javac -classpath ~/rep/YINSmat/java/:. SimulPathTree.java
 *
 */

package lapsolver.lsst;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.UnionFind;

import java.util.Comparator;
import java.util.Random;
import java.util.TreeSet;

public class SlowDelayedSimulPathTree implements SpanningTreeStrategy {
    public double[] edgeCosts;
    public double[] rates;
    public int[] ijvI;
    public int[] ijvJ;
    public double[] ijvV;

    @Override
    public Tree getTree(Graph graph) {
        return new Tree(getTreeEdges(graph));
    }

    public EdgeList getTreeEdges(Graph graph) {
        ijvI = new int[graph.nv - 1];
        ijvJ = new int[graph.nv - 1];
        ijvV = new double[graph.nv - 1];

        Random rand = new Random();

        edgeCosts = new double[graph.ne];
        rates = new double[graph.ne];
        EdgeEvent[] events = new EdgeEvent[graph.ne];

        int edgeNum = 0;
        int[][] edgeNums = new int[graph.nv][];
        for (int u = 0; u < graph.nv; u++) {
            edgeNums[u] = new int[graph.deg[u]];
            for (int i = 0; i < graph.deg[u]; i++) {
                int v = graph.nbrs[u][i];
                if (u < v)
                    edgeNums[u][i] = edgeNum++;
                else
                    edgeNums[u][i] = edgeNums[v][graph.backInd[u][i]];
            }
        }

        double perturbEpsilon = 1. / graph.ne;

        TreeSet<EdgeEvent> pq = new TreeSet<>(new Comparator<EdgeEvent>() {
            public int compare(EdgeEvent X, EdgeEvent Y) {
                if (costFunction(X) == costFunction(Y)) {
                    if (X.u == Y.u) return Integer.compare(X.v, Y.v);
                    return Integer.compare(X.u, Y.u);
                }
                return Double.compare(costFunction(X), costFunction(Y));
            }
        });

        for (int u = 0; u < graph.nv; u++)
            for (int i = 0; i < graph.deg[u]; i++) {
                int v = graph.nbrs[u][i];
                // so, we include each edge just once
                if (u < v) {
                    double rate = rand.nextDouble();
                    int e = edgeNums[u][i];
                    double wt = graph.weights[u][i];
                    edgeCosts[e] = rate * wt;

                    EdgeEvent ev = new EdgeEvent(u, v, 1, rate, wt);
                    events[e] = ev;
                    pq.add(ev);
                }
            }

        int ijvInd = 0;

        UnionFind uf = new UnionFind(graph.nv);

       // for (EdgeEvent it : pq)
         //   System.out.println(it.u + " " + it.v + " " + it.number + " " + it.rate + " " + it.wt);

        while (ijvInd < graph.nv - 1) {
            EdgeEvent ev = pq.pollFirst();

            int u = ev.u;
            int v = ev.v;

            // if in different comps, add that edge
            if (uf.find(u) != uf.find(v)) {
                uf.union(u, v);

                ijvI[ijvInd] = u;
                ijvJ[ijvInd] = v;
                ijvV[ijvInd] = ev.wt;
                ijvInd++;
            }

            ev.number = 0;
            for (int i = 0; i < graph.nv; i++)
                if (uf.find(i) == uf.find(u))
                    ev.number++;

//            System.out.println(ev.number + " " + ev.wt + " " + ev.rate);

            // for each edge attached to u, see if gives a lower time
            for (int i = 0; i < graph.deg[u]; i++) {
                int e = edgeNums[u][i];
                double wt = graph.weights[u][i];
                double c = costFunction(ev);
                double per = rand.nextDouble() * perturbEpsilon;

                if (c < edgeCosts[e]) {
                    edgeCosts[e] = c;

                    // if was on pq, remove it
                    if (events[e] != null)
                        pq.remove(events[e]);

                    EdgeEvent ev2 = new EdgeEvent(u, graph.nbrs[u][i], ev.number, ev.rate + per, wt);
                    pq.add(ev2);
                    events[e] = ev2;
                }
            }

            // for each edge attached to v, see if gives a lower time
            for (int i = 0; i < graph.deg[v]; i++) {
                int e = edgeNums[v][i];
                double wt = graph.weights[v][i];
                double c = costFunction(ev);
                double per = rand.nextDouble() * perturbEpsilon;

                if (c < edgeCosts[e]) {
                    edgeCosts[e] = c;

                    // if was on pq, remove it
                    if (events[e] != null)
                        pq.remove(events[e]);

                    EdgeEvent ev2 = new EdgeEvent(v, graph.nbrs[v][i], ev.number, ev.rate + per, wt);
                    pq.add(ev2);
                    events[e] = ev2;
                }
            }

        }

        return new EdgeList(ijvI, ijvJ, ijvV);
    }

    public double costFunction(SlowDelayedSimulPathTree.EdgeEvent ev) {
        return ev.rate * Math.sqrt(ev.number) * ev.wt;
    }

    public class EdgeEvent {
        int u;
        int v;
        int number;
        double rate;
        double wt;

        public EdgeEvent(int u, int v, int number, double rate, double wt) {
            this.u = u;
            this.v = v;
            this.number = number;
            this.rate = rate;
            this.wt = wt;
        }
    }

}
