/**
 * @file SlowDelayedSimulPathTree.java
 * @author Serban Stan <serban.stan@yale.edu>
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Tue Aug 5 2014
 *
 * Implements slowing down components in Dan's spanning tree algorithm, slow version.
 */

package lapsolver.lsst;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.UnionFind;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;
import java.util.TreeSet;

public class SlowDelayedSimulPathTree implements SpanningTreeStrategy {
    // times[node] is time at which node should fire
    public double[] times;
    public double[] rates;
    // public int[] pArray;
    public int[] ijvI;
    public int[] ijvJ;
    public double[] ijvV;

    @Override
    public Tree getTree(Graph graph) {
        return new Tree(getTreeEdges(graph));
    }

    public EdgeList getTreeEdges(Graph graph) {
//        Logger logger = new Logger();
//        logger.start("SimulPathTree.log");

        ijvI = new int[graph.nv - 1];
        ijvJ = new int[graph.nv - 1];
        ijvV = new double[graph.nv - 1];

        Random rand = new Random();

        times = new double[graph.ne];
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
                if (X.time == Y.time) {
                    if (X.u == Y.u) return Integer.compare(X.v, Y.v);
                    return Integer.compare(X.u, Y.u);
                }
                return Double.compare(X.time, Y.time);
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
                    times[e] = rate * wt;

                    EdgeEvent ev = new EdgeEvent(u, v, rate * wt, rate, wt);
                    events[e] = ev;
                    pq.add(ev);
                }
            }

        int ijvInd = 0;

        UnionFind uf = new UnionFind(graph.nv);
        int[] componentSizes = new int[graph.nv];
        Arrays.fill(componentSizes, 1);

        while (ijvInd < graph.nv - 1) {
            EdgeEvent ev = pq.pollFirst();

            int u = ev.u;
            int v = ev.v;

            // if in different comps, add that edge
            if (uf.find(u) != uf.find(v)) {
                int comp = uf.find(v);
                componentSizes[comp] += componentSizes[uf.find(u)];
                uf.union(u, v);

                for (int i = 0; i < graph.ne; i++) {
                    if (events[i] == null) continue;
                    if ( (uf.find(events[i].u) == comp) || (uf.find(events[i].v) == comp) ) {
                        double t = events[i].time + events[i].rate * componentSizes[comp] * events[i].wt;
                        times[i] = t;
                        pq.remove(events[i]);
                        EdgeEvent ev2 = new EdgeEvent(events[i].u, events[i].v, t, events[i].rate, events[i].wt);
                        pq.add(ev2);
                        events[i] = ev2;
                    }
                }

                ijvI[ijvInd] = u;
                ijvJ[ijvInd] = v;
                ijvV[ijvInd] = ev.wt;
                ijvInd++;

            }

            // for each edge attached to u, see if gives a lower time
            for (int i = 0; i < graph.deg[u]; i++) {
                int e = edgeNums[u][i];
                double wt = graph.weights[u][i];
                double t = ev.time + ev.rate * wt;

                t += rand.nextDouble() * perturbEpsilon;

                if (t < times[e]) {
                    times[e] = t;

                    // if was on pq, remove it
                    if (events[e] != null)
                        pq.remove(events[e]);

                    EdgeEvent ev2 = new EdgeEvent(u, graph.nbrs[u][i], t, ev.rate, wt);
                    pq.add(ev2);
                    events[e] = ev2;
                }
            }

            // for each edge attached to v, see if gives a lower time
            for (int i = 0; i < graph.deg[v]; i++) {
                int e = edgeNums[v][i];
                double wt = graph.weights[v][i];
                double t = ev.time + ev.rate * wt;

                t += rand.nextDouble() * perturbEpsilon;

                if (t < times[e]) {
                    times[e] = t;

                    // if was on pq, remove it
                    if (events[e] != null)
                        pq.remove(events[e]);

                    EdgeEvent ev2 = new EdgeEvent(v, graph.nbrs[v][i], t, ev.rate, wt);
                    pq.add(ev2);
                    events[e] = ev2;
                }
            }

        }

        return new EdgeList(ijvI, ijvJ, ijvV);
    }

    public class NodeEvent {
        int node;
        int from;
        double time;
        double rate;
        double wt;

        public NodeEvent(int node, int from, double time, double rate, double wt) {
            this.node = node;
            this.from = from;
            this.time = time;
            this.rate = rate;
            this.wt = wt;
        }
    }

    public class EdgeEvent {
        int u;
        int v;
        double time;
        double rate;
        double wt;

        public EdgeEvent(int u, int v, double time, double rate, double wt) {
            this.u = u;
            this.v = v;
            this.time = time;
            this.rate = rate;
            this.wt = wt;
        }
    }
}
