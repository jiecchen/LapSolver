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
import lapsolver.util.GraphUtils;
import lapsolver.util.Logger;
import lapsolver.Tree;
import lapsolver.algorithms.UnionFind;

import java.util.Comparator;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.TreeSet;

public class SimulPathTree implements SpanningTreeStrategy {
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

//        for (int u = 0; u < graph.nv; u++)
//            for (int i = 0; i < graph.deg[u]; i++)
//                logger.write("(" + u + ", " + graph.nbrs[u][i] + "), " + edgeNums[u][i]);


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

        while (ijvInd < graph.nv - 1) {
            EdgeEvent ev = pq.pollFirst();

            int u = ev.u;
            int v = ev.v;

            // logger.write(ev.time + ", (" + u + ", " + permutation + "), " + ev.rate);

            // if in different comps, add that edge
            if (uf.find(u) != uf.find(v)) {
                uf.union(u, v);

                ijvI[ijvInd] = u;
                ijvJ[ijvInd] = v;
                ijvV[ijvInd] = ev.wt;
                ijvInd++;

//                logger.write("(" + ev.u + ", " + ev.permutation + ", " + ev.wt + ")");
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

    /*
     * slower, but maybe better
     * so far, much worse on random regular graphs
     */
    /*
    public Graph growTree2() {

	ijvI = new int[graph.nv-1];
	ijvJ = new int[graph.nv-1];
	ijvV = new double[graph.nv-1];
	
	Random rand;
	rand = new Random();

	double[] leastRate = new double[graph.nv];
	rates = new double[graph.nv];

	for (int i = 0; i < graph.nv; i++) {
	    rates[i] = rand.nextDouble();
	    leastRate[i] = Double.POSITIVE_INFINITY;
	}

	int ijvInd = 0;
	
	PriorityQueue<NodeEvent> pq = new PriorityQueue<>(graph.nv, new Comparator<NodeEvent>() {
		public int compare(NodeEvent X, NodeEvent Y) {
		    return (X.time > Y.time ? 1 : -1 );
		} });

	
	// load up the pq with events for every vertex
	for (int u = 0; u < graph.nv; u++) {
	    for (int i = 0; i < graph.deg[u]; i++) {
		int permutation = graph.nbrs[u][i];
		double t = rates[u] * graph.weights[u][i];

		if (rates[u] < leastRate[permutation]) {
		    leastRate[permutation] = rates[u];
		    NodeEvent ev = new NodeEvent(permutation, u, t, rates[u]);
		    pq.add(ev);
		    
		}
	    }
	}

	
	UnionFind uf = new UnionFind(graph.nv);

	while (ijvInd < graph.nv-1) {
	    NodeEvent ev = pq.poll();
	    int u = ev.node;

	    // for each nbr of u, try to place it
	    for (int i = 0; i < graph.deg[u]; i++) {
		int permutation = graph.nbrs[u][i];
		double wt = graph.weights[u][i];

		// if in different comps, add that edge
		if (uf.find(u) != uf.find(permutation)) {
		    uf.union(u,permutation);

		    ijvI[ijvInd] = u;
		    ijvJ[ijvInd] = permutation;
		    ijvV[ijvInd] = wt;
		    ijvInd++;
		}

		double t = ev.time + ev.rate * wt;

		if (rates[u] < leastRate[permutation]) {
		    leastRate[permutation] = rates[u];
		    NodeEvent ev2 = new NodeEvent(permutation, u, t, rates[u]);
		    pq.add(ev2);
		}

	    }
	}
	
	Graph wg = new Graph();
	wg.fromMatlab(ijvI,ijvJ,ijvV);

	return wg;

    }
    */

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
