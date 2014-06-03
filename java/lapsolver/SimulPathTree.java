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

the only good routine in here right now is edgeGrow

example usage from Matlab:

a = del3Graph(10000);
[ai,aj,av] = find(tril(a));
g = WeightedGraph(ai,aj,av);
spt = SimulPathTree(g);
tr = spt.edgeGrow;
trt = tr.treeToTree;
trt.compTotalStretch(g)/length(ai)

 */

package lapsolver;

import java.util.Comparator;
import java.util.PriorityQueue;
import java.util.Random;

public class SimulPathTree {

    public Tree tree;

    public WeightedGraph g;
    // times[node] is time at which node should fire
    public double[] times;
    public double[] rates;
    // public int[] pArray;
    public int[] ijvI;
    public int[] ijvJ;
    public double[] ijvV;

    public SimulPathTree(WeightedGraph g) {
        this.g = g;
    }

    public WeightedGraph growTree() {
        Logger logger = new Logger();
        //	logger.start("SimulPathTree.log");

        ijvI = new int[g.nv - 1];
        ijvJ = new int[g.nv - 1];
        ijvV = new double[g.nv - 1];

        Random rand;
        rand = new Random();

        times = new double[g.nv];
        rates = new double[g.nv];
        NodeEvent[] events = new NodeEvent[g.nv];

        for (int i = 0; i < g.nv; i++) {
            rates[i] = rand.nextDouble();
            times[i] = Double.POSITIVE_INFINITY;
            events[i] = null;
        }

        int ijvInd = 0;

        PriorityQueue<NodeEvent> pq = new PriorityQueue<NodeEvent>(g.nv, new Comparator<NodeEvent>() {
            public int compare(NodeEvent X, NodeEvent Y) {
                return (X.time > Y.time ? 1 : -1);
            }
        });

        // load up the pq with events for every vertex
        for (int u = 0; u < g.nv; u++) {
            for (int i = 0; i < g.deg[u]; i++) {
                int v = g.nbrs[u][i];
                double wt = g.weights[u][i];
                double t = rates[u] / wt;

                if (t < times[v]) {
                    times[v] = t;

                    // if was on pq, remove it
                    if (events[v] != null)
                        pq.remove(events[v]);

                    NodeEvent ev = new NodeEvent(v, u, t, rates[u], wt);
                    pq.add(ev);
                    events[v] = ev;
                }
            }
        }

        UnionFind uf = new UnionFind(g.nv);

        while (ijvInd < g.nv - 1) {
            NodeEvent ev = pq.poll();
            int u = ev.node;
            int v = ev.from;

            logger.write(ev.time + ", (" + u + ", " + v + "), " + ev.rate);

            // if in different comps, add that edge
            if (uf.find(u) != uf.find(v)) {
                uf.union(u, v);

                ijvI[ijvInd] = u;
                ijvJ[ijvInd] = v;
                ijvV[ijvInd] = ev.wt;
                ijvInd++;

                logger.write("(" + u + ", " + v + ", " + ev.wt + ")");
            }

            // for each nbr of u, try to place it
            for (int i = 0; i < g.deg[u]; i++) {
                v = g.nbrs[u][i];
                if (uf.find(v) != uf.find(u)) {
                    double wt = g.weights[u][i];

                    double t = ev.time + ev.rate / wt;
                    if (t < times[v]) {
                        times[v] = t;

                        // if was on pq, remove it
                        // if (events[v] != null)
                        //    pq.remove(events[v]);

                        NodeEvent ev2 = new NodeEvent(v, u, t, ev.rate, wt);
                        pq.add(ev2);
                        events[v] = ev2;
                    }
                }
            }
        }

        WeightedGraph wg = new WeightedGraph();
        wg.fromMatlab(ijvI, ijvJ, ijvV);

        // return wg.treeToTree();
        return wg;
    }

    public WeightedGraph growTree3() {
        Logger logger = new Logger();
        logger.start("SimulPathTree.log");

        ijvI = new int[g.nv - 1];
        ijvJ = new int[g.nv - 1];
        ijvV = new double[g.nv - 1];

        Random rand;
        rand = new Random();

        double[] leastRate = new double[g.nv];
        times = new double[g.nv];
        rates = new double[g.nv];
        NodeEvent[] events = new NodeEvent[g.nv];

        for (int i = 0; i < g.nv; i++) {
            rates[i] = rand.nextDouble();
            times[i] = Double.POSITIVE_INFINITY;
            leastRate[i] = Double.POSITIVE_INFINITY;
            // events[i] = null;
        }

        int ijvInd = 0;

        PriorityQueue<NodeEvent> pq = new PriorityQueue<NodeEvent>(g.nv, new Comparator<NodeEvent>() {
            public int compare(NodeEvent X, NodeEvent Y) {
                return (X.time > Y.time ? 1 : -1);
            }
        });


        // load up the pq with events for every vertex
        for (int u = 0; u < g.nv; u++) {
            for (int i = 0; i < g.deg[u]; i++) {
                int v = g.nbrs[u][i];
                double wt = g.weights[u][i];
                double t = rates[u] / wt;

                if (rates[u] < leastRate[v]) {
                    leastRate[v] = rates[u];
                    NodeEvent ev = new NodeEvent(v, u, t, rates[u], wt);
                    pq.add(ev);

                }
            }
        }

        UnionFind uf = new UnionFind(g.nv);

        while (ijvInd < g.nv - 1) {
            NodeEvent ev = pq.poll();
            int u = ev.node;
            int v = ev.from;

            logger.write(ev.time + ", (" + u + ", " + v + "), " + ev.rate);

            // if in different comps, add that edge
            if (uf.find(u) != uf.find(v)) {
                uf.union(u, v);

                ijvI[ijvInd] = u;
                ijvJ[ijvInd] = v;
                ijvV[ijvInd] = ev.wt;
                ijvInd++;

                logger.write("(" + u + ", " + v + ", " + ev.wt + ")");
            }

            // for each nbr of u, try to place it
            for (int i = 0; i < g.deg[u]; i++) {
                v = g.nbrs[u][i];
                double wt = g.weights[u][i];
                double t = ev.time + ev.rate / wt;

                if (ev.rate < leastRate[v]) {
                    leastRate[v] = rates[u];
                    NodeEvent ev2 = new NodeEvent(v, u, t, rates[u], wt);
                    pq.add(ev2);
                }
            }
        }

        WeightedGraph wg = new WeightedGraph();
        wg.fromMatlab(ijvI, ijvJ, ijvV);

        // return wg.treeToTree();
        return wg;
    }

    public WeightedGraph edgeGrow() {
        Logger logger = new Logger();
        logger.start("SimulPathTree.log");

        ijvI = new int[g.nv - 1];
        ijvJ = new int[g.nv - 1];
        ijvV = new double[g.nv - 1];

        Random rand;
        rand = new Random();

        times = new double[g.ne];
        rates = new double[g.ne];
        EdgeEvent[] events = new EdgeEvent[g.ne];

        g.makeBackEdges();
        int edgeNum = 0;
        int[][] edgeNums = new int[g.nv][];
        for (int u = 0; u < g.nv; u++) {
            edgeNums[u] = new int[g.deg[u]];
            for (int i = 0; i < g.deg[u]; i++) {
                int v = g.nbrs[u][i];
                if (u < v)
                    edgeNums[u][i] = edgeNum++;
                else
                    edgeNums[u][i] = edgeNums[v][g.backInd[u][i]];
            }
        }

        for (int u = 0; u < g.nv; u++)
            for (int i = 0; i < g.deg[u]; i++)
                logger.write("(" + u + ", " + g.nbrs[u][i] + "), " + edgeNums[u][i]);


        PriorityQueue<EdgeEvent> pq = new PriorityQueue<EdgeEvent>(g.ne, new Comparator<EdgeEvent>() {
            public int compare(EdgeEvent X, EdgeEvent Y) {
                return (X.time > Y.time ? 1 : -1);
            }
        });


        for (int u = 0; u < g.nv; u++)
            for (int i = 0; i < g.deg[u]; i++) {
                int v = g.nbrs[u][i];
                // so, we include each edge just once
                if (u < v) {
                    double rate = rand.nextDouble();
                    int e = edgeNums[u][i];
                    double wt = g.weights[u][i];
                    times[e] = rate / wt;

                    EdgeEvent ev = new EdgeEvent(u, v, rate / wt, rate, wt);
                    events[e] = ev;
                    pq.add(ev);
                }
            }

        int ijvInd = 0;

        UnionFind uf = new UnionFind(g.nv);

        while (ijvInd < g.nv - 1) {
            EdgeEvent ev = pq.poll();

            int u = ev.u;
            int v = ev.v;

            // logger.write(ev.time + ", (" + u + ", " + v + "), " + ev.rate);

            // if in different comps, add that edge
            if (uf.find(u) != uf.find(v)) {
                uf.union(u, v);

                ijvI[ijvInd] = u;
                ijvJ[ijvInd] = v;
                ijvV[ijvInd] = ev.wt;
                ijvInd++;

                logger.write("(" + ev.u + ", " + ev.v + ", " + ev.wt + ")");
            }

            // for each edge attached to u, see if gives a lower time
            for (int i = 0; i < g.deg[u]; i++) {
                int e = edgeNums[u][i];
                double wt = g.weights[u][i];
                double t = ev.time + ev.rate / wt;
                if (t < times[e]) {
                    times[e] = t;

                    // if was on pq, remove it
                    if (events[e] != null)
                        pq.remove(events[e]);

                    EdgeEvent ev2 = new EdgeEvent(u, g.nbrs[u][i], t, ev.rate, wt);
                    pq.add(ev2);
                    events[e] = ev2;
                }
            }

            // for each edge attached to v, see if gives a lower time
            for (int i = 0; i < g.deg[v]; i++) {
                int e = edgeNums[v][i];
                double wt = g.weights[v][i];
                double t = ev.time + ev.rate / wt;
                if (t < times[e]) {
                    times[e] = t;

                    // if was on pq, remove it
                    if (events[e] != null)
                        pq.remove(events[e]);

                    EdgeEvent ev2 = new EdgeEvent(v, g.nbrs[v][i], t, ev.rate, wt);
                    pq.add(ev2);
                    events[e] = ev2;
                }
            }

        }

        WeightedGraph wg = new WeightedGraph();
        wg.fromMatlab(ijvI, ijvJ, ijvV);

        // return wg.treeToTree();
        return wg;

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
    public WeightedGraph growTree2() {

	ijvI = new int[g.nv-1];
	ijvJ = new int[g.nv-1];
	ijvV = new double[g.nv-1];
	
	Random rand;
	rand = new Random();

	double[] leastRate = new double[g.nv];
	rates = new double[g.nv];

	for (int i = 0; i < g.nv; i++) {
	    rates[i] = rand.nextDouble();
	    leastRate[i] = Double.POSITIVE_INFINITY;
	}

	int ijvInd = 0;
	
	PriorityQueue<NodeEvent> pq = new PriorityQueue<NodeEvent>(g.nv, new Comparator<NodeEvent>() {
		public int compare(NodeEvent X, NodeEvent Y) {
		    return (X.time > Y.time ? 1 : -1 );
		} });

	
	// load up the pq with events for every vertex
	for (int u = 0; u < g.nv; u++) {
	    for (int i = 0; i < g.deg[u]; i++) {
		int v = g.nbrs[u][i];
		double t = rates[u] * g.weights[u][i];

		if (rates[u] < leastRate[v]) {
		    leastRate[v] = rates[u];
		    NodeEvent ev = new NodeEvent(v, u, t, rates[u]);
		    pq.add(ev);
		    
		}
	    }
	}

	
	UnionFind uf = new UnionFind(g.nv);

	while (ijvInd < g.nv-1) {
	    NodeEvent ev = pq.poll();
	    int u = ev.node;

	    // for each nbr of u, try to place it
	    for (int i = 0; i < g.deg[u]; i++) {
		int v = g.nbrs[u][i];
		double wt = g.weights[u][i];

		// if in different comps, add that edge
		if (uf.find(u) != uf.find(v)) {
		    uf.union(u,v);

		    ijvI[ijvInd] = u;
		    ijvJ[ijvInd] = v;
		    ijvV[ijvInd] = wt;
		    ijvInd++;
		}

		double t = ev.time + ev.rate * wt;

		if (rates[u] < leastRate[v]) {
		    leastRate[v] = rates[u];
		    NodeEvent ev2 = new NodeEvent(v, u, t, rates[u]);
		    pq.add(ev2);
		}

	    }
	}
	
	WeightedGraph wg = new WeightedGraph();
	wg.fromMatlab(ijvI,ijvJ,ijvV);

	// return wg.treeToTree();
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
