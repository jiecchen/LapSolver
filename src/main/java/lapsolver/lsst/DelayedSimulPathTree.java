/**
 * @file DelayedSimulPathTree.java
 * @author Serban Stan <serban.stan@yale.edu>
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Tue Aug 5 2014
 *
 * Implements slowing down components in Dan's spanning tree algorithm
 */

package lapsolver.lsst;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.UnionFind;

import java.util.Comparator;
import java.util.Random;
import java.util.TreeSet;

public class DelayedSimulPathTree implements SpanningTreeStrategy {

    public double[] rates;

    @Override
    public Tree getTree(Graph graph) {
        Random rand = new Random();

        EdgeList edges = new EdgeList(graph.ne);

        // get edge indices
        int edgeIndices[][] = new int[graph.nv][];
        for (int i = 0; i < graph.nv; i++)
            edgeIndices[i] = new int[graph.deg[i]];

        int k = 0;
        for (int i = 0; i < graph.nv; i++)
            for (int j = 0; j < graph.deg[j]; j++) {
                edges.u[k] = i;
                edges.v[k] = graph.nbrs[i][j];
                edges.weight[k] = graph.weights[i][j];

                edgeIndices[i][j] = k++;
            }

        // declare rates
        rates = new double[graph.ne];
        for (int i = 0; i < graph.ne; i++)
            rates[i] = rand.nextDouble();

        // the graph.ne'th treeSet will hold all of the unpicked edges
        TreeSet <EdgeEvent>[] componentEvents = new TreeSet[graph.ne + 1];
        for (int i = 0; i <= graph.ne; i++)
            componentEvents[i] = new TreeSet <EdgeEvent>(EdgeEvent.comparator);
        for (int i = 0; i < graph.ne; i++)
            componentEvents[graph.ne].add(new EdgeEvent(edges.weight[i], i, graph.ne));

        // keep a big treeSet
        TreeSet <EdgeEvent> bestEvent = new TreeSet <EdgeEvent>(EdgeEvent.comparator);
        bestEvent.add(new EdgeEvent(new EdgeEvent(componentEvents[graph.ne].first())));

        // initialize components for each edge
        int[] componentU = new int[graph.ne];
        int[] componentV = new int[graph.ne];

        for (int i = 0; i < graph.ne; i++) {
            componentU[i] = -1;
            componentV[i] = -1;
        }

        // intialize union find
        UnionFind uf = new UnionFind(graph.nv);

        k = 0;
        EdgeList dspt = new EdgeList(graph.nv - 1);

        int compIndex = 0;
        while (k < graph.nv - 1) {
            EdgeEvent currentEdgeEvent = new EdgeEvent(bestEvent.first());
            int index = currentEdgeEvent.edgeIndex;
            int u = edges.u[currentEdgeEvent.edgeIndex];
            int v = edges.v[currentEdgeEvent.edgeIndex];
            double weight = edges.weight[currentEdgeEvent.edgeIndex];

            if (uf.find(u) == uf.find(v)) {
                // the edge adds nothing to the tree
            }
            else {
                // the edge either creates a new component, expands a component, or unites two components
                // case I, the edge creates a new component
                if (componentU[index] == -1 && componentV[index] == -1) {
                    componentU[index] = compIndex;
                    componentV[index] = compIndex;
                    compIndex++;


                } else {
                    // case II, the edge expands a component
                    if ((componentU[index] == -1 || componentV[index] == -1) && (componentU[index] > -1 || componentV[index] > -1)) {
                        if (componentU[index] == -1)
                            componentU[index] = componentV[index];
                        else
                            componentV[index] = componentU[index];


                    } else {
                        // case III, the edge unites two components

                    }
                }
            }

            // remove the edge from all of its components
            if (componentU[index] > -1) {
                currentEdgeEvent.compIndex = componentU[index];
                componentEvents[componentU[index]].remove(currentEdgeEvent);
            }

            if (componentV[index] > -1) {
                currentEdgeEvent.compIndex = componentV[index];
                componentEvents[componentV[index]].remove(currentEdgeEvent);
            }
        }

        return new Tree(dspt);
    }

    public static class EdgeEvent {
        public double cost;
        public int edgeIndex;
        public int compIndex;

        public EdgeEvent(double cost, int edgeIndex, int compIndex) {
            this.cost = cost;
            this.edgeIndex = edgeIndex;
            this.compIndex = compIndex;
        }

        public EdgeEvent(EdgeEvent edgeEvent) {
            this.cost = edgeEvent.cost;
            this.edgeIndex = edgeEvent.edgeIndex;
            this.compIndex = edgeEvent.compIndex;
        }

        public static Comparator<EdgeEvent> comparator = new Comparator<EdgeEvent>() {
            public int compare(EdgeEvent X, EdgeEvent Y) {
                if (X.cost == Y.cost)
                    return Integer.compare(X.edgeIndex, Y.edgeIndex);
                return Double.compare(X.cost, Y.cost);
            }
        };
    }
}
