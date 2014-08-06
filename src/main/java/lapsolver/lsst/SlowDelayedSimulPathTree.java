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

import java.util.Random;

public class SlowDelayedSimulPathTree implements SpanningTreeStrategy {
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

        double perturbEpsilon = 1. / graph.ne;

        EdgeList edges = new EdgeList(graph);

        double[] rates = new double[graph.ne];
        Random rand = new Random();
        for (int i = 0; i < graph.ne; i++)
            rates[i] = rand.nextDouble();

        int[] taken = new int[graph.ne];
        for (int i = 0; i < graph.ne; i++)
            taken[i] = 0;

        int[] comp = new int[graph.nv];
        int[] compSize = new int[graph.nv];
        for (int i = 0; i < graph.nv; i++) {
            comp[i] = i;
            compSize[i] = 1;
        }

        for (int treeCount = 0; treeCount < graph.nv - 1; treeCount++) {
            double vmin = -1;
            int index = 0;
            for (int i = 0; i < graph.ne; i++)
                if (taken[i] == 0 && comp[edges.u[i]] != comp[edges.v[i]]) {
                    int totalCompSize = compSize[comp[edges.u[i]]] + compSize[comp[edges.v[i]]];
                    if (costFunction(rates[i], edges.weight[i], totalCompSize) < vmin || vmin == -1) {
                        vmin = costFunction(rates[i], edges.weight[i], totalCompSize);
                        index = i;
                    }
                }
            taken[index] = 1;

            //System.out.println(index + " " + graph.nv + " " + graph.ne);

            //System.out.println(rates[index] + " " + edges.weight[index] + " " + (compSize[comp[edges.u[index]]] + compSize[comp[edges.v[index]]]));
            //System.out.println(costFunction(rates[index], edges.weight[index], compSize[comp[edges.u[index]]] + compSize[comp[edges.v[index]]]));

            // unite the components determined by the edge
            //System.out.println(comp[edges.u[index]] + " " + comp[edges.v[index]] + " " + compSize[comp[edges.u[index]]]);

            int c = comp[edges.u[index]];
            for (int i = 0; i < graph.nv; i++)
                if (comp[i] == c)
                    comp[i] = comp[edges.v[index]];

            // update the component sizes
            compSize[edges.v[index]] += compSize[c];
            compSize[c] = 0;
            c = comp[edges.v[index]];

            // gets the smallest rate in the component
            double bestRate = -1;
            for (int i = 0; i < graph.ne; i++)
                if (comp[edges.u[i]] == comp[c] && comp[edges.v[i]] == comp[c]) {
                    if (bestRate == -1 || bestRate > rates[i]) {
                        bestRate = rates[i];
                    }
                }

            // update all edges in the component with the smallest rate
            for (int i = 0; i < graph.ne; i++)
                if (comp[edges.u[i]] == comp[c] && comp[edges.v[i]] == comp[c]) {
                    rates[i] = bestRate;
                }

            // update edges surrounding the component
            for (int i = 0; i < graph.ne; i++)
                if ((comp[edges.u[i]] == c || comp[edges.v[i]] == c) && (comp[edges.u[i]] != comp[edges.v[i]])) {
                    double perturbation = rand.nextDouble() * perturbEpsilon;
                    if (rates[i] > bestRate + perturbation)
                        rates[i] = bestRate + perturbation;

                    //if (rates[i] > edges.weight[index] + perturbation)
                    //    rates[i] = edges.weight[index] + perturbation;
                }

            ijvI[treeCount] = edges.u[index];
            ijvJ[treeCount] = edges.v[index];
            ijvV[treeCount] = edges.weight[index];
        }

        return new EdgeList(ijvI, ijvJ, ijvV);
    }


    double costFunction(double rate, double weight, int compSize) {
        return rate * weight * compSize * Math.log(compSize);
    }
}
