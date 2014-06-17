/**
 * @file LDLDecomposition.java
 * @author Serban Stan <serban.stan@yale.edu>
 * @date Tue Jun 10 2014
 *
 * Given a graph and an array x (which represents the diagonal of a matrix), returns the LDL factorization (L and D) after
 * a given number of steps
 * */

package lapsolver.algorithms;

import lapsolver.Graph;
import lapsolver.EdgeList;
import lapsolver.util.GraphUtils;

import java.util.Arrays;

public class LDLDecomposition {
    public int N;
    public Graph graph;
    public double[] X;
    public int Lindex = 0;
    public int Dindex = 0;

    public class returnPair {
        public EdgeList L;
        public EdgeList D;

        public returnPair() {
        }

        public returnPair(EdgeList A, EdgeList B) {
            this.L = A;
            this.D = B;
        }
    }

    // Constructor for making full LDL factorization with graph and addition matrix
    public LDLDecomposition(Graph G, double[] diagValues) {
        this.graph = new Graph(G);
        this.N = graph.nv;
        this.X = new double[N];

        if (N != diagValues.length) {
            throw new Error("Graph and addition matrix have different sizes");
        }

        for (int i = 0; i < N; i++)
            this.X[i] = diagValues[i];
    }

    public returnPair solve(int numSteps) {
        EdgeList L = new EdgeList();
        EdgeList D = new EdgeList();

        // Initialize the answer pair
        initL(L, numSteps);
        initD(D, numSteps);

        // I will first construct the L matrix after removing the degree 1 vertices from the graph
        int[] currentDegree = new int[N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < graph.deg[i]; j++)
                X[i] = X[i] + 1 / graph.weights[i][j];
            currentDegree[i] = graph.deg[i];

            addToEdgeListD(D, i, i, X[i]);
        }

        for (int i = 0; i < numSteps; i++)
            if (currentDegree[i] == 1) {                    // Eliminate the degree one vertices
                for (int j = 0; j < graph.deg[i]; j++) {
                    int u = i;
                    int v = graph.nbrs[i][j];
                    double weight = 1 / graph.weights[i][j];

                    if (u < v) {
                        addToEdgeListL(L, v, u, -weight / X[u]);
                        addToEdgeListD(D, v, v, -weight * weight / X[u]);

                        X[v] = X[v] - weight * weight / X[u];

                        currentDegree[v]--;
                    }
                }
            }
            else {                                          // Eliminate the degree two vertices
                // Will eliminate a chain starting at 'start' and ending at 'stop'
                int start = i;
                int stop = i;
                while (stop < numSteps - 1 && isNeighbor(stop, stop + 1))
                    stop++;

                // outerStart and outerStop represent the higher degree vertices connected to the chain
                // these variables are initialized to 0 so I don't get any warning
                int outerStart = 0, outerStop = 0;

                if (start == stop) {
                    outerStart = graph.nbrs[start][0];
                    outerStop = graph.nbrs[start][1];
                }
                else {
                    for (int j = 0; j < 2; j++) {
                        if (graph.nbrs[start][j] >= numSteps)
                            outerStart = graph.nbrs[start][j];
                        if (graph.nbrs[stop][j] >= numSteps)
                            outerStop = graph.nbrs[stop][j];
                    }
                }

                double outerRatio = -1;
                double newEdgeCost = 0;
                double newValueInLap = 0;

                for (int k = 0; k < 2; k++)
                    if (graph.nbrs[start][k] == outerStart)
                        newValueInLap = 1 / graph.weights[start][k];

                for (int u = start; u <= stop; u++) {
                    for (int k = 0; k < 2; k++) {
                        int v = graph.nbrs[u][k];
                        double weight = 1 / graph.weights[u][k];

                        if ((v > u && currentDegree[v] == 2) || (v == outerStop)) { // vertex v is the next vertex in the chain
                            X[v] = X[v] - weight * weight / X[u];

                            addToEdgeListL(L, v, u, -weight / X[u]);
                            addToEdgeListD(D, v, v, -weight * weight / X[u]);

                            outerRatio = outerRatio * (weight / X[u]);
                            newEdgeCost = outerRatio * (-weight);

                            addToEdgeListD(D, outerStart, outerStart, outerRatio * newValueInLap);
                            newValueInLap = newValueInLap * (weight / X[u]);
                        }
                    }

                    addToEdgeListL(L, outerStart, u, outerRatio);
                }

                addToEdgeListD(D, outerStop, outerStart, -newEdgeCost);
                addToEdgeListD(D, outerStart, outerStop, -newEdgeCost);

                i = stop;
            }

        returnPair answer = new returnPair(L, D);
        return answer;
    }

    public void initL(EdgeList el, int steps) {
        // This number of declared elements is an overestimate for the actual needed memory. It has an extra O(steps) field used.
        int cnt = N;
        for (int i = 0; i < steps; i++)
            cnt += graph.deg[i] + 1;

        el.ne = cnt;
        el.u = new int[cnt];
        el.v = new int[cnt];
        el.weight = new double[cnt];

        for (int i = 0; i < N; i++)
            addToEdgeListL(el, i, i, 1);
    }

    public void initD(EdgeList el, int steps) {
        int cnt = N;    // The number of entries on the diagonal
        cnt += 2 * steps;
        for (int i = steps; i < N; i++)
            cnt += 2 * graph.deg[i];    // The number of edges already existing in L

        el.ne = cnt;
        el.u = new int[cnt];
        el.v = new int[cnt];
        el.weight = new double[cnt];

        for (int i = steps; i < N; i++)
            for (int j = 0; j < graph.deg[i]; j++)
                if (graph.nbrs[i][j] > i) {
                    addToEdgeListD(el, i, graph.nbrs[i][j], -1 / graph.weights[i][j]);
                    addToEdgeListD(el, graph.nbrs[i][j], i, -1 / graph.weights[i][j]);
                }
    }

    public void addToEdgeListL(EdgeList el, int u, int v, double weight) {
        el.u[Lindex] = u;
        el.v[Lindex] = v;
        el.weight[Lindex] = weight;
        Lindex++;
    }

    public void addToEdgeListD(EdgeList el, int u, int v, double weight) {
        el.u[Dindex] = u;
        el.v[Dindex] = v;
        el.weight[Dindex] = weight;
        Dindex++;
    }

    // Check if two degree two vertices are neighbors
    public boolean isNeighbor(int u, int v) {
        if (v == graph.nbrs[u][0] || v == graph.nbrs[u][1])
            return true;
        return false;
    }
}
