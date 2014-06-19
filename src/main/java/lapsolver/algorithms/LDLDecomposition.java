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

public class LDLDecomposition {
    public int N;
    public Graph graph;

    public double[] X;
    public int[] updatedDeg;

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
        this.updatedDeg = new int[N];

        if (N != diagValues.length) {
            throw new Error("Graph and addition matrix have different sizes");
        }

        for (int i = 0; i < N; i++)
            this.X[i] = diagValues[i];
    }

    EdgeList L;
    EdgeList D;

    public returnPair solve(int numSteps) {
        L = new EdgeList();
        D = new EdgeList();

        // Initialize the answer pair
        initL(numSteps);
        initD(numSteps);

        // Update the matrices' diagonal values
        updateDiag();

        for (int i = 0; i < numSteps; i++)
            if (updatedDeg[i] == 1) {                       // Eliminate the degree one vertices
                remDeg1(i);
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

                // Pick outerStart and outerStop
                if (start == stop) {
                    for (int j = 0; j < graph.deg[start]; j++)
                        if (graph.nbrs[start][j] >= numSteps) {
                            outerStart = graph.nbrs[start][j];

                            for (int k = j + 1; k < graph.deg[start]; k++)
                                if (graph.nbrs[start][k] >= numSteps) {
                                    outerStop = graph.nbrs[start][k];
                                    break;
                                }

                            break;
                        }
                }
                else {
                    for (int j = 0; j < graph.deg[start]; j++)
                        if (graph.nbrs[start][j] >= numSteps) {
                            outerStart = graph.nbrs[start][j];
                            break;
                        }

                    for (int j = 0; j < graph.deg[stop]; j++)
                        if (graph.nbrs[stop][j] >= numSteps) {
                            outerStop = graph.nbrs[stop][j];
                            break;
                        }
                }

                removeDeg2Chain(start, stop, outerStart, outerStop);

                i = stop;
            }

        returnPair answer = new returnPair(L, D);
        return answer;
    }

    public void removeDeg2Chain(int start, int stop, int outerStart, int outerStop) {
//        System.out.println(outerStart + " > " + start + " " + stop + " < " + outerStop);

        // outerValue is used to compute the L values of type (outerStart, u)
        double outerValue = -1;
        for (int i = 0; i < graph.deg[start]; i++)
            if (graph.nbrs[start][i] == outerStart)
                outerValue = -1 / graph.weights[start][i];

        // newValueInLap is used to mimic the update of the laplacian matrix
        double newValueInLap = outerValue;

        // Compute the L values adjacent elements in the degree 2 chain
        for (int u = start; u <= stop; u++) {
            outerValue = outerValue / X[u];

            X[outerStart] = X[outerStart] - newValueInLap * outerValue;

            addToEdgeListL(outerStart, u, outerValue);
            addToEdgeListD(outerStart, outerStart, -newValueInLap * outerValue);

            for (int i = 0; i < graph.deg[u]; i++) {
                int v = graph.nbrs[u][i];
                double weight = 1 / graph.weights[u][i];

                if (updatedDeg[v] > 1 && v > u && v != outerStart) {
                    X[v] = X[v] - weight * weight / X[u];
                    addToEdgeListL(v, u, -weight / X[u]);
                    addToEdgeListD(v, v, -weight * weight / X[u]);

                    outerValue = outerValue * weight;
                    newValueInLap = newValueInLap * (weight / X[u]);
                }
            }
        }

        addToEdgeListD(outerStop, outerStart, newValueInLap);
        addToEdgeListD(outerStart, outerStop, newValueInLap);
    }

    public void updateDiag() {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < graph.deg[i]; j++)
                X[i] = X[i] + 1 / graph.weights[i][j];
            updatedDeg[i] = graph.deg[i];

            addToEdgeListD(i, i, X[i]);
        }
    }

    public void remDeg1(int k) {
        for (int i = 0; i < graph.deg[k]; i++) {
            int u = k;
            int v = graph.nbrs[k][i];
            double weight = 1 / graph.weights[k][i];

            if (u < v) {
                addToEdgeListL(v, u, -weight / X[u]);
                addToEdgeListD(v, v, -weight * weight / X[u]);

                X[v] = X[v] - weight * weight / X[u];

                updatedDeg[v]--;
            }
        }
    }

    public void initL(int steps) {
        // This number of declared elements is an overestimate for the actual needed memory. It has an extra O(steps) field used.
        int cnt = N;
        for (int i = 0; i < steps; i++)
            cnt += graph.deg[i] + 1;

        L.ne = cnt;
        L.u = new int[cnt];
        L.v = new int[cnt];
        L.weight = new double[cnt];

        for (int i = 0; i < N; i++)
            addToEdgeListL(i, i, 1);
    }

    public void initD(int steps) {
        int cnt = N;    // The number of entries on the diagonal
        cnt += 2 * steps;
        for (int i = steps; i < N; i++)
            cnt += 2 * graph.deg[i];    // The number of edges already existing in L

        D.ne = cnt;
        D.u = new int[cnt];
        D.v = new int[cnt];
        D.weight = new double[cnt];

        for (int i = steps; i < N; i++)
            for (int j = 0; j < graph.deg[i]; j++)
                if (graph.nbrs[i][j] > i) {
                    addToEdgeListD(i, graph.nbrs[i][j], -1 / graph.weights[i][j]);
                    addToEdgeListD(graph.nbrs[i][j], i, -1 / graph.weights[i][j]);
                }
    }

    public void addToEdgeListL(int u, int v, double weight) {
        L.u[Lindex] = u;
        L.v[Lindex] = v;
        L.weight[Lindex] = weight;
        Lindex++;
    }

    public void addToEdgeListD(int u, int v, double weight) {
        D.u[Dindex] = u;
        D.v[Dindex] = v;
        D.weight[Dindex] = weight;
        Dindex++;
    }

    // Check if two degree two vertices are neighbors
    public boolean isNeighbor(int u, int v) {
        for (int i = 0; i < graph.deg[u]; i++)
            if (graph.nbrs[u][i] == v)
                return true;
        return false;
    }
}
