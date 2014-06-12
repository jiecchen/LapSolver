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
        returnPair answer = new returnPair();

        // Initialize the answer pair
        initAns(answer, numSteps);

        // I will first construct the L matrix after removing the degree 1 vertices from the graph
        int[] currentDegree = new int[N];
        for (int i = 0; i < N; i++) {
            X[i] += graph.deg[i];
            currentDegree[i] = graph.deg[i];
        }

        for (int i = 0; i < numSteps; i++)
            if (currentDegree[i] == 1) {                    // Eliminate the degree one vertices
                for (int j = 0; j < graph.deg[i]; j++) {
                    int u = i;
                    int v = graph.nbrs[i][j];
                    double weight = graph.weights[i][j];

                    if (u < v) {
                        addToEdgeList(answer.L, v, u, -weight / X[u]);
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

                /**
                 * TODO: case when graph is a cycle with vertices of degree 2
                 */

                // outerStart and outerStop represent the higher degree vertices connected to the chain
                // these variables are initialized to 0 so I don't get any warning
                int outerStart = 0, outerStop = 0;
                double outerUpdt;

                if (start == stop) {
                    outerStart = graph.nbrs[start][0];
                    outerStop = graph.nbrs[start][1];
                }
                else {
                    for (int j = 0; j < 2; j++) {
                        if (graph.deg[graph.nbrs[start][j]] > 2)
                            outerStart = graph.nbrs[start][j];
                        if (graph.deg[graph.nbrs[stop][j]] > 2)
                            outerStop = graph.nbrs[stop][j];
                    }
                }
                outerUpdt = -graph.weights[start][0] / X[start];

                for (int u = start; u <= stop; u++) {
                    addToEdgeList(answer.L, outerStart, u, outerUpdt);

                    for (int k = 0; k < 2; k++) {
                        int v = graph.nbrs[u][k];
                        double weight = graph.weights[u][k];

                        if ((v > u && graph.deg[v] == 2) || (v == outerStop)) { // vertex v is the next vertex in the chain
                            addToEdgeList(answer.L, v, u, -weight / X[u]);
                            X[v] = X[v] - weight * weight * X[u];

                            outerUpdt = outerUpdt * (-weight / X[u]);
                        }
                    }
                }

                addToEdgeList(answer.L, outerStop, stop, outerUpdt);

                i = stop;
            }

        return answer;
    }

    public void initAns(returnPair LDLt, int steps) {
        int cnt = N;
        for (int i = 0; i < steps; i++) {
            cnt += graph.deg[i];
        }

        LDLt.L = new EdgeList(cnt);

        for (int i = 0; i < N; i++)
            addToEdgeList(LDLt.L, i, i, 1);
    }

    public void addToEdgeList(EdgeList el, int u, int v, double weight) {
        el.u[Lindex] = u;
        el.v[Lindex] = v;
        el.weight[Lindex] = weight;
        Lindex++;
    }

    // Check if two degree two vertices are neighbors
    public boolean isNeighbor(int u, int v) {
        if (v == graph.nbrs[u][0] || v == graph.nbrs[u][1])
            return true;
        return false;
    }
}
