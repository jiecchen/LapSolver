/**
 * @file GraphVertexRemoval.java
 * @author Serban Stan <serban.stan@yale.edu>
 * @date Mon Jun 9 2014
 *
 * For an input graph graph removes degree one and degree two vertices. Returns a list of the
 * removed vertices and their number.
 */

package lapsolver.algorithms;

import lapsolver.Graph;

public class GraphVertexRemoval {
    public final Graph graph;
    public final int N;

    public final int[] removed;
    public final int[] updatedDeg;

    /**
     * @param g The input graph (weights are either conductances or resistances).
     */
    public GraphVertexRemoval(Graph g) {
        this.graph = g;
        this.N = g.nv;
        this.removed = new int[N];

        this.updatedDeg = graph.deg.clone();
    }

    public static class AnswerPair {
        public int[] permutation;
        public int numRemoved;

        public AnswerPair() {
        }

        public AnswerPair(int[] elementList, int value) {
            this.permutation = elementList.clone();
            this.numRemoved = value;
        }
    }

    public AnswerPair solve() {
        AnswerPair deg1rem = removeDegreeOne();
        if (deg1rem.numRemoved == N) {
            // The given graph is a tree. Will consider one vertex left behind.
            deg1rem.permutation[--deg1rem.numRemoved] = 0;
            return constructPerm(deg1rem);
        }

        AnswerPair deg2rem = removeDegreeTwo();

        AnswerPair finalAnswer = new AnswerPair();
        finalAnswer.numRemoved = deg1rem.numRemoved + deg2rem.numRemoved;
        finalAnswer.permutation = new int[N];

        System.arraycopy(deg1rem.permutation, 0, finalAnswer.permutation, 0, deg1rem.numRemoved);
        System.arraycopy(deg2rem.permutation, 0, finalAnswer.permutation, deg1rem.numRemoved, deg2rem.numRemoved);
        finalAnswer = constructPerm(finalAnswer);

        return finalAnswer;
    }

    /**
     * Remove the degree two vertices.
     */
    public AnswerPair removeDegreeTwo() {
        AnswerPair isDeg2Cycle = checkForDeg2Cycle();
        if (isDeg2Cycle.numRemoved != -1)
            return isDeg2Cycle;

        // Will eliminate successive chains of degree two
        return removeDeg2Chains();
    }

    /**
     * Removes the degree two vertices from the graph as chains.
     * If both ends of the chain are connected to the same vertex with degree greater than two,
     * then return all the elements of the chain but one end.
     *
     * @return the elements of the chain
     */
    public AnswerPair removeDeg2Chains() {
        int[] unRemovable = new int[N];
        for (int i = 0; i < N; i++)
            if (updatedDeg[i] > 2)
                unRemovable[i] = 1;

        int[] chains = new int[N];
        int index = 0;

        for (int i = 0; i < N; i++)
            if (unRemovable[i] == 1) {
                for (int j = 0; j < graph.deg[i]; j++) {
                    int u = graph.nbrs[i][j];

                    int newIndex = index;
                    if (updatedDeg[u] == 2 && removed[u] == 0) {
                        // This is a degree two chain
                        int chainCount = 0;

                        while (updatedDeg[u] == 2) {
                            removed[u] = 1;
                            chains[newIndex++] = u;

                            // I create the chains such that they have lengths at most 10
                            chainCount++;
                            if (chainCount % 3 == 0)
                                chains[newIndex]--;

                            int v = u;
                            for (int k = 0; k < graph.deg[u]; k++)
                                if (updatedDeg[graph.nbrs[u][k]] == 2 && removed[graph.nbrs[u][k]] == 0)
                                    v = graph.nbrs[u][k];

                            if (v == u) {
                                // I must exit the chain
                                for (int k = 0; k < graph.deg[u]; k++)
                                    if (updatedDeg[graph.nbrs[u][k]] > 2 && graph.nbrs[u][k] != i)
                                        v = graph.nbrs[u][k];

                                if (v == u) {
                                    // The chain is connected to the same exit vertex
                                    v = i;
                                }
                            }

                            u = v;
                        }

                        if (u == i) {
                            // If the chain is connected to the same outer vertex at both ends, then its end
                            // should be ignored
                            unRemovable[chains[newIndex - 1]] = 1;
                            chains[--newIndex] = 0;
                        }

                        index = newIndex;
                    }
                }
            }

        return new AnswerPair(chains, index);
    }

    /**
     * Checks if the graph composed by the remaining vertices is a degree 2 cycle.
     * @return all but two of the remaining vertices, if it is
     */
    public AnswerPair checkForDeg2Cycle() {
        AnswerPair ans = new AnswerPair();
        ans.numRemoved = -1;

        int maxDeg = 0;
        for (int i = 0; i < N; i++)
            if (updatedDeg[i] > maxDeg)
                maxDeg = updatedDeg[i];

        if (maxDeg > 2)
            return ans;

        ans = new AnswerPair(new int[N], 0);

        int last = 0;
        for (int i = 0; i < N; i++)
            if (updatedDeg[i] == 2)
                last = i;

        for (int i = last - 1; i >= 0; i--)
            if (updatedDeg[i] == 2) {
                last = i;
                break;
            }

        // Add the remaining vertices to the vertex pool
        for (int i = 0; i < last; i++)
            if (updatedDeg[i] == 2)
                ans.permutation[ans.numRemoved++] = i;

        return ans;
    }

    /**
     * Remove the degree one vertices
     * @return
     */
    public AnswerPair removeDegreeOne() {
        int[] bfsQueue = new int[N];
        int bfsQSize = 0;
        for (int i = 0; i < N; i++)
            if (graph.deg[i] == 1) {
                bfsQueue[bfsQSize++] = i;
                removed[i] = 1;
            }

        return bfsRemoval(bfsQueue, bfsQSize);
    }

    public AnswerPair bfsRemoval(int[] bfsQueue, int right) {
        int left = 0;
        while (left < right) {
            int u = bfsQueue[left++];

            for (int i = 0; i < graph.deg[u]; i++) {
                int v = graph.nbrs[u][i];
                updatedDeg[v]--;

                if (removed[v] == 0 && updatedDeg[v] == 1) {
                    bfsQueue[right++] = v;
                    removed[v] = 1;
                }
            }
        }

        return new AnswerPair(bfsQueue, right);
    }

    public AnswerPair constructPerm(AnswerPair ap) {
        int[] permutation = new int[N];
        int[] use = new int[N];

        for (int i = 0; i < ap.numRemoved; i++) {
            use[ap.permutation[i]] = 1;
            permutation[i] = ap.permutation[i];
        }

        int index = ap.numRemoved;
        for (int i = 0; i < N; i++)
            if (use[i] == 0)
                permutation[index++] = i;

        ap.permutation = permutation.clone();

        return ap;
    }
}