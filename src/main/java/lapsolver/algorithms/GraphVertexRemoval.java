/**
 * @file GraphVertexRemoval.java
 * @author Serban Stan <serban.stan@yale.edu>
 * @date Mon Jun 9 2014
 *
 * For an input graph graph removes degree one and degree two vertices. Returns a list of the
 * removed vertices.
 */

package lapsolver.algorithms;

import lapsolver.Graph;

import java.lang.reflect.Array;
import java.util.Arrays;

public class GraphVertexRemoval {
    public Graph graph;
    public int N;
    public int[] auxiliary;
    public int auxiliarySize;

    public GraphVertexRemoval(Graph inputG) {
        this.graph = new Graph(inputG);
        this.N = graph.nv;
        this.auxiliary = new int[N];
    }

    /*
        A subclass that helps the solve function return an (array,int) pair
     */
    public class AnswerPair {
        public int[] v;
        public int n;

        public AnswerPair(int[] elementList, int value) {
            this.n = value;
            this.v = new int[N];
            System.arraycopy(elementList, 0, this.v, 0, N);
        }
    }

    /*
        This is the class's main function. It returns an array of elements representing the order of the
        removed vertices as a prefix of a N size permutation. It also return the size of this prefix.
     */
    public AnswerPair solve() {
        int[] updatedDegree = new int[N];
        System.arraycopy(graph.deg, 0, updatedDegree, 0, N);

        boolean[] eliminated = new boolean[N];

        /*
            Remove the degree one vertices from the graph. Reorder the vertices
         */
        int[] answer = new int[N];
        int answerCnt = 0;

        removeDegreeOnes(updatedDegree, eliminated);
        int[] dfsQueue = new int[N];
        int dfsQSize = 0;

        int[] canUse = new int[N];
        for (int i = 0; i < auxiliarySize; i++) {
            int v = auxiliary[i];

            int willAdd = 0;
            for (int j = 0; j < graph.deg[v]; j++)
                if (updatedDegree[graph.nbrs[v][j]] >= 2)
                    willAdd = 1;
            if (willAdd == 1)
                dfsQueue[dfsQSize++] = v;

            canUse[v] = 1;
        }
        if (dfsQSize == 0)                              // In case the tree graph is a tree
            dfsQueue[dfsQSize++] = auxiliary[0];

        answer = dfsOrdering(dfsQueue, dfsQSize, canUse);
        answerCnt = auxiliarySize;

        /*
            Remove the degree two vertices from the graph. Reorder these vertices as well. Concatenate with previous answer.
         */
        auxiliary = new int[N];
        auxiliarySize = 0;

        removeDegreeTwos(updatedDegree, eliminated);

        dfsQueue = new int[auxiliarySize];
        canUse = new int[N];
        for (int i = 0; i < auxiliarySize; i++)
            canUse[auxiliary[i]] = 1;

        System.arraycopy(auxiliary, 0, answer, answerCnt, auxiliarySize);
        answerCnt += auxiliarySize;

        return new AnswerPair(buildPermutation(answer, answerCnt), answerCnt);
    }

    public void removeDegreeTwos(int[] updatedDegree, boolean[] eliminated) {
        int[] cantUse = new int[N];
        for (int i = 0; i < N; i++)
            if (updatedDegree[i] > 2)
                cantUse[i] = 1;

        int maxDeg = 0;
        for (int i = 0; i < N; i++)
            if (graph.deg[i] > maxDeg)
                maxDeg = graph.deg[i];

        if (maxDeg > 2) {
            // Case I, the given graph is not a degree two cycle
            for (int i = 0; i < N; i++)
                if (updatedDegree[i] > 2) {
                    // Find the start of a degree two chain
                    for (int j = 0; j < graph.deg[i]; j++) {
                        int neighbor = graph.nbrs[i][j];

                        if (graph.deg[neighbor] == 2 && !eliminated[neighbor]) {
                            // Eliminate the whole chain neighbor is part of
                            bfsDeg2Path(neighbor, updatedDegree, eliminated);
                        }
                    }
                }
        }
        else {
            // Case II, the given graph is a degree two cycle
            for (int i = 0; i < N; i++)
                if (graph.deg[i] == 2)
                    auxiliary[auxiliarySize++] = i;
        }
    }

    void bfsDeg2Path(int v, int[] updatedDegree, boolean[] eliminated) {
        int[] queue = new int[N];
        int left = 0;
        int right = 0;

        queue[right++] = v;
        while (left < right) {
            v = queue[left++];

            eliminated[v] = true;
            auxiliary[auxiliarySize++] = v;

            for (int i = 0; i < graph.deg[v]; i++) {
                int neighbor = graph.nbrs[v][i];

                if (updatedDegree[neighbor] == 2 && !eliminated[neighbor])
                    queue[right++] = neighbor;
            }
        }
    }

    // Remove all vertices of degrees 0 (in case of a tree) or 1
    public void removeDegreeOnes(int[] updatedDegree, boolean[] eliminated) {
        int[] queue = new int[N];
        boolean[] inQ = new boolean[N];
        int left = 0;
        int right = 0;

        for (int i = 0; i < N; i++)
            if (graph.deg[i] < 2) {
                queue[right++] = i;
                inQ[i] = true;
            }

        while (left < right) {
            int v = queue[left++];

            eliminated[v] = true;
            auxiliary[auxiliarySize++] = v;

            for (int i = 0; i < graph.deg[v]; i++) {
                int neighbor = graph.nbrs[v][i];
                updatedDegree[neighbor]--;

                if (updatedDegree[neighbor] < 2 && !inQ[neighbor]) {
                    queue[right++] = neighbor;
                    inQ[neighbor] = true;
                }
            }
        }
    }

    public int[] dfsOrder;
    public int[] dfsVisited;
    public int dfsCount;

    public int[] dfsOrdering(int[] Q, int count, int[] canUse) {
        this.dfsOrder = new int[N];
        this.dfsVisited = new int[N];
        this.dfsOrder = new int[N];

        for (int index = 0; index < count; index++) {
            if (dfsVisited[Q[index]] == 0)
                dfs(Q[index], canUse);
        }

        // Reverse the DFS ordering vector, so leafs come first
        for (int i = 0; i < dfsCount / 2; i++) {
            int aux = dfsOrder[i];
            dfsOrder[i] = dfsOrder[dfsCount - 1 - i];
            dfsOrder[dfsCount - 1 - i] = aux;
        }

        return dfsOrder;
    }

    public void dfs(int vertex, int[] canUse) {
        if (dfsVisited[vertex] == 1)
            return;

        dfsVisited[vertex] = 1;
        dfsOrder[dfsCount++] = vertex;

        for (int i = 0; i < graph.deg[vertex]; i++)
            if (canUse[graph.nbrs[vertex][i]] == 1)
                dfs(graph.nbrs[vertex][i], canUse);
    }

    /*
        Adds to vertexList the rest of the elements in range [1..N]
     */
    public int[] buildPermutation(int[] vertexList, int cnt) {
        int[] ret = new int[N];
        int[] use = new int[N];

        for (int i = 0; i < cnt; i++) {
            use[vertexList[i]] = 1;
            ret[i] = vertexList[i];
        }

        int aux = cnt;
        for (int i = 0; i < N; i++)
            if (use[i] == 0)
                ret[aux++] = i;

        return ret;
    }
}