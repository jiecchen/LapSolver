/**
 * @file GraphVertexRemoval.java
 * @author Serban Stan <serban.stan@yale.edu>
 * @date Mon Jun 9 2014
 *
 * For an input graph G removes degree one and degree two vertices, and possibly degree three vertices. Returns a list of the
 * removed vertices. Can remove vertices of degrees 0...K, where K is given.
 */

package lapsolver.algorithms;

import lapsolver.Graph;

import java.util.Arrays;

public class GraphVertexRemoval {
    public Graph G;
    public int N;

    public GraphVertexRemoval(Graph inputG) {
        this.G = new Graph(inputG);
        this.N = G.nv;
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
            for (int i = 0; i < N; i++)
                this.v[i] = elementList[i];
        }
    }

    /*
        This is the class's main function. It returns an array of elements representing the order of the
        removed vertices as a prefix of a N size permutation. It also return the size of this prefix.
     */
    public AnswerPair solve(int maxRemovedDegree) {
        int[] answer = new int[N];
        int cntAnswer = 0;

        int[][] Q = new int[maxRemovedDegree + 1][N];
        int first[] = new int[maxRemovedDegree + 1];
        int last[] = new int[maxRemovedDegree + 1];

        int eliminated[] = new int[N];

        // A local copy of the vertices degrees. Will come in handy later.
        int updatedDegree[] = new int[N];
        for (int i = 0; i < N; i++)
            updatedDegree[i] = G.deg[i];

        // Introduces all vertices that can be initially removed into separate queues.
        for (int i = 0; i < N; i++) {
            if (G.deg[i] <= maxRemovedDegree) {
                int degree = G.deg[i];

                Q[degree][last[degree]] = i;
                last[G.deg[i]]++;
            }
        }

        // Pick vertices sequentially, and remove them from the graph
        int vertex = canPickVertex(Q, first, last, eliminated, maxRemovedDegree);
        while (vertex != -1) {
            eliminated[vertex] = 1;
            answer[cntAnswer++] = vertex;

            for (int i = 0; i < G.deg[vertex]; i++) {
                int neighbor = G.nbrs[vertex][i];
                updatedDegree[neighbor]--;

                if (updatedDegree[neighbor] <= maxRemovedDegree && eliminated[neighbor] == 0) {
                    int degree = updatedDegree[neighbor];

                    Q[degree][last[degree]] = neighbor;
                    last[degree]++;
                }
            }

            vertex = canPickVertex(Q, first, last, eliminated, maxRemovedDegree);
        }

        // Place the eliminated vertices in a cache-friendly sequencing
        answer = cacheFriendly(G, answer, cntAnswer);

        // Add the rest of the vertices in the graph to the list of vertices I will return,
        // so that answer will become a permutation
        answer = buildPermutation(answer, cntAnswer);

        AnswerPair result = new AnswerPair(answer, cntAnswer);
        return result;
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

    /*
        Returns the answer in an order that has consecutive nodes close-by in DFS/BFS tree search algorithms
     */
    public int[] cacheFriendly(Graph G, int[] retArray, int L) {
        int[] answer = new int[N];
        int[] isInGraph = new int[N];

        for (int i = 0; i < N; i++)
            isInGraph[i] = 1;
        for (int i = 0; i < L; i++)
            isInGraph[retArray[i]] = 0;

        // Will pick as starting vertices for the BFS the vertices that are neighbors
        // of some vertices that have not been removed

        int[] queue = new int[N];
        int count = 0;

        for (int i = 0; i < N; i++)
            if (isInGraph[i] == 0) {
                int hasNeighborInGraph = 0;
                for (int j = 0; j < G.deg[i]; j++) {
                    if (isInGraph[G.nbrs[i][j]] == 1) {
                        hasNeighborInGraph = 1;
                        break;
                    }
                }

                if (hasNeighborInGraph == 1) {
                    queue[count++] = i;
                }
            }

        // Will add the rest of the vertices in the queue
        for (int i = 0; i < N; i++)
            if (isInGraph[i] == 0) {
                int hasNeighborInGraph = 0;
                for (int j = 0; j < G.deg[i]; j++) {
                    if (isInGraph[G.nbrs[i][j]] == 1) {
                        hasNeighborInGraph = 1;
                        break;
                    }
                }

                if (hasNeighborInGraph == 0) {
                    queue[count++] = i;
                }
            }

        // This can be swapped between bsfOrdering/dfsOrdering
        //answer = bfsOrdering(G, queue, count, isInGraph);
        answer = dfsOrdering(G, queue, count, isInGraph);

        return answer;
    }

    public int[] dfsOrder;
    public int[] dfsVisited;
    public int dfsCount;

    public int[] dfsOrdering(Graph G, int[] Q, int count, int[] isInGraph) {
        this.dfsOrder = new int[N];
        this.dfsVisited = new int[N];
        this.dfsOrder = new int[N];

        for (int index = 0; index < count; index++) {
            if (dfsVisited[Q[index]] == 0)
                dfs(Q[index], G, isInGraph);
        }

        // Reverse the DFS ordering vector, so leafs come first
        for (int i = 0; i < dfsCount / 2; i++) {
            int aux = dfsOrder[i];
            dfsOrder[i] = dfsOrder[dfsCount - 1 - i];
            dfsOrder[dfsCount - 1 - i] = aux;
        }

        return dfsOrder;
    }

    public void dfs(int vertex, Graph G, int[] isInGraph) {
        if (dfsVisited[vertex] == 1)
            return;

        dfsVisited[vertex] = 1;
        dfsOrder[dfsCount++] = vertex;

        for (int i = 0; i < G.deg[vertex]; i++)
            if (isInGraph[G.nbrs[vertex][i]] == 0)
                dfs(G.nbrs[vertex][i], G, isInGraph);
    }

    public int[] bfsOrdering(Graph G, int[] Q, int count, int[] isInGraph) {
        int[] bfsOrder = new int[N];
        int[] visited = new int[N];
        int left = 0;
        int right = 0;

        for (int index = 0; index < count; index++) {
            if (visited[Q[index]] == 0) {
                bfsOrder[right++] = Q[index];
                visited[Q[index]] = 1;
            }

            // do Iterative BFS
            while (left < right) {
                int vertex = bfsOrder[left];
                int degree = G.deg[vertex];

                for (int i = 0; i < degree; i++) {
                    int neighbor = G.nbrs[vertex][i];
                    if (visited[neighbor] == 0 && isInGraph[neighbor] == 0) {
                        visited[neighbor] = 1;
                        bfsOrder[right++] = neighbor;
                    }
                }

                left++;
            }
        }

        // Reverse the BFS ordering vector, so leafs come first
        for (int i = 0; i < right / 2; i++) {
            int aux = bfsOrder[i];
            bfsOrder[i] = bfsOrder[right - 1 - i];
            bfsOrder[right - 1 - i] = aux;
        }

        return bfsOrder;
    }

    /*
        A utility function that can pick the vertex with smallest degree
     */
    public int canPickVertex(int[][] Q, int[] first, int[] last, int[] eliminated, int maxDegree) {
        for (int degree = 0; degree <= maxDegree; degree++) {
            while (first[degree] < last[degree]) {
                int ret = Q[degree][first[degree]];
                first[degree]++;

                // If the vertex was already eliminated with a lower degree, I am skipping it at this point
                if (eliminated[ret] == 0)
                    return ret;
            }
        }

        return -1;
    }
}
