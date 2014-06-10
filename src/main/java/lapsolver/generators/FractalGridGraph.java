/**
 * @file FractalTreeGenerator.java
 * @author Serban Stan <serban.stan@yale.edu>
 * @date Fri Jun 6 2014
 *
 * A fractal tree generator that whose result should be a target for other spanning tree generators.
 */

package lapsolver.generators;

import lapsolver.EdgeList;
import lapsolver.Graph;

public class FractalGridGraph {
    public int N;       // the number of vertices
    public int index;   // a global iterator to make recursion easier
    public int[] u;
    public int[] v;

    public FractalGridGraph(int vertexCount) {
        // the number of vertices is 2^N, the number of edges is 4^(N+1)-1

        this.N = (int) Math.pow(2, vertexCount);
        this.u = new int[N * N - 1];
        this.v = new int[N * N - 1];
        this.index = 0;
    }

    public Graph getFractalTree() {
        EdgeList edgeList = new EdgeList(N * N -1);

        cover(0, 0, N - 1, N - 1);

        for (int i = 0; i < u.length; i++) {
            edgeList.u[i] = u[i];
            edgeList.v[i] = v[i];
        }

        Graph answer = new Graph(edgeList);
        return answer;
    }

    public void cover(int x0, int y0, int x1, int y1) {
        if (x0 == x1 && y0 == y1) {
            return;
        }

        int X = x0 + (x1 - x0) / 2;
        int Y = y0 + (y1 - y0) / 2;

        u[index] = getIndex(X, Y);
        v[index] = getIndex(X + 1, Y);
        index++;

        u[index] = getIndex(X + 1, Y);
        v[index] = getIndex(X + 1, Y + 1);
        index++;

        u[index] = getIndex(X + 1, Y + 1);
        v[index] = getIndex(X, Y + 1);
        index++;

        cover(x0, y0, X, Y);
        cover(X + 1, y0, x1, Y);
        cover(x0, Y + 1, X, y1);
        cover(X + 1, Y + 1, x1, y1);
    }

    public int getIndex(int p, int q) {
        return p + q * N;
    }
}
