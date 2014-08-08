/**
 * @file RandomCut.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Fri Aug 8 2014
 *
 * Finds a random s-t cut in a graph, by building a maximal set of edges that leave s and t disconnected.
 */

package lapsolver.algorithms;

import lapsolver.EdgeList;
import lapsolver.Graph;

import java.util.Random;

public class RandomCut {
    // Compute an s-t cut.
    public static EdgeList compute(Graph graph, int s, int t) {
        UnionFind uf = new UnionFind(graph.nv);
        Random rand = new Random();
        EdgeList edges = new EdgeList(graph);
        EdgeList cut = new EdgeList(graph.ne);
        int cutSize = 0;

        // idea: first nGoodEdges elements of goodEdges can be sampled without disconnecting the graph
        int[] goodEdges = new int[graph.ne];
        for (int i = 0; i < graph.ne; i++) {
            goodEdges[i] = i;
        }
        int nGoodEdges = graph.ne;
        int temp; // swap space

        while (true) {
            // add all cut edges that we're forced to add
            // exclude all cut edges that don't matter
            int findS = uf.find(s), findT = uf.find(t);
            for (int i = 0; i < nGoodEdges; i++) {
                int findU = uf.find(edges.u[goodEdges[i]]), findV = uf.find(edges.v[goodEdges[i]]);
                if (findU == findV) {
                    // already in same component with excluded edges, so exclude this edge too
                    temp = goodEdges[i];
                    goodEdges[i] = goodEdges[nGoodEdges - 1];
                    goodEdges[nGoodEdges - 1] = temp;
                    nGoodEdges--; i--;
                }
                else if ((findS == findU && findT == findV) || (findS == findV && findT == findU)) {
                    // bridge between s and t, so add this edge to cut
                    cut.u[cutSize] = edges.u[goodEdges[i]];
                    cut.v[cutSize] = edges.v[goodEdges[i]];
                    cut.weight[cutSize] = edges.weight[goodEdges[i]];
                    cutSize++;

                    // swap with last good edge
                    temp = goodEdges[i];
                    goodEdges[i] = goodEdges[nGoodEdges - 1];
                    goodEdges[nGoodEdges - 1] = temp;
                    nGoodEdges--; i--;
                }
            }

            // check if we have a maximal cut
            if (nGoodEdges == 0) {
                break;
            }

            // pick a random edge to exclude from cut
            int sample = rand.nextInt(nGoodEdges);
            uf.union(edges.u[goodEdges[sample]], edges.v[goodEdges[sample]]);
            temp = goodEdges[sample];
            goodEdges[sample] = goodEdges[nGoodEdges - 1];
            goodEdges[nGoodEdges - 1] = temp;
            nGoodEdges--;
        }

        cut.resize(cutSize);
        return cut;
    }

    // Compute a cut between random vertices.
    public static EdgeList compute(Graph graph) {
        if (graph.nv == 1) {
            // avoid infinite loop in sampling
            return new EdgeList(0);
        }

        Random rand = new Random();
        int s = rand.nextInt(graph.nv);
        int t;

        do {
            t = rand.nextInt(graph.nv);
        } while(t == s);

        return compute(graph, s, t);
    }
}
