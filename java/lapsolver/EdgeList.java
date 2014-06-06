/**
 * @file EdgeList.java
 * @author Serban Stan <serban.stan@gmail.com>
 * @date Thu Jun 5 2014
 *
 * An EdgeList of the form (u, v, weight)
 * From a number N and a three arrays it will create a new element of the class
 * From a tree/graph it will create a new element of the class
 * From an EdgeList it will create a new tree/graph
 */

package lapsolver;

import lapsolver.util.GraphUtils;

public class EdgeList {
    public int ne;
    public int[] u;
    public int[] v;
    public double[] weight;

    // an empty EdgeList
    public EdgeList() {}

    // EdgeList from a number of elements
    public EdgeList (int ne) {
        this.ne = ne;
        u = new int[ne];
        v = new int[ne];
        weight = new double[ne];
    }

    // EdgeList from a three lists of length N
    public EdgeList (int[] u, int[] v, double[] weight) {
        ne = u.length;
        this.u = u;
        this.v = v;
        this.weight = weight;
    }

    // return an EdgeList from a tree
    public EdgeList (Tree tree) {
        ne = tree.nv - 1;
        u = new int[ne];
        v = new int[ne];
        weight = new double[ne];

        int index = 0;
        for (int i = 0; i < tree.nv; i++) {
            // get all parents except root
            if (i != tree.getRoot()) {
                u[index] = i;
                v[index] = tree.getNode(i).getParent().getId();
                weight[index] = tree.getNode(i).getLength();
                index++;
            }
        }
    }

    // EdgeList from a graph
    public EdgeList (Graph graph) {
        ne = graph.ne;

        u = new int[ne];
        v = new int[ne];
        weight = new double[ne];

        int index = 0;
        int vertexCount = graph.nv;
        for (int i = 0; i < vertexCount; i++) {
            for (int j = 0; j < graph.deg[i]; j++) {
                int dest = graph.nbrs[i][j];
                double wt = graph.weights[i][j];

                // only count an edge once
                if (i < dest) {
                    u[index] = i;
                    v[index] = dest;
                    weight[index] = wt;
                    index++;
                }
            }
        }
    }

}
