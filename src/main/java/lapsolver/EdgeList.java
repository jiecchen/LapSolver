/**
 * @file EdgeList.java
 * @author Serban Stan <serban.stan@gmail.com>
 * @date Thu Jun 5 2014
 *
 * An EdgeList of the form (u, permutation, weight)
 * From a number N and a three arrays it will create a new element of the class
 * From a tree/graph it will create a new element of the class
 * From an EdgeList it will create a new tree/graph
 */

package lapsolver;

import java.util.Arrays;

public class EdgeList {
    public int ne;
    public int[] u;
    public int[] v;
    public double[] weight;

    // an empty EdgeList
    public EdgeList() {
        this(0);
    }

    // EdgeList from a number of elements
    public EdgeList(int ne) {
        this.ne = ne;
        u = new int[ne];
        v = new int[ne];
        weight = new double[ne];
    }

    // copy constructor
    public EdgeList(EdgeList other) {
        ne = other.ne;
        u = other.u.clone();
        v = other.v.clone();
        weight = other.weight.clone();
    }

    // EdgeList from two lists of length N (consider the weights = 1)
    public EdgeList(int[] u, int[] v) {
        ne = u.length;
        this.u = u;
        this.v = v;
        this.weight = new double[ne];
        for (int i = 0; i < ne; i++)
            this.weight[i] = 1;
    }

    // EdgeList from three lists of length N
    public EdgeList(int[] u, int[] v, double[] weight) {
        ne = u.length;
        this.u = u;
        this.v = v;
        this.weight = weight;
    }

    // return an EdgeList from a tree
    public EdgeList(Tree tree) {
        ne = tree.nv - 1;
        u = new int[ne];
        v = new int[ne];
        weight = new double[ne];

        int index = 0;
        for (int i = 0; i < tree.nv; i++) {
            // get all parents except root
            if (i != tree.root) {
                u[index] = i;
                v[index] = tree.parent[i];
                weight[index] = tree.weight[i];
                index++;
            }
        }
    }

    // EdgeList from a graph
    public EdgeList(Graph graph) {
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

    // change number of edges, without affecting old elements
    public void resize(int ne) {
        this.ne = ne;
        u = Arrays.copyOf(u, ne);
        v = Arrays.copyOf(v, ne);
        weight = Arrays.copyOf(weight, ne);
    }

}
