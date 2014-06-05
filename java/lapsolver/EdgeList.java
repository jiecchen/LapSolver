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

public class EdgeList {
    public int edgeCount;
    public int[] firstEntry;
    public int[] secondEntry;
    public double[] weight;

    // A null EdgeList
    public EdgeList() {}

    // EdgeList from a number of elements
    public EdgeList (int N) {
        this.edgeCount = N;
        this.firstEntry = new int[N];
        this.secondEntry = new int[N];
        this.weight = new double[N];
    }

    // EdgeList from a three lists of length N
    public EdgeList (int[] U, int V[], double[] W) {
        this.edgeCount = U.length;
        this.firstEntry = U;
        this.secondEntry = V;
        this.weight = W;
    }

    // return an EdgeList from a tree
    public EdgeList (Tree tree) {
        this.edgeCount = tree.nv;
        this.firstEntry = new int[edgeCount - 1];
        this.secondEntry = new int[edgeCount - 1];
        this.weight = new double[edgeCount - 1];

        int index = 0;
        for (int i = 0; i < edgeCount; i++)
            if (i != tree.getRoot()) {
                this.firstEntry[index] = i;
                this.secondEntry[index] = tree.getNode(i).getParent().getId();
                this.weight[index] = tree.getNode(i).getLength();
                index++;
            }
    }

    // EdgeList from a graph
    public EdgeList (Graph graph) {
        this.edgeCount = graph.ne;
        this.firstEntry = new int[edgeCount];
        this.secondEntry = new int[edgeCount];
        this.weight = new double[edgeCount];
        
        int index = 0;
        int vertexCount = graph.nv;
        for (int i = 0; i < vertexCount; i++)
            for (int j = 0; j < graph.deg[i]; j++) {
                int v = graph.nbrs[i][j];
                double w = graph.weights[i][j];

                if (i < v) {
                    this.firstEntry[index] = i;
                    this.secondEntry[index] = v;
                    this.weight[index] = w;
                    index++;
                }
            }
    }

    public Graph toGraph() {
        return new Graph(firstEntry, secondEntry, weight);
    }
    
    public Tree toTree() {
        return toTree(0);
    }

    public Tree toTree(int root) {
        int nV = 0;
        for (int i : firstEntry) if (i > nV) nV = i;
        for (int i : secondEntry) if (i > nV) nV = i;

        if (nV != edgeCount)
            throw new RuntimeException("Error: edge list does not represent a tree.");

        int [] parentList = new int[nV];
        for (int i = 0; i < edgeCount; i++) {

        }
        return new Tree(parentList);
    }
        
}
