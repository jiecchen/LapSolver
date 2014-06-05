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
    public int elementCount;
    public int[] firstEntry;
    public int[] secondEntry;
    public double[] weight;


    // A null EdgeList
    public EdgeList() {}

    // EdgeList from a number of elements
    public EdgeList (int N) {
        this.elementCount = N;
        this.firstEntry = new int[N];
        this.secondEntry = new int[N];
        this.weight = new double[N];
    }

    // EdgeList from a number of elements and a list of triplets
    public EdgeList (int[] U, int V[], double[] W) {
        this.elementCount = U.length;
        this.firstEntry = U;
        this.secondEntry = V;
        this.weight = W;
    }

    // return an EdgeList from a tree
    public EdgeList (Tree tree) {
        this.elementCount = tree.nv;

        int index = 0;
        for (int i = 0; i < elementCount; i++)
            if (i != tree.getRoot()) {
                this.firstEntry[index] = i;
                this.secondEntry[index] = tree.getNode(i).getParent().getId();
                this.weight[index] = tree.getNode(i).getLength();
                index++;
            }
    }


}
