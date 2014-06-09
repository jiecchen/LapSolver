/**
 * @file KelnerFlowTree.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Thu Jun 5 2014
 *
 * An implementation of the cycle update data structure in the paper of Kelner et al.
 */

package lapsolver.solvers.kelner;

import lapsolver.EdgeList;
import lapsolver.Tree;

public class DirectFlowTree extends FlowTree {
    private double[] treeFlows;

    // base constructor
    public DirectFlowTree (Tree tree, EdgeList offEdges) {
        super(tree, offEdges);
    }

    // push `alpha` units of flow along the tree path on edge e
    public void treeUpdate(int e, double alpha) {
        int u = super.offEdges.u[e];
        int v = super.offEdges.v[e];
        int root = super.tree.getRoot();

        do {
            treeFlows[u] += alpha;
            u = super.tree.getNode(u).getParent().getId();
        } while(u != root);

        do {
            treeFlows[v] -= alpha;
            v = super.tree.getNode(v).getParent().getId();
        } while(v != root);
    }

    // find sum of V = IR along the tree path on edge e
    public double treeQuery(int e) {
        int u = super.offEdges.u[e];
        int v = super.offEdges.v[e];
        int root = super.tree.getRoot();

        double total = 0;

        do {
            total += treeFlows[u];
            u = super.tree.getNode(u).getParent().getId();
        } while(u != root);

        do {
            total -= treeFlows[v];
            v = super.tree.getNode(v).getParent().getId();
        } while(v != root);

        return total;
    }

    // initialize the structure with some flows
    public void setTreeFlows(double[] treeFlows) {
        this.treeFlows = treeFlows.clone();
    }

    // retrieve the tree flows (off-tree flows are redundant)
    public double[] getTreeFlows() {
        return treeFlows;
    }
}
