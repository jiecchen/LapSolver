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
        treeFlows = new double[tree.nv];
    }

    // push `alpha` units of flow along the tree path on edge e
    public void treeUpdate(int e, double alpha) {
        int u = super.offEdges.u[e];
        int v = super.offEdges.v[e];
        int root = super.tree.root;

        while (u != root) {
            treeFlows[u] += alpha;
            u = super.tree.parent[u];
        }

        while (v != root) {
            treeFlows[v] -= alpha;
            v = super.tree.parent[v];
        }
    }

    // find sum of V = IR along the tree path on edge e
    public double treeQuery(int e) {
        int u = super.offEdges.u[e];
        int v = super.offEdges.v[e];
        int root = super.tree.root;

        double total = 0;

        while (u != root) {
            double length = super.tree.weight[u];
            total += treeFlows[u] * length;
            u = super.tree.parent[u];
        }

        while (v != root) {
            double length = super.tree.weight[v];
            total -= treeFlows[v] * length;
            v = super.tree.parent[v];
        }

        return total;
    }

    // initialize the structure with some flows
    public void setTreeFlows(double[] treeFlows) {
        super.clearOffFlows();
        this.treeFlows = treeFlows.clone();
    }

    // retrieve the tree flows (off-tree flows are redundant)
    public double[] getTreeFlows() {
        return treeFlows;
    }
}
