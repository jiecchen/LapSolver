/**
 * @file FlowTree.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Thu Jun 5 2014
 *
 * An interface for the cycle update data structure, needed in Kelner's electrical flow solver.
 */

package lapsolver.solvers.kelner;

import lapsolver.EdgeList;
import lapsolver.Tree;

public abstract class FlowTree {
    // common implementation for off-tree edges
    public Tree tree;
    public EdgeList offEdges;
    public double[] offFlow;

    // base constructor
    public FlowTree(Tree tree, EdgeList offEdges) {
        this.tree = tree;
        this.offEdges = offEdges;
        offFlow = new double[offEdges.ne];
    }

    // push `alpha` units of flow along the tree cycle on edge e
    public void update(int e, double alpha) {
        treeUpdate(e, alpha);
        offFlow[e] += alpha;
    }

    // find sum of V = IR along the tree cycle on edge e
    public double query(int e) {
        return treeQuery(e) + offFlow[e] * offEdges.weight[e];
    }

    // push `alpha` units of flow along the tree path on edge e
    public abstract void treeUpdate(int e, double alpha);

    // find sum of V = IR along the tree path on edge e
    public abstract double treeQuery(int e);

    // initialize the structure with some flows
    public abstract void setTreeFlows(double[] treeFlows);

    // retrieve the tree flows (off-tree flows are redundant)
    public abstract double[] getTreeFlows();
}

