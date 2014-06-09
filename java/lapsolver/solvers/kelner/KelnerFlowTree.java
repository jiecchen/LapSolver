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

public class KelnerFlowTree extends FlowTree {
    // base constructor
    public KelnerFlowTree (Tree tree, EdgeList offEdges) {
        super(tree, offEdges);
    }

    // push `alpha` units of flow along the tree path on edge e
    public void treeUpdate(int e, double alpha) {

    }

    // find sum of V = IR along the tree path on edge e
    public double treeQuery(int e) {
        return 0.0;
    }

    // initialize the structure with some flows
    public void setTreeFlows(double[] treeFlows) {

    }

    // retrieve the tree flows (off-tree flows are redundant)
    public double[] getTreeFlows() {
        return null;
    }
}
