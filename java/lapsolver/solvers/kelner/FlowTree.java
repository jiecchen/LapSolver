/**
 * @file FlowTree.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Thu Jun 5 2014
 *
 * An interface for the cycle update data structure, needed in Kelner's electrical flow solver.
 */

package lapsolver.solvers.kelner;

import lapsolver.Tree;

public interface FlowTree {
    // initialize the data structure
    public void init(Tree t);

    // push `alpha` units of flow along the tree path from u to v
    public void update(int u, int v, double alpha);

    // find the total flow along the tree path from u to v
    public double query(int u, int v);
}

