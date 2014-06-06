/**
 * @file PetalDecompositionLSST.java
 * @author Alex Reinking <alexander.reinking@yale.edu>
 * @date Tue Jun 3 2014
 *
 * An LSST solver based on the Petal Decomposition method.
 * Still a stub.
 *
 */
package lapsolver.lsst;

import lapsolver.Tree;
import lapsolver.Graph;

public class PetalDecompositionTree implements SpanningTreeStrategy {
    public final Tree tree;
    public final Graph graph;

    public PetalDecompositionTree(Graph graph) {
        this.graph = graph;
        this.tree = null;
    }

    @Override
    public Tree getTree() {
        return tree;
    }
}
