/**
 * @file OrderedSpanningTreeStrategy.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Thu Aug 7 2014
 *
 * A tree strategy that provides an ordered EdgeList, for MATLAB animation purposes.
 */

package lapsolver.lsst;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;

public abstract class OrderedSpanningTreeStrategy implements SpanningTreeStrategy {
    public abstract EdgeList getTreeEdges (Graph graph);

    @Override
    public Tree getTree (Graph graph) {
        return new Tree(getTreeEdges(graph));
    }
}
