/**
 * @file SpanningTreeStrategy.java
 * @author Alex Reinking <alexander.reinking@yale.edu>
 * @date Tue Jun 3 2014
 *
 * An interface for various SpanningTreeStrategy solvers to inherit from.
 * This is basically an implementation of the strategy pattern.
 *
 */
package lapsolver.lsst;

import lapsolver.Tree;
import lapsolver.WeightedGraph;

public interface SpanningTreeStrategy {
    public Tree solve(WeightedGraph in);
}
