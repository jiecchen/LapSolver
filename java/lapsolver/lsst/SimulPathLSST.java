package lapsolver.lsst;

import lapsolver.SimulPathTree;
import lapsolver.Tree;
import lapsolver.WeightedGraph;

/**
 * Created by alex on 6/3/14.
 */
public class SimulPathLSST implements LSST {
    @Override
    public Tree solve(WeightedGraph in) {
        return new SimulPathTree(in).growTree().treeToTree();
    }
}
