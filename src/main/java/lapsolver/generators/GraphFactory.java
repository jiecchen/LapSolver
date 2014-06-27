package lapsolver.generators;

import lapsolver.Graph;

public interface GraphFactory {
    /**
     * @return The output graph (weights are conductances).
     */
    public Graph generateGraph();
}
