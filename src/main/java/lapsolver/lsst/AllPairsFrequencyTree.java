package lapsolver.lsst;

import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.ShortestPathTree;
import lapsolver.util.GraphUtils;

import java.util.Arrays;

public class AllPairsFrequencyTree implements SpanningTreeStrategy {
    @Override public Tree getTree(Graph graph) {
        GraphUtils.reciprocateWeights(graph);

        Graph freqGraph = new Graph(graph);
        for (int i = 0; i < freqGraph.nv; i++)
            for (int j = 0; j < freqGraph.weights[i].length; j++)
                Arrays.fill(freqGraph.weights[i], 0);

        for (int i = 0; i < graph.nv; i++) {
            Tree spt = new ShortestPathTree(graph, i).getTree();
            for (int u = 0; u < spt.parent.length; u++) {
                int v = spt.parent[u];
                if (u != v) // can probably do better than adding 1 (something proportional to total stretch?)
                    freqGraph.setWeight(u, v, freqGraph.getWeight(u, v) + 1);
            }
        }

        GraphUtils.reciprocateWeights(graph);
        GraphUtils.reciprocateWeights(freqGraph);

        Tree finalTree = new KruskalTree().getTree(freqGraph);
        for (int u = 0; u < finalTree.nv; u++) {
            int v = finalTree.parent[u];
            if (u != v)
                finalTree.weight[u] = graph.getWeight(u, v);
        }

        return finalTree;
    }
}
