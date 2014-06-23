/**
 * @file KMPSolver.java
 * @author Serban Stan <serban.stan@yale.edu>
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Fri Jun 20 2014
 *
 * The recursive preconditioning solver from KMP1.
 */

package lapsolver.solvers;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.DiscreteSampler;
import lapsolver.algorithms.GraphVertexRemoval;
import lapsolver.algorithms.LDLDecomposition;
import lapsolver.algorithms.Stretch;
import lapsolver.lsst.SpanningTreeStrategy;
import lapsolver.util.GraphUtils;
import lapsolver.util.TreeUtils;

import matlabcontrol.MatlabConnectionException;
import matlabcontrol.MatlabInvocationException;
import matlabcontrol.MatlabProxy;
import matlabcontrol.MatlabProxyFactory;

import java.util.ArrayList;

public class KMPSolver {
    public Tree spanningTree;
    public Graph reweightedGraph;
    public Graph sparsifier;
    public Graph reducedSparsifier;
    private SpanningTreeStrategy treeStrategy;

    // initialize solver with a spanning tree strategy
    public KMPSolver(SpanningTreeStrategy treeStrategy) {
        this.treeStrategy = treeStrategy;
    }

    public static void main(String[] args) throws MatlabConnectionException, MatlabInvocationException
    {
        //Create a proxy, which we will use to control MATLAB
        MatlabProxyFactory factory = new MatlabProxyFactory();
        MatlabProxy proxy = factory.getProxy();

        //Display 'hello world' just like when using the demo
        proxy.eval("disp('hello world')");

        //Disconnect the proxy from MATLAB
        proxy.disconnect();
    }

    // initialize solver on a particular graph, and perform preprocessing
    public double[] solve(Graph graph, double b[], double err) {
        if (graph.nv < 500)
            return thirdPartySolver(graph, b, err);

        // compute LSST, cache BFS order
        spanningTree = treeStrategy.getTree(graph);

        // build the preconditioner for the given graph
        buildPreconditioner(graph);


        // create the graph that will be used in the recursion (laplacian given by D matrix)
        GraphVertexRemoval gvrElement = new GraphVertexRemoval(sparsifier);
        GraphVertexRemoval.AnswerPair gvr = gvrElement.solve();

        Graph permSparsifier = GraphUtils.permuteGraph(sparsifier, gvr.v);
        LDLDecomposition ldlElement = new LDLDecomposition(permSparsifier, new double[permSparsifier.nv]);
        LDLDecomposition.ReturnPair ldl = ldlElement.solve(gvr.n);


        int[] inversePerm = new int[gvr.n];
        for (int i = 0; i < graph.nv; i++) {
            inversePerm[gvr.v[i]] = i;
        }

        buildRecursionGraph(graph, gvr, ldl);

        return reconstructSolution(b, err, graph.nv, gvr.n, inversePerm, ldl);
    }

    public double[] thirdPartySolver(Graph graph, double[] b, double err) {
        double[] x = new double[graph.nv];

        return x;
    }

    public double[] reconstructSolution(double[] b, double err, int currentN, int smallN, int[] invPerm, LDLDecomposition.ReturnPair ldl) {
        b = LDLDecomposition.applyInvL(ldl.L, currentN, b);

        double[] x = new double[currentN];
        System.arraycopy(b, 0, x, 0, smallN);

        x = LDLDecomposition.applyLtransInv(ldl.L, currentN, x);

        double[] smallb = new double[currentN - smallN];
        System.arraycopy(b, currentN - smallN, smallb, 0, smallb.length);

        double[] KMPx = solve(reducedSparsifier, smallb, err);
        System.arraycopy(KMPx, 0, x, currentN - smallN, KMPx.length);

        double[] answer = new double[currentN];
        for (int i = 0; i < x.length; i++)
            answer[invPerm[i]] = x[i];

        return answer;
    }

    public void buildPreconditioner(Graph graph) {
        EdgeList offEdges;

        // get off-tree edges, find stretches
        offEdges = TreeUtils.getOffTreeEdges(graph, spanningTree);
        Stretch.StretchResult stretch = Stretch.compute(graph, spanningTree, offEdges);

        // blow up graph by 4 * avgstretch * log(n)
        double k = 4. * stretch.total / offEdges.ne * Math.log(graph.nv);
        reweightedGraph = blowUpTreeEdges(graph, spanningTree, k);

        // expect to grab q = O(m / log(m)) edges
        double q = 4. * graph.ne / Math.log(graph.ne);

        // assign p_e = stretch(e) / (total stretch)
        double[] p = stretch.allStretches.clone();
        for (int i = 0; i < offEdges.ne; i++) {
            p[i] = q * p[i] / stretch.total;
        }

        // sample the edges
        ArrayList<Integer> edgesToAdd = new ArrayList<Integer>();
        for (int i = 0; i < offEdges.ne; i++) {
            if (Math.random() < p[i]) {
                // with probability p[i], take edge i
                edgesToAdd.add(i);
            }
        }

        // generate EdgeList
        EdgeList sparsifierEdges = new EdgeList(graph.nv - 1 + edgesToAdd.size());
        int index = 0;

        // add tree edges
        for (int u = 0; u < graph.nv; u++) {
            if (u == spanningTree.root) continue;
            sparsifierEdges.u[index] = u;
            sparsifierEdges.v[index] = spanningTree.parent[u];
            sparsifierEdges.weight[index] = spanningTree.weight[u];
            index++;
        }

        // add off-tree edges
        for (int i : edgesToAdd) {
            sparsifierEdges.u[index] = offEdges.u[i];
            sparsifierEdges.v[index] = offEdges.v[i];
            sparsifierEdges.weight[index] = offEdges.weight[i];
            index++;
        }

        // build sparsified graph
        sparsifier = new Graph(sparsifierEdges);
    }

    // given graph and tree, blow up edge weights (not lengths!!) of tree edges in G
    public Graph blowUpTreeEdges(Graph graph, Tree spanningTree, double k) {
        Graph auxGraph = new Graph(graph);

        for (int u = 0; u < auxGraph.nv; u++) {
            for (int i = 0; i < auxGraph.deg[u]; i++) {
                int v = auxGraph.nbrs[u][i];

                if (spanningTree.parent[u] == v || spanningTree.parent[v] == u) {
                    auxGraph.weights[u][i] /= k;
                }
            }
        }

        return auxGraph;
    }

    public void buildRecursionGraph(Graph graph, GraphVertexRemoval.AnswerPair gvr, LDLDecomposition.ReturnPair ldl) {
        ldl.D = GraphUtils.sanitizeEdgeList(ldl.D);
        EdgeList reducedSparsifierEdges = new EdgeList(ldl.D.ne);
        int index = 0;
        for (int i = 0; i < ldl.D.ne; i++) {
            if (ldl.D.u[i] >= gvr.n && ldl.D.v[i] >= gvr.n) {
                reducedSparsifierEdges.u[index] = ldl.D.u[i] - gvr.n;
                reducedSparsifierEdges.v[index] = ldl.D.v[i] - gvr.n;
                reducedSparsifierEdges.weight[index] = 1 / ldl.D.weight[i];
                index++;
            }
        }

        reducedSparsifier = new Graph(reducedSparsifierEdges);
    }
}
