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
import lapsolver.algorithms.GraphVertexRemoval;
import lapsolver.algorithms.LDLDecomposition;
import lapsolver.algorithms.Stretch;
import lapsolver.lsst.SpanningTreeStrategy;
import lapsolver.util.GraphUtils;
import lapsolver.util.TreeUtils;

import matlabcontrol.*;
import matlabcontrol.extensions.*;

import java.util.ArrayList;
import java.util.Arrays;

public class KMPSolver {
    private SpanningTreeStrategy treeStrategy;

    public double tol;
    public int maxit;

    //Initialize solver with a spanning tree strategy
    public KMPSolver(SpanningTreeStrategy treeStrategy) {
        this.treeStrategy = treeStrategy;
    }

    public EdgeList LL;
    public EdgeList DD;

    //Initialize solver on a particular graph, and perform preprocessing
    public double[] solve(Graph graph, double b[], double tol, int maxit, double[] addDiag) {
        this.tol = tol;
        this.maxit = maxit;

        System.out.println(graph.nv);

        if (graph.nv < 500)
            return runMain(graph, b, addDiag);

        //Build the preconditioner for the given graph
        Graph sparsifier = buildPreconditioner(graph, treeStrategy.getTree(graph));

        //Create the graph that will be used in the recursion (laplacian given by D matrix)
        GraphVertexRemoval gvrElement = new GraphVertexRemoval(sparsifier);
        GraphVertexRemoval.AnswerPair gvr = gvrElement.solve();

        b = applyPerm(gvr.permutation, b);
        addDiag = applyPerm(gvr.permutation, addDiag);
        Graph permSparsifier = GraphUtils.permuteGraph(sparsifier, gvr.permutation);

        LDLDecomposition ldlElement = new LDLDecomposition(permSparsifier, addDiag);
        LDLDecomposition.ReturnPair ldl = ldlElement.solve(gvr.numRemoved);

        LL = new EdgeList(ldl.L);
        DD = new EdgeList(ldl.D);

        b = LDLDecomposition.applyLInv(ldl.L, b);

        // grab diagonal elements from D matrix
        double[] diagD = new double[graph.nv];
        for (int i = 0; i < ldl.D.ne; i++) {
            if (ldl.D.u[i] == ldl.D.v[i]) {
                diagD[ldl.D.u[i]] += ldl.D.weight[i];
            }
        }

        Graph reducedSparsifier = buildRecursionGraph(graph, gvr, ldl);

        // for eliminated vertices, x[i] = b[i]/D[i,i]
        double[] x = new double[graph.nv];
        for (int i = 0; i < gvr.numRemoved; i++) {
            x[i] = b[i] / diagD[i];
        }

        double[] smallb = new double[graph.nv - gvr.numRemoved];
        System.arraycopy(b, gvr.numRemoved, smallb, 0, smallb.length);

        double[] smallAddDiag = new double[graph.nv - gvr.numRemoved];
        System.arraycopy(addDiag, gvr.numRemoved, smallAddDiag, 0, smallAddDiag.length);

        double[] KMPx = solve(reducedSparsifier, smallb, tol, maxit, smallAddDiag);
        System.arraycopy(KMPx, 0, x, gvr.numRemoved, KMPx.length);

        x = LDLDecomposition.applyLTransInv(ldl.L, x);

        int[] inversePerm = new int[graph.nv];
        for (int i = 0; i < graph.nv; i++)
            inversePerm[gvr.permutation[i]] = i;

        return applyPerm(inversePerm, x);
    }

    public static double[] applyPerm (int[] perm, double[] x) {
        double[] answer = new double[perm.length];

        for (int i = 0; i < x.length; i++)
            answer[i] = x[perm[i]];

        return answer;
    }

    //Construct a preconditioner for graph
    public Graph buildPreconditioner(Graph graph, Tree spanningTree) {
        EdgeList offEdges;

        //Get off-tree edges, find stretches
        offEdges = TreeUtils.getOffTreeEdges(graph, spanningTree);
        Stretch.StretchResult stretch = Stretch.compute(graph, spanningTree, offEdges);

        //Blow up graph by 4 * avgstretch * log(numRemoved)
        double k = 4. * stretch.total / (offEdges.ne + 1) * (Math.log(graph.nv) + 1) + 1;
        Graph blownUpGraph = blowUpTreeEdges(graph, spanningTree, k);

        // find stretches in blown-up graph
        Stretch.StretchResult blownUpStretch = Stretch.compute(blownUpGraph, spanningTree, offEdges);

        // Expect to grab q = O(m / log(m)) edges
        double q = 10. * graph.ne / Math.log(graph.ne);

        //Assign p_e = stretch(e) / (total stretch)
        double[] p = blownUpStretch.allStretches.clone();
        for (int i = 0; i < offEdges.ne; i++) {
            p[i] = q * p[i] / blownUpStretch.total;
            if (p[i] > 1) p[i] = 1;
        }

        //Sample the edges
        ArrayList<Integer> edgesToAdd = new ArrayList<Integer>();
        for (int i = 0; i < offEdges.ne; i++) {
            if (Math.random() < p[i]) {
                // with probability p[i], take edge i
                edgesToAdd.add(i);
            }
        }

        //Generate EdgeList
        EdgeList sparsifierEdges = new EdgeList(graph.nv - 1 + edgesToAdd.size());
        int index = 0;

        //Add tree edges
        for (int u = 0; u < graph.nv; u++) {
            if (u == spanningTree.root) continue;
            sparsifierEdges.u[index] = u;
            sparsifierEdges.v[index] = spanningTree.parent[u];
            sparsifierEdges.weight[index] = spanningTree.weight[u];
            index++;
        }

        //Add off-tree edges
        for (int i : edgesToAdd) {
            sparsifierEdges.u[index] = offEdges.u[i];
            sparsifierEdges.v[index] = offEdges.v[i];
            sparsifierEdges.weight[index] = offEdges.weight[i] / p[i];
            index++;
        }

        //Build sparsified graph
        return new Graph(sparsifierEdges);
    }

    // given graph and tree, blow up edge weights (not lengths!!) of tree edges in G
    public static Graph blowUpTreeEdges(Graph graph, Tree spanningTree, double k) {
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

    //Construct the graph for the next step of the recursion
    public Graph buildRecursionGraph(Graph graph, GraphVertexRemoval.AnswerPair gvr, LDLDecomposition.ReturnPair ldl) {
        ldl.D = GraphUtils.sanitizeEdgeList(ldl.D);
        EdgeList reducedSparsifierEdges = new EdgeList(ldl.D.ne);
        int index = 0;
        for (int i = 0; i < ldl.D.ne; i++) {
            if (ldl.D.u[i] >= gvr.numRemoved && ldl.D.v[i] >= gvr.numRemoved) {
                reducedSparsifierEdges.u[index] = ldl.D.u[i] - gvr.numRemoved;
                reducedSparsifierEdges.v[index] = ldl.D.v[i] - gvr.numRemoved;
                reducedSparsifierEdges.weight[index] = ldl.D.weight[i];
                index++;
            }
        }

        return new Graph(reducedSparsifierEdges);
    }

    //Construct the graph laplacian and call the pcg solver
    public double[] runMain(Graph graph, double[] b, double[] diag) {
        double[][] lap = new double[graph.nv][graph.nv];
        for (int i = 0; i < graph.nv; i++) {
            for (int j = 0; j < graph.deg[i]; j++) {
                lap[i][graph.nbrs[i][j]] -= graph.weights[i][j];
                lap[i][i] += graph.weights[i][j] / 2;
                lap[graph.nbrs[i][j]][graph.nbrs[i][j]] += graph.weights[i][j] / 2;
            }
        }

        for (int i = 0; i < diag.length; i++)
            lap[i][i] += diag[i];

        try {
            return main(lap, b, tol, maxit);
        } catch (MatlabConnectionException e) {
            e.printStackTrace();
        } catch (MatlabInvocationException e) {
            e.printStackTrace();
        }

        return null;
    }

    //Call pcg in matlab from java
    public static double[] main(double[][] lap, double[] b, double tol, int maxit) throws MatlabConnectionException, MatlabInvocationException
    {
        double[] x = new double[b.length];

        //Create a proxy, which we will use to control MATLAB
        MatlabProxyFactory factory = new MatlabProxyFactory();
        MatlabProxy proxy = factory.getProxy();

        //Send the laplacian to MATLAB
        MatlabTypeConverter processor = new MatlabTypeConverter(proxy);
        processor.setNumericArray("internal_Lap", new MatlabNumericArray(lap, null));
        proxy.setVariable("internal_b", b);

//        proxy.setVariable("internal_x", x);
//        proxy.setVariable("internal_tol", tol);
//        proxy.setVariable("internal_maxit", maxit);
//        proxy.eval("internal_x = pcg(internal_Lap, internal_b', internal_tol, internal_maxit);");

        proxy.eval("internal_x = pinv(internal_Lap) * internal_b';");

        x = (double[]) proxy.getVariable("internal_x");

        //Disconnect the proxy from MATLAB
        proxy.disconnect();

        return x;
    }
}
