/**
 * Created by serbanstan on 8/7/14.
 *
 * Finds an approximation to the diameter of the tree, then creates a cut in the tree by growing components (lssts)
 * on each edges of the diameter. Repeats the process recursively for each side of the cut.
 *
 */

package lapsolver.lsst;

import lapsolver.EdgeList;
import lapsolver.Graph;
import lapsolver.Tree;
import lapsolver.algorithms.ShortestPathTree;

import java.util.Arrays;

public class DiameterSplitTree implements SpanningTreeStrategy{
    @Override
    public Tree getTree(Graph graph) {
        return new Tree(getTreeEdges(graph));
    }

    public EdgeList getTreeEdges(Graph graph) {
        if (graph.nv < 5) {
            SpanningTreeStrategy lsst = new StarDecompositionTree();
            EdgeList treeEdges = new EdgeList(lsst.getTree(graph));
            return treeEdges;
        }

        // find diameter in graph
        int v = getFurthestVertex(0, graph);
        int u = getFurthestVertex(v, graph);
        v = getFurthestVertex(u, graph);

        int start = u, stop = v;

        // k should be smaller than 0.5
        double k = 0.4;
        u = getRatioVertex(k, start, stop, graph);
        v = getRatioVertex(k, stop, start, graph);

        // create the split
        ShortestPathTree sptU = new ShortestPathTree(graph, u);
        ShortestPathTree sptV = new ShortestPathTree(graph, v);

        double[] distU = sptU.getDist();
        double[] distV = sptV.getDist();

        int[] side = new int[graph.nv];
        for (int i = 0; i < graph.nv; i++)
            if (distU[i] < distV[i])
                side[i] = 0;
            else
                side[i] = 1;

        int sizeU = 0, sizeV = 0;
        for (int i = 0; i < graph.nv; i++)
            for (int j = 0; j < graph.nbrs[i].length; j++) {
                int p = i;
                int q = graph.nbrs[i][j];

                if (p < q) {
                    if (side[p] == 0 && side[q] == 0)
                        sizeU++;
                    if (side[p] == 1 && side[q] == 1)
                        sizeV++;
                }
            }

        EdgeList sideU = new EdgeList(sizeU);
        EdgeList sideV = new EdgeList(sizeV);

        int[] verticesU = new int[graph.nv];
        int[] verticesV = new int[graph.nv];

        sizeU = 0; sizeV = 0;
        for (int i = 0; i < graph.nv; i++)
            for (int j = 0; j < graph.nbrs[i].length; j++) {
                int p = i;
                int q = graph.nbrs[i][j];

                if (p < q) {
                    if (side[p] == 0 && side[q] == 0) {
                        verticesU[p] = 1;
                        verticesU[q] = 1;

                        sideU.u[sizeU] = p;
                        sideU.v[sizeU] = q;
                        sideU.weight[sizeU++] = graph.weights[i][j];
                    }
                    if (side[p] == 1 && side[q] == 1) {
                        verticesV[p] = 1;
                        verticesV[q] = 1;

                        sideV.u[sizeV] = p;
                        sideV.v[sizeV] = q;
                        sideV.weight[sizeV++] = graph.weights[i][j];
                    }
                }
            }

        int[] indexU = new int[graph.nv];
        int[] indexV = new int[graph.nv];

        int[] backIndexU = new int[graph.nv];
        int[] backIndexV = new int[graph.nv];

        int xU = 0, xV = 0;
        for (int i = 0; i < graph.nv; i++) {
            if (verticesU[i] == 1) {
                backIndexU[xU] = i;
                indexU[i] = xU++;
            }
            if (verticesV[i] == 1) {
                backIndexV[xV] = i;
                indexV[i] = xV++;
            }
        }

        //System.out.println(Arrays.toString(indexU));
        //System.out.println(Arrays.toString(backIndexU));

        //System.out.println(Arrays.toString(indexV));
        //System.out.println(Arrays.toString(backIndexV));

        for (int i = 0; i < sizeU; i++) {
            sideU.u[i] = indexU[sideU.u[i]];
            sideU.v[i] = indexU[sideU.v[i]];
        }

        for (int i = 0; i < sizeV; i++) {
            sideV.u[i] = indexV[sideV.u[i]];
            sideV.v[i] = indexV[sideV.v[i]];
        }

        // do the same algorithm recursively on the split graph
        SpanningTreeStrategy dstU = new DiameterSplitTree();
        SpanningTreeStrategy dstV = new DiameterSplitTree();

        Tree lsstU = dstU.getTree(new Graph(sideU));
        Tree lsstV = dstV.getTree(new Graph(sideV));

        // add the edge of minimum resistance and return the unified graph
        EdgeList treeEdges = new EdgeList(graph.nv);
        int cnt = 0;

        for (int i = 0; i < lsstU.nv; i++) {
            if (i == lsstU.root) continue;

            treeEdges.u[cnt] = backIndexU[i];
            treeEdges.v[cnt] = backIndexU[lsstU.parent[i]];
            treeEdges.weight[cnt++] = lsstU.weight[i];
        }

        for (int i = 0; i < lsstV.nv; i++) {
            if (i == lsstV.root) continue;

            treeEdges.u[cnt] = backIndexV[i];
            treeEdges.v[cnt] = backIndexV[lsstV.parent[i]];
            treeEdges.weight[cnt++] = lsstV.weight[i];
        }

        //System.out.println("***********");

        int p = 0, q = 0;
        double cost = -1;
        for (int i = 0; i < graph.nv; i++)
            for (int j = 0; j < graph.nbrs[i].length; j++) {
                if (i < graph.nbrs[i][j] && side[i] != side[graph.nbrs[i][j]]) {
                    //System.out.println("! " + graph.weights[i][j] + " " + i + " " + graph.nbrs[i][j]);

                    if (cost == -1 || graph.weights[i][j] < cost) {
                        cost = graph.weights[i][j];
                        p = i;
                        q = graph.nbrs[i][j];
                    }
                }
            }

        //System.out.println("**************");

        //System.out.println(p + " " + q);

        treeEdges.u[cnt] = p;
        treeEdges.v[cnt] = q;
        treeEdges.weight[cnt++] = cost;

        //System.out.println(cnt + " " + graph.nv);

        return treeEdges;
    }

    public int getFurthestVertex(int start, Graph graph) {
        ShortestPathTree spt = new ShortestPathTree(graph, start);
        double[] dist = spt.getDist();

        double maxDist = 0;
        int index = start;
        for (int i = 0; i < graph.nv; i++)
            if (maxDist < dist[i]) {
                maxDist = dist[i];
                index = i;
            }

        return index;
    }

    public int getRatioVertex(double k, int start, int stop, Graph graph) {
        //System.out.println(k);

        ShortestPathTree spt = new ShortestPathTree(graph, start);
        double[] distStart = spt.getDist();

        spt = new ShortestPathTree(graph, stop);
        double[] distStop = spt.getDist();

        double distTotal = distStart[stop];

        double vmax = 0;
        int ret = start;
        for (int i = 0; i < graph.nv; i++)
            if (distStart[i] + distStop[i] == distTotal) {
                if (distStart[i] < k * distTotal && distStart[i] > vmax) {
                    vmax = distStart[i];
                    ret = i;
                }
            }

        return ret;
    }
}
