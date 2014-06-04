/**
 * @file TreePath.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Wed Jun 4 2014
 *
 * A data structure for computing tree path lengths. Uses Tarjan's offline LCA.
 */

package lapsolver.algorithms;

import lapsolver.Tree;
import lapsolver.algorithms.TarjanLCA;

public class TreePath {
    private Tree tree;
    private double depth[];

    // constructor: do precomputations
    public TreePath (Tree tree) {
        this.tree = tree;
        depth = new double[tree.nv];

        // compute depths
        depth[tree.root] = 0;
        dfs(tree.root);
    }

    // query path lengths
    public double[] query(int[] a, int[] b) {
        // compute LCAs
        TarjanLCA lca = new TarjanLCA(tree);
        int[] c = lca.solve(a, b);
        double[] answer = new double[c.length];

        // len(a,b) = d(a) + d(b) - 2*d(lca(a,b))
        for (int i = 0; i < c.length; i++) {
            answer[i] = depth[a[i]] + depth[b[i]] - 2*depth[c[i]];
        }

        return answer;
    }

    // populate depth array
    public void dfs(int u) {
        if (u != tree.root) {
            depth[u] = depth[tree.nodes[u].parent] + tree.nodes[u].length;
        }

        for (int v : tree.nodes[u].children) {
            dfs(v);
        }
    }
}
