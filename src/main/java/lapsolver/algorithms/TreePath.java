/**
 * @file TreePath.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Wed Jun 4 2014
 *
 * A data structure for computing tree path lengths. Uses Tarjan's offline LCA.
 *
 * Usage:
 * - query(a, b) takes endpoint arrays, and returns array of path lengths
 * - public field depth[u] gives length of path to root
 *
 * TreePath tp = new TreePath(tree);
 * lens = tp.query(a, b);
 *
 */

package lapsolver.algorithms;

import lapsolver.Tree;
import lapsolver.util.TreeUtils;

public class TreePath {
    private final Tree tree;
    public double depth[];

    /**
     * Initialize the data structure with tree depths.
     * @param tree The input tree (weights are resistances).
     */
    public TreePath (Tree tree) {
        this.tree = tree;
        depth = new double[tree.nv];

        // compute depths
        depth = TreeUtils.getDepths(tree);
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
}
