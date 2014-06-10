/**
 * @file KelnerFlowTree.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Thu Jun 5 2014
 *
 * An implementation of the cycle update data structure in the paper of Kelner et al.
 */

package lapsolver.solvers.kelner;

import lapsolver.EdgeList;
import lapsolver.Tree;
import lapsolver.algorithms.TreeSeparator;
import lapsolver.util.TreeUtils;

import java.lang.reflect.Array;
import java.util.Arrays;

public class KelnerFlowTree extends FlowTree {
    private KelnerStructure rootStructure;

    // base constructor
    public KelnerFlowTree (Tree tree, EdgeList offEdges) {
        super(tree, offEdges);
        rootStructure = new KelnerStructure(tree);
    }

    // push `alpha` units of flow along the tree path on edge e
    public void treeUpdate(int e, double alpha) {
        rootStructure.update(offEdges.u[e], alpha);
        rootStructure.update(offEdges.v[e], -alpha);
    }

    // find sum of V = IR along the tree path on edge e
    public double treeQuery(int e) {
        return rootStructure.query(offEdges.u[e]) - rootStructure.query(offEdges.v[e]);
    }

    // initialize the structure with some flows
    public void setTreeFlows(double[] treeFlows) {
        int[] order = TreeUtils.dfsOrder(tree);
        double[] flowTo = new double[tree.nv];

        for (int i = tree.nv-1; i >= 1; i--) {
            int v = order[i];
            int parent = tree.getNode(v).getParent().getId();
            flowTo[parent] += treeFlows[v];
        }

        for (int i = 1; i < tree.nv; i++) {
            rootStructure.update(i, treeFlows[i] - flowTo[i]);
        }
    }

    // retrieve the tree flows
    public double[] getTreeFlows() {
        double[] voltages = new double[tree.nv];
        double[] flowUp = new double[tree.nv];

        for (int i = 0; i < tree.nv; i++) {
            voltages[i] = rootStructure.query(i);
        }

        for (int i = 0; i < tree.nv; i++) {
            int parent = tree.getNode(i).getParent().getId();
            double resistance = tree.getNode(i).getLength();
            flowUp[i] = (voltages[i] - voltages[parent]) / resistance;
        }

        return flowUp;
    }

    // dump the tree structure
    public void dump () {
        dump(rootStructure,0);
    }

    public void dump (KelnerStructure root, int depth) {
        for(int i=0; i<depth; ++i) System.out.print(" ");
        System.out.println(root.tree.nv + " " + root.drop + " " + root.ext);
        for(int i=0; i<depth; ++i) System.out.print(" ");
        System.out.println(Arrays.toString(root.height));
        if(root.children != null)
        for(KelnerStructure ch : root.children) {
            if (ch != null) dump(ch,depth+1);
        }
    }

    // recursively defined data structure
    private class KelnerStructure {
        // structure
        public Tree tree;
        public int separator;
        public KelnerStructure[] children;
        public int[] component; // component id for each vertex
        public int[] relabel; // names of vertices in substructures
        public double[] height; // length of root -> LCA(v, separator) path

        // contents at each substructure
        public double drop;
        public double ext;

        // recursive constructor
        public KelnerStructure(Tree tree) {
            this.tree = tree;
            drop = 0;
            ext = 0;

            if (tree.nv == 2) {
                // base case
                separator = 1 - tree.getRoot(); // whichever isn't the root

                height = new double[2];
                height[separator] = tree.getNode(separator).getLength();

                children = null;
            }
            else {
                // split into subtrees by separator
                separator = TreeSeparator.find(tree);
                int[] order = TreeUtils.dfsOrder(tree);

                // separator has index 0 in each component
                relabel = new int[tree.nv];
                relabel[separator] = 0;

                // need to flip child-parent relation for edges going from sep to root
                boolean[] sepToRoot = new boolean[tree.nv];
                int sepAncestor = separator;
                while (sepAncestor != tree.getRoot()) {
                    sepToRoot[sepAncestor] = true;
                    sepAncestor = tree.getNode(sepAncestor).getParent().getId();
                }

                // compute heights
                height = new double[tree.nv];
                for (int v : order) {
                    int parent = tree.getNode(v).getParent().getId();

                    if (sepToRoot[v]) {
                        height[v] = height[parent] + tree.getNode(v).getLength();
                    }
                    else {
                        height[v] = height[parent];
                    }
                }

                // colour component ids
                // 0 = parent's subtree
                component = new int[tree.nv];
                int[] componentSize = new int[tree.nv];
                int componentCount = 1;
                component[separator] = -1;

                // dfs to find child components
                for (int v : order) {
                    int parent = tree.getNode(v).getParent().getId();
                    if (v != tree.getRoot() && parent == separator) { // new child component
                        component[v] = componentCount;
                        componentSize[componentCount] = 1;
                        relabel[v] = 1;
                        componentCount++;
                    }
                    else if (component[parent] > 0) { // old child component
                        component[v] = component[parent];
                        componentSize[component[v]]++;
                        relabel[v] = componentSize[component[v]];
                    }
                    else if (v != separator) { // parent component
                        componentSize[0]++;
                        relabel[v] = componentSize[component[v]];
                    }
                }

                // initialize parent arrays
                int[][] parentArrays = new int[componentCount][];
                double[][] weights = new double[componentCount][];

                for (int i = 0; i < componentCount; i++) {
                    parentArrays[i] = new int[componentSize[i]+1];
                    weights[i] = new double[componentSize[i]+1];
                }

                parentArrays[0][relabel[tree.getRoot()]] = relabel[tree.getRoot()];

                // build parent arrays
                for (int i = 0; i < tree.nv; i++) {
                    int parent = tree.getNode(i).getParent().getId();
                    int comp = component[i];

                    if (comp == -1) comp = component[parent];
                    if (comp == -1) continue;

                    parentArrays[comp][relabel[i]] = relabel[parent];
                    weights[comp][relabel[i]] = 1 / tree.getNode(i).getLength();
                }

                // build trees from parent arrays
                Tree[] subtrees = new Tree[componentCount];
                for (int i = 0; i < componentCount; i++) {
                    if (separator == tree.getRoot() && i == 0) {
                        // leave parent subtree empty
                        subtrees[i] = null;
                    }
                    else {
                        subtrees[i] = new Tree(parentArrays[i], weights[i]);
                    }
                }

                // call constructor on child structures
                children = new KelnerStructure[componentCount];
                for (int i = 0; i < componentCount; i++) {
                    if (subtrees[i] != null) {
                        children[i] = new KelnerStructure(subtrees[i]);
                    }
                }
            }
        }

        // recursive update
        public void update(int v, double alpha) {
            drop += alpha * height[v];
            if (tree.nv == 2) return;
            if (component[v] > 0) ext += alpha;
            if (v != separator) {
                children[component[v]].update(relabel[v], alpha);
            }
        }

        // recursive query
        public double query(int v) {
            if (v == separator) return drop;
            if (tree.nv == 2) return 0;
            if (component[v] <= 0) {
                return ext * height[v] + children[0].query(relabel[v]);
            }
            return drop + children[component[v]].query(relabel[v]);
        }

    }
}
