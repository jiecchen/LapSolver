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

import java.util.Arrays;

public class KelnerFlowTree extends FlowTree {
    private KelnerStructure rootStructure;

    // delete me
    public KelnerFlowTree (Tree tree) {
        super(tree, new EdgeList(0));
        rootStructure = new KelnerStructure(tree);
    }

    // base constructor
    public KelnerFlowTree (Tree tree, EdgeList offEdges) {
        super(tree, offEdges);
        rootStructure = new KelnerStructure(tree);
    }

    // push `alpha` units of flow along the tree path on edge e
    public void treeUpdate(int e, double alpha) {

    }

    // find sum of V = IR along the tree path on edge e
    public double treeQuery(int e) {
        return 0.0;
    }

    // initialize the structure with some flows
    public void setTreeFlows(double[] treeFlows) {

    }

    // retrieve the tree flows (off-tree flows are redundant)
    public double[] getTreeFlows() {
        return null;
    }

    // recursively defined data structure
    private class KelnerStructure {
        // structure
        public Tree tree;
        public int separator;
        public KelnerStructure[] children;
        public int[] component; //
        public int[] relabel; // names of vertices in substructures

        // contents at each substructure
        public double drop;
        public double ext;

        // recursive constructor
        public KelnerStructure(Tree tree) {
            this.tree = tree;
            drop = 0;
            ext = 0;

            if (tree.nv <= 2) {
                // base case
                children = null;
            }
            else {
                // split into subtrees by separator
                EdgeList edges = new EdgeList(tree);
                separator = TreeSeparator.find(tree);

                // separator has index 0 in each component
                relabel = new int[tree.nv];
                relabel[separator] = 0;

                // need to flip child-parent relation for edges going from sep to root
                boolean[] flip = new boolean[tree.nv];
                int sepAncestor = separator;
                while (sepAncestor != tree.getRoot()) {
                    flip[sepAncestor] = true;
                    sepAncestor = tree.getNode(sepAncestor).getParent().getId();
                }

                // colour component ids
                // 0 = parent's subtree
                component = new int[tree.nv];
                int[] componentSize = new int[tree.nv];
                int[] order = TreeUtils.dfsOrder(tree);
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

                // build parent arrays
                for (int i = 0; i < tree.nv; i++) {
                    int parent = tree.getNode(i).getParent().getId();
                    int comp = component[i];

                    if (comp == -1) comp = component[parent];
                    if (comp == -1) continue;

                    if (flip[i]) {
                        parentArrays[comp][relabel[parent]] = relabel[i];
                        weights[comp][relabel[parent]] = 1 / tree.getNode(i).getLength();
                    }
                    else {
                        parentArrays[comp][relabel[i]] = relabel[parent];
                        weights[comp][relabel[i]] = 1 / tree.getNode(i).getLength();
                    }
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
    }
}
