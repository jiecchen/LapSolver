/**
 * @file Tree.java
 * @author Daniel Spielman <spielman@cs.yale.edu>
 * @date Wed Jan 8 2014
 *
 * (c) Daniel Spielman, 2014, part of the YINSmat package
 *
 * A static tree data structure.  Intended for use with Matlab.
 * It's main use is for computing the stretch of a tree.
 *
 * It may also compute some simple trees.
 *
 * This is not optimized for performance.
 *
 * Edges of the tree can have lengths, which should be reciprocals of weights.
 * When given as inputs, they will be as weights.
 *
 * Issue: this code does not check if the input it gets is really a tree.
 */

package lapsolver;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Random;

/**
 * a class enabling tree operations
 * that are useful for preconditioning
 * Fat because it is implemented for easy code,
 * not speed
 */
public class Tree {
    static final boolean trace = false;

    public int nv;  // number of verts

    public TreeNode[] nodes;

    public int root;

    public double[] depth;
    public double[] compDepth;
    public int[] order;


    public ArrayDeque<Integer> aux;
    public ArrayDeque<Integer> aux2;
    public int[] comp;
    public int[] nodeSize;
    WeightedGraph stretchGraph;
    private WeightedGraph E;

    /**
     * Make an empty Tree
     *
     * @param nv number of nodes
     * @return an empty Tree of nv nodes
     */
    public Tree(int nv) {
        nodes = new TreeNode[nv];
    }

    /**
     * Make a Tree from parent array
     *
     * @param p the parent array
     * @return a Tree
     */
    public Tree(int[] p) {
        nodes = new TreeNode[p.length];
        this.setFromArray(p);
    }

    /**
     * Builds tree from parent pointers.
     * the root should have a self-loop.
     *
     * @param p pointers to parents
     */
    public void setFromArray(int[] p) {
        nv = p.length;

        nodes = new TreeNode[nv];

        for (int i = 0; i < nv; i++) {
            nodes[i] = new TreeNode(p[i]);
            if (p[i] == i)
                root = i;
        }

        setKidsFromParents();
    }

    /**
     * Assuming just parents are set,
     * fills in all the kids.
     */
    public void setKidsFromParents() {
        for (int i = 0; i < nv; i++)
            if (i != root)
                nodes[nodes[i].parent].addKid(i);
    }


    /**
     * Make a Tree from parent array and weights on edges
     *
     * @param p the parent array
     * @param w the weight array
     * @return a Tree
     */
    public Tree(int[] p, double[] w) {
        nodes = new TreeNode[p.length];
        this.setFromArray(p, w);
    }

    /**
     * Builds tree from parent pointers.
     * the root should have a self-loop.
     *
     * @param p pointers to parents
     * @param w weights on edges (which now turn to 1/length)
     */
    public void setFromArray(int[] p, double[] w) {
        nv = p.length;

        nodes = new TreeNode[nv];

        for (int i = 0; i < nv; i++) {
            nodes[i] = new TreeNode(p[i], 1 / w[i]);
            if (p[i] == i)
                root = i;
        }

        setKidsFromParents();
    }

    /**
     * @param root from which we grow the tree
     * @param G    a Graph
     * @return the Tree
     * <p/>
     * based on Anup Rao's idea for a random tree
     */
    public static Tree growRandTree(Graph G, int root) {
        if (trace) System.out.println("growRandTree");

        boolean[] visited;
        int[] pArray;

        Sampler2 s;

        // Logger log = new Logger();
        // log.start("logx");


        pArray = new int[G.nv];
        visited = new boolean[G.nv];

        int ne = 0;  // number edges
        for (int x = 0; x < G.nv; x++)
            ne += G.deg[x];

        s = new Sampler2(ne, java.lang.System.currentTimeMillis());

        for (int x = 0; x < G.nv; x++)
            visited[x] = false;

        pArray[root] = root;

        // add the edges from root
        visited[root] = true;

        // log.write("root: " + root);

        for (int i = 0; i < G.deg[root]; i++) {
            int nbr = G.nbrs[root][i];
            s.add(root, nbr);
        }

        while (s.last >= 0) {
            int u, v;
            // log.write("s.last : " + s.last);

            int[] out = s.poprand();
            u = out[0];
            v = out[1];
            // log.write("popped : (" + u + " , " + v + ")");

            // only do if bdry edge
            if (!visited[v]) {
                // u is old, visited is new
                pArray[v] = u;

                // log.write("parent of : " + v + " is " + u);

                visited[v] = true;

                for (int i = 0; i < G.deg[v]; i++) {
                    int nbr = G.nbrs[v][i];
                    // log.write("nbr : " + nbr + " of " + v);
                    s.add(v, nbr);
                }
            }
        }
        return new Tree(pArray);
    }

    /**
     * @param root from which we grow the tree
     * @param G    a Graph
     * @return the parent array of the tree
     * <p/>
     * based on Anup Rao's idea for a random tree,
     * with Dan's mod of 2^(-d) sampling
     */
    public static int[] growRandTreeD(Graph G, int root) {
        if (trace) System.out.println("growRandTreeD");

        Random rand = new Random();
        boolean[] visited = new boolean[G.nv];
        int[] pArray = new int[G.nv];
        int[] depth = new int[G.nv];
        int maxdepth, mindepth;
        Sampler2[] s = new Sampler2[G.nv];
        int edgesInQueue = 0;

        depth[root] = 0;

        for (int i = 0; i < G.nv; i++)
            s[i] = new Sampler2(4, i + java.lang.System.currentTimeMillis());

        int ne = 0;  // number edges
        for (int x = 0; x < G.nv; x++)
            ne += G.deg[x];

        for (int x = 0; x < G.nv; x++)
            visited[x] = false;

        pArray[root] = root;

        // add the edges from root
        visited[root] = true;

        for (int i = 0; i < G.deg[root]; i++) {
            int nbr = G.nbrs[root][i];
            s[0].add(root, nbr);
            edgesInQueue++;
        }

        // here, our first depth is 0 for nbrs of root
        mindepth = 0;
        maxdepth = 0;

        // Logger log = new Logger();
        // log.start("logx");

        while (edgesInQueue > 0) {
            int u, v;

            // need to choose
            double[] cums = new double[maxdepth - mindepth + 1];
            double cs = 0;
            double pow2 = 1;
            for (int i = mindepth; i <= maxdepth; i++) {
                cs = cs + pow2 * ((double) s[i].last);
                cums[i - mindepth] = cs;
                pow2 = pow2 / 2;
                // log.write("cums " + i + " val " + cs);
            }

            double r = cs * rand.nextFloat();
            // log.write("r: " + r);

            int lev;
            for (lev = mindepth; cums[lev - mindepth] < r; lev++) ;

            // log.write("lev: " + lev + " last " + s[lev].last);

            int[] out = s[lev].poprand();
            u = out[0];
            v = out[1];

            edgesInQueue--;

            if (lev == mindepth)
                while ((s[mindepth].last < 0) && (mindepth < maxdepth))
                    mindepth++;

            // only do if bdry edge
            if (!visited[v]) {
                // u is old, visited is new
                pArray[v] = u;
                depth[v] = depth[u] + 1;
                if (depth[v] + 1 > maxdepth)
                    maxdepth = depth[v] + 1;

                visited[v] = true;

                for (int i = 0; i < G.deg[v]; i++) {
                    int nbr = G.nbrs[v][i];
                    s[depth[v] + 1].add(v, nbr);
                    edgesInQueue = edgesInQueue + 1;
                }
            }
        }

        return pArray;
    }

    /**
     * Outputs a parent pointer array for the tree
     * The root points to itself
     *
     * @return a parent pointer array
     */
    public int[] parentArray() {
        int[] p = new int[nv];

        for (int i = 0; i < nv; i++) {
            p[i] = nodes[i].parent;
        }

        return p;
    }

    /**
     * Computes the total stretch
     *
     * @param Ein graph of edges to stretch
     * @return the sum of the stretches
     */
    public double compTotalStretch(WeightedGraph Ein) {
        E = Ein;

        this.bfsOrderAndDepth();

        comp = new int[nv];   // the component of a node

        // the algorithm will create a new component
        // each time it encounters a node with no kids

        nodeSize = new int[nv];  // how much is below a node

        int orderPtr; // pointer into the order
        int node;     // and, its node
        int numComps = 0;

        aux  = new ArrayDeque<Integer>(nv);  // will use for traversing
        aux2 = new ArrayDeque<Integer>(nv);  // will use for traversing

        double totStretch = 0;

        for (orderPtr = nv - 1; orderPtr >= 0; orderPtr--) {
            node = order[orderPtr];

            if (nodes[node].numKids() == 0) {
                comp[node] = ++numComps;
                nodeSize[node] = E.deg[node];
                // nothing to merge
            } else {
                // find largest kid, and set this node size
                nodeSize[node] = E.deg[node];

                int maxKid = -1;
                int maxKidSize = -1;

                for (int i = 0; i < nodes[node].numKids(); i++) {
                    int kid = nodes[node].getKid(i);
                    int kidSize = nodeSize[kid];

                    nodeSize[node] += kidSize;

                    if (kidSize > maxKidSize) {
                        maxKidSize = kidSize;
                        maxKid = kid;
                    }
                }

                int curComp = comp[maxKid];

                // go through other kids,
                // re-componentizing,
                // and checking edges for new merge (if other end is this comp)

                for (int i = 0; i < nodes[node].numKids(); i++) {
                    int kid = nodes[node].getKid(i);
                    if (kid != maxKid)
                        totStretch += totStretchTraverse(kid, node, curComp);
                }

                // and, still have to fix up the root!
                // note that this code is easy,
                // but costs a branch in totStretchTraverse
                // that is unnecessary
                // would be faster to handle this node
                // on its own (because do not go beneath)
                totStretch += totStretchTraverse(node, node, curComp);
            }
        }
        return totStretch;
    }

    /**
     * walk in BFS order, setting the fields
     * order and depth
     */
    public void bfsOrderAndDepth() {
        order = new int[nv];
        depth = new double[nv];

        depth[root] = 0;
        order[0] = root;

        int orderPtr = 1;
        int curPtr = 1;
        int curNode = root;

        while (curPtr < nv) {
            for (int i = 0; i < nodes[curNode].numKids(); i++)
                order[orderPtr++] = nodes[curNode].getKid(i);

            curNode = order[curPtr++];

            depth[curNode] = depth[nodes[curNode].parent] + nodes[curNode].length;
        }
    }

    /*
     * the subroutine for compStretch
     *
     * first, account for all edges betwen this subtree
     *  and those in comp curComp
     *
     * then, make the component of all these curComp
     */
    public double totStretchTraverse(int kid, int node, int curComp) {
        double stretch = 0;

        // now, traverse
        aux.clear();
        aux.add(kid);

        while (!aux.isEmpty()) {
            int curNode = aux.poll();

            // note: following test is unnecessary
            // except when we are handling the root of the tree

            if (comp[curNode] != curComp) {
                for (int i = 0; i < nodes[curNode].numKids(); i++)
                    aux.add(nodes[curNode].getKid(i));

                // put node on list to change comp
                aux2.add(curNode);

                // now, go over the edges to check if can now compute stretch
                for (int i = 0; i < E.deg[curNode]; i++) {
                    if (comp[E.nbrs[curNode][i]] == curComp) {
                        int nbr = E.nbrs[curNode][i];
                        double len = 1 / E.weights[curNode][i];
                        stretch += (depth[curNode] + depth[nbr] - 2 * depth[node]) / len;
                        /*
                        System.out.println("(" + curNode + ", " + nbr + ") : " +
                                   + (depth[curNode] + depth[nbr] - 2*depth[node]));
                        */
                    }
                }
            }
        }

        //--------------------------------------------
        //
        // finally, recomp all the nodes we've seen

        while (!aux2.isEmpty()) {
            int curNode = aux2.poll();
            comp[curNode] = curComp;
        }

        return stretch;
    }

    /**
     * Computes the stretch of every edge
     *
     * @param Ein graph of edges to stretch
     * @return a weighted graph where the weight is the stretch of every edge
     */
    public WeightedGraph compStretches(WeightedGraph Ein) {
        stretchGraph = Ein.copy();

        // make all stretches zero, just to find errors
        for (int a = 0; a < Ein.nv; a++) {
            for (int i = 0; i < Ein.deg[a]; i++) {
                stretchGraph.weights[a][i] = 0;
            }
        }

        E = Ein;

        this.bfsOrderAndDepth();

        comp = new int[nv];   // the component of a node
        //
        // the algorithm will create a new component
        // each time it encounters a node with no kids

        nodeSize = new int[nv];  // how much is below a node

        int orderPtr; // pointer into the order
        int node;     // and, its node
        int numComps = 0;

        aux  = new ArrayDeque<Integer>(nv);  // will use for traversing
        aux2 = new ArrayDeque<Integer>(nv);  // will use for traversing

        double totStretch = 0;

        for (orderPtr = nv - 1; orderPtr >= 0; orderPtr--) {
            node = order[orderPtr];

            if (nodes[node].numKids() == 0) {
                comp[node] = ++numComps;
                nodeSize[node] = E.deg[node];
                // nothing to merge
            } else {

                // find largest kid, and set this node size
                nodeSize[node] = E.deg[node];

                int maxKid = -1;
                int maxKidSize = -1;

                for (int i = 0; i < nodes[node].numKids(); i++) {
                    int kid = nodes[node].getKid(i);
                    int kidSize = nodeSize[kid];

                    nodeSize[node] += kidSize;

                    if (kidSize > maxKidSize) {
                        maxKidSize = kidSize;
                        maxKid = kid;
                    }
                }

                int curComp = comp[maxKid];

                // go through other kids,
                // re-componentizing,
                // and checking edges for new merge (if other end is this comp)

                for (int i = 0; i < nodes[node].numKids(); i++) {
                    int kid = nodes[node].getKid(i);
                    if (kid != maxKid)
                        totStretch += allStretchTraverse(kid, node, curComp);
                }

                // and, still have to fix up the root!
                // note that this code is easy,
                // but costs a branch in allStretchTraverse
                // that is unnecessary
                // would be faster to handel this node
                // on its own (because do not go beneath)
                totStretch += allStretchTraverse(node, node, curComp);
            }
        }

        return stretchGraph;
    }

    /*
     * the subroutine for compStretch
     *
     * first, account for all edges betwen this subtree
     *  and those in comp curComp
     *
     * then, make the component of all these curComp
     */
    public double allStretchTraverse(int kid, int node, int curComp) {
        double stretch = 0;

        // now, traverse
        aux.clear();
        aux.add(kid);

        while (!aux.isEmpty()) {
            int curNode = aux.poll();

            // note: following test is unnecessary
            // except when we are handling the root of the tree

            if (comp[curNode] != curComp) {
                for (int i = 0; i < nodes[curNode].numKids(); i++)
                    aux.add(nodes[curNode].getKid(i));

                // put node on list to change comp
                aux2.add(curNode);

                // now, go over the edges to check if can now compute stretch
                for (int i = 0; i < E.deg[curNode]; i++) {
                    if (comp[E.nbrs[curNode][i]] == curComp) {
                        int nbr = E.nbrs[curNode][i];
                        double len = 1 / E.weights[curNode][i];
                        double thisStretch = (depth[curNode] + depth[nbr] - 2 * depth[node]) / len;
                        stretch += thisStretch;
                        stretchGraph.weights[curNode][i] = thisStretch;

                /*
                System.out.println("(" + curNode + ", " + nbr + ") : " +
                           + (depth[curNode] + depth[nbr] - 2*depth[node]));
                */
                    }
                }
            }
        }

        //--------------------------------------------
        //
        // finally, recomp all the nodes we've seen
        while (!aux2.isEmpty()) {
            int curNode = aux2.poll();
            comp[curNode] = curComp;
        }

        return stretch;
    }

    /*
     * A node of the tree
     * default length is 1
     * length is length of edge to the parent
     */
    public class TreeNode {
        public int parent;
        public double length;
        public ArrayList<Integer> kids;

        /**
         * Create a node by specifying its parent
         *
         * @param p parent
         */
        public TreeNode(int p) {
            parent = p;
            length = 1;
            kids = new ArrayList<Integer>();
        }

        public TreeNode(int p, double length) {
            parent = p;
            this.length = length;
            kids = new ArrayList<Integer>();
        }

        public void setParent(int p) {
            parent = p;
        }

        public void addKid(int k) {
            kids.add(k);
        }

        public int numKids() {
            return kids.size();
        }

        public int getKid(int i) {
            return kids.get(i);
        }
    }
}