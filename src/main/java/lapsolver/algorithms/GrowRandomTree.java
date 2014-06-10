/**
 * @file Congestion.java
 * @author Cyril Zhang <cyril.zhang@yale.edu>
 * @date Thu Jun 5 2014
 *
 * Anup and Dan's GrowRandTree algorithms, originally members of Tree.
 */

package lapsolver.algorithms;

import lapsolver.UnweightedGraph;
import lapsolver.Tree;
import java.util.Random;

public class GrowRandomTree {
    /**
     * @param root from which we grow the tree
     * @param G    an UnweightedGraph
     * @return the Tree
     * based on Anup Rao's idea for a random tree
     */
    public static Tree growRandTree(UnweightedGraph G, int root) {
        boolean[] visited;
        int[] pArray;

        PairSampler s;

        // Logger log = new Logger();
        // log.start("logx");

        pArray = new int[G.nv];
        visited = new boolean[G.nv];

        int ne = 0;  // number edges
        for (int x = 0; x < G.nv; x++)
            ne += G.deg[x];

        s = new PairSampler(ne, java.lang.System.currentTimeMillis());

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
     * @param G    an UnweightedGraph
     * @return the parent array of the tree
     * <p/>
     * based on Anup Rao's idea for a random tree,
     * with Dan's mod of 2^(-d) sampling
     */
    public static int[] growRandTreeD(UnweightedGraph G, int root) {
        Random rand = new Random();
        boolean[] visited = new boolean[G.nv];
        int[] pArray = new int[G.nv];
        int[] depth = new int[G.nv];
        int maxdepth, mindepth;
        PairSampler[] s = new PairSampler[G.nv];
        int edgesInQueue = 0;

        depth[root] = 0;

        for (int i = 0; i < G.nv; i++)
            s[i] = new PairSampler(4, i + java.lang.System.currentTimeMillis());

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

            int lev = mindepth;
            while (cums[lev - mindepth] < r)
                lev++;

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
}
