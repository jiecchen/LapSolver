#!/usr/bin/env python
from scipy.sparse import *
from scipy.spatial import Delaunay
import numpy as np
import lapsolver as ls


def del3_graph(n):
    points = np.random.rand(n, 3)
    tri = Delaunay(points).simplices
    (nt, _) = tri.shape

    ones = np.ones((nt, 1))[:, 0]
    a = csr_matrix((ones, (tri[:, 0], tri[:, 1])), shape=(n, n))
    a = a + csr_matrix((ones, (tri[:, 0], tri[:, 2])), shape=(n, n))
    a = a + csr_matrix((ones, (tri[:, 0], tri[:, 3])), shape=(n, n))
    a = a + csr_matrix((ones, (tri[:, 1], tri[:, 2])), shape=(n, n))
    a = a + csr_matrix((ones, (tri[:, 1], tri[:, 3])), shape=(n, n))
    a = a + csr_matrix((ones, (tri[:, 2], tri[:, 3])), shape=(n, n))

    return (a + a.transpose()).sign()


def test_del3():
    """
    Tests the Java -> Python pipeline.
        Generates a 3D Delaunay graph, uses Dan's SimulPathLSST to
        generate a low-stretch spanning tree, it then computes its
        mean stretch.
    :return:
    """
    a = del3_graph(10000)

    [ai, aj] = tril(a).nonzero()
    av = np.asarray(a[(ai, aj)])[0, :]

    g = ls.WeightedGraph(ai, aj, av)

    spt = ls.SimulPathLSST(g)
    return spt.edgeGrow().treeToTree().compTotalStretch(g.toJava()) / len(ai.tolist())

print test_del3()
