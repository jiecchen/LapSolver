# !/usr/bin/env python
from scipy.sparse import *
from scipy.spatial import Delaunay
import numpy as np
import time

import lapsolver
import lapsolver.lsst
import lapsolver.algorithms

# Timer from StackOverflow question "tic, toc functions analog in Python"
# -- http://stackoverflow.com/questions/5849800
class Timer(object):
    def __init__(self, name=None):
        self.name = name
        self.tstart = None

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        if self.name:
            print '[%s]' % self.name,
        print 'Elapsed: %s' % (time.time() - self.tstart)


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
        Generates a 3D Delaunay graph, uses Dan's SimulPathTree to
        generate a low-stretch spanning tree, it then computes its
        mean stretch.
    :return: mean stretch of the random Delaunay tree
    """
    with Timer('generate'):
        a = del3_graph(10000)

    [ai, aj] = tril(a).nonzero()
    av = np.asarray(a[(ai, aj)])[0, :]

    g = lapsolver.Graph(ai, aj, av)
    with Timer('generate the tree'):
        spt = lapsolver.lsst.SimulPathTree()
        tree = spt.getTree(g)
    with Timer('compute stretch'):
        total_stretch = lapsolver.algorithms.Stretch.compute(g, tree).total
    return total_stretch / len(ai.tolist())

with Timer('total'):
    print test_del3()
