__author__ = 'YINS'  # or whoever really owns the rights
import os
from jpype import *

try:
    prefix = os.environ['LAP_PATH']
    if prefix[-1] != '/':
        prefix += '/'
except KeyError:
    prefix = "../"  # Assume we're being called from the python/ directory

startJVM(getDefaultJVMPath(), '-Djava.class.path={}'.format(prefix + "build/LapSolver.jar"))
del prefix

# TODO: implement these
# Graph = JPackage("lapsolver").Graph
# Logger = JPackage("lapsolver").Logger
# Sampler2 = JPackage("lapsolver").Sampler2
# Tree = JPackage("lapsolver").Tree


class WeightedGraph(object):
    def __init__(self, *args):
        self._java_class = JPackage("lapsolver").WeightedGraph
        if len(args) == 0:
            self._instance = self._java_class()
        elif len(args) == 3:
            [ai, aj, av] = args
            self._instance = self._java_class(JArray(JInt, ai.ndim)(ai.tolist()),
                                              JArray(JInt, aj.ndim)(aj.tolist()),
                                              JArray(JDouble, av.ndim)(av.tolist()))
        else:
            raise RuntimeError("[ERROR] lapsolver: WeightedGraph: incorrect arguments")

    def __getattr__(self, name):
        try:
            return object.__getattribute__(self, name)
        except AttributeError:
            return getattr(self._instance, name)

    def toJava(self):
        return self._instance


class SimulPathLSST(object):
    def __init__(self, g):
        self._java_class = (JPackage("lapsolver")).lsst.SimulPathLSST
        self._instance = self._java_class(g.toJava())

    def __getattr__(self, name):
        try:
            return object.__getattribute__(self, name)
        except AttributeError:
            return getattr(self._instance, name)

    def toJava(self):
        return self._instance

