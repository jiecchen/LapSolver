from jpype import *
from lapsolver import _GeneratedObject, _java_typecast


class PetalDecompositonTree(_GeneratedObject):
    def __init__(self, *args, **kwargs):
        super(PetalDecompositonTree, self).__init__(*args, **kwargs)
        try:
            self._java_class = JPackage("lapsolver").lsst.PetalDecompositonTree
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = self._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.lsst: PetalDecompositonTree: incorrect arguments")


class SimulPathTree(_GeneratedObject):
    def __init__(self, *args, **kwargs):
        super(SimulPathTree, self).__init__(*args, **kwargs)
        try:
            self._java_class = JPackage("lapsolver").lsst.SimulPathTree
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = self._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.lsst: SimulPathTree: incorrect arguments")


class KruskalTree(_GeneratedObject):
    def __init__(self, *args, **kwargs):
        super(KruskalTree, self).__init__(*args, **kwargs)
        try:
            self._java_class = JPackage("lapsolver").lsst.KruskalTree
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = self._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.lsst: KruskalTree: incorrect arguments")


