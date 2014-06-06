from jpype import *
from lapsolver import _GeneratedObject, _java_typecast, _from_java


class SimulPathTree(_GeneratedObject):
    _java_class = JPackage("lapsolver").lsst.SimulPathTree
    def __init__(self, *args, **kwargs):
        super(SimulPathTree, self).__init__(*args, **kwargs)
        try:
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = SimulPathTree._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.lsst: SimulPathTree: incorrect arguments")


class PetalDecompositionTree(_GeneratedObject):
    _java_class = JPackage("lapsolver").lsst.PetalDecompositionTree
    def __init__(self, *args, **kwargs):
        super(PetalDecompositionTree, self).__init__(*args, **kwargs)
        try:
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = PetalDecompositionTree._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.lsst: PetalDecompositionTree: incorrect arguments")


class KruskalTree(_GeneratedObject):
    _java_class = JPackage("lapsolver").lsst.KruskalTree
    def __init__(self, *args, **kwargs):
        super(KruskalTree, self).__init__(*args, **kwargs)
        try:
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = KruskalTree._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.lsst: KruskalTree: incorrect arguments")


