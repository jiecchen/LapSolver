from jpype import *
from lapsolver import _GeneratedObject, _java_typecast


class PetalDecompositonLSST(_GeneratedObject):
    def __init__(self, *args, **kwargs):
        super(PetalDecompositonLSST, self).__init__(*args, **kwargs)
        try:
            self._java_class = JPackage("lapsolver").lsst.PetalDecompositonLSST
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = self._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.lsst: PetalDecompositonLSST: incorrect arguments")


class PrimsLSST(_GeneratedObject):
    def __init__(self, *args, **kwargs):
        super(PrimsLSST, self).__init__(*args, **kwargs)
        try:
            self._java_class = JPackage("lapsolver").lsst.PrimsLSST
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = self._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.lsst: PrimsLSST: incorrect arguments")


class SimulPathLSST(_GeneratedObject):
    def __init__(self, *args, **kwargs):
        super(SimulPathLSST, self).__init__(*args, **kwargs)
        try:
            self._java_class = JPackage("lapsolver").lsst.SimulPathLSST
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = self._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.lsst: SimulPathLSST: incorrect arguments")


