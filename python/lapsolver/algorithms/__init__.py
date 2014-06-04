from jpype import *
from lapsolver import _GeneratedObject, _java_typecast


class Congestion(_GeneratedObject):
    def __init__(self, *args, **kwargs):
        super(Congestion, self).__init__(*args, **kwargs)
        try:
            self._java_class = JPackage("lapsolver").algorithms.Congestion
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = self._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.algorithms: Congestion: incorrect arguments")


class TarjanLCA(_GeneratedObject):
    def __init__(self, *args, **kwargs):
        super(TarjanLCA, self).__init__(*args, **kwargs)
        try:
            self._java_class = JPackage("lapsolver").algorithms.TarjanLCA
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = self._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.algorithms: TarjanLCA: incorrect arguments")


class Stretch(_GeneratedObject):
    def __init__(self, *args, **kwargs):
        super(Stretch, self).__init__(*args, **kwargs)
        try:
            self._java_class = JPackage("lapsolver").algorithms.Stretch
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = self._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.algorithms: Stretch: incorrect arguments")


class UnionFind(_GeneratedObject):
    def __init__(self, *args, **kwargs):
        super(UnionFind, self).__init__(*args, **kwargs)
        try:
            self._java_class = JPackage("lapsolver").algorithms.UnionFind
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = self._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.algorithms: UnionFind: incorrect arguments")


