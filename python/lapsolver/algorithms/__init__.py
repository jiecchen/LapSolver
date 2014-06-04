from jpype import *
from lapsolver import _GeneratedObject, _java_typecast


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


