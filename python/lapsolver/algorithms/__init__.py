from jpype import *
from lapsolver import _GeneratedObject, _java_typecast, _from_java


class GrowRandomTree(_GeneratedObject):
    _java_class = JPackage("lapsolver").algorithms.GrowRandomTree
    def __init__(self, *args, **kwargs):
        super(GrowRandomTree, self).__init__(*args, **kwargs)
        try:
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = GrowRandomTree._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.algorithms: GrowRandomTree: incorrect arguments")

    @staticmethod
    def growRandTree(*args):
        return _from_java(GrowRandomTree._java_class.growRandTree(*_java_typecast(*args)))

    @staticmethod
    def growRandTreeD(*args):
        return _from_java(GrowRandomTree._java_class.growRandTreeD(*_java_typecast(*args)))


class TreePath(_GeneratedObject):
    _java_class = JPackage("lapsolver").algorithms.TreePath
    def __init__(self, *args, **kwargs):
        super(TreePath, self).__init__(*args, **kwargs)
        try:
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = TreePath._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.algorithms: TreePath: incorrect arguments")


class Congestion(_GeneratedObject):
    _java_class = JPackage("lapsolver").algorithms.Congestion
    def __init__(self, *args, **kwargs):
        super(Congestion, self).__init__(*args, **kwargs)
        try:
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = Congestion._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.algorithms: Congestion: incorrect arguments")

    @staticmethod
    def compute(*args):
        return _from_java(Congestion._java_class.compute(*_java_typecast(*args)))


class UnionFind(_GeneratedObject):
    _java_class = JPackage("lapsolver").algorithms.UnionFind
    def __init__(self, *args, **kwargs):
        super(UnionFind, self).__init__(*args, **kwargs)
        try:
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = UnionFind._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.algorithms: UnionFind: incorrect arguments")


class Stretch(_GeneratedObject):
    _java_class = JPackage("lapsolver").algorithms.Stretch
    def __init__(self, *args, **kwargs):
        super(Stretch, self).__init__(*args, **kwargs)
        try:
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = Stretch._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.algorithms: Stretch: incorrect arguments")

    @staticmethod
    def compute(*args):
        return _from_java(Stretch._java_class.compute(*_java_typecast(*args)))


class StretchDan(_GeneratedObject):
    _java_class = JPackage("lapsolver").algorithms.StretchDan
    def __init__(self, *args, **kwargs):
        super(StretchDan, self).__init__(*args, **kwargs)
        try:
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = StretchDan._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.algorithms: StretchDan: incorrect arguments")


class PairSampler(_GeneratedObject):
    _java_class = JPackage("lapsolver").algorithms.PairSampler
    def __init__(self, *args, **kwargs):
        super(PairSampler, self).__init__(*args, **kwargs)
        try:
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = PairSampler._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.algorithms: PairSampler: incorrect arguments")


class DiscreteSampler(_GeneratedObject):
    _java_class = JPackage("lapsolver").algorithms.DiscreteSampler
    def __init__(self, *args, **kwargs):
        super(DiscreteSampler, self).__init__(*args, **kwargs)
        try:
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = DiscreteSampler._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.algorithms: DiscreteSampler: incorrect arguments")


class TarjanLCA(_GeneratedObject):
    _java_class = JPackage("lapsolver").algorithms.TarjanLCA
    def __init__(self, *args, **kwargs):
        super(TarjanLCA, self).__init__(*args, **kwargs)
        try:
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = TarjanLCA._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.algorithms: TarjanLCA: incorrect arguments")


