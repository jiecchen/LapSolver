from jpype import *
from lapsolver import _GeneratedObject, _java_typecast, _from_java


class TreeUtils(_GeneratedObject):
    _java_class = JPackage("lapsolver").util.TreeUtils

    def __init__(self, *args, **kwargs):
        super(TreeUtils, self).__init__(*args, **kwargs)
        try:
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = TreeUtils._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.util: TreeUtils: incorrect arguments")

    @staticmethod
    def bfsOrder(*args):
        return _from_java(TreeUtils._java_class.bfsOrder(*_java_typecast(*args)))

    @staticmethod
    def depthFromTreeOrder(*args):
        return _from_java(TreeUtils._java_class.depthFromTreeOrder(*_java_typecast(*args)))

    @staticmethod
    def dfsOrder(*args):
        return _from_java(TreeUtils._java_class.dfsOrder(*_java_typecast(*args)))

    @staticmethod
    def getDepths(*args):
        return _from_java(TreeUtils._java_class.getDepths(*_java_typecast(*args)))

    @staticmethod
    def dumpBFSTree(*args):
        return _from_java(TreeUtils._java_class.dumpBFSTree(*_java_typecast(*args)))

    @staticmethod
    def getOffTreeEdges(*args):
        return _from_java(TreeUtils._java_class.getOffTreeEdges(*_java_typecast(*args)))


class Logger(_GeneratedObject):
    _java_class = JPackage("lapsolver").util.Logger

    def __init__(self, *args, **kwargs):
        super(Logger, self).__init__(*args, **kwargs)
        try:
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = Logger._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.util: Logger: incorrect arguments")


class GraphUtils(_GeneratedObject):
    _java_class = JPackage("lapsolver").util.GraphUtils

    def __init__(self, *args, **kwargs):
        super(GraphUtils, self).__init__(*args, **kwargs)
        try:
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = GraphUtils._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.util: GraphUtils: incorrect arguments")

    @staticmethod
    def dump(*args):
        return _from_java(GraphUtils._java_class.dump(*_java_typecast(*args)))

    @staticmethod
    def toTree(*args):
        return _from_java(GraphUtils._java_class.toTree(*_java_typecast(*args)))

    @staticmethod
    def toTree(*args):
        return _from_java(GraphUtils._java_class.toTree(*_java_typecast(*args)))


