from jpype import *
from lapsolver import _GeneratedObject, _java_typecast


class TreeUtils(_GeneratedObject):
    def __init__(self, *args, **kwargs):
        super(TreeUtils, self).__init__(*args, **kwargs)
        try:
            self._java_class = JPackage("lapsolver").util.TreeUtils
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = self._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.util: TreeUtils: incorrect arguments")


class Logger(_GeneratedObject):
    def __init__(self, *args, **kwargs):
        super(Logger, self).__init__(*args, **kwargs)
        try:
            self._java_class = JPackage("lapsolver").util.Logger
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = self._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.util: Logger: incorrect arguments")


class GraphUtils(_GeneratedObject):
    def __init__(self, *args, **kwargs):
        super(GraphUtils, self).__init__(*args, **kwargs)
        try:
            self._java_class = JPackage("lapsolver").util.GraphUtils
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = self._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.util: GraphUtils: incorrect arguments")


