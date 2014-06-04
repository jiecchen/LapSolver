import numpy as np
import sys, os
from jpype import *

prefix = os.path.dirname(os.path.realpath(__file__))
startJVM(getDefaultJVMPath(), '-Djava.class.path={}'.format(prefix + "/_src/LapSolver.jar"))
del prefix


def _from_java(val):
    try:
        className = type(val).__name__.split('.')
        module = __import__('.'.join(className[0:-1]))
        class_ = getattr(module, className[-1])
        return class_(val, fromJVM=True)
    except:
        return val


def _java_typecast(*args):
    values = []
    for arg in args:
        value = arg  # hope for the best
        if isinstance(arg, _GeneratedObject):
            value = arg.toJava()
        elif isinstance(arg, np.ndarray):
            if "int8" in str(arg.dtype):
                jtype = JByte
            elif "int" in str(arg.dtype):
                jtype = JInt
            elif "bool" in str(arg.dtype):
                jtype = JBoolean
            elif "float64" in str(arg.dtype):
                jtype = JDouble
            elif "float" in str(arg.dtype):
                jtype = JFloat
            else:
                print "error: incompatible type in ndarray `{}`".format(arg.dtype)
                sys.exit(0)
            value = JArray(jtype, arg.ndim)(arg.tolist())
        values.append(value)
    return values


class _GeneratedObject(object):
    def __init__(self, *args, **kwargs):
        self._instance = None

    def __getattr__(self, name):
        try:
            return object.__getattribute__(self, name)
        except AttributeError:
            attr = getattr(self._instance, name)
            if type(attr).__name__ == "JavaBoundMethod":
                def wrap(*args):
                    return _from_java(attr(*_java_typecast(*args)))
                return wrap
            return attr

    def toJava(self):
        return self._instance


class Logger(_GeneratedObject):
    def __init__(self, *args, **kwargs):
        super(Logger, self).__init__(*args, **kwargs)
        try:
            self._java_class = JPackage("lapsolver").Logger
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = self._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver: Logger: incorrect arguments")


class WeightedGraph(_GeneratedObject):
    def __init__(self, *args, **kwargs):
        super(WeightedGraph, self).__init__(*args, **kwargs)
        try:
            self._java_class = JPackage("lapsolver").WeightedGraph
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = self._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver: WeightedGraph: incorrect arguments")


class Tree(_GeneratedObject):
    def __init__(self, *args, **kwargs):
        super(Tree, self).__init__(*args, **kwargs)
        try:
            self._java_class = JPackage("lapsolver").Tree
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = self._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver: Tree: incorrect arguments")


class Sampler2(_GeneratedObject):
    def __init__(self, *args, **kwargs):
        super(Sampler2, self).__init__(*args, **kwargs)
        try:
            self._java_class = JPackage("lapsolver").Sampler2
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = self._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver: Sampler2: incorrect arguments")


class Graph(_GeneratedObject):
    def __init__(self, *args, **kwargs):
        super(Graph, self).__init__(*args, **kwargs)
        try:
            self._java_class = JPackage("lapsolver").Graph
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = self._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver: Graph: incorrect arguments")


