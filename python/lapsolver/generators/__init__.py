from jpype import *
from lapsolver import _GeneratedObject, _java_typecast, _from_java


class Grid2(_GeneratedObject):
    _java_class = JPackage("lapsolver").generators.Grid2
    def __init__(self, *args, **kwargs):
        super(Grid2, self).__init__(*args, **kwargs)
        try:
            if "fromJVM" in kwargs and kwargs["fromJVM"]:
                self._instance = args[0]
            else:
                self._instance = Grid2._java_class(*_java_typecast(*args))
        except:
            raise RuntimeError("[ERROR] lapsolver.generators: Grid2: incorrect arguments")


