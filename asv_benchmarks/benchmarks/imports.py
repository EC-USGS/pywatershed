from . import _is_pws


class Import:
    """Benchmark importing pywatershed"""

    def timeraw_import_pywatershed(self):
        if _is_pws:
            return "import pywatershed"
        else:
            return "import pynhm"

    def timeraw_import_pywatershed_only(self):
        if _is_pws:
            return "import pywatershed", "import numpy"
        else:
            return "import pynhm", "import numpy"
