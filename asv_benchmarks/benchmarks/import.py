class Import:
    """Benchmark importing pywatershed"""

    def timeraw_import_pywatershed(self):
        return "import pynhm"

    def timeraw_import_pywatershed_only(self):
        return "import pynhm", "import numpy"
