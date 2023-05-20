class Import:
    """Benchmark importing pywatershed"""

    def timeraw_import_pywatershed(self):
        return "import pywatershed"

    def timeraw_import_pywatershed_only(self):
        return "import pywatershed", "import numpy"
