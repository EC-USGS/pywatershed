from enum import Enum
from io import StringIO

import pandas as pd


class CbhPosition(Enum):
    first = 0
    comments = 1
    dimensions = 2
    data = 3


class CbhInput:
    def __init__(self, fpth, verbose=False):
        self.fpth = fpth
        self.verbose = verbose
        self._df = self._read_cbh()

    @property
    def get_dataframe(self):
        return self._df

    @property
    def get_dataframe_columns(self):
        return tuple(self._df.columns)

    @property
    def get_dataframe_column(self, column):
        if column in self._df.columns:
            return self._df[column]
        else:
            raise KeyError(f"{column} not in {self.fpth}")

    def _read_cbh(self):
        if self.verbose:
            print(f"Loading {self.fpth}")
        with open(self.fpth, "rb+") as cbh_file:
            str_io = "date"
            max_dimensions = 0
            nitems = 0
            dimensions = {}
            cbh_position = CbhPosition.first
            file_position0 = 0
            for line in cbh_file:
                # skip first line
                if cbh_position == CbhPosition.first:
                    # update position
                    cbh_position = CbhPosition.comments
                # skip comments
                elif cbh_position == CbhPosition.comments:
                    if not line.startswith(b"//"):
                        # update position
                        cbh_position = CbhPosition.dimensions
                        cbh_file.seek(file_position0)
                elif cbh_position == CbhPosition.dimensions:
                    if line.startswith(b"#"):
                        # create header
                        pad_size = len(f"{max_dimensions}")
                        fmt = "{:0" + f"{pad_size}" + "d}"
                        for key, value in dimensions.items():
                            for n in range(value):
                                column_name = f"{key}_hru" + fmt.format(n + 1)
                                str_io += f",{column_name}"
                        str_io += "\n"
                        # update position
                        cbh_position = CbhPosition.data
                    else:
                        line_split = line.strip().split()
                        dimension = int(line_split[1])
                        nitems += dimension
                        max_dimensions = max(max_dimensions, dimension)
                        dimensions[line_split[0].decode("ascii")] = dimension
                elif cbh_position == CbhPosition.data:
                    ll = line.strip().split()
                    yr = int(ll[0])
                    mo = int(ll[1])
                    da = int(ll[2])
                    mn = int(ll[3])
                    sec = int(ll[4])
                    data = [
                        str(vv)
                        for vv in [float(v) for v in ll[5 : 5 + nitems]]
                    ]
                    str_io += (
                        f"{da:02d}/{mo:02d}/{yr:04d} {mn:02d}:{sec:02d},"
                        + f"{','.join(data)}\n"
                    )

                # set previous file position to the start of the next line
                file_position0 = cbh_file.tell()

        # create dataframe from str_io
        df = pd.read_csv(
            StringIO(str_io),
            parse_dates=["date"],
            index_col=["date"],
            dtype=float,
            float_precision="high",
            na_values=(-999,),
        )

        return df
