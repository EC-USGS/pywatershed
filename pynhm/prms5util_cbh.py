from enum import Enum
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
        self._date_columns = [
            0,
            1,
            2,
            3,
            4,
            5,
        ]

        # create dataframe
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
        skip_lines = 0
        column_names = []
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
                                column_names.append(column_name)
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
                    break

                # set previous file position to the start of the next line
                file_position0 = cbh_file.tell()

                # increment skip lines
                skip_lines += 1

        # read the data in file directly using pandas
        incl_cols = [n for n in range(nitems + len(self._date_columns))]

        df = pd.read_csv(
            self.fpth,
            sep=" ",
            skipinitialspace=True,
            usecols=incl_cols,
            skiprows=skip_lines,
            engine="c",
            memory_map=True,
            parse_dates=[self._date_columns],
            index_col=0,
            header=None,
            na_values=(-99.0, -999.0),
        )

        df.index = pd.to_datetime(df.index, format="%Y %m %d %H %M %S")
        df.index.name = "date"
        df.columns = column_names

        return df
