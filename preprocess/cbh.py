import math
import pandas as pd

class CBH:
    def __init__(
            self,
            files,
            convert_units=False,
            adjust=False,
            output_file=None,
            verbosity=0):

        # Data -----
        self.files = files
        self.convert_units = convert_units
        self.adjust = adjust
        self.output_file = output_file
        self.verbosity = verbosity
        self.df = None

        # Methods -----
        # self.to_csv = _df_to_csv  # pre-deprecated?
        self.to_netcdf = _df_to_netcdf
        self.convert_units_cbh = _convert_units_cbh
        self.adjust_cbh = _adjust_cbh
        self.check_cbh = _check_cbh

        # How to get a data frame depends on the files argument
        if isinstance(self.files, (str)):
            self.to_df = _cbh_file_to_df
        elif isinstance(self.files, (dict)):
            self.to_df = _cbh_files_to_df
        else:
            raise ValueError(f'"files" argument of type {type(files)} not accepted.')

        # Special action -----
        # If the output_file is specified: get it all done
        # (if you dont want the object and you just want the written file side-effect)
        if output_file is not None:
            self.df = self.to_df(self.files)
            if convert_units:
                self.df = self.convert_units(self.df)
            if adjust:
                self.df = self._
                # print the path that was written



        # self.checks

        # accept output file to just write the converted files?

        return None


def _df_to_netcdf():
    pass


def _convert_units_cbh():
    pass


def _adjust_cbh():
    pass


def _check_cbh():
    pass


created_line = 'Created '
written_line = 'Written '
hash_line = '##############'
hash_line_official = '########################################'


def _cbh_file_to_df(the_file):
        # This is trying to handle as many input formats for these kinds of files
        # as we can find. It may not be comprehensive. See tests for what is currently
        # handled
        meta_lines = []
        with open(the_file, 'r') as file_open:
            file_lines = file_open.readlines()
            wh_hash_line = 0
            while hash_line not in file_lines[wh_hash_line]:
                the_line = file_lines[wh_hash_line]
                if (
                        ('//' not in the_line)
                        and (created_line not in the_line)
                        and (written_line not in the_line)):
                    meta_lines += [the_line.strip()]
                wh_hash_line += 1

        col_names = []
        for posn in range(len(meta_lines)):
            key, count = meta_lines[posn].split(' ')
            count = int(count)
            zs = math.ceil(math.log(count, 10))
            col_names += [f"{key}{str(ii).zfill(zs)}" for ii in range(count)]

        # Can we get the hru info? is that standarized at all?

        dtypes = (['str'] * 6) + (['float64'] * len(col_names))
        date_cols = ["Y", "m", "d", "H", "M", "S"]
        col_names = date_cols + col_names
        assert len(dtypes) == len(col_names)
        dtype_dict = dict(zip(col_names, dtypes))
        assert len(dtype_dict) == len(col_names)

        data = pd.read_csv(
            the_file, names=col_names,
            skiprows=wh_hash_line+1,
            delim_whitespace=True, dtype=dtype_dict)

        data['date'] = pd.to_datetime(
            data.Y + '-'
            + data.m.str.zfill(2) + '-'
            + data.d.str.zfill(2))
        data = data.drop(columns=set(date_cols))
        data = data.set_index('date')

        return data


def _cbh_files_to_df():
        pass
