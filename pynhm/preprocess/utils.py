import math
import pandas as pd

created_line = 'Created '
written_line = 'Written '
hash_line = '##############'
hash_line_official = '########################################'


def convert_units_cbh():
    pass


def adjust_cbh():
    pass


def check_cbh():
    pass


def cbh_file_to_df(the_file):
    # This is trying to handle as many input formats for these kinds of files
    # as we can find. It may not be comprehensive. See tests for what is currently
    # handled
    meta_lines = []
    with open(the_file, 'r') as file_open:
        wh_hash_line = -1
        the_line = ''
        while hash_line not in the_line:
            wh_hash_line += 1
            the_line = file_open.readline()
            if (
                    ('//' not in the_line)
                    and (created_line not in the_line)
                    and (written_line not in the_line)
                    and (hash_line not in the_line)):
                meta_lines += [the_line.strip()]

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

    # Some files erroneously specify the number of columns
    # catch that here
    data = pd.read_csv(
        the_file,
        nrows=1,
        # names=col_names,  # dont specify names
        header=None,  # dont use default header from first line
        index_col=False,
        skiprows=wh_hash_line+1,
        delim_whitespace=True, dtype=dtype_dict)
    msg = f"Number of actual data columns does not match metadata info: {meta_lines}"
    assert len(data.columns) == len(col_names), msg
    # JLM: is the above sufficient?

    data = pd.read_csv(
        the_file,
        names=col_names,
        index_col=False,
        skiprows=wh_hash_line+1,
        delim_whitespace=True, dtype=dtype_dict)

    data['date'] = pd.to_datetime(
        data.Y + '-'
        + data.m.str.zfill(2) + '-'
        + data.d.str.zfill(2))
    data = data.drop(columns=set(date_cols))
    data = data.set_index('date')

    return data


def cbh_files_to_df(file_dict):
    dfs = [cbh_file_to_df(val) for val in file_dict.values()]
    return pd.concat(dfs, axis=1)
