import math
import pandas as pd

created_line = 'Created '
hash_line = '##############'
hash_line_official = '########################################'


def cbh_to_df(the_file):
    meta_lines = []
    with open(the_file, 'r') as file_open:
    # file_open =  open(the_file, 'r')
        file_lines = file_open.readlines()
        wh_hash_line = 0
        while not hash_line in file_lines[wh_hash_line]:
            the_line = file_lines[wh_hash_line]
            if (not '//' in the_line) and (not created_line in the_line):
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
