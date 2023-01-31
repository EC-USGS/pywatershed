from ..constants import fileish

import numpy as np
import pandas as pd
import xarray as xr


skip_blocks = [
    "periods",
]


def mf6_file_writer(
    selfish,
    file_struct: dict,
    required: list,
    output_file: fileish,
    idt=" " * 2,
):
    """A common writer for MF6 block-structured text files.

    Args:
      selfish: pass self from another object to be queried on its available
        attributes for writing to file.
      file_struct: a dictionary that defines the block struture of the MF6
        output file.
      require: a list of keys in file_struct that are required to be in the
        output file.
      idt: the indentation definition, default is 2 spaces.
    """

    if not output_file:
        msg = "no output file specified"
        raise ValueError(msg)

    with open(output_file, "w") as outfile:
        for block, fields in file_struct.items():
            wrote_block = False

            if block in required:
                wrote_block = True
                outfile.write(f"BEGIN {block.upper()}\n")

            for field, ftype in fields.items():

                if hasattr(selfish, field):
                    if not wrote_block and (block not in skip_blocks):
                        wrote_block = True
                        outfile.write(f"BEGIN {block.upper()}\n")

                    # Some times the value is not the type requested
                    # for the field output
                    # I suppose thse could be handled in transformations
                    # on self
                    val = getattr(selfish, field)
                    if val is None:
                        val_type = None
                    elif isinstance(val, np.ndarray):
                        if len(val.shape):
                            val_type = "vector"
                        else:
                            val_type = "scalar"
                    elif isinstance(val, xr.DataArray):
                        val_type = "DataArray"
                    else:
                        val_type = "scalar"

                    if ftype is None:
                        outfile.write(f"{idt}{field.upper()}\n")

                    elif (ftype == "scalar") and (val_type == "scalar"):
                        outfile.write(f"{idt}{field.upper()} {val}\n")

                    elif (ftype == "vector") and (val_type == "scalar"):
                        outfile.write(
                            f"{idt}{field.upper()}\n"
                            f"{idt * 2}CONSTANT {val}\n"
                        )

                    elif (ftype == "vector") and (val_type == "vector"):
                        val_repr = (
                            str(val)
                            .replace("\n", "")
                            .replace("[", "")
                            .replace("]", "")
                            .lstrip()
                        )
                        outfile.write(
                            f"{idt}{field.upper()}\n" f"{idt * 2}{val_repr}\n"
                        )

                    elif ftype == "DataArray" and block == "periods":
                        # Period data
                        loc_ind1 = (
                            np.arange(len(val.nhm_id), dtype="int64") + 1
                        )  # 1-based
                        for ii, tt in enumerate(val.time):
                            period_slice = val[ii, :].values
                            df = pd.DataFrame(
                                {"loc_ind": loc_ind1, "val": period_slice}
                            )
                            loc_vals_str = df.to_csv(
                                index=False, header=False, sep=" "
                            )
                            time = str(val.time.values[ii])
                            outfile.write(f"BEGIN period {ii+1}  # {time}\n")
                            outfile.write(
                                idt
                                + loc_vals_str.replace("\n", f"\n{idt}")[
                                    0 : -(len(idt))
                                ]
                            )
                            outfile.write(f"END period {ii+1}\n\n")

                    else:
                        msg = "unhandled combo of expected and actual types"
                        raise ValueError(msg)

                elif field in required:
                    msg = "Required disu field not found: {field}"
                    raise ValueError(msg)

            if wrote_block:
                outfile.write(f"END {block.upper()}\n\n")

    # print a message that file was written? verbose?
    return
