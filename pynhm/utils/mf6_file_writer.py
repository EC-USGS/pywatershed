from ..constants import fileish

import numpy as np


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
            for field, type in fields.items():

                if hasattr(selfish, field):
                    if not wrote_block:
                        wrote_block = True
                        outfile.write(f"BEGIN {block.upper()}\n")

                    val = getattr(selfish, field)
                    # handle None?
                    if isinstance(val, np.ndarray):
                        if len(val.shape):
                            val_type = "vector"
                        else:
                            val_type = "scalar"
                    else:
                        val_type = "scalar"

                    if (type == "scalar") and (val_type == "scalar"):
                        outfile.write(f"{idt}{field.upper()} {val}\n")

                    elif (type == "vector") and (val_type == "scalar"):
                        outfile.write(
                            f"{idt}{field.upper()}\n"
                            f"{idt * 2}CONSTANT {val}\n"
                        )

                    elif (type == "vector") and (val_type == "vector"):
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
