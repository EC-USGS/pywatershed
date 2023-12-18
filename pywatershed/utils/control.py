import pathlib as pl
from typing import Union
from warnings import warn

import numpy as np

from .prms5_file_util import PrmsFile

fileish = Union[str, pl.PosixPath, dict]
listish = Union[str, list, tuple]


class ControlVariables:
    """PRMS control file class

    Args:
        control_dict: control variable dictionary
    """

    def __init__(self, control_dict: dict) -> "ControlVariables":
        self.control = control_dict["control"]

    def get_variables(self, keys: listish) -> "ControlVariables":
        """Get a subset of keys in the control variable dictionary

        Args:
            keys: keys to retrieve from the full PRMS control variable object

        Returns:
            ControlVariables : subset of full control variable dictionary
                Passed keys that do not exist in the full control variable
                dictionary are skipped.

        """
        if isinstance(keys, str):
            keys = [keys]

        return ControlVariables(
            {
                key: self.control.get(key)
                for key in keys
                if key in self.control.keys()
            }
        )

    @staticmethod
    def load(control_file: fileish) -> "ControlVariables":
        """Load variables from a PRMS control file

        Args:
            control_file: control file path

        Returns:
            ControlVariables: full PRMS control variable dictionary

        """
        return ControlVariables(
            PrmsFile(control_file, file_type="control").get_data()
        )


def compare_control_files(file0, file1, silent=False):
    """Compare the contents of two control files.

    Args:
        file0: str or pathlib.Path for the first PRMS control file
        file1: str or pathlib.Path for the second PRMS control file
        silent: suppress informative warnings (if you only want the return
            value)

    Returns:
        A bash-style integer where zero indicates the files are identical and
        all other values indicate the files are not identical and the total
        number of keys with differences.

    """

    return_code = 0

    ctl0 = ControlVariables.load(file0).control
    ctl1 = ControlVariables.load(file1).control

    ctl0_keys_only = set(ctl0.keys()) - set(ctl1.keys())
    ctl1_keys_only = set(ctl1.keys()) - set(ctl0.keys())

    if len(ctl0_keys_only):
        return_code = len(ctl0_keys_only)
        if not silent:
            warn(f"Keys only in control file0: {sorted(ctl0_keys_only)}")
    if len(ctl1_keys_only):
        return_code = len(ctl1_keys_only)
        if not silent:
            warn(f"Keys only in control file1: {sorted(ctl1_keys_only)}")

    common_keys = set(ctl0.keys()).intersection(set(ctl1.keys()))

    for ck in common_keys:
        try:
            np.testing.assert_equal(ctl0[ck], ctl1[ck])
        except AssertionError:
            return_code += 1
            if not silent:
                msg = (
                    f"Keys '{ck}' differ\n"
                    f"file0: {ctl0[ck]}\n"
                    f"file1: {ctl1[ck]}\n"
                )
                warn(msg)

    return return_code
