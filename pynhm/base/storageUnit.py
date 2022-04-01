import os
import pandas as pd

from ..atmosphere.NHMBoundaryLayer import NHMBoundaryLayer
from ..utils.parameters import PrmsParameters


class StorageUnit:
    def __init__(
        self,
        storage_type,
        id: list,
        params: PrmsParameters,
        atm: NHMBoundaryLayer,
        verbose: bool,
    ):

        self.storage_type = storage_type
        self.id = id
        self.params = params
        self.atm = atm
        self.verbose = verbose

        # Go through list of parameters for this process and assign them
        # to self with a value of None
        for param in self.get_required_parameters():
            setattr(self, param, None)

        # Go through the parameters for this process and see if self
        # has a variable with that name.  If so, then assign the parameter
        # value to the self variable.
        for key in params.parameters:
            if hasattr(self, key):
                setattr(self, key, params.parameters[key])

        # if any of the required parameters are still none,
        # then we should terminate with an error
        for key in self.get_required_parameters():
            value = getattr(self, key)
            if value is None:
                print(
                    f"{storage_type} storage unit requires {key} but it was not found in parameters."
                )

        return

    @staticmethod
    def get_required_parameters() -> list:
        raise Exception("This must be overridden")

    def initialize_output_data(self):
        self.output_column_names = ["date"]
        if "nhm_id" in self.params.parameters:
            self.output_column_names += [
                f"nhru_hru_intcpstor_{nhmid}"
                for nhmid in self.params.parameters["nhm_id"]
            ]
        else:
            self.output_column_names += [
                f"intcpstor_{id}" for id in range(self.nhru)
            ]

        # Initialize the output_data dictionary for each entry with an empty list
        self.output_data = {}
        for output_name in self.output_data_names:
            self.output_data[output_name] = []
        return

    def get_output_dataframes(self):
        """
        Return a dictionary of data frames of output data from this process

        """
        output_data = {}
        for key in self.output_data:
            value = self.output_data[key]
            df = pd.DataFrame(value, columns=self.output_column_names)
            df.set_index("date", inplace=True)
            output_data[key] = df
        return output_data

    def output_to_csv(self, pth):
        """
        Save each output variable to separate csv file in specified path

        """
        output_data = self.get_output_dataframes()
        for key in output_data:
            df = output_data[key]
            fname = os.path.join(pth, f"{key}.csv")
            df.to_csv(fname)
        return