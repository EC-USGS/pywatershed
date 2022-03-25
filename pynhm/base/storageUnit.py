from ..atmosphere.NHMBoundaryLayer import NHMBoundaryLayer
from ..utils.parameters import PrmsParameters


class StorageUnit:
    @staticmethod
    def get_required_parameters() -> list:
        raise Exception("This must be overridden")

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

        # self.area = area
        self.atm = atm
        self.verbose = verbose
        return
