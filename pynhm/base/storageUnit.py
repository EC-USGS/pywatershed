from ..atmosphere.NHMBoundaryLayer import NHMBoundaryLayer
from ..utils.parameters import PrmsParameters


class StorageUnit:
    @staticmethod
    def get_required_parameters():
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

        for param in self.get_required_parameters():
            setattr(self, param, None)

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
        self.residual_old = None
        self.residual_new = None
        self.inflow_volumes = None
        self.outflow_volumes = None
        self.recipients = []
        self.depencencies = []
        self.output_data = []
        self.output_column_names = []
        # self.advance(0)
        return

    def register_recipient(self, recipient, process_name):
        recipient_info = (recipient, process_name)
        self.recipients.append(recipient_info)
        return

    def register_dependency(self, dependency, process_name):
        dependency_info = (dependency, process_name)
        self.depencencies.append(dependency_info)

    def advance(self, itime_step):
        self.residual_old = self.residual_new
        self.inflow_volumes = {}
        self.outflow_volumes = {}
        return

    def add_water(self, inflow_name, inflow_rate, time_length):
        self.inflow_volumes[inflow_name] = inflow_rate * time_length
        return

    def remove_water(self, outflow_name, outflow_rate, time_length):
        self.outflow_volumes[outflow_name] = outflow_rate * time_length
        return

    def calculate(self, time_length):
        residual_new = 0.0
        for inflow_name in self.inflow_volumes:
            vi = self.inflow_volumes[inflow_name]
            residual_new += vi
        for outflow_name in self.outflow_volumes:
            vo = self.outflow_volumes[outflow_name]
            residual_new += vo
        self.residual_new = residual_new
        return

    def get_budget_summary_str(self):
        s = self.get_name()
        s += (
            f" inflow {self.inflow_volumes} outflow {self.outflow_volumes} "
            f"residual {self.residual_new} "
        )
        return s

    def get_name(self):
        return f"{self.storage_type}{self.id}"

    def finalize(self):
        return
