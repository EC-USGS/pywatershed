import numpy as np

from ..atmosphere.NHMBoundaryLayer import NHMBoundaryLayer
from ..base.storageUnit import StorageUnit
from ..utils.parameters import PrmsParameters


class PRMSGroundwater(StorageUnit):
    def __init__(
        self,
        params: PrmsParameters,
        atm: NHMBoundaryLayer,
    ) -> "PRMSGroundwater":

        verbose = True
        super().__init__("gwflow", id, params, atm, verbose)

        # define self variables that will be used for the calculation
        self.gw_stor = self.gwstor_init.copy()
        self.gw_stor_old = self.gwstor_init.copy()

        for name in PRMSGroundwater.get_input_variables():
            setattr(self, name, np.zeros(self.nhru, dtype=float))

        # self.soil_to_gw = np.zeros(self.nhru, dtype=float)
        # self.ssr_to_gw = np.zeros(self.nhru, dtype=float)
        # self.dprst_seep = np.zeros(self.nhru, dtype=float)
        # self.gwres_in = np.zeros(self.nhru, dtype=float)

        # define information on the output data that will be created
        # output_variables = (
        #     "gwres_flow",
        #     "gwres_in",
        #     "gwres_sink",
        #     "gwres_stor",
        # )
        self.output_column_names = ["date"]
        self.output_data = []
        for name in PRMSGroundwater.get_output_variables():
            setattr(self, name, np.zeros(self.nhru, dtype=float))
            if "nhm_id" in params.parameters:
                self.output_column_names += [
                    f"nhru_{name}_{nhmid}"
                    for nhmid in params.parameters["nhm_id"]
                ]
            else:
                self.output_column_names += [
                    f"{name}_{id}" for id in range(self.nhru)
                ]
        return

    @staticmethod
    def get_required_parameters() -> tuple:
        """
        Return a tuple of the parameters required for this process

        """
        return (
            "nhru",
            "ngw",
            "hru_area",
            "gwflow_coef",
            "gwsink_coef",
            "gwstor_init",
            "gwstor_min",
        )

    @staticmethod
    def get_input_variables() -> tuple:
        """

        Returns:

        """
        return (
            "soil_to_gw",
            "ssr_to_gw",
            "dprst_seep",
        )

    @staticmethod
    def get_output_variables() -> tuple:
        return (
            "gwres_flow",
            "gwres_in",
            "gwres_sink",
            "gwres_stor",
        )

    def advance(self, itime_step):
        self.gw_stor_old = self.gw_stor

        return

    def calculate(self, time_length):

        gwarea = self.hru_area

        # calculate volume terms
        gwstor = self.gw_stor * gwarea
        soil_to_gw_vol = self.soil_to_gw * gwarea
        ssr_to_gw_vol = self.ssr_to_gw * gwarea
        dprst_seep_vol = self.dprst_seep * gwarea

        # initialize calculation variables
        gwres_in = soil_to_gw_vol + ssr_to_gw_vol + dprst_seep_vol

        # todo: what about route order

        gwstor += gwres_in

        gwflow = gwstor * self.gwflow_coef

        gwstor -= gwflow

        # output variables
        self.gw_stor = gwstor / gwarea
        self.gwres_in = gwres_in / gwarea
        self.gw_flow = gwflow / gwarea

        return
