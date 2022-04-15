import pathlib as pl

from pynhm.atmosphere.NHMBoundaryLayer import NHMBoundaryLayer
from pynhm.base.control import Control
from pynhm.hydrology.PRMSGroundwater import PRMSGroundwater
from pynhm.utils.netcdf_utils import NetCdfCompare
from pynhm.utils.parameters import PrmsParameters


class TestPRMSGroundwaterDomain:
    def test_init(self, domain, tmp_path):
        tmp_path = pl.Path(tmp_path)
        prms_params = PrmsParameters.load(domain["param_file"])

        # Set information from the control file
        control = Control.load(domain["control_file"])

        # load csv files into dataframes
        output_files = domain["prms_outputs"]
        input_variables = {}
        for key in PRMSGroundwater.get_input_variables():
            output_path = output_files[key]
            nc_path = output_path.with_suffix(".nc")
            input_variables[key] = nc_path

        gw = PRMSGroundwater(control, prms_params, **input_variables)
        nc_parent = tmp_path / domain["domain_name"]
        gw.initialize_netcdf(nc_parent)

        output_compare = {}
        for key in PRMSGroundwater.get_variables():
            output_path = output_files[key]
            base_nc_path = output_path.with_suffix(".nc")
            compare_nc_path = tmp_path / domain["domain_name"] / f"{key}.nc"
            output_compare[key] = (base_nc_path, compare_nc_path)

        print(f"base_nc_path: {base_nc_path}")
        print(f"compare_nc_path: {compare_nc_path}")

        for istep in range(control.n_times):

            gw.advance()

            gw.calculate(float(istep))

            gw.output()

        gw.finalize()

        assert_error = False
        for key, (base, compare) in output_compare.items():
            success, diff = NetCdfCompare(base, compare).compare()
            if not success:
                print(f"comparison for {key} failed: maximum error {diff}")
                assert_error = True
        assert not assert_error, "comparison failed"

        # # create data frame of prms interception storage (from nhru_hru_intcpstor.csv)
        # output_files = domain["prms_outputs"]
        # fname = output_files["intcpstor"]
        # print(f"loading {fname}")
        # intcpstor_prms = pynhm.preprocess.CsvFile(fname)
        # intcpstor_prms = intcpstor_prms.to_dataframe()
        #
        # # create data frame of pynhm interception storage
        # intcpstor_pynhm = pd.DataFrame(
        #     self.cnp.output_data, columns=self.cnp.output_column_names
        # )
        # intcpstor_pynhm.set_index("date", inplace=True)
        #
        # print("intcp_stor  min  max")
        # a1 = intcpstor_prms.to_numpy()
        # a2 = intcpstor_pynhm.to_numpy()
        # print(f"prms  {a1.min()}    {a1.max()}")
        # print(f"pynhm  {a2.min()}    {a2.max()}")
        #
        # makeplot = False
        # if makeplot:
        #     import matplotlib.pyplot as plt
        #
        #     f = plt.figure()
        #     ax = plt.subplot(1, 1, 1, aspect="equal")
        #     xmin = 1e30
        #     xmax = -1e30
        #     for colname in intcpstor_prms.columns:
        #         x1 = intcpstor_prms.loc[:, colname].values
        #         x2 = intcpstor_pynhm.loc[:, colname].values
        #         ax.scatter(
        #             x1, x2, facecolors="none", edgecolors="k", linewidth=0.25
        #         )
        #         xmin = min(xmin, x1.min(), x2.min())
        #         xmax = max(xmax, x1.max(), x2.max())
        #     ax.set_title("Interception Storage")
        #     ax.set_xlabel("PRMS Interception Storage, in inches")
        #     ax.set_ylabel("PYNHM Interception Storage, in inches")
        #     ax.set_xlim(xmin, xmax)
        #     ax.set_ylim(xmin, xmax)
        #     pname = domain["domain_name"] + "_1-1.png"
        #     print(f"Creating plot {pname}")
        #     plt.savefig(pname, dpi=300)
        #
        #     f = plt.figure()
        #     ax = plt.subplot(1, 1, 1)
        #     pc = plt.imshow(a1 - a2, cmap="jet")
        #     plt.colorbar(pc, shrink=0.5)
        #     ax.set_title("Difference Between PRMS and PYNHM")
        #     ax.set_xlabel("HRU number")
        #     ax.set_ylabel("Time step")
        #     pname = domain["domain_name"] + "_2d.png"
        #     print(f"Creating plot {pname}")
        #     plt.savefig(pname, dpi=300)

        return
