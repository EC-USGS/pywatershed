import pathlib as pl

from pynhm.base.control import Control
from pynhm.hydrology.PRMSChannel import PRMSChannel
from pynhm.utils.netcdf_utils import NetCdfCompare
from pynhm.utils.parameters import PrmsParameters


class TestPRMSChannelDomain:
    def test_init(self, domain, tmp_path):
        tmp_path = pl.Path(tmp_path)
        prms_params = PrmsParameters.load(domain["param_file"])

        # Set information from the control file
        control = Control.load(domain["control_file"])

        # load csv files into dataframes
        output_dir = domain["prms_output_dir"]
        input_variables = {}
        for key in PRMSChannel.get_inputs():
            nc_path = output_dir / f"{key}.nc"
            input_variables[key] = nc_path

        channel = PRMSChannel(control, prms_params, **input_variables)
        nc_parent = tmp_path / domain["domain_name"]
        channel.initialize_netcdf(nc_parent)

        output_compare = {}
        for key in PRMSChannel.get_variables():
            base_nc_path = output_dir / f"{key}.nc"
            compare_nc_path = tmp_path / domain["domain_name"] / f"{key}.nc"
            output_compare[key] = (base_nc_path, compare_nc_path)

        print(f"base_nc_path: {base_nc_path}")
        print(f"compare_nc_path: {compare_nc_path}")

        for istep in range(control.n_times):
            control.advance()

            channel.advance()

            channel.calculate(float(istep))

            channel.output()

        channel.finalize()

        assert_error = False
        for key, (base, compare) in output_compare.items():
            success, diff = NetCdfCompare(base, compare).compare()
            if not success:
                print(f"comparison for {key} failed: maximum error {diff}")
                assert_error = True
        assert not assert_error, "comparison failed"

        return
