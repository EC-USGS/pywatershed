import pathlib as pl

from pynhm.base.control import Control
from pynhm.hydrology.PRMSChannel import PRMSChannel
from pynhm.utils.netcdf_utils import NetCdfCompare
from pynhm.utils.parameters import PRMSParameters


class TestPRMSChannelDomain:
    def test_init(self, domain, tmp_path):
        tmp_path = pl.Path(tmp_path)
        params = PRMSParameters.load(domain["param_file"])

        # Set information from the control file
        control = Control.load(domain["control_file"], params=params)

        # load csv files into dataframes
        output_dir = domain["prms_output_dir"]
        input_variables = {}
        for key in PRMSChannel.get_inputs():
            nc_path = output_dir / f"{key}.nc"
            input_variables[key] = nc_path

        channel = PRMSChannel(control, **input_variables, budget_type="error")
        nc_parent = tmp_path / domain["domain_name"]
        channel.initialize_netcdf(nc_parent)

        for istep in range(control.n_times):
            control.advance()

            channel.advance()

            channel.calculate(float(istep))

            channel.output()

        channel.finalize()

        output_compare = {}
        for key in PRMSChannel.get_variables():
            base_nc_path = output_dir / f"{key}.nc"
            compare_nc_path = tmp_path / domain["domain_name"] / f"{key}.nc"
            # PRMS does not output the storage change in the channel
            if not base_nc_path.exists():
                continue
            output_compare[key] = (base_nc_path, compare_nc_path)

        print(f"base_nc_path: {base_nc_path}")
        print(f"compare_nc_path: {compare_nc_path}")

        assert_error = False
        for key, (base, compare) in output_compare.items():
            success, diff = NetCdfCompare(base, compare).compare()
            if not success:
                print(
                    f"comparison for {key} failed: "
                    + f"maximum error {diff[key][0]} "
                    + f"(maximum allowed error {diff[key][1]}) "
                    + f"in column {diff[key][2]}"
                )
                assert_error = True
            else:
                print(f"comparison for {key} passed")

        assert not assert_error, "comparison failed"

        return
