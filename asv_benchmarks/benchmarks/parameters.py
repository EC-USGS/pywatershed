import pynhm as pws

from . import parameterized, test_data_dir


class Parameters:
    """Benchmark Parameters of pywatershed"""

    def setup(self, *args, **kwargs):
        self.PrmsParameters = pws.PrmsParameters
        return

    @parameterized(
        ["domain"],
        (["drb_2yr", "ucb_2yr"]),
    )
    def time_parameter_read_prms(self, domain):
        parameter_file = test_data_dir / f"{domain}/myparam.param"
        _ = self.PrmsParameters.load(parameter_file)
        return
