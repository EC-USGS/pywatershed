import pynhm as pws
from . import parameterized, test_data_dir


class Parameters:
    """Benchmark Control object of pywatershed"""

    def setup(self, *args, **kwargs):
        self.Control = pws.Control
        return

    @parameterized(
        ["domain"],
        (["drb_2yr", "ucb_2yr"]),
    )
    def time_control_read(self, domain):
        control_file = test_data_dir / f"{domain}/control.test"
        _ = self.Control.load(control_file)
        return
