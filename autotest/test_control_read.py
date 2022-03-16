import pytest

from pynhm.utils import get_all_variables
from utils import assert_or_print


@pytest.fixture
def control_keys():
    return tuple(("start_time", "end_time", "initial_deltat"))


def test_control_read(domain, control_keys):
    control_file = domain["control_file"]
    print(f"parsing...'{control_file}'")

    variables = get_all_variables(control_file)

    # check dimensions
    answers = domain["test_ans"]["control_read"]
