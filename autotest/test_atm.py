import pathlib as pl
from copy import deepcopy
from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from pynhm.atmosphere.AtmBoundaryLayer import AtmBoundaryLayer
from pynhm.atmosphere.NHMBoundaryLayer import NHMBoundaryLayer
from pynhm.atmosphere.NHMSolarGeometry import NHMSolarGeometry
from pynhm.utils.parameters import PrmsParameters
from pynhm.utils.prms5util import load_prms_statscsv, load_soltab_debug

test_time = np.arange(
    datetime(1979, 1, 1), datetime(1979, 1, 7), timedelta(days=1)
).astype(np.datetime64)

atm_init_test_dict = {
    "start_time": np.datetime64("1979-01-03T00:00:00.00"),
    "time_step": np.timedelta64(1, "D"),
    "verbosity": 3,
    "height_m": 5,
}

test_time_steps = [atm_init_test_dict["time_step"], np.timedelta64(1, "h")]


class TestNHMSolarGeometry:
    def test_init(self, domain):
        params = PrmsParameters(domain["param_file"])
        solar_geom = NHMSolarGeometry(params)
        potential_sw_rad_ans, sun_hrs_ans = load_soltab_debug(
            domain["prms_outputs"]["soltab"]
        )

        # check the shapes
        assert potential_sw_rad_ans.shape == solar_geom.potential_sw_rad.shape
        assert sun_hrs_ans.shape == solar_geom.sun_hrs.shape

        # check the values
        assert np.isclose(
            potential_sw_rad_ans, solar_geom.potential_sw_rad, atol=1e-04
        ).all()
        assert np.isclose(sun_hrs_ans, solar_geom.sun_hrs, atol=1e-03).all()

        # debug/examine
        # ans = potential_sw_rad_ans[:, 0]
        # res = solar_geom.potential_sw_rad[:, 0]

        # ans = sun_hrs_ans[:, 0]
        # res = solar_geom.sun_hrs[:, 0]

        return


@pytest.fixture(
    scope="function", params=[None, test_time], ids=["markov", "datetime"]
)
def atm_init(request):
    init_dict = deepcopy(atm_init_test_dict)
    init_dict["datetime"] = request.param
    atm = AtmBoundaryLayer(**init_dict)
    return atm


# There really is no state, just testing datetime here
class TestAtmBoundaryLayer:
    def test_init(self, atm_init):
        for key, val in atm_init_test_dict.items():
            key_check = key
            if key not in list(atm_init.__dict__.keys()):
                key_check = f"_{key}"
            assert getattr(atm_init, key_check) == val
        return

    @pytest.mark.parametrize(
        "the_state", ["datetime", "foo"], ids=["valid", "invalid"]
    )
    def test_state_roundtrip(self, atm_init, the_state):
        # Has to be a round trip because there's no defined time or state
        # get and set state via __setitem__ and __getitem__
        # Could test other errors of set and get
        try:
            if the_state == "datetime":
                atm_init.datetime = test_time  # skip setitem
            else:
                atm_init[the_state] = test_time
            if the_state in ["foo"]:
                assert False
        except KeyError:
            if the_state in ["foo"]:
                assert True
            else:
                assert False

        try:
            assert (atm_init[the_state] == test_time).all()
            if the_state in ["foo"]:
                assert False
        except KeyError:
            if the_state in ["foo"]:
                assert True
            else:
                assert False

        return

    # Most time methods tested in Time, except current_state

    def test_current_state_roundtrip(self, atm_init):
        # Must be a round trip since there's no defined time/state
        n_adv = 5
        for aa in range(n_adv):
            atm_init.advance()

        assert atm_init.current_time == atm_init._start_time + (
            n_adv * atm_init.time_step
        )
        if atm_init["datetime"] is None:
            assert atm_init.current_time_index == 1
        else:
            assert atm_init.current_time_index == n_adv

        return


@pytest.fixture(scope="function")
def atm_nhm_init(domain):
    atm = NHMBoundaryLayer(domain["cbh_nc"], **atm_init_test_dict)
    return atm


state_dict = {
    "datetime": np.array(
        [
            datetime(1979, 1, 3, 0, 0),
            datetime(1979, 1, 4, 0, 0),
        ]
    ),
    "spatial_id": np.array([5307, 5308]),
    "tmin": np.array(
        [
            [46.1209, 45.76805],
            [37.7609, 37.4881],
        ]
    ),
    "rhavg": np.array(
        [
            [82.45999908447266, 82.5999984741211],
            [81.98999786376953, 82.3499984741211],
        ]
    ),
    "tmax": np.array(
        [
            [57.41188049316406, 56.47270965576172],
            [55.511878967285156, 55.032711029052734],
        ]
    ),
    "snowfall": np.array([[0.0, 0.0], [0.0, 0.0]]),
    "prcp": np.array(
        [
            [0.31392958760261536, 0.24780480563640594],
            [0.6605601906776428, 0.5214226245880127],
        ]
    ),
    "rainfall": np.array(
        [
            [0.31392958760261536, 0.24780480563640594],
            [0.6605601906776428, 0.5214226245880127],
        ]
    ),
}


# JLM: Should try to set a state with time length less than 2.
@pytest.fixture(scope="function")
def atm_nhm_init_dict():
    atm = NHMBoundaryLayer(state_dict, **atm_init_test_dict)
    return atm


class TestNHMBoundaryLayer:
    def test_init_dict(self, atm_nhm_init_dict):

        # local atm init with an invalid start time
        test_dict = deepcopy(atm_init_test_dict)
        test_dict["start_time"] = np.datetime64(datetime(1979, 1, 12))
        try:
            atm = NHMBoundaryLayer(state_dict, **test_dict)
            assert False
        except ValueError:
            pass

        # local atm init with an invalid start time
        test_dict = deepcopy(atm_init_test_dict)
        current_time_ind = 1
        test_dict["start_time"] = state_dict["datetime"][current_time_ind]
        atm = NHMBoundaryLayer(state_dict, **test_dict)
        assert atm.current_time_index == current_time_ind
        assert atm.current_time == test_dict["start_time"]
        try:
            atm.advance()
            assert False
        except ValueError:
            pass

        # Check all the init arguments
        # back to atm_nhm_init_test_dict
        current_time_ind = 0
        assert atm_nhm_init_dict.current_time_index == current_time_ind
        assert (
            atm_nhm_init_dict.current_time == atm_init_test_dict["start_time"]
        )

        for key, val in atm_init_test_dict.items():
            key_check = key
            if key not in list(atm_nhm_init_dict.__dict__.keys()):
                key_check = f"_{key}"
            assert getattr(atm_nhm_init_dict, key_check) == val
        # Check state - this tests get_state.
        for vv in atm_nhm_init_dict.variables:
            assert (atm_nhm_init_dict[vv] == state_dict[vv]).all()

        return

    def test_init(self, domain, atm_nhm_init):
        init_test_dict = deepcopy(atm_init_test_dict)
        init_test_dict["nc_file"] = domain["cbh_nc"]
        # Check all the init arguments
        for key, val in init_test_dict.items():
            key_check = key
            if key not in list(atm_nhm_init.__dict__.keys()):
                key_check = f"_{key}"
            assert getattr(atm_nhm_init, key_check) == val
        # Check state - this tests get_state.
        results = {
            var: atm_nhm_init[var].mean() for var in atm_nhm_init.variables
        }
        ds = xr.open_dataset(domain["cbh_nc"])
        # preference for adjusted over raw fields in the code
        file_vars = [
            f"{var}_adj" if f"{var}_adj" in ds.variables else var
            for var in atm_nhm_init.variables
        ]
        answers = {var: ds[var].mean().values.tolist() for var in file_vars}
        for key_res in results.keys():
            key_ans = key_res
            if key_ans not in answers:
                key_ans = f"{key_res}_adj"
            assert np.isclose(results[key_res], answers[key_ans])
        return

    def test_state_roundtrip(self, atm_nhm_init):
        # Roundtrip and a half just to use existing values: get, set, get
        mult = 2
        values = atm_nhm_init["prcp"] * mult
        assert (values / atm_nhm_init.prcp == mult).all()  # while we are here
        atm_nhm_init["prcp"] = values
        # One better, these are actually the same object, could test "is"
        assert (atm_nhm_init["prcp"] == values).all()
        return

    def test_current_state(self, atm_nhm_init):
        n_adv = 5
        start_time_index = atm_nhm_init.current_time_index
        answer = deepcopy(atm_nhm_init["prcp"])[start_time_index + n_adv, :]
        for aa in range(n_adv):
            atm_nhm_init.advance()
        assert (atm_nhm_init.get_current_state("prcp") == answer).all()
        return

    @pytest.mark.parametrize(
        "read_vars",
        [["prcp", "tmax", "tmin"], ["foo"], ["tmin_adj"]],
        ids=["toadj", "invalid", "preadj"],
    )
    def test_nc_read_vars(self, atm_nhm_init, read_vars):
        atm_read_dict = deepcopy(atm_init_test_dict)
        atm_read_dict["nc_read_vars"] = read_vars
        try:
            atm_read = NHMBoundaryLayer(atm_nhm_init._nc_file, **atm_read_dict)
            if read_vars[0] in ["foo"]:
                assert False
            else:
                assert set(atm_read.variables) == set(
                    [
                        var.split("_")[0]
                        for var in atm_read_dict["nc_read_vars"]
                    ]
                )
        except ValueError:
            if read_vars[0] in ["foo"]:
                assert True
            else:
                assert False
        return

    # JLM: def test_nc_read_vars

    def test_state_adj(self, domain, atm_nhm_init):
        atm_adj_dict = deepcopy(atm_init_test_dict)
        atm_adj_dict["nc_read_vars"] = ["prcp", "tmax", "tmin"]
        atm_to_adj = NHMBoundaryLayer(atm_nhm_init._nc_file, **atm_adj_dict)
        params = PrmsParameters(domain["param_file"])
        atm_to_adj.param_adjust(params)

        for vv in atm_to_adj.variables:
            # JLM: we loose rhavg in the adjustments... it's a bit of a mystery
            # if that should also be adjusted. checking on that w parker
            if vv in atm_nhm_init.variables:
                assert np.isclose(
                    atm_to_adj[vv],
                    atm_nhm_init[vv],
                    atol=1e-06,
                ).all()

        return

    # def test_solar_rad_deg_day(self, domain, atm_nhm_init):
    #     params = PrmsParameters(domain["param_file"])
    #     # assert output matches output on file.
    #     prms_output_file = domain["prms_outputs"]["swrad"]
    #     prms_output = load_prms_statscsv(prms_output_file)
    #     swrad_ans_dates = prms_output.index.values
    #     swrad_ans_array = prms_output.to_numpy()
    #     wh_dates = np.where(np.isin(atm_nhm_init["datetime"], swrad_ans_dates))

    #     swrad = atm_nhm_init.calculate_sw_rad_degree_day(params)

    #     result = np.isclose(
    #         swrad_ans_array,
    #         swrad[wh_dates, :],
    #         rtol=1e-121,
    #         atol=1e-04,  # Only the atol matters here, if atol < 1e-4 fails
    #     )
    #     assert result.all()

    #     return
