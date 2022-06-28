from copy import deepcopy
from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from pynhm.atmosphere.AtmBoundaryLayer import AtmBoundaryLayer
from pynhm.atmosphere.NHMBoundaryLayer import NHMBoundaryLayer
from pynhm.atmosphere.NHMSolarGeometry import NHMSolarGeometry
from pynhm.utils.parameters import PrmsParameters
from pynhm.utils.prms5util import load_soltab_debug

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
        params = PrmsParameters.load(domain["param_file"])
        solar_geom = NHMSolarGeometry(params)
        (
            potential_sw_rad_ans,
            potential_sw_rad_flat_ans,
            sun_hrs_ans,
        ) = load_soltab_debug(domain["prms_run_dir"] / "soltab_debug")

        # check the shapes
        assert potential_sw_rad_ans.shape == solar_geom.potential_sw_rad.shape
        assert (
            potential_sw_rad_flat_ans.shape
            == solar_geom.potential_sw_rad.shape
        )
        assert sun_hrs_ans.shape == solar_geom.sun_hrs.shape

        # check the values
        atol = 1e-4
        assert np.isclose(solar_geom.sun_hrs, sun_hrs_ans, atol=atol).all()
        assert np.isclose(
            potential_sw_rad_ans, solar_geom.potential_sw_rad, atol=atol
        ).all()
        assert np.isclose(
            solar_geom.potential_sw_rad_flat,
            potential_sw_rad_flat_ans,
            atol=atol,
        ).all()

        data_size = np.product(potential_sw_rad_ans.shape)
        # for vv in range(3):
        #     if vv == 0:
        #         ans = potential_sw_rad_ans
        #         res = solar_geom.potential_sw_rad
        #         var = "potential_sw_rad"
        #     elif vv == 1:
        #         ans = potential_sw_rad_flat_ans
        #         res = solar_geom.potential_sw_rad_flat
        #         var = "potential_sw_rad_flat"
        #     else:
        #         ans = sun_hrs_ans
        #         res = solar_geom.sun_hrs
        #         var = "sun_hrs"

        #     print(f"\nvariable: {var}")
        #     for exp in range(9):
        #         zz = abs(res - ans)
        #         nn = len(np.where(zz > 10 ** (-1 * exp))[0])

        #         if not (zz < 10 ** (-1 * exp)).all():
        #             if vv == 1:
        #                 asdf
        #             zz_absmax = zz.max()
        #             wh_absmax = np.where(zz == zz_absmax)
        #             wh_absmax
        #             zz[wh_absmax]
        #             ans[wh_absmax]
        #             res[wh_absmax]

        #         max_pct_err = (zz / potential_sw_rad_flat_ans).max()
        #         print(
        #             f"exp: {exp}, nn: {nn}, frac: {nn/data_size}, "
        #             f"max_pct_err: {max_pct_err}"
        #         )

        # zz = potential_sw_rad_ans - solar_geom.potential_sw_rad
        # zz_absmax = abs(zz).max()
        # wh_absmax = np.where(zz == zz_absmax)
        # # (array([48]), array([227]))
        # zz[wh_absmax]

        # debug/examine
        # ans = potential_sw_rad_ans[:, 0]
        # res = solar_geom.potential_sw_rad[:, 0]

        # ans = sun_hrs_ans[:, 0]
        # res = solar_geom.sun_hrs[:, 0]

        return


@pytest.fixture(
    scope="function", params=[None, test_time], ids=["markov", "timeseries"]
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


state_dict = {
    "datetime": np.array(
        [
            datetime(1979, 1, 3, 0, 0),
            datetime(1979, 1, 4, 0, 0),
        ]
    ),
    "spatial_id": np.array([5307, 5308]),
    "tminf": np.array(
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
    "tmaxf": np.array(
        [
            [57.41188049316406, 56.47270965576172],
            [55.511878967285156, 55.032711029052734],
        ]
    ),
    "hru_snow": np.array([[0.0, 0.0], [0.0, 0.0]]),
    "hru_ppt": np.array(
        [
            [0.31392958760261536, 0.24780480563640594],
            [0.6605601906776428, 0.5214226245880127],
        ]
    ),
    "hru_rain": np.array(
        [
            [0.31392958760261536, 0.24780480563640594],
            [0.6605601906776428, 0.5214226245880127],
        ]
    ),
}


@pytest.fixture(scope="function")
def params(domain):
    return PrmsParameters.load(domain["param_file"])


@pytest.fixture(scope="function")
def atm_nhm_init_prms_output(domain):
    vars = [
        "hru_ppt",
        "hru_rain",
        "hru_snow",
        "tmaxf",
        "tminf",
        "swrad",
        "potet",
    ]
    var_file_dict = {
        var: domain["prms_output_dir"] / f"{var}.csv" for var in vars
    }
    atm = NHMBoundaryLayer.load_prms_output(
        **var_file_dict, **atm_init_test_dict
    )
    return atm


@pytest.fixture(scope="function")
def atm_nhm_init_nc(domain, params):
    atm = NHMBoundaryLayer.load_netcdf(
        domain["cbh_nc"], parameters=params, **atm_init_test_dict
    )
    return atm


# JLM: Should try to set a state with time length less than 2.
@pytest.fixture(scope="function")
def atm_nhm_init_dict(domain, params):
    atm = NHMBoundaryLayer(state_dict, parameters=params, **atm_init_test_dict)
    return atm


class TestNHMBoundaryLayer:
    def test_init_dict(self, atm_nhm_init_dict, params):
        # local atm init with an invalid start time
        test_dict = deepcopy(atm_init_test_dict)
        test_dict["start_time"] = np.datetime64(datetime(1979, 1, 12))
        try:
            atm = NHMBoundaryLayer(state_dict, parameters=params, **test_dict)
            assert False
        except ValueError:
            pass

        # local atm init with an invalid start time
        test_dict = deepcopy(atm_init_test_dict)
        current_time_ind = 1
        test_dict["start_time"] = state_dict["datetime"][current_time_ind]
        atm = NHMBoundaryLayer(state_dict, parameters=params, **test_dict)
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

    def test_init_prms_output(self, atm_nhm_init_prms_output):
        # Tests thats the atm_nhm_init_prms_output fixture is properly created.
        assert True
        return

    def test_init_nc(self, domain, atm_nhm_init_nc, params):
        init_test_dict = deepcopy(atm_init_test_dict)
        init_test_dict["nc_file"] = domain["cbh_nc"]
        # Check all the init arguments
        for key, val in init_test_dict.items():
            key_check = key
            if key not in list(atm_nhm_init_nc.__dict__.keys()):
                key_check = f"_{key}"
            assert getattr(atm_nhm_init_nc, key_check) == val
        # Check state - this tests get_state.
        results = {
            var: atm_nhm_init_nc[var].mean()
            for var in atm_nhm_init_nc.variables
        }
        ds = xr.open_dataset(domain["cbh_nc"])
        # preference for adjusted over raw fields in the code
        file_vars = [
            f"{var}_adj" if f"{var}_adj" in ds.variables else var
            for var in atm_nhm_init_nc.variables
        ]
        answers = {var: ds[var].mean().values.tolist() for var in file_vars}
        for key_res in results.keys():
            key_ans = key_res
            if key_ans not in answers:
                key_ans = f"{key_res}_adj"
            assert np.isclose(results[key_res], answers[key_ans])
        return

    def test_state_roundtrip(self, atm_nhm_init_nc):
        # Roundtrip and a half just to use existing values: get, set, get
        mult = 2
        values = atm_nhm_init_nc["prcp"] * mult
        assert (values / mult == atm_nhm_init_nc.prcp).all()
        atm_nhm_init_nc["prcp"] = values
        # One better, these are actually the same object, could test "is"
        assert (atm_nhm_init_nc["prcp"] == values).all()
        return

    def test_current_state(self, atm_nhm_init_nc):
        n_adv = 5
        start_time_index = atm_nhm_init_nc.current_time_index
        answer = deepcopy(atm_nhm_init_nc["prcp"])[start_time_index + n_adv, :]
        for aa in range(n_adv):
            atm_nhm_init_nc.advance()
        assert (atm_nhm_init_nc.get_current_state("prcp") == answer).all()
        return

    @pytest.mark.parametrize(
        "read_vars",
        [["prcp", "tmax", "tmin"], ["notavar"], ["tmin"]],
        ids=["toadj", "invalid", "preadj"],
    )
    def test_nc_read_vars(self, atm_nhm_init_nc, params, read_vars):
        atm_read_dict = deepcopy(atm_init_test_dict)
        atm_read_dict["nc_read_vars"] = read_vars
        if read_vars[0] in ["notavar"]:
            with pytest.raises(KeyError):
                atm_read = NHMBoundaryLayer.load_netcdf(
                    atm_nhm_init_nc._nc_file,
                    parameters=params,
                    **atm_read_dict,
                )
        else:
            atm_read = NHMBoundaryLayer.load_netcdf(
                atm_nhm_init_nc._nc_file, parameters=params, **atm_read_dict
            )
            assert set(atm_read.variables) == set(
                [var.split("_")[0] for var in atm_read_dict["nc_read_vars"]]
            )
        return

    def test_state_adj(self, domain, atm_nhm_init_nc, params):
        atm_adj_dict = deepcopy(atm_init_test_dict)
        atm_adj_dict["nc_read_vars"] = ["prcp", "tmax", "tmin"]
        atm_to_adj = NHMBoundaryLayer.load_netcdf(
            atm_nhm_init_nc._nc_file, parameters=params, **atm_adj_dict
        )

        atm_to_adj.param_adjust()
        for vv in atm_to_adj.variables:
            # JLM: we loose rhavg in the adjustments... it's a bit of a mystery
            # if that should also be adjusted. checking on that w parker
            if vv in atm_nhm_init_nc.variables:
                assert np.isclose(
                    atm_to_adj[vv],
                    atm_nhm_init_nc[vv],
                    atol=1e-06,
                ).all()

        return

    def test_calculations_vs_prms_output(
        self, domain, atm_nhm_init_nc, atm_nhm_init_prms_output
    ):
        atm_nhm_init_nc.calculate_sw_rad_degree_day()
        atm_nhm_init_nc.calculate_potential_et_jh()

        var_tol = {
            "hru_ppt": 1e-5,
            "hru_rain": 1e-5,
            "hru_snow": 1e-5,
            "tmaxf": 1e-4,
            "tminf": 1e-4,
            "swrad": 1e-3,
            "potet": 1e-5,
        }

        for var, tol in var_tol.items():
            result = np.isclose(
                atm_nhm_init_nc[var],
                atm_nhm_init_prms_output[var],
                rtol=1e-121,
                atol=tol,
            )  # Only the atol matters here

            # dd = np.abs(
            #     atm_nhm_init_prms_output[var] - atm_nhm_init_nc[var]
            # ).max()
            # print(dd, var)

            assert result.all(), var

        var = "hru_ppt"
        for tt in range(3):
            result = np.isclose(
                atm_nhm_init_nc.get_current_state(var),
                atm_nhm_init_prms_output.get_current_state(var),
                rtol=1e-121,
                atol=var_tol[var],
            )  # Only the atol matters here
            assert result.all()

            atm_nhm_init_nc.advance()
            result = np.isclose(
                atm_nhm_init_nc.get_current_state(var),
                atm_nhm_init_prms_output.get_current_state(var),
                rtol=1e-121,
                atol=var_tol[var],
            )  # Only the atol matters here
            assert not result.all()

            atm_nhm_init_prms_output.advance()

        return
