import numpy as np
import pytest
from pynhm.base.DataAccess import DataAccess


class TestDataAccess:
    def test_init(self):
        da = DataAccess()
        assert da.name == "DataAccess"
        assert len(da.coords) == 0
        assert len(da.variables) == 0
        assert len(da._potential_variables) == 0
        return

    @pytest.mark.parametrize(
        "data", [np.arange(4), list(range(4))], ids=["valid", "invalid"]
    )
    def test_coord(self, data):
        da = DataAccess()
        da._coords = ["foo"]  # strictly verboten!
        try:
            da["foo"] = data
            assert isinstance(data, np.ndarray)
            assert np.isclose(da["foo"], data).all()
            # can not delete a coord
        except TypeError:
            assert not isinstance(data, np.ndarray)
            return

        # can not set or delete coords directly
        try:
            da.coords = ["foo"]
            assert False
        except AttributeError:
            assert True

        try:
            da.coords["foo"] = data
            assert False
        except TypeError:
            assert True

        # can not delete a coord
        try:
            del da["foo"]
            assert False
        except KeyError:
            assert True

        return

    @pytest.mark.parametrize(
        "data", [np.arange(4), list(range(4))], ids=["valid", "invalid"]
    )
    def test_variable(self, data):
        da = DataAccess()
        da._potential_variables = ["foo"]  # strictly verboten!
        try:
            da["foo"] = data
            assert "foo" in da.variables
            assert isinstance(data, np.ndarray)
            assert np.isclose(da["foo"], data).all()
            del da["foo"]
            assert len(da.variables) == 0
        except TypeError:
            assert not isinstance(data, np.ndarray)
            return

        # can not set or delete coords directly
        try:
            da.variables = ["foo"]
            assert False
        except AttributeError:
            assert True

        return
