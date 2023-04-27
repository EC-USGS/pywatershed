import pytest

from pywatershed.base.accessor import Accessor


class TestAccessor:
    def test_accessor(self):
        da = Accessor()
        assert isinstance(da, Accessor)
        key, value = ("a", 0)
        da[key] = value
        assert da[key] is value
        del da[key]
        with pytest.raises(AttributeError):
            da[key]
        return
