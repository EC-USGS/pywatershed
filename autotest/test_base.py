from datetime import datetime, timedelta

import numpy as np
import pytest

from pynhm.base.accessor import Accessor
from pynhm.base.control import Control
from pynhm.base.storageUnit import StorageUnit
from pynhm.base.Time import Time
from pynhm.utils.parameters import PrmsParameters


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
