from collections.abc import Mapping


class DictionaryAsProperties(dict):
    """
    Quick and dirty implementation of a dot-able dict, which allows access and
    assignment via object properties rather than dict indexing.
    source: https://stackoverflow.com/questions/13520421/recursive-dotdict
    modified OrderdDict to standard dict
    """

    def __init__(self, *args, **kwargs):
        od = dict(*args, **kwargs)
        for key, val in od.items():
            if isinstance(val, Mapping):
                value = DictionaryAsProperties(val)
            else:
                value = val
            self[key] = value

    def __delattr__(self, name):
        try:
            del self[name]
        except KeyError as ex:
            raise AttributeError(f"No attribute called: {name}") from ex

    def __getattr__(self, k):
        try:
            return self[k]
        except KeyError as ex:
            raise AttributeError(f"No attribute called: {k}") from ex

    __setattr__ = dict.__setitem__
