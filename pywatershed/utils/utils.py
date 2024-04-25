import functools
from time import time

import numpy as np


def timer(func):
    """Use as a decorator to print the execution time of the passed function"""

    @functools.wraps(func)
    def wrap_func(*args, **kwargs):
        t1 = time()
        result = func(*args, **kwargs)
        t2 = time()
        print(f"Function {func.__name__!r} executed in {(t2-t1):.4f}s")
        return result

    return wrap_func


def diff_dicts(dict_a: dict, dict_b: dict, ignore_keys: list = []):
    """Diff two dictionaries

    Args:
        dict_a: the first dictionary
        dict_b: the second dictionary
    """
    keys_a = dict_a.keys()
    keys_b = dict_b.keys()
    if keys_a_not_b := keys_a - keys_b:
        print(f"keys in a but not b: {keys_a_not_b}")
        print("")

    if keys_b_not_a := keys_b - keys_a:
        print(f"keys in b but not a: {keys_b_not_a}")
        print("")

    keys_both = keys_a & keys_b
    for kk in keys_both:
        if kk in ignore_keys:
            continue
        val_a = dict_a[kk]
        val_b = dict_b[kk]
        try:
            np.testing.assert_equal(val_a, val_b)
        except AssertionError:
            print(f'key "{kk}" does not match')
            print("value for a: ")
            print(f"    {val_a}")
            print("value for b: ")
            print(f"    {val_b}")
            print("")
