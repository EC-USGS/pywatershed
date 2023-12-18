import functools
from time import time


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
