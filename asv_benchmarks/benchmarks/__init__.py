import pywatershed as pws

pws_root = pws.constants.__pywatershed_root__
test_data_dir = pws_root / "../test_data"


def parameterized(names, params):
    def decorator(func):
        func.param_names = names
        func.params = params
        return func

    return decorator
