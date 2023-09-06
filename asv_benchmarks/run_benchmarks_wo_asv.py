from itertools import product
from benchmarks.prms import domains, outputs, model_tests, PRMSModels

# Just running the PRMSModels benchmarks

all_tests = list(
    product(
        *[
            domains,
            model_tests.values(),
            outputs,
        ]
    )
)


all_tests = all_tests

for args in all_tests:
    print("\n======================")
    print(args)
    mm = PRMSModels()
    mm.setup(*args)
    mm.time_prms_run(*args)
    mm.teardown(*args)
