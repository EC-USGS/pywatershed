import numpy as np

print_ans = False


def assert_or_print(results, answers, test_name=None, print_ans=False):
    print("\n")
    if test_name is not None:
        test_name = f"{test_name}:"
    n_space = 4
    if print_ans and (test_name is not None):
        sp = "".join(n_space * [" "])
        n_space = n_space + 2
        print(f"{sp}{test_name}")

    # Always check every test in the answers
    for key in answers.keys():
        if print_ans:
            sp = "".join(n_space * [" "])
            print(f"{sp}{key}: {results[key]}")
        else:
            msg = f"{test_name}{key}"
            assert np.isclose(results[key], answers[key]), msg

    # Always fail if printing answers
    assert not print_ans
    return
