import numpy as np
import pytest

from pynhm.base.budget import Budget


def test_budget():
    nhru = 5
    inputs = {
        "in1": np.ones([nhru]),
        "in2": np.ones([nhru]),
    }
    outputs = {
        "out1": np.zeros([nhru]),
        "out2": np.ones([nhru]),
    }
    storage_changes = {
        "stor1": np.ones([nhru]),
        "stor2": np.zeros([nhru]),
    }

    argz = {
        "inputs": inputs,
        "outputs": outputs,
        "storage_changes": storage_changes,
    }

    # budget = Budget(**argz)
    argz2 = {key: list(val.keys()) for key, val in argz.items()}

    budget = Budget(**argz2)
    budget.set(argz)

    assert isinstance(budget, Budget)

    for component in budget.components:
        assert budget[component] == argz[component]

    budget.calculate()

    assert (budget._inputs_sum == 2).all()
    assert (budget._outputs_sum == 1).all()
    assert (budget._storage_changes_sum == 1).all()

    for component in budget.components:
        for var in budget[component].keys():
            assert (
                budget.accumulations[component][var]
                == argz[component][var] * 1
            ).all()

    budget.calculate()

    assert (budget._inputs_sum == 2).all()
    assert (budget._outputs_sum == 1).all()
    assert (budget._storage_changes_sum == 1).all()
    assert (budget.balance == 1).all()

    for component in budget.components:
        for var in budget[component].keys():
            assert (
                budget.accumulations[component][var]
                == argz[component][var] * 2
            ).all()

    # change inputs
    inputs["in1"][:] = -1 * np.ones([nhru])
    storage_changes["stor2"][:] = -2 * np.ones([nhru])

    for comp_name in budget.components:
        for var_key, var_val in budget[comp_name].items():
            assert budget[comp_name][var_key] is var_val

    budget.calculate()
    assert (budget._inputs_sum == 0).all()
    assert (budget._outputs_sum == 1).all()
    assert (budget._storage_changes_sum == -1).all()
    assert (budget.balance == -1).all()

    accum_answers = {
        "inputs": {"in1": np.ones([nhru]), "in2": np.ones([nhru]) * 3},
        "outputs": {"out1": np.zeros([nhru]), "out2": np.ones([nhru]) * 3},
        "storage_changes": {
            "stor1": np.ones([nhru]) * 3,
            "stor2": np.ones([nhru]) * -2,
        },
    }
    for component in budget.components:
        for var in budget[component].keys():
            assert (
                budget.accumulations[component][var]
                == accum_answers[component][var]
            ).all()

    with pytest.raises(ValueError):
        storage_changes["stor2"][:] = -23 * np.ones([nhru])
        budget.calculate()

    return
