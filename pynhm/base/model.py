from copy import deepcopy
from pprint import pprint

import numpy as np

from pynhm.base.adapter import adapter_factory
from pynhm.base.control import Control
from pynhm.utils.parameters import PrmsParameters


class Model:
    def __init__(
        self,
        *component_classes,
        control: Control,
        params: PrmsParameters,
        input_dir: str = None,
        budget_type: str = "strict",  # also pass dict
        verbose: bool = False,
    ):

        self.control = control
        self.params = params
        self.input_dir = input_dir
        self.verbose = verbose

        class_dict = {comp.__name__: comp for comp in component_classes}
        class_inputs = {kk: vv.get_inputs() for kk, vv in class_dict.items()}
        class_vars = {kk: vv.get_variables() for kk, vv in class_dict.items()}

        # Solve where inputs come from
        inputs_from = {}
        for comp in class_dict.keys():
            inputs = deepcopy(class_inputs)
            vars = deepcopy(class_vars)
            c_inputs = inputs.pop(comp)
            c_vars = vars.pop(comp)
            inputs_from[comp] = {}
            # inputs
            for input in c_inputs:
                inputs_from[comp][input] = []  # could use None
                for other in inputs.keys():
                    if input in vars[other]:
                        inputs_from[comp][input] += [other]
                        # this should be a list of length one
                        # check?

        # determine component order
        class_deps = {}
        for k0, v0 in inputs_from.items():
            class_deps[k0] = set([dep[0] for dep in list(v0.values()) if dep])

        dep_lens = {key: len(val) for key, val in class_deps.items()}
        # this will need to be more sophisticated...
        self.component_order = [None for dd in dep_lens.keys()]
        for key, val in dep_lens.items():
            self.component_order[val] = key

        # If inputs dont come from other components, assume they come from
        # file in input_dir
        file_input_names = set([])
        for k0, v0 in inputs_from.items():
            for k1, v1 in v0.items():
                if not v1:
                    file_input_names = file_input_names.union([k1])

        # initiate the file inputs here rather than in the components
        file_inputs = {}
        for name in file_input_names:
            nc_path = input_dir / f"{name}.nc"
            file_inputs[name] = adapter_factory(nc_path, name, control=control)

        # instance dict
        self.components = {}
        for component in self.component_order:
            component_inputs = {
                input: None for input in class_dict[component].get_inputs()
            }
            self.components[component] = class_dict[component](
                control=control,
                params=params,
                budget_type=budget_type,
                **component_inputs,
            )

        for component in self.component_order:
            print(f"{component}:")
            for input, frm in inputs_from[component].items():
                if not frm:
                    print(f"    {input}: from file")
                    self.components[component].set_input_to_adapter(
                        input, file_inputs[input]
                    )
                else:
                    print(f"    {input}: from {frm}")
                    self.components[component].set_input_to_adapter(
                        input,
                        adapter_factory(
                            self.components[frm[0]][input], control=control
                        ),  # drop list above
                    )

    def advance(self):
        self.control.advance()
        for cls in self.component_order:
            if self.verbose:
                print(f"advancing component: {cls}")
            self.components[cls].advance()
        return

    def calculate(self):
        for cls in self.component_order:
            if self.verbose:
                print(f"calculating component: {cls}")
            self.components[cls].calculate(1.0)
        return
