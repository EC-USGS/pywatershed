import re
from copy import deepcopy
from pprint import pprint

import numpy as np

from pynhm.base.adapter import adapter_factory
from pynhm.base.control import Control


class Model:
    def __init__(
        self,
        *component_classes,
        control: Control,
        input_dir: str = None,
        budget_type: str = "strict",  # also pass dict
        verbose: bool = False,
    ):

        self.control = control
        self.input_dir = input_dir
        self.verbose = verbose

        class_dict = {comp.__name__: comp for comp in component_classes}
        class_inputs = {kk: vv.get_inputs() for kk, vv in class_dict.items()}
        class_vars = {kk: vv.get_variables() for kk, vv in class_dict.items()}

        # Solve where inputs come from
        inputs_from = {}
        # inputs_from_prev = {}
        for comp in class_dict.keys():
            inputs = deepcopy(class_inputs)
            vars = deepcopy(class_vars)
            c_inputs = inputs.pop(comp)
            _ = vars.pop(comp)
            inputs_from[comp] = {}
            # inputs_from_prev[comp] = {}
            # inputs
            for input in c_inputs:
                inputs_ptr = inputs_from
                # if re.compile(".*_(ante|prev|old)").match(input):
                #    inputs_ptr = inputs_from_prev
                inputs_ptr[comp][input] = []  # could use None
                for other in inputs.keys():
                    if input in vars[other]:
                        inputs_ptr[comp][input] += [other]
                        # this should be a list of length one
                        # check?

        # determine component order
        # for now this is sufficient. when more classes exist, we can analyze
        # dependencies, might inputs_from_prev above.
        component_order_prms = [
            "PRMSSolarGeometry",
            "PRMSBoundaryLayer",
            "PRMSCanopy",
            "PRMSSnow",
            "PRMSRunoff",
            "PRMSSoilzone",
            "PRMSEt",
            "PRMSGroundwater",
            "PRMSChannel",
        ]
        self.component_order = []
        for comp in component_order_prms:
            if comp in class_dict.keys():
                self.component_order += [comp]

        missing_components = set(self.component_order).difference(
            set(class_dict.keys())
        )
        if missing_components:
            raise ValueError(
                f"Component {missing_components} not currently "
                f"handled by Model class."
            )

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

        # instantiate components: instance dict
        self.components = {}
        for component in self.component_order:
            component_inputs = {
                input: None for input in class_dict[component].get_inputs()
            }
            # This is a hack. Need to make StorageUnit a subclass of a more
            # general model/element/component class
            args = {
                "control": control,
                **component_inputs,
            }
            if component not in ["PRMSSolarGeometry"]:
                args["budget_type"] = budget_type
            self.components[component] = class_dict[component](**args)

        # Wire it up
        self.component_input_from = {}
        for component in self.component_order:
            self.component_input_from[component] = {}
            for input, frm in inputs_from[component].items():
                if not frm:
                    fname = file_inputs[input]._fname
                    self.component_input_from[component][input] = fname
                    self.components[component].set_input_to_adapter(
                        input, file_inputs[input]
                    )
                else:
                    self.component_input_from[component][input] = frm[0]
                    self.components[component].set_input_to_adapter(
                        input,
                        adapter_factory(
                            self.components[frm[0]][input], control=control
                        ),  # drop list above
                    )

        return

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
