import pathlib as pl

import yaml

import pywatershed

yaml_path = "../ucb_2yr.yaml"
dd_file = pl.Path(yaml_path)
with dd_file.open("r") as yaml_file:
    domain_dict = yaml.safe_load(yaml_file)

# Construct/derive some convenience quantities
domain_dict["file"] = dd_file
domain_dict["dir"] = dd_file.parent

for ff in (
    "param_file",
    "control_file",
    "cbh_nc",
):
    domain_dict[ff] = domain_dict["dir"] / domain_dict[ff]

for fd_key in ("cbh_inputs",):
    domain_dict[fd_key] = {
        key: (domain_dict["dir"] / val)
        for key, val in domain_dict[fd_key].items()
    }

# Construct a dictionary that gets used in CBH
domain_dict["input_files_dict"] = {
    key: val for key, val in domain_dict["cbh_inputs"].items()
}

parameter_file = domain_dict["param_file"]
parameters = pywatershed.PrmsParameters.load(domain_dict["param_file"])


cbh = pywatershed.CBH(
    domain_dict["input_files_dict"],
    parameters=parameters,
    adjust=True,
)

global_atts = {"domain_name": domain_dict["domain_name"]}
_ = cbh.to_netcdf(domain_dict["cbh_nc"], global_atts=global_atts)

print(f"created...'{domain_dict['cbh_nc']}'")
