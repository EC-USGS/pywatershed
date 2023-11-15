import pathlib as pl
import re
from collections import OrderedDict

import xmltodict
import yaml

# This requires prms6 to be checkout from the prms repo.

xml_dir = pl.Path("../../../prms/src/xml")
xml_files = {
    "control": xml_dir / "control.xml",
    "dimensions": xml_dir / "dimensions.xml",
    "parameters": xml_dir / "parameters.xml",
    "variables": xml_dir / "variables.xml",
}

metadata_dir = pl.Path("../metadata")


def guess_type(str_in):
    if str_in == "none":
        return None
    elif len(re.findall("-", str_in)) > 1 or (
        re.search("[a-zA-Z]", str_in) is not None
    ):
        return str(str_in).replace("'", "")
    elif "." in str_in:
        return float(str_in)
    else:
        return int(str_in)


for xml_key, xml_file in xml_files.items():
    yaml_file = metadata_dir / (xml_key + ".yaml")

    with xml_file.open() as xml_stream:
        xml_dict = xmltodict.parse(xml_stream.read())

    if xml_key == "dimensions":
        # Dimensions
        xml_list = xml_dict[xml_key]["dimension"]
        xml_dict = {}
        for xx in xml_list:
            name = xx["@name"]
            _ = xx.pop("@name")
            xx["default"] = int(xx["default"])
            xx["size"] = int(xx["size"])
            xml_dict[name] = dict(xx)

    elif xml_key == "control":
        # Control
        xml_list = xml_dict[xml_key]["control_param"]
        # xml_list = xml_list[0:2]  # just testing
        xml_dict = {}
        for xx in xml_list:
            name = xx["@name"]
            _ = xx.pop("@name")
            xx["default"] = guess_type(xx["default"])
            xx["type"] = guess_type(xx["type"])
            xx["numvals"] = int(xx["numvals"])
            if "force_default" in xx.keys():
                xx["force_default"] = int(xx["force_default"])

            xml_dict[name] = dict(xx)

            data = xml_dict[name]

            # values

            if "values" in data.keys() and "value" in data.keys():
                print(name)

            if "values" in data.keys():
                vals = dict(data["values"])
                del data["values"]
                data["values_type"] = vals["@type"]
                if "value" in vals.keys():
                    xml_dict[name]["values"] = {}
                    for vv in vals["value"]:
                        key, val = tuple(vv.values())
                        xml_dict[name]["values"][int(key)] = val
                        # print(guess_type(val))
            elif "value" in data.keys():
                # print(f"value: {name}")
                pass
            else:
                # print(f"neither values nor value: {name}")
                pass

            # related variables
            if "related_variables" in data.keys():
                data["related_variables"] = [
                    dict(od)["@name"]
                    for od in data["related_variables"]["variable"]
                ]

    elif xml_key == "parameters":
        xml_list = xml_dict[xml_key]["parameter"]
        # xml_list = xml_list[325:328]  # just testing
        xml_dict = {}
        for xx in xml_list:
            name = xx["@name"]
            _ = xx.pop("@name")

            if "default" in xx.keys():
                xx["default"] = guess_type(xx["default"])
            xx["maximum"] = guess_type(xx["maximum"])
            xx["minimum"] = guess_type(xx["minimum"])

            if "requires" in xx.keys():
                xx["requires"] = dict(xx["requires"])

            xml_dict[name] = dict(xx)
            data = xml_dict[name]

            # dimensions
            dims = dict(data["dimensions"])["dimension"]
            if len(dims) > 1:
                # print(f"dims: {name}")
                xml_dict[name]["dimensions"] = {}
            if isinstance(dims, OrderedDict):
                dims = [dims]
            for dd in dims:
                pos = int(dd["position"]) - 1
                xml_dict[name]["dimensions"][pos] = dd["@name"]

            # modules
            mods = dict(data["modules"])
            mod_val = list(mods.values())[0]  # too convoluted
            if isinstance(mod_val, str):
                xml_dict[name]["modules"] = [mod_val]
            elif isinstance(mod_val, list):
                xml_dict[name]["modules"] = mod_val
            else:
                raise ValueError("more than one module hasnt happened before")

    elif xml_key == "variables":
        xml_list = xml_dict[xml_key]["variable"]
        # xml_list = xml_list[325:328]  # just testing
        xml_dict = {}
        for xx in xml_list:
            name = xx["@name"]
            _ = xx.pop("@name")
            xml_dict[name] = dict(xx)
            data = xml_dict[name]

            # dimensions
            dims = dict(data["dimensions"])["dimension"]
            if len(dims) > 1:
                # print(f"dims: {name}")
                xml_dict[name]["dimensions"] = {}
            if isinstance(dims, OrderedDict):
                dims = [dims]
            for dd in dims:
                pos = int(dd["position"]) - 1
                xml_dict[name]["dimensions"][pos] = dd["@name"]

            # modules
            mods = dict(data["modules"])
            mod_val = list(mods.values())[0]  # too convoluted
            if isinstance(mod_val, str):
                xml_dict[name]["modules"] = [mod_val]
            elif isinstance(mod_val, list):
                xml_dict[name]["modules"] = mod_val
            else:
                raise ValueError("more than one module hasnt happened before")

    with yaml_file.open("w") as yaml_stream:
        yaml.dump(xml_dict, yaml_stream)
