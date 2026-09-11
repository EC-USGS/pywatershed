"""Snapshot pywatershed's public API surface as sorted text.

Writes (or prints) one fact per line, sorted, so that any change to the
public surface appears as a line-local diff. The committed baseline is
autotest/api_surface.txt (beside this script; a repo path, not a package
path, so it works for non-editable installs); autotest/test_api_surface.py
fails when the installed package no longer matches it.

The snapshot records *resolved* values, not source text: the declared
name-set methods (get_inputs, get_variables, ...) are called, so a
removed override that falls through to a base class implementation
changes the snapshot even though the subclass source shows no method.

Usage (from the repo root):
    python autotest/api_surface.py           # print to stdout
    python autotest/api_surface.py --write   # update the baseline
"""

import inspect
import pathlib as pl
import sys

import yaml

import pywatershed as pws

BASELINE = pl.Path(__file__).parent / "api_surface.txt"

# Class methods declaring the API's name sets. Called (not parsed), so
# inherited implementations are recorded as the user sees them.
DECLARED_METHODS = [
    "get_dimensions",
    "get_inputs",
    "get_init_values",
    "get_mass_budget_terms",
    "get_energy_budget_terms",
    "get_parameters",
    "get_restart_variables",
    "get_variables",
]


def _format_signature(func) -> str:
    """A signature as text: argument names and defaults, no annotations.

    Annotations are excluded deliberately: annotation-only corrections
    (e.g. a wrong return type fixed) are not API changes and should not
    churn the baseline.
    """
    parts = []
    for p in inspect.signature(func).parameters.values():
        name = p.name
        if p.kind is inspect.Parameter.VAR_POSITIONAL:
            name = "*" + name
        elif p.kind is inspect.Parameter.VAR_KEYWORD:
            name = "**" + name
        if p.default is not inspect.Parameter.empty:
            if isinstance(p.default, pl.PurePath):
                # repr is PosixPath/WindowsPath, which varies by
                # platform; the baseline must not.
                name += f"=Path('{p.default.as_posix()}')"
            else:
                name += f"={p.default!r}"
        parts.append(name)
    return "(" + ", ".join(parts) + ")"


def _declared_lines(cls_name: str, cls) -> list:
    lines = []
    for meth_name in DECLARED_METHODS:
        meth = getattr(cls, meth_name, None)
        if meth is None or not callable(meth):
            continue
        try:
            result = meth()
        except Exception as ee:
            lines.append(
                f"set: {cls_name}.{meth_name} raises {type(ee).__name__}"
            )
            continue
        if isinstance(result, dict):
            if meth_name in (
                "get_mass_budget_terms",
                "get_energy_budget_terms",
            ):
                for term, names in sorted(result.items()):
                    names = " ".join(sorted(names))
                    lines.append(
                        f"set: {cls_name}.{meth_name} {term}: {names}"
                    )
            else:
                keys = " ".join(sorted(result.keys()))
                lines.append(f"set: {cls_name}.{meth_name} dict-keys: {keys}")
        else:
            container = type(result).__name__
            names = " ".join(sorted(str(vv) for vv in result))
            lines.append(f"set: {cls_name}.{meth_name} {container}: {names}")
    return lines


def _metadata_lines() -> list:
    lines = []
    meta_dir = pl.Path(pws.__file__).parent / "static" / "metadata"
    for kind, fname in [
        ("meta-var", "variables.yaml"),
        ("meta-param", "parameters.yaml"),
    ]:
        data = yaml.safe_load((meta_dir / fname).read_text())
        for name, entry in data.items():
            units = entry.get("units", "-")
            dtype = entry.get("type", "-")
            lines.append(f"{kind}: {name} type={dtype} units={units}")
    return lines


def _control_option_lines() -> list:
    from pywatershed.base.control import (
        prms_legacy_options_avail,
        pws_control_options_avail,
    )

    return [
        f"control-option: pws {oo}" for oo in pws_control_options_avail
    ] + [
        f"control-option: prms-legacy {oo}" for oo in prms_legacy_options_avail
    ]


def generate() -> str:
    """The full API surface, one sorted fact per line."""
    lines = []
    for name in pws.__all__:
        lines.append(f"export: {name}")
        obj = getattr(pws, name, None)
        if obj is None:
            continue
        if inspect.isclass(obj):
            lines.append(f"init: {name}{_format_signature(obj.__init__)}")
            lines += _declared_lines(name, obj)
        elif inspect.isfunction(obj):
            lines.append(f"func: {name}{_format_signature(obj)}")

    lines += _control_option_lines()
    lines += _metadata_lines()

    return "\n".join(sorted(lines)) + "\n"


if __name__ == "__main__":
    text = generate()
    if "--write" in sys.argv[1:]:
        BASELINE.write_text(text)
        print(f"wrote {BASELINE}")
    else:
        sys.stdout.write(text)
