import os
import pathlib as pl
import sys
from fnmatch import fnmatch
from platform import processor

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--force",
        required=False,
        action="store_true",
        default=False,
        help=("Force run without scheduler present"),
    )

    parser.addoption(
        "--domain",
        required=False,
        action="append",
        default=[],
        help=(
            "Domain(s) to run (name of domain dir and NOT path to it). "
            "You can pass multiples of this argument. If not used, "
            "all self-contained domains in test_data/ will be run."
        ),
    )


@pytest.fixture()
def exe():
    platform = sys.platform.lower()
    if platform == "win32":
        exe_name = "prms_win_ifort_mixed_prec.exe"
    elif platform == "darwin":
        if processor() == "arm":
            exe_name = "prms_mac_m1_ifort_dbl_prec"
        else:
            exe_name = "prms_mac_intel_gfort_dbl_prec"
    elif platform == "linux":
        exe_name = "prms_linux_gfort_dbl_prec"
    exe_pth = pl.Path(f"../../bin/{exe_name}").resolve()
    return exe_pth


rootdir = ".."
test_dirs = sorted(
    [path for path in pl.Path(rootdir).iterdir() if path.is_dir()]
)


# This would change to handle other/additional schedulers
domain_globs_schedule = ["*conus*"]

# For generating timeseries of previous states
previous_vars = [
    "dprst_stor_hru",
    "freeh2o",
    "hru_impervstor",
    "pk_ice",
    "pref_flow_stor",
    "slow_stor",
    "soil_lower",
    "soil_moist",
    "soil_rechr",
    "ssres_stor",
]

misc_nc_file_vars = [
    "infil",
    "sroff",
    "ssres_flow",
    "gwres_flow",
]


final_nc_file_vars = [
    "through_rain",
    "pk_ice",
    "freeh2o",
    "hru_impervstor",
    "dprst_stor_hru",
]


def scheduler_active():
    slurm = os.getenv("SLURM_JOB_ID") is not None
    pbs = os.getenv("PBS_JOBID") is not None
    # add more scheduler checks
    return slurm or pbs


def enforce_scheduler(test_dir):
    if scheduler_active():
        return None
    glob_match = list(
        fnmatch(str(test_dir), gg) for gg in domain_globs_schedule
    )
    if any(glob_match):
        raise RuntimeError(
            f"Domain '{test_dir}' must be scheduled (use --force to override)"
        )
    return None


def collect_simulations(domain_list: list, force: bool):
    simulations = {}
    for test_dir in test_dirs:
        for pth in test_dir.iterdir():
            # checking for prcp.cbh ensure this is a self-contained run (all
            # files in repo)
            if (
                (test_dir / "prcp.cbh").exists()
                and pth.is_file()
                and pth.name == "control.test"
            ):
                if len(domain_list) and (test_dir.name not in domain_list):
                    continue

                if not force:
                    enforce_scheduler(test_dir)

                # add simulation
                simulations[str(test_dir)] = pth.name

    if len(domain_list) and (len(simulations) < len(domain_list)):
        requested = set(domain_list)
        found = [pl.Path(dd).name for dd in simulations.keys()]
        requested_not_found = requested.difference(found)
        msg = (
            f"The following requested domains were not found: "
            f"{requested_not_found}"
        )
        pytest.exit(msg)

    print("\nrun_domains.py found the following domains to run:\n")
    print(f"{list(simulations.keys())}")
    return simulations


def collect_csv_files(domain_list: list, force: bool):
    simulations = collect_simulations(domain_list, force)
    csv_files = []
    for key, value in simulations.items():
        output_pth = pl.Path(key) / "output"
        csv_files_dom = sorted(output_pth.glob("*.csv"))
        csv_files += [ff for ff in csv_files_dom if ff.name != "stats.csv"]
    return csv_files


def collect_misc_nc_files(domain_list: list, var_list: list, force: bool):
    simulations = collect_simulations(domain_list, force)
    sim_dirs = list(simulations.keys())
    misc_nc_files = []
    for var in var_list:
        for sim in sim_dirs:
            the_file = pl.Path(sim) / f"output/{var}.nc"
            # assert the_file.exists()
            misc_nc_files += [the_file.with_suffix("")]

    return misc_nc_files


def pytest_generate_tests(metafunc):
    domain_list = metafunc.config.getoption("domain")
    force = metafunc.config.getoption("force")

    if "simulations" in metafunc.fixturenames:
        simulations = collect_simulations(domain_list, force)
        sim_list = [
            {"ws": key, "control_file": val}
            for key, val in simulations.items()
        ]
        ids = [pl.Path(ss).name for ss in simulations.keys()]
        metafunc.parametrize("simulations", sim_list, ids=ids)

    if "csv_files" in metafunc.fixturenames:
        csv_files = collect_csv_files(domain_list, force)
        ids = [ff.parent.parent.name + ":" + ff.name for ff in csv_files]
        metafunc.parametrize("csv_files", csv_files, ids=ids)

    if "csv_files_prev" in metafunc.fixturenames:
        csv_files = collect_csv_files(domain_list, force)
        csv_files = [
            ff for ff in csv_files if ff.with_suffix("").name in previous_vars
        ]
        ids = [ff.parent.parent.name + ":" + ff.name for ff in csv_files]
        metafunc.parametrize("csv_files_prev", csv_files, ids=ids)

    if "misc_nc_files_input" in metafunc.fixturenames:
        misc_nc_files = collect_misc_nc_files(
            domain_list, misc_nc_file_vars, force
        )
        ids = [ff.parent.parent.name + ":" + ff.name for ff in misc_nc_files]
        metafunc.parametrize("misc_nc_files_input", misc_nc_files, ids=ids)

    if "misc_nc_final_input" in metafunc.fixturenames:
        misc_nc_files = collect_misc_nc_files(
            domain_list, final_nc_file_vars, force
        )
        ids = [ff.parent.parent.name + ":" + ff.name for ff in misc_nc_files]
        metafunc.parametrize("misc_nc_final_input", misc_nc_files, ids=ids)

    if "soltab_file" in metafunc.fixturenames:
        simulations = collect_simulations(domain_list, force)
        soltab_files = [
            pl.Path(kk) / "soltab_debug" for kk in simulations.keys()
        ]
        ids = [ff.parent.name + ":" + ff.name for ff in soltab_files]
        metafunc.parametrize("soltab_file", soltab_files, ids=ids)
