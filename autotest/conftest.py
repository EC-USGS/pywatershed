import pathlib
import yaml
from yaml import Loader

def pytest_addoption(parser):
    parser.addoption(
        '--domain_yaml',
        required=False,
        action='append',
        default=['../test_data/drb_2yr/drb_2yr.yaml'],
        help=(
            'YAML file(s) for indiv domain tests '
            '(can pass multiples of this argument)'))


def pytest_generate_tests(metafunc):
    if "domain" in metafunc.fixturenames:
        domain_file_list = metafunc.config.getoption("domain_yaml")

        # open and read in the yaml and
        domain_list = []
        for dd in domain_file_list:
            dd_file = pathlib.Path(dd)
            with dd_file.open('r') as yaml_file:
                domain_dict = yaml.safe_load(yaml_file)
                # Construct/derive some convenience quantities
                # JLM probably just want to transform all paths to the rel path
                domain_dict['file'] = dd_file
                domain_dict['dir'] = dd_file.parent
                domain_dict['input_files_dict'] = {
                    key: domain_dict['dir'] / val
                    for key, val in domain_dict['cbh_inputs'].items()}
                # append to the list of all domains
                domain_list += [domain_dict]

        metafunc.parametrize("domain", domain_list)
