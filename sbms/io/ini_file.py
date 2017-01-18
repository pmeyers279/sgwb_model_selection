import ConfigParser
from collections import OrderedDict

def config_section_map(c, section):
    dict1 = OrderedDict()
    options = c.options(section)
    for option in options:
        dict1[option] = c.get(section, option).strip()
    return dict1


def config_list(c):
    dict1 = {}
    for section in c.sections():
        dict1[section] = config_section_map(c, section)
    return dict1


def config_pipeline_ini(c):
    dict1 = {}
    for section in c.sections():
        dict1[section] = config_section_map(c, section)
    return dict1


def read_ini(file):
    """
    Read ini file for analysis and put
    it into a dict.

    Parameters
    ----------
    file : str
        ini file for coherence pipeline

    Returns
    -------
    dict1 : dict
        dictionary with parameters for pipeline
    """
    c = ConfigParser.ConfigParser()
    c.optionxform=str
    c.read(file)
    dict1 = config_pipeline_ini(c)
    return dict1

