"""
Utility functions for configuring and running ebtel++ simulations
"""
import os
import subprocess
import warnings
from collections import OrderedDict
import tempfile
import xml.etree.ElementTree as ET
import xml.dom.minidom as xdm

import numpy as np

__all__ = ['run_ebtel', 'read_xml', 'write_xml']


class EbtelPlusPlusError(Exception):
    """
    Raise this exception when there's an ebtel++ error
    """
    pass


def run_ebtel(config, ebtel_dir):
    """
    Run an ebtel++ simulation

    Parameters
    ----------
    config: `dict`
        Dictionary of configuration options
    ebtel_dir: `str`
        Path to directory containing ebtel++ source code.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        config_filename = os.path.join(tmpdir, 'ebtelplusplus.tmp.xml')
        results_filename = os.path.join(tmpdir, 'ebtelplusplus.tmp')
        config['output_filename'] = results_filename
        write_xml(config, config_filename)
        cmd = subprocess.run(
            [os.path.join(ebtel_dir, 'bin/ebtel++.run'), '-c', config_filename],
            shell=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if cmd.stderr:
            raise EbtelPlusPlusError(f"{cmd.stderr.decode('utf-8')}")
        data = np.loadtxt(results_filename)
    
    results = {
        'time': data[:, 0],
        'electron_temperature': data[:, 1],
        'ion_temperature': data[:, 2],
        'density': data[:, 3],
        'electron_pressure': data[:, 4],
        'ion_pressure': data[:, 5],
        'heat': data[:, 7],
    }
    
    results_dem = {}
    if config['calculate_dem']:
        results_dem['dem_tr'] = np.loadtxt(
            config['output_filename'] + '.dem_tr')
        results_dem['dem_corona'] = np.loadtxt(
            config['output_filename'] + '.dem_corona')
        # The first row of both is the temperature bins
        results_dem['dem_temperature'] = results['dem_tr'][0, :]
        results_dem['dem_tr'] = results_dem['dem_tr'][1:, :]
        results_dem['dem_corona'] = results_dem['dem_corona'][1:, :]
    
    return {**results, **results_dem}



def read_xml(input_filename,):
    """
    For all input variables, find them in the XML tree and return them to a
    dictionary

    Parameters
    ----------
    input_filename : `str`
    """
    tree = ET.parse(input_filename)
    root = tree.getroot()
    var_list = [child.tag for child in root]
    input_dict = {}
    for var in var_list:
        # Find node
        node = root.find(var)
        # Check if found
        if node is None:
            warnings.warn(f'No value found for input {var}. Returning None.')
            input_dict[var] = None
        else:
            input_dict[node.tag] = read_node(node)

    return input_dict


def read_node(node):
    """
    Read in node values for different configurations
    """
    if node.getchildren():
        _child_tags = [child.tag for child in node.getchildren()]
        if len(_child_tags) != len(set(_child_tags)):
            tmp = []
            for child in node.getchildren():
                tmp.append({child.tag: read_node(child)})
        else:
            tmp = OrderedDict()
            for child in node.getchildren():
                tmp[child.tag] = read_node(child)
        return tmp
    else:
        if node.text:
            return type_checker(node.text)
        elif node.attrib:
            return {key: type_checker(node.attrib[key]) for key in node.attrib}
        else:
            warnings.warn(
                f'Unrecognized node format for {node.tag}. Returning None.')
            return None


def bool_filter(val):
    """
    Convert true/false string to Python bool. Otherwise, return string.
    """
    trues = ['True', 'TRUE', 'true', 'yes', 'Yes']
    falses = ['False', 'FALSE', 'false', 'no', 'No']
    if any([val == t for t in trues]):
        return True
    elif any([val == f for f in falses]):
        return False
    else:
        return val


def type_checker(val):
    """
    Convert to int or float if possible
    """
    try:
        return int(val)
    except ValueError:
        pass
    try:
        return float(val)
    except ValueError:
        pass
    return bool_filter(val)


def write_xml(output_dict, output_filename):
    """
    Print dictionary to XML file

    Parameters
    ----------
    output_dict : `dict`
        structure to print to file
    output_filename : `str`
        filename to print to
    """
    root = ET.Element('root')
    for key in output_dict:
        set_element_recursive(root, output_dict[key], key)
    with open(output_filename, 'w') as f:
        f.write(pretty_print_xml(root))


def set_element_recursive(root, node, keyname):
    """
    Set element tags, values, and attributes. Recursive for arrays/lists.
    """
    element = ET.SubElement(root, keyname)
    if type(node) is list:
        for item in node:
            sub_keyname = [k for k in item][0]
            set_element_recursive(element, item[sub_keyname], sub_keyname)
    elif type(node).__name__ == 'OrderedDict':
        for key in node:
            set_element_recursive(element, node[key], key)
    elif type(node) is dict:
        for key in node:
            element.set(key, str(node[key]))
    else:
        element.text = str(node)


def pretty_print_xml(element):
    """
    Formatted XML output for writing to file.
    """
    unformatted = ET.tostring(element)
    xdmparse = xdm.parseString(unformatted)
    return xdmparse.toprettyxml(indent="    ")
