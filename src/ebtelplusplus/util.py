"""
Utility functions for configuring and running ebtel++ simulations
"""
import warnings
import xml.dom.minidom as xdm
import xml.etree.ElementTree as ET
from collections import OrderedDict

__all__ = ['read_xml', 'write_xml']


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
    if list(node):
        _child_tags = [child.tag for child in node]
        if len(_child_tags) != len(set(_child_tags)):
            tmp = []
            for child in node:
                tmp.append({child.tag: read_node(child)})
        else:
            tmp = OrderedDict()
            for child in node:
                tmp[child.tag] = read_node(child)
        return tmp

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
