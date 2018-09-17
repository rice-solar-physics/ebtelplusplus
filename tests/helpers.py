"""
Helper functions for tests
"""
import os
import sys
import subprocess
import json

import numpy as np
import hissw

TOPDIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.join(TOPDIR, 'rsp_toolkit/python'))
from xml_io import OutputHandler


def run_ebtelplusplus(config, config_filename):
    """
    Given a config `dict`, run ebtel++ and return results
    """
    oh = OutputHandler(config_filename, config)
    oh.print_to_xml()
    subprocess.call([os.path.join(TOPDIR, 'bin', 'ebtel++.run'), '-c', config_filename])
    data = np.loadtxt(config['output_filename'])
    return {
        'time': data[:, 0],
        'electron_temperature': data[:, 1],
        'ion_temperature': data[:, 2],
        'density': data[:, 3],
        'electron_pressure': data[:, 4],
        'ion_pressure': data[:, 5], 
        'velocity': data[:, 6]
    }


def generate_idl_test_data(ebtel_idl_path, config):
    flags = []
    if 'dem' not in config or not config['dem']['use_new_method']:
        flags += ['dem_old']
    if not config['use_flux_limiting']:
        flags += ['classical']
    time = np.arange(0, config['total_time']-config['tau'], config['tau'])
    heat = np.ones(time.shape) * config['heating']['background']
    for _e in config['heating']['events']:
        e = _e['event']
        # Rise
        i = np.where(np.logical_and(time >= e['rise_start'], time < e['rise_end']))
        heat[i] += e['magnitude'] * (time[i] - e['rise_start']) / (e['rise_end'] - e['rise_start'])
        # Plateau
        i = np.where(np.logical_and(time >= e['rise_end'], time < e['decay_start']))
        heat[i] += e['magnitude']
        # Decay
        i = np.where(np.logical_and(time >= e['decay_start'], time <= e['decay_end']))
        heat[i] += e['magnitude'] * (e['decay_end'] - time[i])/(e['decay_end'] - e['decay_start'])

    args = {
        'time': time.tolist(),
        'loop_length': config['loop_length'],
        'heat': heat.tolist(),
        'flags': flags,
    }
    idl = hissw.ScriptMaker(extra_paths=[ebtel_idl_path])
    script = """time={{ time }}
heat = {{ heat }}
loop_length = {{ loop_length }}
ebtel2,time,heat,loop_length,temperature,density,pressure,velocity{% if flags %}, /{{ flags | join(', /') }}{% endif %}
    """
    return idl.run(script, args=args)


def read_idl_test_data(data_filename, ebtel_idl_path, config):
    varnames = ['time', 'temperature', 'density', 'pressure', 'velocity']
    # Generate and save if it does not exist
    if not os.path.isfile(data_filename):
        data = generate_idl_test_data(ebtel_idl_path, config)
        data_array = np.zeros(data['time'].shape+(len(data),))
        for i, v in enumerate(varnames):
            data_array[:, i] = data[v]
        np.savetxt(data_filename, data_array)
    # Load data into a dictionary
    data = np.loadtxt(data_filename)
    return {v: data[:, i] for i, v in enumerate(varnames)}
