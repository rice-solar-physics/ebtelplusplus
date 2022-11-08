"""
Helper functions for tests
"""
import os
import sys

import numpy as np

TOPDIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.join(TOPDIR, 'examples'))
from util import run_ebtel


def run_ebtelplusplus(config):
    return run_ebtel(config, TOPDIR)


def generate_idl_test_data(ebtel_idl_path, config):
    # Import here to avoid making this a hard dependency for this whole
    # script.
    import hissw
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
    idl = hissw.Environment(extra_paths=[ebtel_idl_path])
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
