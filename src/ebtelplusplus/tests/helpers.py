"""
Helper functions for tests
"""
import pathlib

import astropy.units as u
import h5py
import numpy as np
from astropy.utils.data import get_pkg_data_path


def generate_idl_test_data(ebtel_idl_path, config):
    # Import here to avoid making this a hard dependency for this whole
    # script.
    import hissw
    flags = []
    if not config['dem_use_new_tr_method']:
        flags += ['dem_old']
    if not config['use_flux_limiting']:
        flags += ['classical']
    # Set expansions keywords
    # NOTE: It is assumed that A_TR is always 1 here which we can do wlog
    keywords = {
        'a_tr': 1.0,
        'a_c': 1/config['area_ratio_tr_corona'],
        'a_0': config['area_ratio_0_corona']/config['area_ratio_tr_corona'],
        'l_fact': 1.0 - config['loop_length_ratio_tr_total'],
    }
    # keywords = [f'{k}={v:.6f}' for k,v in keywords.items()]
    time = np.arange(0, config['total_time']-config['tau'], config['tau'])
    heat = np.ones(time.shape) * config['heating']['background']
    for e in config['heating']['events']:
        # Rise
        i = np.where(np.logical_and(time >= e['rise_start'], time < e['rise_end']))
        heat[i] += e['rate'] * (time[i] - e['rise_start']) / (e['rise_end'] - e['rise_start'])
        # Plateau
        i = np.where(np.logical_and(time >= e['rise_end'], time < e['decay_start']))
        heat[i] += e['rate']
        # Decay
        i = np.where(np.logical_and(time >= e['decay_start'], time <= e['decay_end']))
        heat[i] += e['rate'] * (e['decay_end'] - time[i])/(e['decay_end'] - e['decay_start'])

    args = {
        'time': time.tolist(),
        'loop_length': config['loop_length'],
        'heat': heat.tolist(),
        'flags': flags,
        **keywords,
    }
    idl = hissw.Environment(extra_paths=[ebtel_idl_path], idl_only=True)
    script = """time = {{ time | force_double_precision }}
heat = {{ heat | force_double_precision }}
loop_length = {{ loop_length | force_double_precision }}
ebtel2,time,heat,loop_length,temperature,density,pressure,velocity,$
    a_tr={{ a_tr | force_double_precision }},$
    a_c={{ a_c | force_double_precision }},$
    a_0={{ a_0 | force_double_precision }},$
    l_fact={{ l_fact | force_double_precision }}{% if flags %}, /{{ flags | join(', /') }}{% endif %}
    """
    return idl.run(script, args=args)


def read_idl_test_data(data_filename, ebtel_idl_path, config):
    varnames = ['time', 'temperature', 'density', 'pressure', 'velocity']
    varunits = ['s', 'K', 'cm-3', 'dyne cm-2', 'cm s-1']
    # Generate and save if it does not exist
    data_dir = pathlib.Path(get_pkg_data_path('data', package='ebtelplusplus.tests'))
    if ebtel_idl_path is not None:
        data = generate_idl_test_data(ebtel_idl_path, config)
        data_array = np.zeros(data['time'].shape+(len(data),))
        for i, v in enumerate(varnames):
            data_array[:, i] = data[v]
        np.savetxt(data_dir / data_filename, data_array)
    # Load data into a dictionary
    data = np.loadtxt(data_dir / data_filename)
    return {v: u.Quantity(data[:, i], vu) for i, (v,vu) in enumerate(zip(varnames, varunits))}


def read_hydrad_test_data(data_filename, tau, heating):
    data_dir = pathlib.Path(get_pkg_data_path('data', package='ebtelplusplus.tests'))
    data = {}
    with h5py.File(data_dir / data_filename, 'r') as hf:
        grp = hf[f'/{heating}/tau{tau:.0f}']
        data['time'] = u.Quantity(np.asarray(hf['time']), hf['time'].attrs['unit'])
        data['electron_temperature'] = u.Quantity(np.asarray(grp['electron_temperature']),
                                                  grp['electron_temperature'].attrs['unit'])
        data['ion_temperature'] = u.Quantity(np.asarray(grp['ion_temperature']),
                                             grp['ion_temperature'].attrs['unit'])
        data['density'] = u.Quantity(np.asarray(grp['density']),
                                     grp['density'].attrs['unit'])
    return data


def plot_comparison(r_cpp, r_idl):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(8,9),layout='constrained')
    ax = fig.add_subplot(311)
    ax.plot(r_cpp.time, r_cpp.electron_temperature, label='ebtel++')
    ax.plot(r_idl['time'], r_idl['temperature'], label='IDL')
    ax.set_xlabel('$t$ [s]')
    ax.set_ylabel('$T$ [K]')
    ax.legend()
    ax = fig.add_subplot(312,sharex=ax)
    ax.plot(r_cpp.time, r_cpp.density)
    ax.plot(r_idl['time'], r_idl['density'])
    ax.set_xlabel('$t$ [s]')
    ax.set_ylabel('$n$ [cm$^{-3}$]')
    ax = fig.add_subplot(313)
    ax.plot(r_cpp.time, (r_cpp.electron_temperature-r_idl['temperature'])/r_idl['temperature'],
             label='$T$')
    ax.plot(r_cpp.time, (r_cpp.density-r_idl['density'])/r_idl['density'],
             label='$n$')
    ax.axhline(y=0, color='k', ls='--')
    ax.axhline(y=0.1, color='k', ls=':')
    ax.axhline(y=-0.1, color='k', ls=':')
    ax.set_ylim(-0.5,0.5)
    ax.set_xlabel('$t$ [s]')
    ax.set_ylabel(r'$(q_\mathrm{C++} - q_\mathrm{IDL}) / q_\mathrm{IDL}$')
    ax.legend()
    plt.show()
