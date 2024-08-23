"""
Copyright (c) 2024 Will Barnes. All rights reserved.

ebtelplusplus: Zero-dimensional hydrodynamics of coronal loops
"""
import pathlib
import tempfile
from collections import namedtuple

import astropy.units as u

from astropy.utils.data import get_pkg_data_path
from ebtelplusplus.util import write_xml

from ._core import run as _run_cpp

__all__ = ["run"]


_UNITS_MAPPING = {
    'time': 's',
    'electron_temperature': 'K',
    'ion_temperature': 'K',
    'density': 'cm-3',
    'electron_pressure': 'dyne cm-2',
    'ion_pressure': 'dyne cm-2',
    'velocity': 'cm s-1',
    'heat': 'erg cm-3 s-1',
    'dem_temperature': 'K',
    'dem_tr': 'cm-5 K-1',
    'dem_corona': 'cm-5 K-1',
    'electron_thermal_conduction': 'erg cm-2 s-1',
    'ion_thermal_conduction': 'erg cm-2 s-1',
    'radiative_loss': 'erg cm-3 s-1',
    'tr_corona_radiative_loss_ratio': '',
}

EbtelResult = namedtuple('EbtelResult',
                         _UNITS_MAPPING.keys(),
                         defaults=[None,]*len(_UNITS_MAPPING))


def run(config):
    """
    Run an ebtel++ simulation

    Parameters
    ----------
    config: `dict`
        Dictionary of configuration options

    Returns
    -------
    results: `dict`
        Dictionary of ebtel results
    """
    # TODO: refactor to accept inputs directly and avoid roundtripping config to disk
    config['radiation_data_dir'] = get_pkg_data_path('radiation', package='ebtelplusplus.data')
    with tempfile.TemporaryDirectory() as tmpdir:
        config_filename = pathlib.Path(tmpdir) / 'ebtelplusplus.tmp.xml'
        write_xml(config, config_filename)
        results = _run_cpp(str(config_filename))

    return EbtelResult(**{k: u.Quantity(v, _UNITS_MAPPING[k]) for k,v in results.items()})
