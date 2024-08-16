"""
Tests for the different kinds of radiative loss functions
"""
from collections import OrderedDict

import pytest

from .helpers import run_ebtelplusplus


@pytest.fixture(scope='module')
def base_config():
    base_config = {
        'total_time': 5e3,
        'tau': 1.0,
        'tau_max': 10.0,
        'loop_length': 4e9,
        'loop_length_ratio_tr_total': 0.0,
        'area_ratio_tr_corona': 1.0,
        'area_ratio_0_corona': 1.0,
        'saturation_limit': 1/6,
        'force_single_fluid': False,
        'use_c1_loss_correction': True,
        'use_c1_grav_correction': True,
        'use_flux_limiting': True,
        'calculate_dem': False,
        'save_terms': False,
        'use_adaptive_solver': True,
        'adaptive_solver_error': 1e-6,
        'adaptive_solver_safety': 0.5,
        'c1_cond0': 2.0,
        'c1_rad0': 0.6,
        'helium_to_hydrogen_ratio': 0.075,
        'surface_gravity': 1.0,
        'heating': OrderedDict({
            'partition': 1.0,
            'background': 1e-6,
            'events': [
                {'event': {'rise_start': 0.0, 'rise_end': 100.0, 'decay_start': 100.0,
                           'decay_end': 200.0, 'magnitude': 0.1}}
            ]
        }),
    }
    return base_config

@pytest.mark.parametrize('radiation', ['power_law', 'variable', 'coronal', 'photospheric'])
def test_rad_loss_options(base_config, radiation):
    # Just a smoke test to make sure the new radiative loss options work
    base_config['radiation'] = radiation
    results = run_ebtelplusplus(base_config)
    quantities = [
        'electron_temperature',
        'ion_temperature',
        'density',
        'electron_pressure',
        'ion_pressure',
        'velocity'
    ]
    for q in quantities:
        assert q in results
