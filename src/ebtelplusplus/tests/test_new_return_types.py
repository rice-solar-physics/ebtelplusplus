import os
from collections import OrderedDict
import copy
import pytest
import tempfile

import ebtelplusplus

from .util import write_xml


@pytest.fixture
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
        'force_single_fluid': True,
        'use_c1_loss_correction': True,
        'use_c1_grav_correction': True,
        'use_flux_limiting': True,
        'calculate_dem': False,
        'save_terms': False,
        'use_adaptive_solver': True,
        'radiation': 'power_law',
        'adaptive_solver_error': 1e-6,
        'adaptive_solver_safety': 0.5,
        'c1_cond0': 2.0,
        'c1_rad0': 0.6,
        'helium_to_hydrogen_ratio': 0.075,
        'surface_gravity': 1.0,
        'heating': OrderedDict({
            'partition': 0.5,
            'background': 1e-6,
            'events': [],
        }),
    }
    return base_config


def test_return_data_struct(base_config):
    config = copy.deepcopy(base_config)
    with tempfile.TemporaryDirectory() as tmpdir:
        config_filename = os.path.join(tmpdir, 'ebtelplusplus.tmp.xml')
        results_filename = os.path.join(tmpdir, 'ebtelplusplus.tmp')
        config['output_filename'] = results_filename
        write_xml(config, config_filename)
        res = ebtelplusplus.run(config_filename)

    assert False
    