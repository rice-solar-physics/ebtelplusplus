"""
Compare output of EBTEL IDL and ebtel++
"""
import astropy.units as u
import pytest

from scipy.interpolate import interp1d

import ebtelplusplus

from ebtelplusplus.models import (
    HeatingModel,
    PhysicsModel,
    SolverModel,
    TriangularHeatingEvent,
)

from .helpers import read_hydrad_test_data


@pytest.fixture()
def solver_model():
    return SolverModel(tau=0.1*u.s,
                       tau_max=10*u.s,
                       adaptive_solver_error=1e-8)


@pytest.mark.parametrize('tau', [200.0*u.s, 500.0*u.s])
@pytest.mark.parametrize('heating_type', ['single', 'electron', 'ion'])
def test_compare_hydrad_single_event_peak_values(tau, heating_type, solver_model):
    if heating_type == 'single':
        force_single_fluid = True
        partition = 0.5
    elif heating_type == 'electron':
        force_single_fluid = False
        partition = 1.0
    elif heating_type == 'ion':
        force_single_fluid = False
        partition = 0.0
    total_energy = 10.0 * u.erg / u.cm**3
    heating_model = HeatingModel(
        background=3.5e-5*u.Unit('erg cm-3 s-1'),
        partition=partition,
        events=[TriangularHeatingEvent(0*u.s, tau, 2 * total_energy/tau)]
    )
    physics_model = PhysicsModel(saturation_limit=1,
                                 force_single_fluid=force_single_fluid)
    r_ebtel = ebtelplusplus.run(5e3*u.s,
                                40*u.Mm,
                                heating_model,
                                physics=physics_model,
                                solver=solver_model,)
    r_hydrad = read_hydrad_test_data(tau.to_value('s'), heating_type)
    # Require all quantities at the peak to be <20% in accordance with the comparisons
    # done in Barnes et al. (2016a)
    for name in ['electron_temperature', 'ion_temperature', 'density']:
        f_interp = interp1d(r_ebtel.time.to_value('s'),
                            getattr(r_ebtel, name).to_value(r_hydrad[name].unit),
                            fill_value='extrapolate')
        x_ebtel_interp = u.Quantity(f_interp(r_hydrad['time'].to_value('s')), r_hydrad[name].unit)
        i_peak = r_hydrad[name].argmax()
        assert u.allclose(x_ebtel_interp[i_peak], r_hydrad[name][i_peak], rtol=0.20)
