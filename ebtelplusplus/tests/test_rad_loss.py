"""
Tests for the different kinds of radiative loss functions
"""
import astropy.units as u
import pytest

import ebtelplusplus

from ebtelplusplus.models import HeatingModel, PhysicsModel, TriangularHeatingEvent


@pytest.mark.parametrize('radiation', ['power_law', 'variable', 'coronal', 'photospheric'])
def test_rad_loss_options(radiation):
    # Just a smoke test to make sure the new radiative loss options work
    physics = PhysicsModel(radiative_loss=radiation)
    heating = HeatingModel(
        background=1e-6*u.Unit('erg cm-3 s-1'),
        events=[TriangularHeatingEvent(0.0*u.s, 200*u.s, 0.1*u.Unit('erg cm-3 s-1'))]
    )
    results = ebtelplusplus.run(5e3*u.s, 40*u.Mm, heating, physics=physics)
    quantities = [
        'electron_temperature',
        'ion_temperature',
        'density',
        'electron_pressure',
        'ion_pressure',
        'velocity'
    ]
    for q in quantities:
        assert getattr(results, q) is not None
