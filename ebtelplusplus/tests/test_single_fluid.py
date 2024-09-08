"""
Test whether single-fluid option keeps electron and ion temperatures equal
"""
import astropy.units as u
import pytest

import ebtelplusplus

from ebtelplusplus.models import HeatingEvent, HeatingModel, PhysicsModel


@pytest.fixture()
def physics_model():
    return PhysicsModel(force_single_fluid=True,
                        saturation_limit=1/6,
                        c1_conduction=2.0,)


@pytest.mark.parametrize('event', [
    HeatingEvent(0*u.s, 250*u.s, 1000*u.s, 2000*u.s, 0.005*u.Unit('erg cm-3 s-1')),  # gentle case
    HeatingEvent(0*u.s, 100*u.s, 100*u.s, 200*u.s, 0.1*u.Unit('erg cm-3 s-1')),  # more impulsive case
])
def test_single_fluid_gentle(physics_model, event):
    heating = HeatingModel(background=1e-6*u.Unit('erg cm-3 s-1'),
                           partition=0.5,
                           events=[event])
    results = ebtelplusplus.run(5e3*u.s, 40*u.Mm, heating, physics=physics_model)
    assert u.allclose(results.electron_temperature, results.ion_temperature, rtol=1e-10)
    assert u.allclose(results.electron_pressure, results.ion_pressure, rtol=1e-10)
