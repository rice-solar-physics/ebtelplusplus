"""
Compare results of adaptive and static solvers
"""

import astropy.units as u
import numpy as np
import pytest

import ebtelplusplus

from ebtelplusplus.models import (
    HeatingModel,
    PhysicsModel,
    SolverModel,
    TriangularHeatingEvent,
)


@pytest.fixture()
def physics_model():
    return PhysicsModel(
        saturation_limit=1/6,
        c1_conduction=2.0,
    )


@pytest.fixture()
def adaptive_solver():
    return SolverModel(tau_max=10*u.s,
                       use_adaptive_solver=True)


@pytest.fixture()
def static_solver():
    return SolverModel(tau_max=10*u.s,
                       use_adaptive_solver=False)


@pytest.fixture()
def heating_model():
    return HeatingModel(
        background=1e-6*u.Unit('erg cm-3 s-1'),
        partition=1.0,
        events=[TriangularHeatingEvent(rise_start=0.0*u.s,
                                       duration=200*u.s,
                                       rate=0.1*u.Unit('erg cm-3 s-1'))]
    )


@pytest.fixture()
def adaptive_results(physics_model, adaptive_solver, heating_model):
    return ebtelplusplus.run(5e3*u.s,
                             40*u.Mm,
                             heating_model,
                             physics=physics_model,
                             solver=adaptive_solver)


@pytest.fixture()
def static_results(physics_model, static_solver, heating_model):
    return ebtelplusplus.run(5e3*u.s,
                             40*u.Mm,
                             heating_model,
                             physics=physics_model,
                             solver=static_solver)


@pytest.mark.parametrize(('name', 'atol'), [
    ('electron_temperature', 1e4),
    ('ion_temperature', 1e4),
    ('density', 1e7),
    ('electron_pressure', 1e-4),
    ('ion_pressure', 1e-4),
    ('velocity', 1e4),
])
def test_quantities_equal_adaptive_static(adaptive_results, static_results, name, atol):
    t_static = static_results.time.to_value('s')
    t_adaptive = adaptive_results.time.to_value('s')
    adapt_interp = np.interp(t_static, t_adaptive, getattr(adaptive_results, name).value)
    # NOTE: Skip the first 5 steps b/c there is always one anomalous point that gives
    # an error > 10%; due to static case not rising fast enough
    assert np.allclose(adapt_interp[5:],
                       getattr(static_results, name)[5:].to_value(getattr(adaptive_results,name).unit),
                       rtol=1e-2,
                       atol=atol)


@pytest.mark.parametrize('bad_heating', [
    HeatingModel(background=-1e-5*u.Unit('erg cm-3 s-1')),
    HeatingModel(background=0*u.Unit('erg cm-3 s-1')),
    HeatingModel(background=1e-15*u.Unit('erg cm-3 s-1')),
])
def test_insufficient_heating(bad_heating, static_solver,):
    with pytest.raises(RuntimeError):
        _ = ebtelplusplus.run(5e3*u.s, 40*u.Mm, bad_heating, solver=static_solver)


@pytest.mark.parametrize('use_adaptive_solver', [True, False])
def test_NaNs_in_solver(use_adaptive_solver):
    solver = SolverModel(use_adaptive_solver=use_adaptive_solver)
    heating = HeatingModel(
        partition=1.0,
        background=1e-6*u.Unit('erg cm-3 s-1'),
        events=[TriangularHeatingEvent(0.0*u.s, 200*u.s, -10*u.Unit('erg cm-3 s-1'))]
    )
    with pytest.raises(RuntimeError):
        _ = ebtelplusplus.run(5e3*u.s, 40*u.Mm, heating, solver=solver)


@pytest.mark.parametrize(('A_c', 'A_0', 'A_tr'), [
    (3, 1, 1),
    (3, 2, 1),
    (1, 1, 1),
])
def test_area_expansion(A_c, A_0, A_tr, adaptive_solver, heating_model):
    # This is just a smoke test for the area expansion functionality
    physics = PhysicsModel(
        saturation_limit=1/6,
        c1_conduction=2.0,
        loop_length_ratio_tr_total=0.15,
        area_ratio_tr_corona=A_tr/A_c,
        area_ratio_0_corona=A_0/A_c,
    )
    results = ebtelplusplus.run(5e3*u.s,
                                40*u.Mm,
                                heating_model,
                                physics=physics,
                                solver=adaptive_solver)
    vars = [
        'electron_temperature',
        'ion_temperature',
        'density',
        'electron_pressure',
        'ion_pressure',
        'velocity',
        'time',
    ]
    for v in vars:
        assert getattr(results, v) is not None
        assert not np.any(np.isnan(getattr(results, v)))
