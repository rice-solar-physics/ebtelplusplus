"""
Compare output of EBTEL IDL and ebtel++
"""
import astropy.units as u
import pytest

import ebtelplusplus

from ebtelplusplus.models import (
    HeatingModel,
    PhysicsModel,
    SolverModel,
    TriangularHeatingEvent,
)

from .helpers import plot_comparison, read_idl_test_data

# Tolerated error between IDL and C++ results
RTOL = 0.01


@pytest.fixture()
def physics_model():
    return PhysicsModel(
        force_single_fluid=True,
    )


@pytest.fixture()
def solver_model():
    return SolverModel(
        tau=0.5*u.s,
        use_adaptive_solver=False,
    )


@pytest.fixture()
def heating_model():
    return HeatingModel(
        background=1e-6*u.erg/(u.cm**3*u.s),
        partition=0.5,
        events=[TriangularHeatingEvent(0.0*u.s, 200*u.s, 0.1*u.Unit('erg cm-3 s-1'))]
    )


@pytest.mark.xfail()
def test_compare_idl_single_event(physics_model,
                                  solver_model,
                                  heating_model,
                                  ebtel_idl_path,
                                  plot_idl_comparisons):
    total_time = 5e3 * u.s
    loop_length = 40 * u.Mm
    r_cpp = ebtelplusplus.run(total_time,
                              loop_length,
                              heating_model,
                              physics=physics_model,
                              solver=solver_model)

    r_idl = read_idl_test_data(
        'idl_single_event.txt',
        ebtel_idl_path,
        r_cpp.inputs,
    )
    if plot_idl_comparisons:
        plot_comparison(r_cpp, r_idl)
    assert u.allclose(r_cpp.electron_temperature, r_idl['temperature'], rtol=RTOL)
    assert u.allclose(r_cpp.ion_temperature, r_idl['temperature'], rtol=RTOL)
    assert u.allclose(r_cpp.density, r_idl['density'], rtol=RTOL)
    assert u.allclose(r_cpp.electron_pressure+r_cpp.ion_pressure,
                      r_idl['pressure'],
                      rtol=RTOL)
    assert u.allclose(r_cpp.velocity, r_idl['velocity'], rtol=RTOL)


@pytest.mark.parametrize(('A_c', 'A_0', 'A_tr'), [
    pytest.param(3, 1, 1, marks=pytest.mark.xfail),
    pytest.param(3, 2, 1, marks=pytest.mark.xfail),
    pytest.param(1, 1, 1, marks=pytest.mark.xfail),
])
def test_compare_idl_area_expansion(
        A_c, A_0, A_tr, solver_model, heating_model, ebtel_idl_path, plot_idl_comparisons
    ):
    physics_model = PhysicsModel(force_single_fluid=True,
                                 loop_length_ratio_tr_total=0.15,
                                 area_ratio_tr_corona=A_tr/A_c,
                                 area_ratio_0_corona=A_0/A_c)
    total_time = 5e3 * u.s
    loop_length = 40 * u.Mm
    r_cpp = ebtelplusplus.run(total_time,
                              loop_length,
                              heating_model,
                              physics=physics_model,
                              solver=solver_model)
    r_idl = read_idl_test_data(f'idl_area_expansion_{A_c=}_{A_0=}_{A_tr=}.txt', ebtel_idl_path, r_cpp.inputs)
    if plot_idl_comparisons:
        plot_comparison(r_cpp, r_idl)
    assert u.allclose(r_cpp.electron_temperature, r_idl['temperature'], rtol=RTOL)
    assert u.allclose(r_cpp.ion_temperature, r_idl['temperature'], rtol=RTOL)
    assert u.allclose(r_cpp.density, r_idl['density'], rtol=RTOL)
    assert u.allclose(r_cpp.electrong_pressure+r_cpp.ion_pressure, r_idl['pressure'], rtol=RTOL)
    assert u.allclose(r_cpp.velocity, r_idl['velocity'], rtol=RTOL)
