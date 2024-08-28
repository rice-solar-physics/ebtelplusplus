"""
High level interface to running an ebtel++ simulation
"""
import astropy.units as u

from astropy.utils.data import get_pkg_data_path
from collections import namedtuple

import ebtelplusplus._low_level

from ebtelplusplus.models import DemModel, PhysicsModel, SolverModel

__all__ = ["run", "build_configuration"]


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


@u.quantity_input
def run(total_time: u.s, loop_length: u.cm, heating, physics=None, solver=None, dem=None):
    """
    Run an ebtel++ simulation

    Parameters
    ----------
    total_time: `~astropy.units.Quantity`
        Total duration of the simulation
    loop_length: `~astropy.units.Quantity`
        Loop half length
    heating_model: `~ebtelplusplus.models.HeatingModel`
        Configuration of heating model
    physics_model: `~ebtelplusplus.models.PhysicsModel`, optional
        Configuration parameters related to the physics of the simulation
    solver_model: `~ebtelplusplus.models.SolverModel`, optional
        Configuration parameters related to the numerical solver
    dem_model: `~ebtelplusplus.models.DemModel`, optional
        Configuration parameters related to the DEM calculation

    Returns
    -------
    results: `EbtelResult`
        Data structure holding the results of the ebtel++ calculation
    """
    if solver is None:
        solver = SolverModel()
    if physics is None:
        physics = PhysicsModel()
    if dem is None:
        dem = DemModel()
    config = build_configuration(total_time, loop_length, solver, physics, dem, heating)
    results = ebtelplusplus._low_level.run(config)
    return EbtelResult(**{k: u.Quantity(v, _UNITS_MAPPING[k]) for k,v in results.items()})


def build_configuration(total_time, loop_length, solver, physics, dem, heating):
    """
    Helper function for building a dictionary of ebtel++ configuration options
    """
    return {
        'total_time': total_time.to_value('s'),
        'loop_length': loop_length.to_value('cm'),
        'radiation_data_dir': get_pkg_data_path('radiation', package='ebtelplusplus.data'),
        **solver.to_dict(),
        **physics.to_dict(),
        **dem.to_dict(),
        'heating': heating.to_dict(),
    }
