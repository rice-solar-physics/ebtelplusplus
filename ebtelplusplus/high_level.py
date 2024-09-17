"""
High level interface to running an ebtel++ simulation
"""
import astropy.units as u
import dataclasses
import typing

from astropy.utils.data import get_pkg_data_path

import ebtelplusplus._low_level

from ebtelplusplus.models import DemModel, HeatingModel, PhysicsModel, SolverModel

__all__ = ["run", "EbtelResult"]


@u.quantity_input
def run(total_time: u.s, loop_length: u.cm, heating=None, physics=None, solver=None, dem=None):
    """
    Run an ebtelplusplus simulation

    To run an ebtelplusplus simulation, at a minimum, you must specify a total simulation time
    and a loop half length. Additional input parameters are configured through high-level objects
    as listed below. Most users will likely want to modify the ``heating`` input by constructing
    a `~ebtelplusplus.models.HeatingModel` that specifies the amount of energy injected into the
    simulation as a function of time. Additional inputs are documented in the docstring of each
    respective input model.

    Parameters
    ----------
    total_time: `~astropy.units.Quantity`
        Total duration of the simulation
    loop_length: `~astropy.units.Quantity`
        Loop half length
    heating: `~ebtelplusplus.models.HeatingModel`, optional
        Configuration of heating model
    physics: `~ebtelplusplus.models.PhysicsModel`, optional
        Configuration parameters related to the physics of the simulation
    solver: `~ebtelplusplus.models.SolverModel`, optional
        Configuration parameters related to the numerical solver
    dem: `~ebtelplusplus.models.DemModel`, optional
        Configuration parameters related to the DEM calculation

    Returns
    -------
    results: `~ebtelplusplus.EbtelResult`
        Data structure holding the results of the ebtelplusplus calculation
    """
    if heating is None:
        heating = HeatingModel()
    if solver is None:
        solver = SolverModel()
    if physics is None:
        physics = PhysicsModel()
    if dem is None:
        dem = DemModel()
    config = _build_configuration(total_time, loop_length, heating, physics, solver, dem)
    results = ebtelplusplus._low_level.run(config)
    return EbtelResult(results, config)


@u.quantity_input
def _build_configuration(total_time:u.s, loop_length:u.cm, heating, physics, solver, dem):
    """
    Helper function for building a dictionary of ebtelplusplus configuration options

    .. note::
        It is not necessary to use this function if you are running a simulation using
        `~ebtelplusplus.run`. However, it may be useful to use this function directly
        if you want to store the resulting inputs a dictionary.

    Parameters
    ----------
    total_time: `~astropy.units.Quantity`
        Total duration of the simulation
    loop_length: `~astropy.units.Quantity`
        Loop half length
    heating_model: `~ebtelplusplus.models.HeatingModel`
        Configuration of heating model
    physics_model: `~ebtelplusplus.models.PhysicsModel`
        Configuration parameters related to the physics of the simulation
    solver_model: `~ebtelplusplus.models.SolverModel`
        Configuration parameters related to the numerical solver
    dem_model: `~ebtelplusplus.models.DemModel`
        Configuration parameters related to the DEM calculation
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


@dataclasses.dataclass
class EbtelResult:
    """
    Result of an ebtelplusplus simulation

    .. note::

        This class is not meant to be instantiated directly.
        Rather, it is meant to be returned by `~ebtelplusplus.run`.

    Parameters
    ----------
    time: `~astropy.units.Quantity`
        Simulation time at each time step. Has shape ``(N,)``
    electron_temperature: `~astropy.units.Quantity`
        The temperature of the electrons at each time step
    ion_temperature: `~astropy.units.Quantity`
        The temperature of the ions at each time step
    density: `~astropy.units.Quantity`
        The density of the electrons and ions at each time step
    electron_pressure: `~astropy.units.Quantity`
        The pressure of the electrons at each time step
    ion_pressure: `~astropy.units.Quantity`
        The pressure of the ions at each time step
    total_pressure: `~astropy.units.Quantity`
        The total pressure at each time step. This is the sum of
        the electron and ion pressures and is equivalent to the
        pressure returned for single-fluid models.
    velocity: `~astropy.units.Quantity`
        The velocity at each timestep. This parameter should be
        used cautiously as the EBTEL does not model the actual
        transport of material through the spatial extent of the
        loop.
    heat: `~astropy.units.Quantity`
        The energy released into the loop at each time step.
        This is the amount of heating as specified in the
        input `~ebtelplusplus.models.HeatingModel`.
    electron_thermal_conduction: `~astropy.units.Quantity`
        Electron thermal conductive flux at the TR-corona interface
        at each time step.
    ion_thermal_conduction: `~astropy.units.Quantity`
        Ion thermal conductive flux at the TR-corona interface
        at each time step.
    radiative_loss: `~astropy.units.Quantity`
        Energy lost due to radiation at each time step
    tr_corona_radiative_loss_ratio: `~astropy.units.Quantity`
        The ratio of the average radiative losses in the TR and the
        corona. In the EBTEL model, this is often referred to as :math:`c_1`.
    dem_temperature: `~astropy.units.Quantity`
        The temperature grid on which the differential emission measure
        distribution (DEM) is computed. Has shape ``(M,)``. These can be
        interpreted as the centers of the temperature bins.
    dem_tr: `~astropy.units.Quantity`
        The TR DEM. Has shape ``(N,M)``. For more details of how this is
        calculated, see Section 3 and the Appendix of :cite:t:`klimchuk_highly_2008`.
    dem_corona: `~astropy.units.Quantity`
        The coronal DEM. Has shape ``(N,M)``.
    inputs: `dict`
        All model inputs used to run the simulation.
    """
    inputs: dict
    time: u.Quantity[u.s]
    electron_temperature: u.Quantity[u.K]
    ion_temperature: u.Quantity[u.K]
    density: u.Quantity[u.cm**(-3)]
    electron_pressure: u.Quantity[u.dyne*u.cm**(-2)]
    ion_pressure: u.Quantity[u.dyne*u.cm**(-2)]
    velocity: u.Quantity[u.cm/u.s]
    heat: u.Quantity[u.erg*u.cm**(-3)*u.s**(-1)]
    electron_thermal_conduction: u.Quantity[u.erg*u.cm**(-2)*u.s**(-1)]
    ion_thermal_conduction: u.Quantity[u.erg*u.cm**(-2)*u.s**(-1)]
    radiative_loss: u.Quantity[u.erg*u.cm**(-3)*u.s**(-1)]
    tr_corona_radiative_loss_ratio: u.Quantity[u.dimensionless_unscaled]
    dem_temperature: u.Quantity[u.K] = None
    dem_tr: u.Quantity[u.K**(-1)*u.cm**(-5)] = None
    dem_corona: u.Quantity[u.K**(-1)*u.cm**(-5)] = None
    # derived fields
    total_pressure: u.Quantity[u.dyne*u.cm**(-2)] = dataclasses.field(init=None)

    def __init__(self, results, inputs):
        # Use type hinting above to assign units to outputs; this is so
        # we only have to specify the units of the outputs in one place
        type_hints = type(self).__annotations__
        for k, v in results.items():
            _, unit = typing.get_args(type_hints[k])
            setattr(self, k, u.Quantity(v, unit))
        self.inputs = inputs

    def __post_init__(self):
        self.total_pressure = self.electron_pressure + self.ion_pressure
