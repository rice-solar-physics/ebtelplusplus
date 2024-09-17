"""
Classes for configuring inputs to `ebtelplusplus.run`.
All model inputs have default values such that each model only need be instantiated
using the values of the parameters a user wants to change.
The default values of each input are listed below.
"""
import astropy.units as u
import dataclasses

__all__ = [
    'PhysicsModel',
    'SolverModel',
    'DemModel',
    'HeatingModel',
    'HeatingEvent',
    'TriangularHeatingEvent',
    'SquareHeatingEvent',
]


@dataclasses.dataclass
class PhysicsModel:
    """
    ebtelplusplus input parameters related to the physics of the simulation

    Parameters
    ----------
    saturation_limit: `float`
        Flux limiter constant. See Eq. 21 of :cite:t:`klimchuk_highly_2008`
    force_single_fluid: `bool`
        If true, electron and ion populations forced into equilibrium.
    c1_conduction: `float`
        Nominal value of :math:`c_1` during the conductive cooling phase
        See Appendix A of :cite:t:`barnes_inference_2016`.
    c1_radiation: `float`
        Nominal value of :math:`c_1` during radiative phase. See Eq. 16 of
        :cite:t:`cargill_enthalpy-based_2012`.
    use_c1_gravity_correction: `bool`
        Use correction in Eq. 12 of :cite:t:`cargill_enthalpy-based_2012`.
    use_c1_radiation_correction: `bool`
        Use correction in Eq. 16 of :cite:t:`cargill_enthalpy-based_2012`.
    surface_gravity: `float`
        Surface gravity in units of solar surface gravity. Should be set
        to 1.0 unless using for extra-solar cases.
    helium_to_hydrogren_ratio: `float`
        Ratio of helium to hydrogen abundance; used in correction to
        ion mass and ion equation of state.
    radiative_loss: `str`
        The kind of radiative loss function to use. Must be either "power_law" (to use radiative losses of :cite:t:`klimchuk_highly_2008`,
        "coronal" (to use radiative losses computed with coronal abundances),
        "photospheric" (to use radiative losses computed with photospheric
        abundances), or "variable" (to vary the radiative loss function from
        coronal to photospheric as a function of density and temperature).
    loop_length_ratio_tr_total: `float`
        Ratio between the length of the transition region and the total loop
        length. For a transition region of finite length, a typical value of
        0.15 is used :cite:p:`cargill_static_2022`.
    area_ratio_tr_corona: `float`
        Ratio between the cross-sectional area averaged over the transition
        region and averaged over the corona
    area_ratio_0_corona: `float`
        Ratio between the cross-sectional area at the TR-corona boundary and
        the cross-sectional area averaged over the corona
    """
    saturation_limit: float = None
    force_single_fluid: bool = False
    c1_conduction: float = 6.0
    c1_radiation: float = 0.6
    use_c1_gravity_correction: bool = True
    use_c1_radiation_correction: bool = True
    surface_gravity: float = 1.0
    helium_to_hydrogen_ratio: float = 0.075
    radiative_loss: str = 'power_law'
    loop_length_ratio_tr_total: float = 0.0
    area_ratio_tr_corona: float = 1.0
    area_ratio_0_corona: float = 1.0
    # derived fields
    use_flux_limiting: bool = dataclasses.field(init=False)

    def __post_init__(self):
        self.use_flux_limiting = self.saturation_limit is not None

    def to_dict(self):
        config = dataclasses.asdict(self)
        if self.saturation_limit is None:
            config['saturation_limit'] = 1.0
        return config


@dataclasses.dataclass
class SolverModel:
    """
    ebtelplusplus input parameters related to the numerical solver

    Parameters
    ----------
    tau: `~astropy.units.Quantity`
        Time step if using adaptive solver, the initial timestep
    tau_max: `~astropy.units.Quantity`
        Maximum allowed time step when using adaptive solver
    use_adaptive_solver: `bool`
        If true, use an adaptive timestepping routine
    adaptive_solver_error: `float`
        Allowed truncation error in adaptive timestep routine
    adaptive_solver_safety: `float`
        Refinement factor, between 0 and 1, used if timestep becomes too large
        and solution contains NaNs. Especially important for short,
        infrequently heated loops. Also controls decreases in timestep due to
        thermal conduction timestep.
    """
    tau: u.Quantity[u.s] = 1.0*u.s
    tau_max: u.Quantity[u.s] = 10*u.s
    use_adaptive_solver: bool = True
    adaptive_solver_error: float = 1e-6
    adaptive_solver_safety: float = 0.5

    def to_dict(self):
        config = dataclasses.asdict(self)
        config['tau'] = config['tau'].to_value('s')
        config['tau_max'] = config['tau_max'].to_value('s')
        return config


@dataclasses.dataclass
class DemModel:
    """
    ebtelplusplus input parameters related to differential emission measure (DEM) calculation.

    Optionally, ebtelplusplus can can also calculate the differential emission
    measure (DEM) in both the transition region and the corona. See sections 2.2
    and 3 of :cite:t:`klimchuk_highly_2008` for the details of this
    calculation. Note that this will result in much longer computation times.

    Parameters
    ----------
    calculate_dem: `bool`
        If true, calculate the coronal and transition region DEM
    use_new_tr_method: `bool`
        If true, the transition region DEM is calculated using the method
        outlined in section 3 (the appendix) of :cite:t:`klimchuk_highly_2008`.
    temperature_bins: `int`
        Number of bins to use when calculating the DEM
    temperature_min: `~astropy.units.Quantity`
        Lower bound on the temperature range for the DEM calculation
    temperature_max: `~astropy.units.Quantity`
        Upper bound on the temperature range for the DEM calculation
    """
    calculate_dem: bool = False
    use_new_tr_method: bool = True
    temperature_bins: int = 451
    temperature_min: u.Quantity[u.K] = 10**4*u.K
    temperature_max: u.Quantity[u.K] = 10**8.5*u.K

    def to_dict(self):
        return {
            'calculate_dem': self.calculate_dem,
            'dem_use_new_tr_method': self.use_new_tr_method,
            'dem_temperature_bins': self.temperature_bins,
            'dem_temperature_min': self.temperature_min.to_value('K'),
            'dem_temperature_max': self.temperature_max.to_value('K'),
        }


class HeatingEvent:
    """
    Single heating event

    Each event has a linear rise phase, a constant phase, and
    a linear decay phase. Using this format, it is easy to specify either
    symmetric or asymmetric events of many different shapes.

    Parameters
    ----------
    rise_start: `~astropy.units.Quantity`
        Time at which the heating event starts
    rise_end: `~astropy.units.Quantity`
        Time at which the rise phase stops (and the constant phase starts)
    decay_start: `~astropy.units.Quantity`
        Time at which the decay phase starts (and the constant phase stops)
    decay_end: `~astropy.units.Quantity`
        Time at which the decay phase and the event ends
    rate: `~astropy.units.Quantity`
        The maximum heating rate of the event

    See also
    --------
    TriangularHeatingEvent
    SquareHeatingEvent
    """

    @u.quantity_input
    def __init__(self, rise_start: u.s, rise_end: u.s, decay_start: u.s, decay_end: u.s, rate: u.Unit('erg cm-3 s-1')):
        self.rise_start = rise_start
        self.rise_end = rise_end
        self.decay_start = decay_start
        self.decay_end = decay_end
        self.rate = rate

    def to_dict(self):
        return {
            'rise_start': self.rise_start.to_value('s'),
            'rise_end': self.rise_end.to_value('s'),
            'decay_start': self.decay_start.to_value('s'),
            'decay_end': self.decay_end.to_value('s'),
            'rate': self.rate.to_value('erg cm-3 s-1'),
        }


class TriangularHeatingEvent(HeatingEvent):
    """
    A single event with a linear rise phase immediately followed by a
    linear decay phase of equal duration

    Parameters
    ----------
    rise_start: `~astropy.units.Quantity`
        Time at which the heating event starts
    duration: `~astropy.units.Quantity`
        Total duration of the event
    rate: `~astropy.units.Quantity`
        The maximum heating rate of the event
    """

    def __init__(self, rise_start, duration, rate):
        super().__init__(
            rise_start,
            rise_start+duration/2,
            rise_start+duration/2,
            rise_start+duration,
            rate,
        )


class SquareHeatingEvent(HeatingEvent):
    """
    A single event with no rise or decay phase and only a constant phase

    Parameters
    ----------
    rise_start: `~astropy.units.Quantity`
        Time at which the heating event starts
    duration: `~astropy.units.Quantity`
        Total duration of the event
    rate: `~astropy.units.Quantity`
        The maximum heating rate of the event
    """

    def __init__(self, rise_start, duration, rate):
        super().__init__(
            rise_start,
            rise_start,
            rise_start+duration,
            rise_start+duration,
            rate,
        )


@dataclasses.dataclass
class HeatingModel:
    """
    ebtelplusplus input parameters for time-dependent heating

    The ebtelplusplus time-dependent heating model is parameterized by a
    series of discrete events (`HeatingEvent`) combined with a constant
    background heating rate. Furthermore, this energy can be injected into
    either the electrons or the ions or some admixture of the two.

    Parameters
    ----------
    background: `~astropy.units.Quantity`
        Constant background heating rate; primarily used to keep the loop
        from reaching unphysical temperatures
    partition: `float`
        Partition of heating between electrons and ions, between 0 and 1;
        1 is pure electron heating, 0 pure ion heating
    events: `list`
        List of `HeatingEvent` objects that parameterize the energy injected
        into the loop by a series of discrete heating events.
    """
    background: u.Quantity[u.erg/(u.cm**3*u.s)] = 1e-6 * u.Unit('erg cm-3 s-1')
    partition: float = 0.5
    events: list[HeatingEvent] = None

    def __post_init__(self):
        if self.events is None:
            self.events = []

    def to_dict(self):
        return {
            'background': self.background.to_value('erg cm-3 s-1'),
            'partition': self.partition,
            'events': [event.to_dict() for event in self.events],
        }
