"""
Models for configuring ebtel++ runs
"""
import dataclasses

import astropy.units as u

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
    ebtel++ input parameters related to the physics of the simulation

    Parameters
    ----------
    saturation_limit
    force_single_fluid
    c1_conduction
    c1_radiation
    surface_gravity
    helium_to_hydrogren_ratio
    radiative_loss
    loop_length_ratio_tr_total
    area_ratio_tr_corona
    area_ratio_0_corona
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
    ebtel++ input parameters related to the numerical solver

    Parameters
    ----------
    tau
    tau_max
    use_adaptive_solver
    adaptive_solver_error
    adaptive_solver_safety
    """
    tau: u.Quantity[u.s] = 1.0*u.s
    tau_max: u.Quantity[u.s] = 1e300*u.s
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
    ebtel++ input parameters related to DEM calculation

    Parameters
    ----------
    calculate_dem
    use_new_method
    temperature_bins
    log_temperature_min
    log_temperature_max
    """
    calculate_dem: bool = False
    use_new_method: bool = True
    temperature_bins: int = 451
    temperature_min: u.Quantity[u.K] = 10**4*u.K 
    temperature_max: u.Quantity[u.K] = 10**8.5*u.K

    def to_dict(self):
        config = dataclasses.asdict(self)
        config['temperature_min'] = config['temperature_min'].to_value('K')
        config['temperature_max'] = config['temperature_max'].to_value('K')
        return config


class HeatingEvent:
    """
    Single heating event

    Each event has a linear rise phase, a constant phase, and
    a linear decay phase

    Parameters
    ----------
    rise_start
    rise_end
    decay_start
    decay_end
    rate
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
    rise_start
    duration
    rate
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
    rise_start
    duration
    rate
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
    ebtel++ input parameters for the time-dependent heating

    Parameters
    ----------
    background
    partition
    events
    """
    background: u.Quantity[u.erg/(u.cm**3*u.s)]
    partition: float = 1.0
    events: list[HeatingEvent] = None

    def __post_init__(self):
        if self.events is None:
            self.events = []

    def to_dict(self):
        return {
            'background': self.background.to_value('erg cm-3 s-1'),
            'partition': self.partition,
            'events': [event.to_dict() for event in self.events]
        }


