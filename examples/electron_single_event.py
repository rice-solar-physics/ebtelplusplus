"""
Heating only the Electrons with One Triangular Event
====================================================
In this example, only the electrons are heated by a single triangular pulse lasting 500 seconds
and injecting 10 ergs per cubic centimeter into the loop plasma.
"""
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

from astropy.visualization import quantity_support

import ebtelplusplus

from ebtelplusplus.models import DemModel, HeatingModel, TriangularHeatingEvent

quantity_support()


##################################################
# Set up a single heating event that is triangular
# in shape such that it rises linearly for 250 s
# and then falls linearly for 250 s starting at
# the beginning of the simulation (0 s).
event = TriangularHeatingEvent(0.0*u.s,
                               500*u.s,
                               0.04*u.Unit('erg cm-3 s-1'))

##################################################
# Let's add this to a heating model in which we
# partition all of the energy into the electrons
heating = HeatingModel(background=3.5e-5*u.Unit('erg cm-3 s-1'),
                       partition=1.0,
                       events=[event])

##################################################
# Now run the simulation for a 40 Mm loop lasting
# a total of 5000 s. We'll also specify that we
# want to compute the DEM
result = ebtelplusplus.run(5e3*u.s,
                           40*u.Mm,
                           heating,
                           dem=DemModel(calculate_dem=True))

##################################################
# Let's visualize the heating profile, temperature,
# and density as a function of time.
fig, axes = plt.subplots(3, 1, sharex=True)
axes[0].plot(result.time, result.heat)
axes[1].plot(result.time, result.electron_temperature, label='electron')
axes[1].plot(result.time, result.ion_temperature, label='ion')
axes[2].plot(result.time, result.density)
axes[1].legend()

##################################################
# Finally, let's visualize the DEM distribution.
# We'll first time-average each component over the
# duration of the simulation.
delta_t = np.gradient(result.time)
dem_avg_total = np.average(result.dem_tr+result.dem_corona,
                           axis=0,
                           weights=delta_t)
dem_avg_tr = np.average(result.dem_tr,
                        axis=0,
                        weights=delta_t)
dem_avg_corona = np.average(result.dem_corona,
                            axis=0,
                            weights=delta_t)

##################################################
# And now we can plot each component
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(result.dem_temperature, dem_avg_total, label='Total')
ax.plot(result.dem_temperature, dem_avg_tr, label='TR')
ax.plot(result.dem_temperature, dem_avg_corona, label='Corona')
ax.set_xlim([10**(4.5), 10**(7.5)]*u.K)
ax.set_ylim([10**(20.0), 10**(23.5)]*u.Unit('cm-5 K-1'))
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend()

plt.show()
