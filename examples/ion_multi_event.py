"""
Ion heating: Multiple square events
===================================
In this example, only the ions are heated by multiple square pulses each lasting
200 seconds with heating rates chosen from a uniform distribution.
"""
import astropy.units as u
from astropy.visualization import quantity_support
import matplotlib.pyplot as plt
import numpy as np

import ebtelplusplus
from ebtelplusplus.models import DemModel, HeatingModel, SquareHeatingEvent

quantity_support()


##################################################
# Set up a series of 5 square heating events each
# lasting 200 s with the first starting at 1500 s.
# The rates are chosen from a uniform distribution
# between 0.001 and 0.1 erg per cubic centimeter
# per second.
events = []
time_start = 1500 * u.s
for i in range(5):
    events.append(SquareHeatingEvent(
        (i+1)*time_start,
        200*u.s,
        np.random.uniform(0.001, 0.1)*u.Unit('erg cm-3 s-1')
    ))

##################################################
# Let's add this to a heating model in which we
# partition all of the energy into the ions
heating = HeatingModel(background=3.5e-5*u.Unit('erg cm-3 s-1'),
                       partition=0,
                       events=events)

##################################################
# Now run the simulation for a 40 Mm loop lasting
# a total of 3 h. We'll also specify that we
# want to compute the DEM
result = ebtelplusplus.run(3*u.h,
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
