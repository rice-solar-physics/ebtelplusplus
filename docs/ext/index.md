#<cldoc:index>

ebtel++

# ebtel++
[![Build Status](https://travis-ci.org/rice-solar-physics/ebtelPlusPlus.svg?branch=master)](https://travis-ci.org/rice-solar-physics/ebtelPlusPlus)

[ebtel++](https://github.com/rice-solar-physics/ebtelPlusPlus) is a two-fluid version of the original enthalpy-based thermal evolution of loops (EBTEL) model implemented in C++. This code provides an enhanced description of plasma behavior above roughly 5 MK. Further generic details about EBTEL can be found in the [repository for the original IDL code](https://github.com/rice-solar-physics/EBTEL) and in the references listed below.

The EBTEL model, originally developed by [Klimchuk et al. (2008)][klimchuk_2008] efficiently computes spatially-averaged, time-dependent plasma parameters ( e.g. temperature,  pressure, density) of dynamically-heated coronal loops. It is often desirable to compute solutions for a large number of coronal loops, but the spatial and temporal scales needed to solve the full 1D-hydrodynamic equations lead to long compute times for even 1D hydrodynamic codes. EBTEL computes quick and accurate solutions for spatially-averaged quantities, allowing efficient insight into how these monolithic structures evolve. [Barnes et al. (2016)][barnes_2016] improved upon this model by extending the treatment to the two-fluid hydrodynamic equations, allowing for differential heating between electrons and ions. Experienced EBTEL users should note the slightly different model for calculating the the c1 parameter during conductive cooling (see Appendix A of [Barnes et al., 2016][barnes_2016]).

EBTEL also calculates the differential emission measure (DEM) for both the transition region and the corona. Details regarding this formulation can be found in [Klimchuk et al. (2008)][klimchuk_2008]

## Citation
If you use ebtel++ in any published work, please include the following citations and mention the use of this code.

* [Klimchuk et al. (2008)][klimchuk_2008]
* [Cargill et al. (2012a)][cargill_2012a]
* [Cargill et al. (2012b)](cargill_2012b)
* [Barnes et al. (2016)][barnes_2016]

The first three papers detail the original single-fluid EBTEL model while the last paper gives the details of the two-fluid model. In particular, the details of how the two-fluid EBTEL equations are derived can be found in the appendix of [Barnes et al. (2016)][barnes_2016].

## Dependencies
To compile ebtel++, first install the following dependencies,

* [git](https://git-scm.com/) (included with OS X, most Linux distros; [cygwin](https://www.cygwin.com/) on Windows)
* [gcc](https://gcc.gnu.org/) (at least v4.7; included with OS X, most Linux distros; [cygwin](https://www.cygwin.com/) on Windows)
* [scons](http://scons.org/) (`pip install scons` via [PyPI](https://pypi.python.org/pypi) or `conda install scons` with [Anaconda](https://www.continuum.io/downloads); requires Python 2.7)
* [boost](http://www.boost.org/) (at least v1.53; `sudo apt-get install libboost-dev` on Debian Linux, `sudo port install boost` via [Macports](https://www.macports.org/) on OS X, [from source](https://github.com/rice-solar-physics/IonPopSolver#installing-boost-from-source) on Windows)

The above list encapsulates all of the _required_ dependencies. If you'd like to use the full radiative loss function,

* [apolloDB](https://github.com/rice-solar-physics/apolloDB)

Additionally, if you'd like to run the included tests and examples, you'll need the following Python dependencies, all easily installed with [anaconda](https://www.continuum.io/downloads),

* [Python](https://www.python.org/) (2.7)
* [numpy](http://www.numpy.org/)
* [matplotlib](http://matplotlib.org/)
* [seaborn](https://stanford.edu/~mwaskom/software/seaborn/index.html)

## Installation
To download the code from GitHub and compile the code,
```Shell
$ git clone --recursive https://github.com/rice-solar-physics/ebtelPlusPlus.git
$ cd ebtelPlusPlus
$ scons
```
If the compile step fails because the compiler cannot find the appropriate headers and/or libraries using the default locations, you can use the `--includepath` and/or `--libpath` flags, respectively. For example, if you've installed the Boost headers and libraries in `/usr/local/include` and `/usr/local/lib`,
```
$ scons --includepath=/usr/local/include --libpath=/usr/local/lib
```
For more information about the available flags that can be passed to `scons`, you can run `scons -h`.

This will create an executable `bin/ebtel++.run`. To see the available command line parameters,
```Shell
$ bin/ebtel++.run --help
```
and to run the executable with the default configuration file `config/ebtel.example.cfg.xml`,
```Shell
$ bin/ebtel++.run
```

If you've installed the above Python dependencies, you can also run the tests using,
```Shell
$ scons --test
```
or run any of the three included examples,
```Shell
$ python examples/ex1.py
$ python examples/ex2.py
$ python examples/ex3.py
```

## Configuration
An ebtel++ run is configured by a single XML configuration file. The table below lists all of the input parameters needed for a run. An [example configuration file](https://github.com/rice-solar-physics/ebtelPlusPlus/blob/master/config/ebtel.example.cfg.xml) is included in the GitHub repository. The configuration file can be written manually or configured via Python as is done in the [included examples](https://github.com/rice-solar-physics/ebtelPlusPlus/tree/master/examples).

| Parameter | Type | Description |
|:---------:|:----:|:------------|
| `total_time` | float | duration of the simulation (in s)|
| `tau` | float | timestep (in s); if using adaptive solver, the initial timestep |
| `loop_length` | float | Loop half-length (in cm) |
| `saturation_limit` | float | Flux limiter, _f_ in section 2.1 of [Barnes et al. (2016)][barnes_2016] |
| `force_single_fluid` | bool | if True, electron and ion populations forced into equilibrium |
| `use_c1_loss_correction` | bool | use correction in Eq. 16 of [Cargill et al. (2012a)][cargill_2012a] |
| `use_c1_grav_correction` | bool | use correction in Eq. 12 of [Cargill et al. (2012a)][cargill_2012a] |
| `use_power_law_radiative_losses` | bool | use Eq. 3 of [Klimchuk et al. (2008)][klimchuk_2008] for radiative loss function |
| `use_flux_limiting` | bool | impose a flux limiter according to Eq. 22 of [Klimchuk et al. (2008)][klimchuk_2008] |
| `calculate_dem` | bool | if True, do the TR and coronal DEM calculation; increases compute time significantly |
| `save_terms` | bool | if True, save heat flux, c1 parameter, and radiative loss to a separate file `output_filename`+`.terms` |
| `use_adaptive_solver` | bool | if True, use adaptive timestep; significantly smaller compute times. In both cases, a Runge-Kutta Cash-Karp integration method is used (see section 16.2 of [Press et al. (1992)][press_num_recipes])  |
| `output_filename` | string | path to output file |
| `adaptive_solver_error` | float | Allowed truncation error in adaptive timestep routine |
| `adaptive_solver_safety` | float | Refinement factor, between 0 and 1, used if timestep becomes too large and solution contains NaNs. Especially important for short, infrequently heated loops. Also controls decreases in timestep due to thermal conduction timestep. Suggested value is 0.5 |
| `c1_cond0` | float | Nominal value of c1 during the conduction phase; see Appendix A of [Barnes et al. (2016)][barnes_2016] |
| `c1_rad0` | float | Nominal value of c1 during radiative phase; see Eq. 16 of [Cargill et al. (2012a)][cargill_2012a] |
| `helium_to_hydrogen_ratio` | float | Ratio of helium to hydrogen abundance; used in correction to ion mass, ion equation of state |

### Heating
The time dependent heating is configured in a separate node. It includes the following parameters,

| Parameter | Type | Description |
|:-------:|:------:|:-----------|
| `partition` | float | partition of heating between electrons and ions, between 0 and 1; 1 is pure electron heating, 0 pure ion heating |
| `background` | float | constant background heating (in ergs cm^-3 s^-1) |

The heating function is constructed by a list of discrete events and should be specified in the following way,
```XML
<events>
  <event magnitude="0.1" rise_start="0.0" rise_end="50.0" decay_start="50.0" decay_end="100.0"/>
  <event magnitude="0.05" rise_start="1000.0" rise_end="1250.0" decay_start="1350.0" decay_end="1450.0"/>
</events>
```
Here, we've configured two separate heating events. The first starts at 0 seconds, rises linearly to a maximum heating rate of 0.1 ergs per cubic centimeter per second in 50 seconds, and then immediately falls off with the event concluding at 100 seconds, i.e. a triangular heating profile. The second starts at 1000 seconds, rises to a maximum heating rate of 0.05 ergs per cubic centimeter per second in 250 seconds, is sustained at 0.05 for 100 seconds and then the event concludes at 1450 seconds.

Using this format, it is easy to specify either symmetric or asymmetric events of many different shapes. For more examples, see the [example configuration file](https://github.com/rice-solar-physics/ebtelPlusPlus/blob/master/config/ebtel.example.cfg.xml) or the included [examples](https://github.com/rice-solar-physics/ebtelPlusPlus/tree/master/examples).

### Differential Emission Measure
Optionally, ebtel++ can can also calculate the differential emission measure (DEM) in both the transition region and the corona. See sections 2.2 and 3 of [Klimchuk et al. (2008)][klimchuk_2008] for the details of this calculation. To enable this calculation, set `calculate_dem` to `True` in the configuration file (as described above). Note that this will result in much longer computation times.

If `calculate_dem` is set to `True`, ebtel++ will read in the `dem` node from the configuration file which should be structured as follows,
```XML
<dem>
  <use_new_method>True</use_new_method>
  <temperature bins="451" log_min="4" log_max="8.5"/>
</dem>
```

If `use_new_method` is set to True (False), the transition region DEM is calculated using the method outlined in section 3 (the appendix) of [Klimchuk et al. (2008)][klimchuk_2008]. The `temperature` node configures the range and number of bins used when calculating the DEM. Here, for example, there are 450 bins of equal width between 10^4 and 10^8.5 K.

If you do not need to calculate the DEM, set the `calculate_dem` parameter to False and this section of the configuration file need not be included.

## Output
Once the EBTEL run has finished, the results are printed to the file specified in `output_filename` in the configuration file (as described above). Several examples of how to parse the results in Python can be found [here](https://github.com/rice-solar-physics/ebtelPlusPlus/tree/master/examples). In general, the results file follows the structure,

| | | | | | | | |
|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|
| _t0_ | _Te_(_t0_) | _Ti_(_t0_) | _n_(_t0_) | _pe_(_t0_) | _pi_(_t0_) | _v_(_t0_) | _h_(_t0_) |
| ... | ... | ... | ... | ... | ... | ... | ... |
| _ti_ | _Te_(_ti_) | _Ti_(_ti_) | _n_(_ti_) | _pe_(_ti_) | _pi_(_ti_) | _v_(_t0_) | _h_(_ti_) |
| ... | ... | ... | ... | ... | ... | ... | ... |
| _tN-1_ | _Te_(_tN-1_) | _Ti_(_tN-1_) | _n_(_tN-1_) | _pe_(_tN-1_) | _pi_(_tN-1_) | _v_(_tN-1_) | _h_(_tN-1_) |

Here _t_ is the time, _Te_ is the electron temperature, _Ti_ is the ion temperature, _pe_ is the electron pressure, _pi_ is the ion pressure, _n_ is the density, _v_ is the velocity, and _h_ is the heating rate.

If `calculate_dem` is set to True, the TR and coronal DEM results are printed to `output_filename`+`.dem_tr` and `output_filename`+`.dem_corona`, respectively. These output files are structured in the following way,

| | | | | |
|:----:|:----:|:----:|:----:|:----:|
| _T0_ | ... | _Tj_ | ... | _TM-1_ |
| DEM(_t0_,_T0_) | ... | DEM(_t0_,_Tj_) | ... | DEM(_t0_,_TM-1_) |
| ... | ... | ... | ... | ... |
| DEM(_ti_,_T0_) | ... | DEM(_ti_,_Tj_) | ... | DEM(_ti_,_TM-1_) |
| ... | ... | ... | ... | ... |
| DEM(_tN-1_,_T0_) | ... | DEM(_tN-1_,_Tj_) | ... | DEM(_tN-1_,_TM-1_) |

## Full Radiative Loss Function
__TODO__: This feature is not yet implemented. If you are interested in using or contributing this feature, create an [issue](https://github.com/rice-solar-physics/ebtelPlusPlus/issues) or [pull request](https://github.com/rice-solar-physics/ebtelPlusPlus/pulls). For more information about configuring the radiation model, see the [Radiation_Model docs](http://rice-solar-physics.github.io/Radiation_Model/).

## Examples
Included below are a few examples of how to run `ebtel++` and plot the results using Python. The results files can be parsed very easily by any language (e.g. IDL,Matlab). Configuration files can also be easily manipulated by hand, but in practice, particularly when computing many EBTEL runs, it is better to script this. The [xml_io.py](https://github.com/rice-solar-physics/rsp_toolkit/blob/master/python/xml_io.py) script provides several convenient utilities for manipulating these configuration files.

### Example 1
In this example, only the electrons are heated by a single triangular pulse lasting 500 seconds and injecting 10 ergs per cubic centimeter into the loop plasma. The code used to make this figure can be found [here](https://github.com/rice-solar-physics/ebtelPlusPlus/blob/master/examples/ex1.py).

![Example 1](ex1.png)

### Example 2
In this second example, only the ions are heated. In this case, we heat the ions with five distinct pulses, each lasting 200 seconds and separated by about 2000 seconds. The energies are uniformly distributed between 0.001 and 0.1 ergs per cubic centimeter per second. The example code can be found [here](https://github.com/rice-solar-physics/ebtelPlusPlus/blob/master/examples/ex2.py). Note how easy it is to programmatically generate the configurations for these different heating profiles.

![Example 2](ex2.png)

### Example 3
Lastly, we show an example where the electron and ion populations are forced into equilibrium at all times, i.e. the single-fluid case. The loop plasma is heated by a single pulse lasting 2000 seconds that rises quickly (in 250 seconds), is sustained at 0.005 ergs per cubic centimeter per second for 750 seconds, and then decays back to the background value over 1000 seconds. The example code can be found [here](https://github.com/rice-solar-physics/ebtelPlusPlus/blob/master/examples/ex3.py).

![Example 3](ex3.png)

[klimchuk_2008]: http://adsabs.harvard.edu/abs/2008ApJ...682.1351K "Klimchuk et al. (2008)"
[cargill_2012a]: http://adsabs.harvard.edu/abs/2012ApJ...752..161C "Cargill et al. (2012a)"
[cargill_2012b]: http://adsabs.harvard.edu/abs/2012ApJ...758....5C "Cargill et al. (2012b)"
[barnes_2016]: http://adsabs.harvard.edu/abs/2016ApJ...829...31B "Barnes et al. (2016)"
[press_num_recipes]: http://dl.acm.org/citation.cfm?id=148286 "Press et al. (1992)"
