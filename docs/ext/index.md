#<cldoc:index>

ebtel++

# ebtel++
[![Build Status](https://travis-ci.org/rice-solar-physics/ebtelPlusPlus.svg?branch=master)](https://travis-ci.org/rice-solar-physics/ebtelPlusPlus)

[ebtel++](https://github.com/rice-solar-physics/ebtelPlusPlus) is a two-fluid version of the original enthalpy-based thermal evolution of loops (EBTEL) model implemented in C++. It replaces the now-deprecated [EBTEL-C code](https://github.com/rice-solar-physics/deprecated_EBTEL). The original IDL implementation can be found [here](https://github.com/rice-solar-physics/EBTEL).

The EBTEL model, originally developed by [Klimchuk et al. (2008)](http://adsabs.harvard.edu/abs/2008ApJ...682.1351K) allows one to efficiently compute spatially-averaged, time-dependent plasma parameters, e.g. _T_,_p_,_n_. It is often desirable to compute solutions for a large number of coronal loops, but the spatial and temporal scales needed to solve the full 1D-hydrodynamic equations lead to long computation times for even 1D hydrodynamic codes. EBTEL computes quick and accurate solutions for spatially-averaged quantities, allowing efficient insight into how these monolithic structures are heated and cool. [Barnes et al. (2016)](http://adsabs.harvard.edu/abs/2016ApJ...829...31B) improved upon this model by extending the treatment to the two-fluid hydrodynamic equations, allowing for differential heating between electrons and ions.

EBTEL also calculates the differential emission measure (DEM) for both the transition region and the corona. Details regarding this formulation can be found in [Klimchuk et al. (2008)](http://adsabs.harvard.edu/abs/2008ApJ...682.1351K).

## Citation
If you use ebtel++ in any published work, please include the following citations and mention the use of this code.

* [Klimchuk et al. (2008)](http://adsabs.harvard.edu/abs/2008ApJ...682.1351K)
* [Cargill et al. (2012a)](http://adsabs.harvard.edu/abs/2012ApJ...752..161C)
* [Cargill et al. (2012b)](http://adsabs.harvard.edu/abs/2012ApJ...758....5C)
* [Barnes et al. (2016)](http://adsabs.harvard.edu/abs/2016ApJ...829...31B)

The first three papers detail the original single-fluid EBTEL model while the last paper gives the details of the two-fluid model. In particular, the details of how the two-fluid EBTEL equations are derived can be found in the appendix of Barnes et al. (2016).

## Dependencies
To compile ebtel++, first install the following dependencies,

* [git](https://git-scm.com/) (included with OS X, most Linux distros)
* [scons](http://scons.org/) ([pip](https://pypi.python.org/pypi) or [anaconda](https://www.continuum.io/downloads))
* [boost](http://www.boost.org/) (`apt-get` on Linux, [Macports](https://www.macports.org/) for OS X, [from source](https://github.com/rice-solar-physics/IonPopSolver#installing-boost-from-source) on Windows)

The above list encapsulates all of the _required_ dependencies. If you'd like to use the full radiative loss function,

* [apolloDB](https://github.com/rice-solar-physics/apolloDB)

Additionally, if you'd like to run the included tests and examples, you'll need the following Python dependencies, all easily installed with [anaconda](https://www.continuum.io/downloads),

* [Python](https://www.python.org/) (2.7)
* [numpy](http://www.numpy.org/)
* [matplotlib](http://matplotlib.org/)
* [seaborn](https://stanford.edu/~mwaskom/software/seaborn/index.html)

## Installation
To download the code from GitHub and build the executable,
```Shell
$ git clone --recursive https://github.com/rice-solar-physics/ebtelPlusPlus.git
$ cd ebtelPlusPlus
$ scons
```
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
__TODO__: Add table with descriptions of input parameters. Name, description, and type

## Output
__TODO__: Add table including description of output file

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
Lastly, we show an example where the electron and ion populations are forced into equilibrium at all times, i.e. the single-fluid case. The loop plasma is heated by a single pulse lasting 1500 seconds that rises very quickly (in 100 seconds), is sustained at 0.005 ergs per cubic centimeter per second for 400 seconds, and then decays back to the background value over 1000 seconds. The example code can be found [here](https://github.com/rice-solar-physics/ebtelPlusPlus/blob/master/examples/ex3.py).

![Example 3](ex3.png)
