# ebtel++

[ebtel++](https://github.com/rice-solar-physics/ebtelPlusPlus) is a two-fluid version of the original enthalpy-based thermal evolution of loops (EBTEL) model implemented in C++. This code provides an enhanced description of plasma behavior above roughly 5 MK. Further generic details about EBTEL can be found in the [repository for the original IDL code](https://github.com/rice-solar-physics/EBTEL) and in the references listed below.

The EBTEL model, originally developed by [Klimchuk et al. (2008)][klimchuk_2008] efficiently computes spatially-averaged, time-dependent plasma parameters ( e.g. temperature,  pressure, density) of dynamically-heated coronal loops. It is often desirable to compute solutions for a large number of coronal loops, but the spatial and temporal scales needed to solve the full 1D-hydrodynamic equations lead to long compute times for even 1D hydrodynamic codes. EBTEL computes quick and accurate solutions for spatially-averaged quantities, allowing efficient insight into how these monolithic structures evolve. [Barnes et al. (2016)][barnes_2016] improved upon this model by extending the treatment to the two-fluid hydrodynamic equations, allowing for differential heating between electrons and ions. Experienced EBTEL users should note the slightly different model for calculating the the c1 parameter during conductive cooling (see Appendix A of [Barnes et al., 2016][barnes_2016]).

EBTEL also calculates the differential emission measure (DEM) for both the transition region and the corona. Details regarding this formulation can be found in [Klimchuk et al. (2008)][klimchuk_2008]

## Dependencies

To compile ebtel++, first install the following dependencies,

* [gcc](https://gcc.gnu.org/) (at least v4.7; included with OS X, most Linux distros; [cygwin](https://www.cygwin.com/) on Windows)
* [scons](http://scons.org/) (`pip install scons` via [PyPI](https://pypi.python.org/pypi) or `conda install scons` with [Anaconda](https://www.continuum.io/downloads))
* [boost](http://www.boost.org/) (at least v1.53; `sudo apt-get install libboost-dev` on Debian Linux, `sudo port install boost` via [Macports](https://www.macports.org/) on OS X, [from source](https://github.com/rice-solar-physics/IonPopSolver#installing-boost-from-source) on Windows; can also be installed via [Anaconda](https://www.continuum.io/downloads))

Additionally, if you'd like to run the included tests and examples, you'll need the following Python dependencies, all easily installed with [anaconda](https://www.continuum.io/downloads),

* [Python](https://www.python.org/)
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

```Shell
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

## Testing

If you've installed the above Python dependencies, you can also run the tests using,

```Shell
$ pip install -r requirements/requirements-test.txt
$ pytest
```

or run any of the three included examples,

```Shell
$ python examples/ex1.py
$ python examples/ex2.py
$ python examples/ex3.py
```

## Citation

If you use ebtel++ in any published work, please include the following citations and mention the use of this code.

* [Klimchuk et al. (2008)][klimchuk_2008]
* [Cargill et al. (2012a)][cargill_2012a]
* [Cargill et al. (2012b)][cargill_2012b]
* [Barnes et al. (2016)][barnes_2016]

The first three papers detail the original single-fluid EBTEL model while the last paper gives the details of the two-fluid model. In particular, the details of how the two-fluid EBTEL equations are derived can be found in the appendix of [Barnes et al. (2016)][barnes_2016].

[klimchuk_2008]: http://adsabs.harvard.edu/abs/2008ApJ...682.1351K "Klimchuk et al. (2008)"
[cargill_2012a]: http://adsabs.harvard.edu/abs/2012ApJ...752..161C "Cargill et al. (2012a)"
[cargill_2012b]: http://adsabs.harvard.edu/abs/2012ApJ...758....5C "Cargill et al. (2012b)"
[barnes_2016]: http://adsabs.harvard.edu/abs/2016ApJ...829...31B "Barnes et al. (2016)"
