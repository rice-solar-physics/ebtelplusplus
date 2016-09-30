# ebtel++
[![Build Status](https://travis-ci.org/rice-solar-physics/ebtelPlusPlus.svg?branch=master)](https://travis-ci.org/rice-solar-physics/ebtelPlusPlus)

## About
ebtel++ is a C++ implementation of the two-fluid EBTEL model, as detailed in [Barnes et al. (2016)](http://adsabs.harvard.edu/abs/2016ApJ...829...31B), for doing efficient hydrodynamics of dynamically-heated solar coronal loops. It is the most current and mainted version of the EBTEL model. The EBTEL model was originally developed by [Klimchuk et al. (2008)](http://adsabs.harvard.edu/abs/2008ApJ...682.1351K) and [Cargill et al. (2012)](http://adsabs.harvard.edu/abs/2012ApJ...752..161C). The IDL code for the original model can be found [here](https://github.com/rice-solar-physics/EBTEL).

## Citation
If you use ebtel++ in any published work, please cite the following papers:

* [Klimchuk et al. (2008)](http://adsabs.harvard.edu/abs/2008ApJ...682.1351K)
* [Cargill et al. (2012a)](http://adsabs.harvard.edu/abs/2012ApJ...752..161C)
* [Cargill et al. (2012b)](http://adsabs.harvard.edu/abs/2012ApJ...758....5C)
* [Barnes et al. (2016)](http://adsabs.harvard.edu/abs/2016ApJ...829...31B)

## Installation
Download the dependencies,
* [git](https://git-scm.com/)
* [scons](http://scons.org/)
* [boost](http://www.boost.org/)

Optionally, if you'd like to use the full radiative loss function option,
* [apolloDB](https://github.com/rice-solar-physics/apolloDB)

Then to install, compile, and run the code using the example configuration file in `config/ebtel.example.cfg.xml`,
```Shell
$ git clone --recursive https://github.com/rice-solar-physics/ebtel++.git
$ cd ebtel++
$ scons
$ bin/ebtel++.run
```
This will create a results file `ebtel++_results_file.txt` in the current directory. For more info about which parameters can be passed to `ebtel++.run`,
```Shell
$ bin/ebtel++.run --help
```

For more information about how to setup the configuration file, see the [documentation](http://rice-solar-physics.github.io/ebtelPlusPlus/).

## Help
* [Documentation](http://rice-solar-physics.github.io/ebtelPlusPlus/)
* [Report a bug](https://github.com/rice-solar-physics/ebtelPlusPlus/issues)
* [Contribute code](https://github.com/rice-solar-physics/ebtelPlusPlus/pulls)
