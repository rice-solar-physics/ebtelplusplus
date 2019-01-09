# ebtel++
[![Build Status](https://travis-ci.org/rice-solar-physics/ebtelPlusPlus.svg?branch=master)](https://travis-ci.org/rice-solar-physics/ebtelPlusPlus)

ebtel++ is a C++ implementation of the two-fluid EBTEL model, as detailed in [Barnes et al. (2016)](http://adsabs.harvard.edu/abs/2016ApJ...829...31B), for doing efficient hydrodynamics of dynamically-heated solar coronal loops. The EBTEL model was originally developed by [Klimchuk et al. (2008)](http://adsabs.harvard.edu/abs/2008ApJ...682.1351K) and [Cargill et al. (2012)](http://adsabs.harvard.edu/abs/2012ApJ...752..161C). This code provides an enhanced description of plasma behavior above roughly 5 MK. Further generic details about EBTEL can be found in the [repository for the original IDL code](https://github.com/rice-solar-physics/EBTEL) and in the references listed below.

## Citation
If you use ebtel++ in any published work, please cite the following papers:

* [Klimchuk et al. (2008)](http://adsabs.harvard.edu/abs/2008ApJ...682.1351K)
* [Cargill et al. (2012a)](http://adsabs.harvard.edu/abs/2012ApJ...752..161C)
* [Cargill et al. (2012b)](http://adsabs.harvard.edu/abs/2012ApJ...758....5C)
* [Barnes et al. (2016)](http://adsabs.harvard.edu/abs/2016ApJ...829...31B)

## Installation
Download the dependencies,
* [git](https://git-scm.com/)
* [gcc](https://gcc.gnu.org/) (v4.7 or later)
* [scons](http://scons.org/)
* [boost](http://www.boost.org/) (v1.53 or later)

Then to install, compile, and run the code using the example configuration file in `config/ebtel.example.cfg.xml`,
```Shell
$ git clone --recursive https://github.com/rice-solar-physics/ebtelPlusPlus.git
$ cd ebtelPlusPlus
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
