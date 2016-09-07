#ebtel++
ebtel++ is a C++ implementation of the two-fluid EBTEL model, as detailed in [Barnes et al. (2016)](http://adsabs.harvard.edu/abs/2016arXiv160804776B), for doing efficient hydrodynamics of dynamically-heated solar coronal loops. The EBTEL model was originally developed by [Klimchuk et al. (2008)](http://adsabs.harvard.edu/abs/2008ApJ...682.1351K) and [Cargill et al. (2012)](http://adsabs.harvard.edu/abs/2012ApJ...752..161C). The IDL code for the original model can be found [here](https://github.com/rice-solar-physics/EBTEL).

## Citation
If you use ebtel++ in any published work, please cite the following papers:

* [Klimchuk et al. (2008)](http://adsabs.harvard.edu/abs/2008ApJ...682.1351K)
* [Cargill et al. (2012a)](http://adsabs.harvard.edu/abs/2012ApJ...752..161C)
* [Cargill et al. (2012b)](http://adsabs.harvard.edu/abs/2012ApJ...758....5C)
* [Barnes et al. (2016)](http://adsabs.harvard.edu/abs/2016arXiv160804776B)

## Installation
The following dependencies are required to install ebtel++,

* [git](#)
* [scons](#)
* [boost](#)
* [apolloDB](#)

To install and compile the code,
```Shell
git clone --recursive https://github.com/rice-solar-physics/ebtel++.git
cd ebtel++
scons
```
The recursive flag ensures that the [`Radiation_Model`]() and [`rsp_toolkit`]() submodules are pulled down as well.

### Move all this to the docs...

Before running the program, you'll need to point `Radiation_Model` at the atomic data installed in `apolloDB`. First, make a copy of the configuration file,
```Shell
cp config/radiation.example.cfg.xml config/radiation.local.cfg.xml
```
and then change the `atomicDB` parameter in ``,
```XML
<atomicDB>/my/custom/path/to/apolloDB/</atomicDB>
```
to the location where you installed `apolloDB`.

## Example

## Help
