#<cldoc:index>

ebtel++

# ebtel++
TRAVIS BADGE HERE

ebtel++ is a two-fluid version of the original enthalpy-based thermal evolution of loops (EBTEL) model implemented in C++. It replaces the now-deprecated [EBTEL-C code](#). The original IDL implementation can be found [here](#).

The EBTEL model, originally developed by [Klimchuk et al. (2008)](http://adsabs.harvard.edu/abs/2008ApJ...682.1351K) allows one to efficiently compute spatially-averaged, time-dependent plasma parameters, e.g. _T_,_p_,_n_. It is often desirable to compute solutions for a large number of coronal loops, but the spatial and temporal scales needed to solve the full 1D-hydrodynamic equations lead to long computation times for even 1D hydrodynamic codes. EBTEL computes quick and accurate solutions for spatially-averaged quantities, allowing efficient insight into how these monolithic structures are heated and cool. [Barnes et al. (2016)](http://adsabs.harvard.edu/abs/2016arXiv160804776B) improved upon this model by extending the treatment to the two-fluid hydrodynamic equations, allowing for differential heating between electrons and ions.

EBTEL also calculates the differential emission measure (DEM) for both the transition region and the corona. Details regarding this formulation can be found in [Klimchuk et al. (2008)](http://adsabs.harvard.edu/abs/2008ApJ...682.1351K).

## Citation
If you use ebtel++ in any published work, please include the following citations and mention the use of this code.

* [Klimchuk et al. (2008)](http://adsabs.harvard.edu/abs/2008ApJ...682.1351K)
* [Cargill et al. (2012a)](http://adsabs.harvard.edu/abs/2012ApJ...752..161C)
* [Cargill et al. (2012b)](http://adsabs.harvard.edu/abs/2012ApJ...758....5C)
* [Barnes et al. (2016)](http://adsabs.harvard.edu/abs/2016arXiv160804776B)

The first three papers detail the original single-fluid EBTEL model while the last paper gives the details of the two-fluid model. In particular, the details of how the two-fluid EBTEL equations are derived can be found in the appendix of Barnes et al. (2016).

## Dependencies and Installation

## Configuration

## Output

## Full Radiative Loss Function
__TODO__: This feature is not yet implemented. If you are interested in using or contributing this feature, create an [issue](#) or [pull request](#).

## Examples
Included below are a few examples of how to run `ebtel++` and plot the results using Python. The results files can be parsed very easily by any language (e.g. IDL,Matlab). Configuration files can also be easily manipulated by hand, but in practice, particularly when computing many EBTEL runs, it is better to script this. The [xml_io.py](https://github.com/rice-solar-physics/rsp_toolkit/blob/master/python/xml_io.py) script provides several convenient utilities for manipulating these configuration files.

### Example 1
In this example, only the electrons are heated by a single triangular pulse lasting 500 seconds and injecting 10 ergs per cubic centimeter into the loop plasma. The code used to make this figure can be found [here](#).
![Example 1](ex1.png)

### Example 2
![Example 2](ex2.png)

### Example 3
![Example 3](ex3.png)
