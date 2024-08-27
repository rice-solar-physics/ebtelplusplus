ebtelplusplus
=============

.. toctree::
   :maxdepth: 1
   :hidden:

   generated/gallery/index
   development
   reference
   bibliography

``ebtelplusplus`` (or ``ebtel++``) is a two-fluid version of the original enthalpy-based thermal evolution of loops (EBTEL) model implemented in C++.
This code provides an enhanced description of plasma behavior above roughly 5 MK.
If you are looking for the original IDL implementation, `repository for the original IDL code can be found here <https://github.com/rice-solar-physics/EBTEL>`__.

The EBTEL model, originally developed by :cite:t:`klimchuk_highly_2008`, efficiently computes spatially-averaged, time-dependent plasma parameters ( e.g. temperature,  pressure, density) of dynamically-heated coronal loops.
It is often desirable to compute solutions for a large number of coronal loops, but the spatial and temporal scales needed to solve the full 1D-hydrodynamic equations lead to long compute times for even 1D hydrodynamic codes.
EBTEL computes quick and accurate solutions for spatially-averaged quantities, allowing efficient insight into how these monolithic structures evolve.
:cite:t:`barnes_inference_2016` improved upon this model by extending the treatment to the two-fluid hydrodynamic equations, allowing for differential heating between electrons and ions.
This also introduced a slightly different model for calculating the the :math:`c_1`` parameter during conductive cooling :cite:p:`{see Appendix A of}barnes_inference_2016`.
EBTEL also calculates the differential emission measure (DEM) for both the transition region and the corona.
Details regarding this formulation can be found in :cite:t:`klimchuk_highly_2008`.

Installation
------------

The easiest way to install ``ebtelplusplus`` is through ``pip``,

.. code::

   pip install ebtelplusplus

This will install all of the needed Python dependencies as well as a precompiled binary for your particular operating system.
Note that it is not necessary to compile the code locally in order to use ``ebtelplusplus``.
If you are interested in modifying or contributing to ``ebtelplusplus``, see the :ref:`ebtelplusplus-development` page.
Once you've successfully installed ``ebtelplusplus``, see the Example Gallery for examples of how to run a simulation.

Why ``ebtelplusplus``?
----------------------

Citation
--------

If you use ``ebtelplusplus`` in any published work, please include the following citations:

* :cite:t:`klimchuk_highly_2008`
* :cite:t:`cargill_enthalpy-based_2012`
* :cite:t:`cargill_enthalpy-based_2012-1`
* :cite:t:`barnes_inference_2016`

Additionally, please add the following line within your methods, conclusion or acknowledgements sections:

*This research used version X.Y.Z (software citation) of the ebtelplusplus open source software package.*

The software citation should be the specific `Zenodo DOI`_ for the version used within your work.
The citation for most current version on Zenodo is,

.. code:: bibtex

   @software{Barnes2024,
      author       = {Barnes, Will and
                      Klimchuk, James and
                      Cargill, Peter and
                      Bradshaw, Stephen and
                      Reep, Jeffrey and
                      Schonfeld, Samuel and
                      Collazzo, Rowan},
      title        = {rice-solar-physics/ebtelPlusPlus: v0.2},
      month        = jul,
      year         = 2024,
      publisher    = {Zenodo},
      version      = {v0.2},
      doi          = {10.5281/zenodo.12675386},
      url          = {https://doi.org/10.5281/zenodo.12675386}
   }


.. _Zenodo DOI: https://zenodo.org/records/12675374