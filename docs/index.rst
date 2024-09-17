ebtelplusplus
=============

.. toctree::
   :maxdepth: 1
   :hidden:

   generated/gallery/index
   topic_guides/index
   development
   reference
   bibliography

``ebtelplusplus`` (or ``ebtel++``) is a two-fluid version of the original enthalpy-based thermal evolution of loops (EBTEL) model implemented in C++ and wrapped in Python.

.. grid:: 1 2 2 2
    :gutter: 2

    .. grid-item-card:: Examples
        :link: generated/gallery
        :text-align: center

        :material-outlined:`palette;8em;sd-text-secondary`

        Examples of how to set up a simulation

    .. grid-item-card:: Topic Guides
        :link: ebtelplusplus-topic-guide
        :link-type: ref
        :text-align: center

        :material-outlined:`school;8em;sd-text-secondary`

        Detailed explanations about the EBTEL model

    .. grid-item-card:: Contributing
        :link: ebtelplusplus-development
        :link-type: ref
        :text-align: center

        :material-outlined:`code;8em;sd-text-secondary`

        Instructions for how to contribute to `ebtelplusplus`

    .. grid-item-card:: Reference
        :link: ebtelplusplus-reference
        :link-type: ref
        :text-align: center

        :material-outlined:`menu_book;8em;sd-text-secondary`

        Technical description of the inputs, outputs, and behavior of each component of ebteplusplus

.. note::

   If you are looking for the original IDL implementation,
   the repository for that code can be found `here <https://github.com/rice-solar-physics/EBTEL>`__.

The EBTEL model, originally developed by :cite:t:`klimchuk_highly_2008`, efficiently computes spatially-averaged, time-dependent plasma parameters ( e.g. temperature,  pressure, density) of dynamically-heated coronal loops.
It is often desirable to compute solutions for a large number of coronal loops, but the spatial and temporal scales needed to solve the full 1D-hydrodynamic equations lead to long compute times for even 1D hydrodynamic codes.
EBTEL computes quick and accurate solutions for spatially-averaged quantities, allowing efficient insight into how these monolithic structures evolve.
EBTEL also calculates the differential emission measure (DEM) for both the transition region and the corona.
Details regarding this formulation can be found in :cite:t:`klimchuk_highly_2008`.

Installation
------------

The easiest way to install ``ebtelplusplus`` is through ``pip``,

.. code:: shell

   pip install ebtelplusplus

This will install all of the needed Python dependencies as well as a precompiled binary for your particular operating system.
Note that it is not necessary to compile the code locally in order to use ``ebtelplusplus``.
If you are interested in modifying or contributing to ``ebtelplusplus``, see the :ref:`ebtelplusplus-development` page.
Once you've successfully installed ``ebtelplusplus``, see the Example Gallery for examples of how to run a simulation.

Citation
--------

If you use ``ebtelplusplus`` in any published work, we ask that you please include citations to the following papers:

* :cite:t:`klimchuk_highly_2008`
* :cite:t:`cargill_enthalpy-based_2012`
* :cite:t:`cargill_enthalpy-based_2012-1`
* :cite:t:`barnes_inference_2016`

If you make use of the variable abundance feature in the radiative losses, please also cite :cite:t:`reep_modeling_2024`.
Additionally, please add the following line within your methods, conclusion or acknowledgements sections:

*This research used version X.Y.Z (software citation) of the ebtelplusplus open source software package.*

The software citation should be the specific `Zenodo DOI`_ for the version used within your work.

.. _Zenodo DOI: https://zenodo.org/records/12675374
