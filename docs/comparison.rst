.. _ebtelplusplus-comparison:

Why ebtelplusplus?
======================

For various historical reasons, there are multiple software implementations of the EBTEL model.
There are currently *three* maintained implementations of the model as described below.
Hereafter, "EBTEL" refers to the model generally while each software implementation is given
its own unique name (e.g. ``ebtelplusplus``).
This page briefly explains the differences between each implementation and the advantages of
the ``ebtelplusplus`` implementation.

The table below summarizes the comparison of the different EBTEL software implementations and the different features provided by each.

.. list-table::
   :header-rows: 1
   :align: center

   * - Feature
     - Citation
     - |EBTEL-IDL|_
     - ``ebtel++``
     - |EBTEL3-IDL|_
   * - Decouple electrons and ions
     - :cite:t:`barnes_inference_2016`
     - :octicon:`x;1em;sd-text-danger`
     - :octicon:`check;1em;sd-text-success`
     - :octicon:`x;1em;sd-text-danger`
   * - Adaptive time step
     - :cite:t:`barnes_inference_2016`
     - :octicon:`x;1em;sd-text-danger`
     - :octicon:`check;1em;sd-text-success`
     - :octicon:`check;1em;sd-text-success`
   * - Area expansion
     - :cite:t:`cargill_static_2022`
     - :octicon:`check;1em;sd-text-success`
     - :octicon:`check;1em;sd-text-success`
     - :octicon:`x;1em;sd-text-danger`
   * - Supersonic flows
     - :cite:t:`rajhans_flows_2022`
     - :octicon:`x;1em;sd-text-danger`
     - :octicon:`x;1em;sd-text-danger`
     - :octicon:`check;1em;sd-text-success`
   * - Time-variable abundances
     - :cite:t:`reep_modeling_2024`
     - :octicon:`x;1em;sd-text-danger`
     - :octicon:`check;1em;sd-text-success`
     - :octicon:`x;1em;sd-text-danger`


``EBTEL-IDL``
-------------

The original software implementation of EBTEL was in the Interactive Data Language (IDL) and
was based off of the model presented in :cite:t:`klimchuk_highly_2008`.
Subsequent improvements in :cite:t:`cargill_enthalpy-based_2012` were made to give better
agreement with spatially-resolved hydrodynamic models.
These improvements were also implemented in IDL.
This version is sometimes referred to as "EBTEL2".
:cite:t:`cargill_static_2022` extended the EBTEL model to include effects due to cross-sectional
area expansion and implemented this in IDL as well.
The IDL software implementation which includes all of these features is referred to as |EBTEL-IDL|_.

``ebtelplusplus``
-----------------

:cite:t:`barnes_inference_2016` improved upon the implementation of :cite:t:`cargill_enthalpy-based_2012`
by extending the treatment to the two-fluid hydrodynamic equations, allowing for differential heating
between electrons and ions.
They also introduced a slightly modified approach for calculating the the :math:`c_1` parameter during
the conductive cooling phase :cite:p:`{see Appendix A of}barnes_inference_2016`.
Modifications to include area expansion in the manner of :cite:t:`cargill_static_2022` were subsequently added
as well as the ability to vary the abundance model for the radiative losses as a function of time :cite:p:`reep_modeling_2024`.
Furthermore, the resulting equations are solved using a Runge-Kutta Cash-Karp integration method
:cite:p:`{see section 16.2 of}press_numerical_1992` and an (optional) adaptive time-stepping technique
to ensure the principal physical timescales are resolved at each phase of the loop evolution.
The software implementation of this version of the model is referred to as ``ebtelplusplus`` (or ``ebtel++``).
The majority of the software is implemented in C++ for computational efficiency and is wrapped in Python
to enable easy installation and a user-friendly API.
This is the implementation provided by this software package.
``ebtelplusplus`` has been benchmarked against both |EBTEL-IDL|_ as well as more advanced field-aligned
hydrodynamic models :cite:p:`barnes_inference_2016`.

``EBTEL3-IDL``
--------------

:cite:t:`rajhans_flows_2022` built upon the model of :cite:t:`cargill_enthalpy-based_2012` and relaxed
the assumption of subsonic flows in EBTEL.
Additionally the Mach numbers and velocities produced are in better agreement with field-aligned
hydrodynamic simulations.
The IDL software implementation of this model is referred to as |EBTEL3-IDL|_.
|EBTEL3-IDL|_ uses an adaptive time grid to ensure the appropriate timescales are resolved in the
impulsive phase.


.. |EBTEL-IDL| replace:: ``EBTEL-IDL``
.. _EBTEL-IDL: https://github.com/rice-solar-physics/EBTEL
.. |EBTEL3-IDL| replace:: ``EBTEL3-IDL``
.. _EBTEL3-IDL: https://github.com/rice-solar-physics/EBTEL3
