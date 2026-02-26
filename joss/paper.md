---
title: "ebtelplusplus: Efficient Hydrodynamic Modeling of Coronal Loops"
tags:
- Python
- hydrodynamics
- astronomy
- astrophysics
- solar physics
authors:
- name: W. T. Barnes
  orcid: 0000-0001-9642-6089
  affiliation: "1,2"
- name: J. W. Reep
  orcid: 0000-0003-4739-1152
  affiliation: "3"
- name: N. Freij
  orcid: 0000-0002-6253-082X
  affiliation: "4,5"
- name: S. J. Schonfeld
  orcid: 0000-0002-5476-2794
  affiliation: "6"
- name: R. Collazzo
  orcid:
  affiliation: "1"
- name: P. J. Cargill
  orcid:
  affiliation: "7,8"
- name: S. J. Bradshaw
  orcid: 0000-0002-3300-6041
  affiliation: "9"
- name: J. A. Klimchuk
  orcid: 0000-0003-2255-0305
  affiliation: "2"
affiliations:
- name: Department of Physics, American University, Washington, D.C. 20016, USA
  index: 1
- name: NASA Goddard Space Flight Center, Greenbelt, MD 20771, USA
  index: 2
- name: Institute for Astronomy, University of Hawai'i at Manoa, Pukalani, HI 96768, USA
  index: 3
- name: SETI Institute, Mountain View, CA 94043, USA
  index: 4
- name: Lockheed Martin Solar and Astrophysics Laboratory, Palo Alto, CA 94304, USA
  index: 5
- name: Air Force Research Laboratory, Space Vehicles Directorate, Kirtland AFB, NM 87117, USA
  index: 6
- name: Space and Atmospheric Physics, The Blackett Laboratory, Imperial College, London SW7 2BW, UK
  index: 7
- name: School of Mathematics and Statistics, University of St Andrews, St Andrews KY16 9SS, UK
  index: 8
- name: Department of Physics and Astronomy, Rice University, Houston, TX 77005, USA
  index: 9
date: 12 February 2026
bibliography: paper.bib
---

# Summary

Understanding the response of the plasma in the solar corona, the outermost layer of the Sun's atmosphere, to heating is of considerable importance in understanding flares and the heating of the quiescent corona.
This requires detailed numerical modeling of the coupled system of the coronal plasma and the solar magnetic field[^stellar].
While solving the full set of three-dimensional *magnetohydrodynamic* (MHD) equations is feasible for small regions of the corona, the needed computational resources and physical complexity of such models means they are not always amenable to simple or even correct interpretation.
Field-aligned hydrodynamic models [e.g. HYDRAD @bradshaw_self-consistent_2003] exploit the fact that the magnetic pressure in the corona is much greater than the gas pressure such that the corona can be considered as a series of mini atmospheres, or *coronal loops*, where the plasma responds hydrodynamically to the heating and traces out the complex coronal magnetic field.
As such, the explicit dependence on the magnetic field can be neglected and the relevant hydrodynamic equations can be reduced to a single-dimension in space: the coordinate along the coronal loop.
However, the large range of spatial and temporal scales necessary to resolve this system, in particular the severe time step limitations imposed by thermal conduction, mean that even field-aligned models are computationally expensive enough to make large parameter explorations prohibitive.
The enthalpy-based thermal evolution of loops (EBTEL) model [@klimchuk_highly_2008;@cargill_enthalpy-based_2012] was originally developed in order to provide a simple and efficient way to study the coronal plasma response to time-dependent plasma heating.
EBTEL accomplishes this by computing spatial integrals of the aforementioned field-aligned hydrodynamic equations.
The `ebtelplusplus` software package is a Python and C++ implementation of the EBTEL model and includes effects due to cross-sectional coronal loop expansion, electron-ion coupling, and variable elemental abundances.

# Statement of Need

The key to EBTEL lies in its treatment of enthalpy between the corona and the transition region (TR), the narrow layer of the solar atmosphere that connects the corona with the denser chromosphere below.
This enables EBTEL to split the coronal loop into two regions that match across the corona/TR boundary, leading to a set of coupled, non-linear ordinary differential equations that model the spatially-averaged, time-dependent behavior of the relevant thermodynamic quantities.
EBTEL equates an enthalpy flux with an imbalance between the heat flux out of the corona and the energy lost due to radiation in the TR.
If TR radiation cannot balance the downward heat flux, this drives an upflow of material into the corona and if the TR is radiating away more energy than the coronal heat flux can supply this drives a downflow.
This approximation is valid for bulk velocities below the local sound speed [@klimchuk_highly_2008].

The EBTEL model was originally developed by @klimchuk_highly_2008.
Subsequent improvements to the gravitational stratification and radiative losses by @cargill_enthalpy-based_2012 gave better agreement with field-aligned hydrodynamic models[^ebtel2].
@barnes_inference_2016 modified the EBTEL model to relax the single-fluid assumption and treat electrons and ions separately, allowing for differential heating between the two species.
@cargill_static_2021 later extended EBTEL to include effects due to cross-sectional area expansion and @reep_modeling_2024 added the ability to vary the abundance model for the radiative losses as a function of time.
`ebtelplusplus` unifies all of the aforementioned features into a single set of equations.
In particular, `ebtelplusplus` solves the following equations for the spatially-averaged electron pressure ($p_e$), ion pressure ($p_i$), and number density ($n$) of a semi-circular coronal loop of half-length $L$,

\begin{eqnarray*}
\frac{1}{\gamma-1}\frac{dp_e}{dt} &=& Q_e + \frac{\psi_c}{L_*}\left(1+\frac{A_{TR}\psi_{TR}}{A_c\psi_c}\right) - \frac{R_c}{L_*}\left(1+c_1\frac{A_{TR}}{A_c}\right), \\
\frac{1}{\gamma-1}\frac{dp_i}{dt} &=& Q_i - \frac{\psi_c}{L_*}\left(1 + \frac{A_{TR}\psi_{TR}}{A_c\psi_c}\right), \\
\frac{dn}{dt} &=& -\frac{(\gamma - 1)\xi c_2}{(\xi + 1)\gamma c_3 k_B L_c T_e}\left(\frac{A_{TR}L_c}{A_cL_*}R_c\left(c_1 - \frac{L_{TR}}{L_c}\right) + \frac{A_0}{A_c}(F_{e,0} + F_{i,0})\right),
\end{eqnarray*}

where $Q_{e,i}$ are the user-specified heating terms for the electrons and ions, $\psi_{c,TR}$ denote integrals of the electron-ion coupling terms over the TR and corona, respectively, $A_{c,TR,0}$ are the cross-sectional area averaged over the corona and TR and at the TR/corona boundary, respectively, $R_c$ is the energy lost to radiation in the corona, $c_1$ is the ratio of energy lost to radiation in the TR and corona, $\xi=T_e/T_i$ is the ratio between the electron and ion temperatures, $F_{e,0;i,0}$ are the conductive heat fluxes at the TR/corona boundary for the electrons and ions, $L_{c,TR}$ are the lengths of the corona and TR such that $L=L_c+L_{TR}$, and $L_*=L_c + (A_{TR}/A_c)L_{TR}$.
The remaining terms are fixed constants.

This set of equations is closed by an ideal gas law for the electrons and ions: $p_e=k_BnT_e,p_i=k_BnT_i$ and $n_e=n_i=n$ due to the assumption of a fully-ionized hydrogen plasma.
These equations and their derivations are explained more fully in the aforementioned publications and the [`ebtelplusplus` documentation](https://ebtelplusplus.readthedocs.io/en/stable/topic_guides/derivation.html).
`ebtelplusplus` solves the above equations using a Runge-Kutta Cash-Karp integration method [see section 16.2 of @press_numerical_1992] and an (optional) adaptive time-stepping scheme to ensure the principal physical timescales are resolved at each phase of the loop evolution[^boost].
\autoref{fig:figure1} shows example output from `ebtelplusplus` with different model parameters for the same time-dependent heating function.

![Temperature (top right), density (bottom left), and temperature-density phase space (bottom right) of a coronal loop with half-length $L=40$ Mm for five different cases with the same heating input (top left panel). In the nominal case (blue), the electron and ion populations are kept in equilibrium, the cross-sectional area of the loop is constant, and the radiative losses are determined by a power-law function. If the electrons (solid) and ions (dashed) are allowed to evolve separately, heating only the electrons (orange) causes the ions to take about 250 s to fully equilibrate with the electrons while heating only the ions (green) causes the ions to become over three times hotter than the electrons due to the relative inefficiency of ion thermal conduction. Incorporating area expansion through the corona (red) leads to a higher peak temperature and a more delayed peak in the density while calculating the radiative losses using a time-varying abundance (purple) leads to a slightly higher peak density.\label{fig:figure1}](figure.pdf)

# State of the Field

There are currently three separate implementations of the EBTEL model.
The original software implementation of the EBTEL model was in the proprietary Interactive Data Language (IDL).
This implementation, which includes features described in @cargill_enthalpy-based_2012 and @cargill_static_2021, is referred to as [`EBTEL-IDL`](https://github.com/rice-solar-physics/EBTEL) [@cargill_2024_13351770].
Comparisons between `EBTEL-IDL` and spatially-averaged results from field-aligned hydrodynamic models show very good agreement [@cargill_enthalpy-based_2012].
@rajhans_flows_2022 relaxed the assumption of subsonic flows in EBTEL such that the Mach numbers and velocities produced are in better agreement with field-aligned hydrodynamic simulations for some heating scenarios.
The IDL software implementation of this model is referred to as `EBTEL3-IDL`.
The initial C++ implementation of `ebtelplusplus` was developed by @barnes_inference_2016 with modifications later made by @reep_modeling_2024 and the Python interface added later.
This is the software implementation described in this paper.
The table below summarizes the features included in each implementation.

| Feature                     | Citation               | `EBTEL-IDL` | `EBTEL3-IDL` | `ebtelplusplus` |
|:----------------------------|:-----------------------|:-----------:|:------------:|:---------------:|
| Decouple electrons and ions | @barnes_inference_2016 | no          | no           | yes             |
| Adaptive time step          | @barnes_inference_2016 | no          | yes          | yes             |
| Area expansion              | @cargill_static_2021   | yes         | no           | yes             |
| Supersonic flows            | @rajhans_flows_2022    | no          | yes          | no              |
| Time-variable abundances    | @reep_modeling_2024    | no          | no           | yes             |

# Software Design

The design of the `ebtelplusplus` software is motivated by two primary needs: 1. computational efficiency and 2. a high-level, intuitive interface.
Both are essential for exploratory analysis of time-dependent heating of the coronal plasma, including large-scale parameter explorations.
`ebtelplusplus` is implemented in C++ for computational efficiency and is wrapped in Python using `pybind11` [@wenzel_jakob_2025_16929811] to enable easier installation and a user-friendly API.
As a result, `ebtelplusplus` is very fast (a single run modeling $10^4$ seconds of simulation time takes only a few milliseconds) and nearly two orders of magnitude faster than previous IDL implementations.
Where appropriate, all inputs and outputs are expressed as `astropy.units.Quantity` objects [@astropy_collaboration_astropy_2022] to maximize flexibility and avoid ambiguity.
As an example, two of the primary inputs for configuring an `ebtelplusplus` simulation are the total simulation time and the loop length.
These inputs can be expressed in any units provided they can be converted to seconds and centimeters, respectively.
High-level Python objects are provided for configuring additional inputs, including the time-dependent heating, and include default values to avoid overly-verbose input configurations.

To make the installation process easier for users, precompiled binary wheels are built using [`cibuildwheel`](https://cibuildwheel.pypa.io/en/stable/) run on GitHub Actions[^oatemplates] and distributed via [PyPI](https://pypi.org/project/ebtelplusplus/) at every release.
These wheels are provided for all major operating systems and the versions of Python recommended by SPEC 0[^spec0].
This alleviates the need to compile the C++ code locally and allows new users to start using the software more quickly.
`ebtelplusplus` is openly-developed on [GitHub](https://github.com/rice-solar-physics/ebtelplusplus).
Documentation, including an example gallery and a guide to contributing to the package, is hosted online on [Read the Docs](https://ebtelplusplus.readthedocs.io).
`ebtelplusplus` also includes a comprehensive test suite built on the [`pytest` testing framework](https://docs.pytest.org/) that is run on GitHub Actions at each check-in.
Test coverage is assessed using [Codecov](https://about.codecov.io/).

# Research Impact Statement

Because of its relative simplicity and computational efficiency, EBTEL has been widely used since its initial development [e.g. @qiu_heating_2012;@ugarte-urra_determining_2014] with @klimchuk_highly_2008 and @cargill_enthalpy-based_2012 having nearly 400 citations combined according to the Astrophysics Data System (ADS).
@barnes_inference_2016, which describes the original implementation of the `ebtelplusplus` model, has over 50 citations per ADS.
In particular, `ebtelplusplus` has been used to model thousands of impulsively-heated loops in coronal active regions [@barnes_understanding_2019] and constrain properties of microflares via comparisons with hard x-ray observations [@marsh_hard_2018].

# AI Usage Disclosure

No generative AI tools were used in the development of this software, the writing of this manuscript, or the preparation of supporting materials.

# References

[^stellar]: This also applies to the coronae of F, G, K, and M dwarf stars.
[^ebtel2]: This version is sometimes referred to as “EBTEL2”.
[^boost]: The Runge-Kutta Cash-Karp integrator is provided by the [Boost Odeint library](https://www.boost.org/library/latest/numericodeint/).
[^spec0]: [Scientific Python Ecosystem Coordination (SPEC) 0](https://scientific-python.org/specs/spec-0000/) recommends a set of minimum supported dependencies for packages commonly used across the scientific Python ecosystem, including Python.
[^oatemplates]: `ebtelplusplus` uses the [reusable GitHub Actions workflows](https://github.com/OpenAstronomy/github-actions-workflows) provided by the OpenAstronomy project.
