.. _ebtelplusplus-topic-guide-derivation:

Deriving the EBTEL Equations
============================

This page will briefly summarize the derivation of the version of the EBTEL
model used in `ebtelplusplus` which includes cross-sectional area expansion
and two-fluid effects.
For a more detailed explanation of the EBTEL model, including the physical
motivation behind the assumptions, see :cite:t:`klimchuk_highly_2008`.
Additionally, :cite:t:`barnes_inference_2016` provides a detailed discussion
of the two-fluid implementation while :cite:t:`cargill_static_2022` discusses
the inclusion of cross-sectional area expansion.

The 1D Hydrodynamic Equations
-----------------------------

Deriving the EBTEL model involves taking spatial integrals,
over both the transition region (TR) and corona, of the one-dimensional (1D)
hydrodynamic equations that govern the evolution of the dynamically-heated
solar atmosphere confined along a magnetic field line.
This yields a set of so-called "zero-dimensional" (0D) equations which
describe the evolution in time of the coronally-averaged temperature, density,
and pressure.
We begin by listing the relevant 1D hydrodynamic equations.
The 1D density, electron energy, and ion energy equations are given by,

.. math::
    :label: density-1d

    \frac{\partial n}{\partial t} + \frac{1}{A}\frac{\partial}{\partial s}(Anv) = 0,

.. math::
    :label: electron-energy-1d

    \frac{\partial E_e}{\partial t} + \frac{1}{A}\frac{\partial}{\partial s}(A(E_e + p_e)v) = &-\frac{1}{A}\frac{\partial}{\partial s}(AF_e) + v\frac{\partial p_e}{\partial s} - n^2\Lambda \\
    &+ \frac{k_Bn}{\gamma-1}\nu_{ie}(T_i - T_e) + Q_e,

.. math::
    :label: ion-energy-1d

    \frac{\partial E_i}{\partial t} + \frac{1}{A}\frac{\partial}{\partial s}(A(E_i + p_i)v) = &-\frac{1}{A}\frac{\partial}{\partial s}(AF_i) - v\frac{\partial p_e}{\partial s} + m_i n v g_{\parallel} + \frac{1}{A}\left(\frac{4A\eta v}{3}\frac{\partial v}{\partial s}\right) \\
    &+ \frac{k_Bn}{\gamma-1}\nu_{ie}(T_e - T_i) + Q_i,

where :math:`s` is the field-aligned spatial coordinate, :math:`n` is the density,
:math:`v` is the velocity, :math:`A` is the cross-sectional area, :math:`F_e` and
:math:`F_i` are the electron and ion thermal conductive fluxes, respectively,
:math:`T_e` and :math:`T_i` are the electron and ion temperatures respectively,
:math:`\Lambda` describes the energy lost to radiation, and :math:`Q_e,Q_i` describe
the *ad hoc* heating applied to the electrons and ions, respectively.These equations
are closed by set of equations of state for the electron and ion energies, :math:`E_e,E_i`,
and pressures, :math:`p_e,p_i`,

.. math::
    :label: eqs-of-state

    E_e &= \frac{p_e}{\gamma - 1}, \\
    E_i &= \frac{p_i}{\gamma - 1} + \frac{1}{2}m_inv^2, \\
    p_e &= nk_BT_e, \\
    p_i &= nk_BT_i.

For additional details about these 1D hydrodynamic equations, see :cite:t:`reep_geometric_2022`.

Using the above equations of state coupled with the assumptions that (1) all flows are
subsonic and (2) gravity is negligible for loops of length :math:`<150` Mm
:cite:p:`{see Section 2 of}klimchuk_highly_2008`, we can simplify :eq:`electron-energy-1d`
and :eq:`ion-energy-1d` to,

.. math::
    :label: electron-energy-1d-simple

    \frac{A}{\gamma-1}\frac{\partial p_e}{\partial t} + \frac{\gamma}{\gamma-1}\frac{\partial}{\partial s}(Ap_ev) = &-\frac{\partial}{\partial s}(AF_e) + Av\frac{\partial p_e}{\partial s} - An^2\Lambda \\
    &+ Ak_Bn\nu_{ie}(T_i - T_e) + AQ_e,

.. math::
    :label: ion-energy-1d-simple

    \frac{A}{\gamma-1}\frac{\partial p_i}{\partial t} + \frac{\gamma}{\gamma-1}\frac{\partial}{\partial s}(Ap_iv) = -\frac{\partial}{\partial s}(AF_i) - Av\frac{\partial p_e}{\partial s} + Ak_Bn\nu_{ie}(T_e - T_i) + AQ_i,

We can now apply the combined methodology of both :cite:t:`barnes_inference_2016` and
:cite:t:`cargill_static_2022` to Eqs. :eq:`density-1d`, :eq:`electron-energy-1d-simple`,
and :eq:`ion-energy-1d-simple` to derive the EBTEL equations including both two-fluid
effects and cross-sectional expansion.

The EBTEL Electron Pressure Equation
------------------------------------

To derive the EBTEL electron pressure equation, we begin by taking a spatial integral
over Eq. :eq:`electron-energy-1d-simple` from the base of the corona to the apex of a
semi-circular loop that is symmetric about the apex,

.. math::
    :label: e-energy-coronal-integral

    \frac{A_cL_c}{\gamma-1}\frac{d p_{e,c}}{dt} - \frac{\gamma}{\gamma-1}(Ap_ev)_0 = (AF_e)_0 + A_cL_cQ_{e,c} + A_c\psi_c - A_cR_c,

where :math:`c` denotes an average taken over the coronal portion of the loop, :math:`0`
denotes evaluation at the TR-corona interface, :math:`L_c` is the coronal portion of the
loop half-length, and,

.. math::
    :label: psi-corona

    \psi_c = \frac{1}{A_c}\int_c\mathrm{d}s\,Av\frac{\partial p_e}{\partial s} + \frac{1}{A_c}\int_c\mathrm{d}s\,\frac{Ak_Bn\nu_{ie}}{\gamma-1}(T_i - T_e),

.. math::
    :label: losses-corona

    R_c = \frac{1}{A_c}\int_c\mathrm{d}s\,An^2\Lambda.

Note that the coronal integral :math:`\int_c\mathrm{d}s=\int_{s=L_{TR}}^{s=L}\mathrm{d}s`,
where :math:`L_{TR}` is the length of the TR and :math:`L=L_{TR}+L_c` is the total loop
half-length from the bottom of the TR to the apex of the loop. Additionally, because the
loop is assumed symmetric about the apex and isolated from the lower atmosphere, the
velocity and heat flux terms vanish at those locations.

Similarly, we can integrate Eq. :eq:`electron-energy-1d-simple` over the TR,

.. math::
    :label: e-energy-tr-integral

    \frac{A_{TR}L_{TR}}{\gamma-1}\frac{d p_{e,TR}}{dt} + \frac{\gamma}{\gamma-1}(Ap_ev)_0 = &-(AF_e)_0 + A_{TR}L_{TR}Q_{e,TR} \\
    &+ A_{TR}\psi_{TR} - A_{TR}R_{TR},

where :math:`TR` denotes an average taken over the TR portion of the loop, :math:`\psi_{TR}`
and :math:`R_{TR}` have the same form as Eqs. :eq:`psi-corona` and :eq:`losses-corona`, respectively,
and the TR integral :math:`\int_{TR}\mathrm{d}s=\int_{s=0}^{s=L_{TR}}`.

Following the approach of :cite:t:`cargill_static_2022`, we add Eqs.
:eq:`e-energy-coronal-integral` and :eq:`e-energy-tr-integral` and let :math:`p_e=p_{e,c}=p_{e,TR}`
and :math:`Q_e=Q_{e,c}=Q_{e,TR}` to get the EBTEL electron pressure equation,

.. math::
    :label: ebtel-electron-pressure

    \boxed{\frac{1}{\gamma-1}\frac{dp_e}{dt} = Q_e + \frac{\psi_c}{L_*}\left(1+\frac{A_{TR}\psi_{TR}}{A_c\psi_c}\right) - \frac{R_c}{L_*}\left(1+c_1\frac{A_{TR}}{A_c}\right)},

where :math:`L_* = L_c + (A_{TR}/A_c)L_{TR}` and :math:`c_1=R_{TR}/R_c`. Eq. :eq:`ebtel-electron-pressure` describes the time-evolution of the spatially-averaged electron pressure.

The EBTEL Ion Pressure Equation
------------------------------------

To derive the EBTEL ion pressure equation, we apply the same procedure as above to Eq.
:eq:`ion-energy-1d-simple`. The spatial integral of Eq. :eq:`ion-energy-1d-simple` over
the coronal portion of the loop is,

.. math::
    :label: i-energy-c-integral

    \frac{A_cL_c}{\gamma-1}\frac{d p_{i,c}}{dt} - \frac{\gamma}{\gamma-1}(Ap_iv)_0 = (AF_i)_0 + A_cL_cQ_{i,c} - A_c\psi_c,

and for the TR portion of the loop is,

.. math::
    :label: i-energy-tr-integral

    \frac{A_{TR}L_{TR}}{\gamma-1}\frac{d p_{i,TR}}{dt} + \frac{\gamma}{\gamma-1}(Ap_iv)_0 = -(AF_i)_0 + A_{TR}L_{TR}Q_{i,TR} - A_{TR}\psi_{TR}.

As above, we can add Eqs. :eq:`i-energy-c-integral` and :eq:`i-energy-tr-integral` together
and let :math:`p_i=p_{i,c}=p_{i,TR}` and :math:`Q_i=Q_{i,c}=Q_{i,TR}` to obtain the EBTEL
ion pressure equation,

.. math::
    :label: ebtel-ion-pressure

    \boxed{\frac{1}{\gamma-1}\frac{dp_i}{dt} = Q_i - \frac{\psi_c}{L_*}\left(1 + \frac{A_{TR}\psi_{TR}}{A_c\psi_c}\right)}.

Eq. :eq:`ebtel-ion-pressure` describes the time-evolution of the spatially-averaged ion pressure.

The EBTEL Density Equation
------------------------------------

Lastly, we derive the EBTEL density equation.
We begin by taking a spatial integral of Eq. :eq:`density-1d` over the coronal portion of the loop,

.. math::

    A_cL_c\frac{dn}{dt} - (Anv)_0 = 0

and using the equation of state for :math:`p_e` from Eq. :eq:`eqs-of-state`,

.. math::
    :label: density-coronal-integral

    A_cL_c\frac{dn}{dt} = \frac{(Ap_ev)_0}{k_BT_{e,0}}.

The quantity :math:`(Ap_ev)_0` is the area-weighted electron enthalpy flux at the TR-corona interface.
We can derive an expression for this term by substituting Eq. :eq:`ebtel-electron-pressure` into
Eq. :eq:`e-energy-tr-integral` and doing a lot of algebra,

.. math::
    :label: electron-enthalpy-flux

    \frac{\gamma}{\gamma-1}(Ap_ev)_0 = -\frac{A_{TR}L_cR_c}{L_*}\left(c_1 - \frac{L_{TR}}{L_c}\right) + \frac{A_{TR}L_c\psi_c}{L_*}\left(\frac{\psi_{TR}}{\psi_c} - \frac{L_{TR}}{L_c}\right) - (AF_e)_0.

Similarly, we can derive an expression for the area-weighted ion enthalpy flux by substituting
Eq. :eq:`ebtel-ion-pressure` into Eq. :eq:`i-energy-tr-integral`,

.. math::
    :label: ion-enthalpy-flux

    \frac{\gamma}{\gamma-1}(Ap_iv)_0 = -\frac{A_{TR}L_c\psi_c}{L_*}\left(\frac{\psi_{TR}}{\psi_c} - \frac{L_{TR}}{L_c}\right) - (AF_i)_0.

Adding Eqs. :eq:`electron-enthalpy-flux` and :eq:`ion-enthalpy-flux`,

.. math::
    :label: total-enthalpy-flux

    \frac{\gamma}{\gamma-1}(Ap_ev)_0 + \frac{\gamma}{\gamma-1}(Ap_iv)_0 = -\frac{A_{TR}L_cR_c}{L_*}\left(c_1 - \frac{L_{TR}}{L_c}\right) - (AF_e)_0 - (AF_i)_0.

Additionally, we define :math:`\xi\equiv T_e/T_i` and again use the electron and ion equations
of state from Eq. :eq:`eqs-of-state` to find,

.. math::
    :label: temperature-ratio

    \xi = \frac{T_e}{T_i} = \frac{T_{e,0}}{T_{i,0}} = \frac{A_0n_0k_BT_{e,0}v_0}{A_0n_0k_BT_{i,0}v_0} = \frac{(Ap_ev)_0}{(Ap_iv)_0}.

Substituting this into Eq. :eq:`total-enthalpy-flux` gives,

.. math::
    :label: electron-enthalpy-flux-simple

    (Ap_ev)_0 = -\frac{\xi(\gamma - 1)}{\gamma(\xi + 1)}\left(\frac{A_{TR}L_cR_c}{L_*}\left(c_1 - \frac{L_{TR}}{L_c}\right) + (AF_e)_0 - (AF_i)_0\right).

Finally, substituting Eq. :eq:`electron-enthalpy-flux-simple` into Eq. :eq:`density-coronal-integral`
yields the EBTEL density equation,

.. math::
    :label: ebtel-density

    \boxed{\frac{dn}{dt} = -\frac{(\gamma - 1)\xi c_2}{(\xi + 1)\gamma c_3 k_B L_c T_e}\left(\frac{A_{TR}L_c}{A_cL_*}R_c\left(c_1 - \frac{L_{TR}}{L_c}\right) + \frac{A_0}{A_c}(F_{e,0} + F_{i,0})\right)},

where :math:`c_2,c_3` are constants relating the base temperature to the spatially-averaged coronal
temperature. Eq. :eq:`ebtel-density` describes the time-evolution of the spatial averaged density.
In summary, Eqs. :eq:`ebtel-electron-pressure`, :eq:`ebtel-ion-pressure`, and :eq:`ebtel-density`
comprise the two-fluid equations including cross-sectional area expansion.

Limiting Behavior
-----------------

Below, we briefly describe how Eqs. :eq:`ebtel-electron-pressure`, :eq:`ebtel-ion-pressure`, and
:eq:`ebtel-density` reduce to the other versions of the EBTEL model.

Constant Cross-section
++++++++++++++++++++++

Under the assumption of constant cross-section, :math:`A_c=A_{TR}=A_0` and :math:`L_{TR}\ll L_c`
such that :math:`L_*\approx L_c = L`. As such, the EBTEL equations simplify to,

.. math::

    \frac{dp_e}{dt} &= (\gamma-1)Q_e + \frac{\gamma-1}{L}(\psi_c+\psi_{TR} - R_c(1+c_1)), \\
    \frac{dp_i}{dt} &= (\gamma-1)Q_i - \frac{\gamma-1}{L}(\psi_c + \psi_{TR}), \\
    \frac{dn}{dt} &= -\frac{(\gamma - 1)\xi c_2}{(\xi + 1)\gamma c_3 k_B L_c T_e}(c_1R_c + F_{e,0} + F_{i,0}).

Note that these expressions are equivalent to the two-fluid EBTEL equations as given in :cite:t:`barnes_inference_2016`.

Single-fluid
++++++++++++

Under the single-fluid assumption, :math:`T_e=T_i` at all times.
Using this and adding Eqs. :eq:`ebtel-electron-pressure` and :eq:`ebtel-ion-pressure`,

.. math::

    \frac{dp}{dt} &= (\gamma - 1)\left(Q - \frac{R_c}{L_*}\left(1 + \frac{A_{TR}}{A_c}c_1\right)\right), \\
    \frac{dn}{dt} &= -\frac{(\gamma - 1)c_2}{2\gamma c_3 k_B L_c T}\left(\frac{A_{TR}L_c}{A_cL_*}R_c\left(c_1 - \frac{L_{TR}}{L_c}\right) + \frac{A_0}{A_c}F_0\right).

Note that these expressions are equivalent to the expanding cross-section EBTEL equations given
in :cite:t:`cargill_static_2022`.
