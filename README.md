# ebtelplusplus

[![CI Status](https://github.com/rice-solar-physics/ebtelplusplus/actions/workflows/ci.yml/badge.svg)](https://github.com/rice-solar-physics/ebtelPlusPlus/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/ebtelplusplus/badge/?version=latest)](https://ebtelplusplus.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/rice-solar-physics/ebtelplusplus/branch/main/graph/badge.svg?token=8G5H9T5AAH)](https://codecov.io/gh/rice-solar-physics/ebtelplusplus)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12675386.svg)](https://doi.org/10.5281/zenodo.12675386)

`ebtelplusplus` is an implementation of the enthalpy-based thermal evolution of loops (EBTEL) model for doing
efficient hydrodynamics of dynamically-heated solar coronal loops.
`ebtelplusplus` decouples the electron and ion energy equations such that the two populations can evolve separately.
This implementation also includes effects to due to cross-sectional area expansion.

If you are looking for the original EBTEL implementation, the you can find the [repository for the IDL code here](https://github.com/rice-solar-physics/EBTEL).

## Installation

The easiest way to install `ebtelplusplus` is through `pip`,

```shell
pip install ebtelplusplus
```

If you would like to compile and build the package from source, see [the instructions here](https://ebtelplusplus.readthedocs.org/en/latest/development.html).

## Usage

The code snippet below shows how to set up a simulation for a 40 Mm loop, lasting 2 hours, heated by a single
heating event lasting 200 s in which all of the energy is injected into the electrons,

```python
import astropy.units as u
import ebtelplusplus
from ebtelplusplus.models import HeatingModel, TriangularHeatingEvent

heating = HeatingModel(
    background=1e-6*u.Unit('erg cm-3 s-1'),
    partition=1,
    events=[TriangularHeatingEvent(0*u.s, 200*u.s, 0.1*u.Unit('erg cm-3 s-1'))]
)
results = ebtelplusplus.run(2*u.h, 40*u.Mm, heating)
```

## Citation

If you use `ebtelplusplus` in any published work, it is greatly appreciated if you follow the [citation instructions here](https://ebtelplusplus.readthedocs.io/en/latest/index.html#citation).
