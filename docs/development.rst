.. _ebtelplusplus-development:

Contributing
============

If you would like to contribute to ``ebtelplusplus``, either by fixing a bug or adding a feature, you will need to download and build the package locally.
If you are new to contributing to open source, we recommend reading the `SunPy newcomer's guide <https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html>`__.
The section on `setting up a development environment <https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html#setting-up-a-development-environment>`__ will be particularly helpful.
While these instructions are particular to SunPy, they are applicable to any open-source project, particularly for Python.

First, create a fork of the `ebtelplusplus`_ repository.
If you're not sure what this is or how to create a fork, see the aforementioned development guide.
Next, clone your fork,

.. code:: shell

    git clone https://github.com/<your-user-name>/ebtelPlusPlus.git

Because ``ebtelplusplus`` implements the actual simulation code in C++ and thus is not a pure Python package, installation of this package requires first compiling the C++ and building the needed binaries.
To do this, you will first need to have the `Boost <http://www.boost.org/>`__ package installed (at least v1.53).
There are various ways to install this package, including:

.. list-table::
    :header-rows: 1

    * - Source
      - OS
      - Command
    * - `conda-forge <https://github.com/conda-forge/boost-feedstock>`__
      - any
      - ``conda install -c conda-forge libboost-devel``
    * - `Homebrew <https://formulae.brew.sh/formula/boost>`__
      - macOS
      - ``brew install boost``
    * - `Macports <https://ports.macports.org/port/boost/>`__
      - macOS
      - ``sudo port install boost``
    * - `Chocolatey <https://community.chocolatey.org/packages/boost-msvc-14.3>`__
      - Windows
      - ``choco install boost-msvc-14.3``

The reason that Boost is required is that ``ebtelplusplus`` uses the ``Boost.Numeric.Odeint`` package to solve the EBTEL equations.
Once you have successfully installed Boost, you can install ``ebtelplusplus``,

.. code:: shell

    cd ebteplusplus
    pip install -e .[dev]

This will compile the C++ code (using `pybind11 <https://pybind11.readthedocs.io/en/stable/index.html>`__ and
`scikit-build-core <https://scikit-build-core.readthedocs.io/en/latest/>`__) and install the ``ebtelplusplus``
package, including all of the dependencies for testing and developing the package.

To verify that you have a working installation, you can run the test suite,

.. code:: shell

    pytest ebtelplusplus

You are now ready to start contributing to or modifying ``ebtelplusplus``.
Note that if you are modifying any of the simulation code, that is contained in ``ebtelplusplus/extern`` and is implemented in C++.
After each change, you will need to recompile using the above installation command.

For more information about proposing new changes back to the main repository, see `this section of the SunPy newcomer's guide <https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html#send-it-back-to-us>`__.

.. _ebtelplusplus: https://github.com/rice-solar-physics/ebtelPlusPlus
