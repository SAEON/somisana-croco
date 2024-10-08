Presentation
-------------

**croco_pytools** are available on gitlab: https://gitlab.inria.fr/croco-ocean/croco_pytools

.. warning:: 

    These tools bring together the methods of several
    users of the CROCO community but do not constitute a definitive toolbox.
    All the features present in croco_tools
    matlab are not available in this version.

    A new, more complete toolbox is being created by the team
    development of the CROCO group.

**croco_pytools/prepro** remains in the footsteps of the matlab **croco_tools** by separating all the steps involveed in creating a simulation. To avoid having a long namelist file, and each routines being independant, parameters need to be specified in the header of each routine.

Steps for creating a configuration (with these routines) are:

* Build the grid
* Build the lateral boundary conditions (3D currents, temperature and salinity, barotropic currents, surface elevation)
* Build the initial conditions

And eventually:

* Build tidal forcing
* Build river forcing

Contact
    mathieu.le.corre@shom.fr

Structure of directories
^^^^^^^^^^^^^^^^^^^^^^^^

* ``env``

  Contains ``yml`` files to install python environment decidated to ``croco_pytools/prepro``

* ``Modules``

  Contains all python routines to run the ``croco_pytools/prepro``
  Fortran routines are in the sub-drectory ``tools_fort_routines``

* ``Readers``

  This contains ``readers`` to decode the input datasets

Description of the routines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The main directory contains the following files:

.. list-table::
   
   * - __init__.py
     - File containing informations to build python environment and compile fortran routines
   * - install.py
     - Script to facilitate the installation of python environnement with conda and the compilation of fortran routines
   * - make_grid.py
     - Script to build CROCO grid and associated nests
   * - make_bry.py
     - Script to build the lateral boundary conditions (surface elevation, 3D currents, barotropic currents, temperature and salinity, other tracers)
   * - make_ini.py
     - Script to build the initial 3D conditions
   * - make_zoom_cond.py
     - Script to build initial and lateral boundary conditions for nested grid (offline nest or AGRIF nest)
   * - make_tides.py
     - Script to build the tidal forcing (amplitude and phase) for elevation and barotropic current
   * - make_rivers.py
     - Script to create netcdf file containing runoff flows

Dependencies
------------

Python
^^^^^^

**croco_pytools/prepro** is using the following package:

.. list-table::
     :widths: 10 90
  
     * - `matplotlib <https://matplotlib.org/>`_
       - Matplotlib is a comprehensive library for creating static, animated,
         and interactive visualizations in Python.
     * - `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_
       - Cartopy is a Python package designed for geospatial data processing 
         in order to produce maps and other geospatial data analyses.   
     * - `wxpython <https://wxpython.org/>`_
       - wxPython is a cross-platform GUI toolkit for the Python 
         programming language. 
     * - `geopandas <https://geopandas.org/en/stable/>`_
       - GeoPandas is an open source project to make working with geospatial 
         data in python easier.
     * - `regionmask <https://regionmask.readthedocs.io/en/stable/>`_
       - Create masks of geographical regions
     * - `pyinterp <https://pangeo-pyinterp.readthedocs.io/en/latest/>`_
       - Python library for optimized geo-referenced interpolation.
     * - `pandas <https://pandas.pydata.org/>`_
       - pandas is a fast, powerful, flexible and easy to use open source
         data analysis and manipulation tool, built on top of the
         Python programming language.
     * - `scipy <https://www.scipy.org/scipylib/index.html>`_
       - Scipy provides many user-friendly and efficient numerical routines,
         such as routines for numerical integration, interpolation,
         optimization, linear algebra, and statistics.
     * - `xarray <http://xarray.pydata.org/en/stable/>`_
       - xarray is an open source project and Python package that makes working
         with labelled multi-dimensional arrays simple, efficient, and fun

Fortran
^^^^^^^

The scripts used to build the CROCO grid and the initial/boundary 
conditions for nests are using fortran routines which have been 
interfaced with python through `f2py <https://numpy.org/doc/stable/f2py/>`_
and must therefore be compiled.

Before compiling make sure that you have the following:

* Open MP-capable Fortran compiler, Ifort or Gfortran, may be others.
* NetCDF library capable of handling netCDF-4/hdf5 format.

