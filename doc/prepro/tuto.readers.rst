Readers for input data
-----------------------

Depending on the input data you are using, you might have to add a reader to read them. 
Several types of data are already managed, but you can any other input data by adding it to the readers in the ``Readers`` directory.

::

    cd Readers

Three sets of readers exist:

.. list-table::

  * - topo_reader.py
    - Reader for topographic/bathymetric data (etopo5, etopo2, etopo1, srtm30, homonim, gebco...)
  * - ibc_reader.py
    - Reader for initial and boundary conditions datasets (mercator, eccov4, soda...)
  * - tides_reader.py
    - Reader for tides datasets (tpxo,fes2014...)

They basically contain a routine, ``lookvar``, to associate the input variable names to commonly used variable names in the croco_pytools with a dictionnary.

In the following we describe the three type of reader:


.. note::
  
 To access variable name in a netcdf, ncdump is your friend


Topography
^^^^^^^^^^

To build bathymetry from a dataset, the dictionnary needs at least four keys:

.. list-table::

  * - lon
    - Longitude 
  * - lat
    - Latitude
  * - topo
    - Name of the 2-D bathymetric field
  * - zaxis
    - Direction of z axis (up or down)

A supplementary key ``srtm`` needs to be added and set to ``True`` if the input 
dataset is `STRM30plus <https://topex.ucsd.edu/WWW_html/srtm30_plus.html>`_

Initial and boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Reader for initial and boundary conditions (ibc) handles 2D and 3D fields and their evolution in time. It requires their coordinates:

.. list-table::

    * - depth
      - Z-coordinates variables 
    * - lonr, lonu, lonv
      - Longitude name for Rho,U,V grid. If only one grid put the same name 
        everywhere
    * - latr, latu, latv
      - Latitude name for Rho,U,V grid. If only one grid put the same name 
        everywhere
    * - time, time_dim
      - Time variable and the associated time dimension

And the following fields:

.. list-table::

    * - ssh
      - Sea Surface Elevation variable name
    * - u
      - Eastward velocity (m/s) variable name
    * - v
      - Northward velocity (m/s) variable name
    * - temp
      - temperature (C) variable name
    * - salt
      - salinity (PSU) variable name

You can also add other tracers if needed following the same pattern::

    'tracer':'tracer_name_in_dataset'

Tides
^^^^^

As for ``ibc_reader``, you need to specify the keys for geographical coordinates

.. list-table:: 

   * - lonr, lonu, lonv
     - Longitude name for Rho, U, V grid. If only one grid put the same name 
       everywhere
   * - latr, latu, latv
     - Latitude name for Rho, U, V grid. If only one grid put the same name 
       everywhere

Tide data consist in three fields (elevation, easwtward velocity, northward velocity), and each field has two components. Depending on the input dataset, those two components can be either Amplitude/Phase or Real/Imaginary part.

If your dataset follows an Amplitude/Phase, keys to declare are:

.. list-table::

   * - H_a
     - Elevation amplitude
   * - H_p
     - Elevation phase
   * - U_a
     - Eastward velocity amplitude
   * - U_p
     - Eastward velocity phase
   * - V_a
     - Northward velocity amplitude
   * - V_p
     -  Northward velocity phase

Or if it follows a Real/Imaginary format, keys to declare are:

.. list-table::

   * - H_r
     - Elevation real part
   * - H_i
     - Elevation imaginary part
   * - U_r
     - Eastward velocity real part
   * - U_i
     - Eastward velocity imaginary part
   * - V_r
     - Northward velocity real part
   * - V_i
     - Northward velocity imaginary part

