Build CROCO tidal forcing
-------------------------

CROCO is able to propagate the different tidal constituents from its lateral 
boundaries by forcing tidal signal both in elevation and velocity. 

.. note:: 
  
  To get a clean signal you need to provide harmonic components from both tide 
  elevation and tide velocity. In case you don’t have velocity harmonics 
  (``undef UV_TIDES`` in ``cppdefs.h``) a set of reduced equations is 
  available to compute velocity from SSH (``OBC_REDUCED_PHYSICS``)
  
Description of make_tides
^^^^^^^^^^^^^^^^^^^^^^^^^

``make_tides`` creates a forcing file containing the amplitude and phase 
of each of the desired tidal components for tide elevation and currents. 
These values are interpolated over the entire CROCO grid, but the simulation 
only uses the values at boundaries for forcing.  

In most cases, the simulations are only forced by the main tidal components. 
To take disturbance waves into account, a nodal correction is applied to these 
main waves. This nodal correction is of low amplitude and only slightly 
variable (considered constant over a year) but can have a significant impact 
on the tidal solution in coastal areas. 

CROCO's common practice is to set the nodal correction value at the simulation 
start date and keep it constant thereafter. In this case, the correction is 
directly included in the forcing file. However, CROCO also has the ability to 
calculate these corrections at each time step in the simulation with the cppkey 
``TIDES_MAS``. In this case, the nodal correction is not included in the 
forcing file. ``make_tides`` manages these different cases.

.. note::

  The ``TIDES_MAS`` key applies nodal correction only on tidal elevation. 
  In the current state, it is therefore necessary to keep the nodal correction 
  on currents in the forcing file.


Fill the Reader
^^^^^^^^^^^^^^^

The associated reader, ``Readers/tides_reader.py``, needs to be filled.
Two kind of input data exist, either providing Real/Imaginary 
part for each wave (like `TPXO <https://www.tpxo.net/global>`_) or providing 
Amplitude/Phase (like 
`FES <https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes.html>`_).

As for the initial conditions, the reader contains information on the 
spatial coordinates at the rho,u,v points, and on the different variables.
For each variable, the first letter indicates the variable, H for elevation, (U,V) for respectively eastward and northward        
currents; and the indices r,i,a,p indicate whether the data are real/imaginary                
or amplitude/phase type.  

For TPXO, which gives elevation and transport data in   
real/imaginary format, the reader writes:

::

  if input == 'tpxo':
        dico={ 'lonr':'lon_z','lonu':'lon_u','lonv':'lon_v',\
               'latr':'lat_z','latu':'lat_u','latv':'lat_v',\
               'H_r':'hRe',\
               'H_i':'hIm',\
               'U_r':'uRe',\
               'U_i':'uIm',\
               'V_r':'vRe',\
               'V_i':'vIm'\
               'topor':'h','topou':'hu','topov':'hv'\
              }


As TPXO does not provide meridional and zonal current transport, an 
extra line is added to read the TPXO topography to be able to compute the transport. 


Using make_tides
^^^^^^^^^^^^^^^^

Let's build a forcing file from TPXO7. This dataset is 
available in the `DATASETS_CROCO <https://data-croco.ifremer.fr/DATASETS/TPXO7.tar.gz>`_  
. The reader for this file is:

::
  
    elif input == 'tpxo7':
        dico={ 'lonr':'lon_r','lonu':'lon_u','lonv':'lon_r',\
               'latr':'lat_r','latu':'lat_r','latv':'lat_v',\
               'H_r':'ssh_r',\
               'H_i':'ssh_i',\
               'U_r':'u_r',\
               'U_i':'u_i',\
               'V_r':'v_r',\
               'V_i':'v_i',\
               'topor':'h','topou':'h','topov':'h'\
              }

Before running ``make_tides.py``, **USER CHANGES** section needs to be filled. 
There are several parts in it:

* **Dates**:
  ::

    # Dates
    # Origin year
    Yorig = 2000 # 1900 if TIDES_MAS defined in cppdef.h
    # Initial date
    Yini, Mini, Dini = 2013, 1, 1
  
* **Input data informations**
  ::
  
    # Input data information and formating
    # Note: if you are using a tpxo dataset please be sure to have somewhere in 
    #       inputdata 'tpxo'. This will allow the code to use the OTIS (TPXO is obtained with it)
    #       convention a-b*i for complex.
    #       Also, if you have already preprocess TPXO, make sure that you have correct units and 
    #       u,v are in m/s and not m**2/s
    inputdata = 'tpxo7_croco' # Input data dictionnary as defined in the Readers/tides_reader.py
    input_dir = '../../DATASETS_CROCOTOOLS/TPXO7/'
    input_file = 'TPXO7.nc' # Leave empty if you have multiple files
    input_type = 'Re_Im' # Format of the input data 'Amp_phase'or 'Re_Im'
    multi_files  = False # Set to True if several input files
    if multi_files:
        waves_separated = True # Set to True if input files waves are separated
        elev_file = 'h_<tides>_tpxo9_atlas_30_v5.nc' # elevation file names. if wave_separated put <tides> where wave name is found
        u_file = 'u_<tides>_tpxo9_atlas_30_v5.nc' # eastward currents file names. if wave_separated put <tides> where wave name is found
        v_file = 'u_<tides>_tpxo9_atlas_30_v5.nc' # northward currents file names. if wave_separated put <tides> where wave name is found
 

Here we select the reader for the ``tpxo7`` data. We also select input file
location and data format (Re/Im or Amp/Pha). Elevation and current data may
not be in the same file. The ``multi_files`` option is then useful for
specifying each of them. If Eastward and Northward components are in the same 
file put the same name in ``u_file`` and ``v_file``. It is also possible to 
have waves that are in different files using the ``<tides>`` key, which will 
be replaced by the wave list specified below. 

.. note:: 

  TPXO follows complex convention a-b*i. You therefore need to pay attention 
  to the sign when calculating the phase. For the scripts to correctly 
  take this convention into account, you need to ensure that ``inputdata`` 
  contains the ``tpxo`` characters. 
  
* **CROCO grid informations**
  ::

    # CROCO grid informations
    croco_dir = '../../CROCO_FILES/'
    croco_grd = 'croco_grd.nc'

Informations about your CROCO grid. Indicate the path (``croco_dir``) and
the input grid to use (``croco_grd``).

* **Tide file and settings**:
  ::
  
    # Tide file informations
    croco_filename = 'croco_frc.nc'
    tides = ['M2','S2','N2','K2','K1','O1','P1','Q1','Mf','Mm']
  
    cur = True # Set to True if you to compute currents
    pot = False # Set to True if you to compute potiential tides
  
    # Nodal correction
    Correction_ssh = True
    Correction_uv = True
 
Contains informations about the output file name(``croco_filename``) and which
tidal components it will contain.

You can also choose whether you want to calculate
currents (``cur``) and/or the generating potential (``pot``).

If you want to apply nodal correction on elevation 
and/or current. 

.. note::

  As previously said, CROCO can compute nodal correction for elevation. In this 
  case you must set ``Correction_ssh = False``, choose (``Yorig = 1900`` and ``define TIDES_MAS`` in ``cppdefs.h``)


To use ``make_tides.py`` do:
::

  python make_tides.py
