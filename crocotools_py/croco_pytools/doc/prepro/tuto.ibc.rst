Build CROCO initial and boundary conditions
--------------------------------------------

This section presents the use of ``make_ini.py`` and ``make_bry.py``
to create initial and boundary condtions (ibc) for
realistic regional configurations with CROCO. 

The different steps to build ibc are:

#. Extrapolating values on mask
#. Interpolating each input for each Z-level on the CROCO horizontal grid 
#. Interpolating on the vertical sigma-coordinates

Fill the Reader
^^^^^^^^^^^^^^^

In this tutorial, we will use as initial and boundary conditions the 
`Mercator product GLORYS12 <https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description>`_, but already processed 
in a CROCO-like format by tools available for downloading Mercator data in croco_tools.  
An extraction from January 2013 to March 2013 is available in the tutorial section of CROCO website: https://www.croco-ocean.org/trainings/, in the section 'Regional configuration (#define REGIONAL in cppdefs.h)' under 'Global datasets needed for preprocessing'. Here is the link to download them: https://data-croco.ifremer.fr/CONFIGS_EXAMPLES/DATASETS_GLOB_INTER/CROCO/MERCATOR_GLOB_201301.tar.gz
The following 3 files are provided:
``mercator_Y2013M1.cdf, mercator_Y2013M2.cdf, mercator_Y2013M3.cdf``.
These files are not "original" mercator files, and we will use the ``mercator_croco`` keyname to differentiate their formatting. 

You can check the variables name from these input data by doing:
::

  ncdump -h mercator_Y2013M1.cdf

You thus, first need to check/edit the ibc reader, ``Readers/ibc_reader.py``, for these data. 
Here we use the ``mercator_croco`` style of formatting:

::

    if input == 'mercator_croco':
        dico={ 'depth':'depth',\
               'lonr':'lonT','lonu':'lonU','lonv':'lonV',\
               'latr':'latT','latu':'latU','latv':'latV',\
               'ssh':'ssh',\
               'temp':'temp',\
               'salt':'salt',\
               'u': 'u',\
               'v': 'v',\
               'time': 'time',\
               'time_dim':'time'\
             }


.. note:: 
    
    If you have downloaded mercator data on your own, with a more "mercator-native" format, you may rather use the ``mercator`` dictionnary:
    ::
 
        if input == 'mercator':
            dico={ 'depth':'depth',\
                'lonr':'longitude','lonu':'longitude','lonv':'longitude',\
                'latr':'latitude','latu':'latitude','latv':'latitude',\
                'ssh':'zos',\
                'temp':'thetao',\
                'salt':'so',\
                'u': 'uo',\
                'v': 'vo',\
                'time': 'time',\
                'time_dim':'time'\
                 }

Using make_ini
^^^^^^^^^^^^^^

Before running ``make_ini.py``, **USER CHANGES** section needs to be filled. 
There are several parts in it:

* **Dates**:

::

  # Dates
  # starting date
  Yini, Mini, Dini  = '2013','01','01' # Month and days need to be 2-digits format
  # reference time (default = ini time)
  Yorig, Morig, Dorig = Yini, Mini, Dini # Month and days need to be 2-digits format

.. note:: 

    Origin time and initial time can be different.

* **Input data informations** contain the general informations about the inputs (reader keyword, directory, prefix)
  ::
  
    # Input data information and formating
    inputdata = 'mercator_croco' # Input data dictionnary as defined in the Readers/ibc_reader.py
    input_dir = '../../MERCATOR_GLOB_2013/'
    input_prefix='mercator_'
  
    input_file  = f'{input_dir}{input_prefix}Y2013M1.cdf'
    multi_files=False # If variables are in different netcdf
    if multi_files: # Mutiple files
        input_file = { 'ssh'  : input_dir + input_prefix + 'ETAN.Y2013M01.nc',\
                       'temp' : input_dir + input_prefix + 'THETA.Y2013M01.nc',\
                       'salt' : input_dir + input_prefix + 'SALT.Y2013M01.nc',\
                       'u'    : input_dir + input_prefix + 'EVEL.Y2013M01.nc',\
                       'v'    : input_dir + input_prefix + 'NVEL.Y2013M01.nc'\
                    }
  
  # time index to use in the file
  tndx = 0

``inputdata`` needs to be set to the keyword to use in the reader, as defined in ``Readers/ibc_reader.py``.

For input source in which variables are in different netcdf files, you can set ``multi_files`` to ``True``, then you
can define one netcdf file for each variable in the following section.

``input_file`` is the input filename (full path + filename). In the example, it is defined with the ``input_dir``, ``input_prefix``, and
simulation start date (``Yini,Mini,Dini``).

If input data files contains several times, the user must select the time index for starting the simulations with
``tndx`` to select one (0 means first index).

``Nzgoodmin`` defines a threshold above which a z-level in the input file is valid. If at a certain depth this threshold is not reached, the level will be considered as not intersecting the ocean.

* **tracers**
  ::

    # tracers
    tracers = ['temp','salt']
  
Set tracer names here, if any. Name defined here must also be defined in your
reader.

* **CROCO grid informations**
  ::

    # CROCO grid informations
    croco_dir = '../../CROCO_FILES/'
    croco_grd = 'croco_grd.nc'
    sigma_params = dict(theta_s=7, theta_b=2, N=32, hc=200) # Vertical streching, sig_surf/sig_bot/ nb level/critical depth

Informations about your CROCO grid. Indicate the path (``croco_dir``),
the input grid to use (``croco_grd``), the parameters for the
`sigma-coordinates <https://croco-ocean.gitlabpages.inria.fr/croco_doc/model/model.grid.html>`_
  
* **Ini filename to be generated** (it will follow the pattern indicated + the date):

::

  # Ini file informations
  ini_filename = 'croco_ini.nc' # output will be put in croco_dir by default

Name of your output file. This file will be written the previoulsy set ``croco_dir``.

To use ``make_ini.py``, do:
::

  python make_ini.py

You will find the ouput in ``croco_dir``::

    croco_ini_mercator_croco_Y2013M01.nc

Using make_bry
^^^^^^^^^^^^^^

``make_bry`` is quite similar to ``make_ini``, but you will have serveral time frames in the output bry file. 
As for ``make_ini``, **USER CHANGES** section of ``make_bry`` needs to be filled.

::

  # Dates
  Yorig = 2013                    # year defining the origin of time as: days since Yorig-01-01
  Ystart, Mstart = '2013', '01'   # Starting month
  Yend, Mend  = '2013','03'       # Ending month 
  
  # Input data information and formating
  inputdata = 'mercator_croco'    # Input data dictionnary as defined in the Readers/ibc_reader.py
  input_dir = '../../MERCATOR_GLOB_2013/'
  input_prefix = 'mercator_*'  # Please use * to include all files
  multi_files = False
  if multi_files: # Multiple data files. Time is read in ssh file
      input_file = {'ssh':sorted(glob.glob(input_dir+input_prefix+'ETAN.*.nc')),\
                    'temp':sorted(glob.glob(input_dir+input_prefix+'THETA.*.nc')),\
                    'salt':sorted(glob.glob(input_dir+input_prefix+'SALT.*.nc')),\
                    'u':sorted(glob.glob(input_dir+input_prefix+'EVEL.*.nc')),\
                    'v':sorted(glob.glob(input_dir+input_prefix+'NVEL.*.nc'))\
                  }
  else:  # glob all files
      input_file  = sorted(glob.glob(input_dir + input_prefix))
  
  # default value to consider a z-level fine to be used
  Nzgoodmin = 4
  
  # Tracers
  tracers = ['temp', 'salt']
  
  # CROCO grid informations
  croco_dir = '../../CROCO_FILES/'
  croco_grd = 'croco_grd.nc'
  sigma_params = dict(theta_s=7, theta_b=2, N=32, hc=200) # Vertical streching, sig_surf/sig_bot/ nb level/critical depth
  
  # Bry file informations
  bry_filename = 'croco_bry.nc' # output will be put in croco_dir by default
  obc_dict = dict(south=1, west=1, east=1, north=1) # open boundaries (1=open , [S W E N])
  output_file_format = "MONTHLY" # How outputs are spit (MONTHLY,YEARLY,FULL)
  cycle_bry = 0.

In bry case , several days/months need to be given. To facilitate input 
selection we use python `glob module <https://docs.python.org/3/library/glob.html>`_ 
which finds all the pathnames matching a specified pattern according to the 
rules used by the Unix shell. 
In ``input_prefix``, select a specified 
pattern and use Unix shell rules (\*, \$) to select all files following it. 

In the last part, 'Bry file informations', output file name (``bry_filename``) is defined. Select 
which boundary to open (by putting 1 to the corresponding boundary in 
``obc_dict``). Several formats (``output_file_format``) exist to gather data 
by month, year or put them all in the same file.

``cycle`` is used if you want to create a cycle in number of days on bry conditions. Typically, to have a repetition of a year of boundary forcing ``cycle`` would be set to 365

Origin year is also defined in this part, along with starting and ending date.

To use ``make_bry.py`` do:
::

  python make_bry.py


