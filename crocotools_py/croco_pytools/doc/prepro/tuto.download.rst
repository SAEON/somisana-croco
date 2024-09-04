Download some datasets
----------------------

.. note::

  Python API used here are included in conda environnement ``croco_pyenv``

Download Glorys dataset
^^^^^^^^^^^^^^^^^^^^^^^

``download_glorys_data_copernicus_cli.py`` allows to download dataset from 
CMEMS website using the new Copernicus Marine API. Further informations on the 
API is given here:
https://help.marine.copernicus.eu/en/articles/7949409-copernicus-marine-toolbox-introduction

A CMEMS account is needed in order to retreive the data. Your login and 
password are asked each time you want to use this script. 

This script takes advantage of "subset" function of copernicusmarine API which 
allows to download a subset of the original files.

.. warning::

  When using ``subset`` functionnality, dataset is converted to Analysis-Ready 
  Cloud-Optimized (ARCO) format. Data are the same, except for time variable 
  where date is centered at midnight instead of noon for original files. For 
  more intel on the difference between netcdf and ARCO format go to this 
  `page <https://help.marine.copernicus.eu/en/articles/8656000-differences-between-netcdf-and-arco-formats>`_

  To use ``download_glorys_data_copernicus_cli.py``, **USER CHANGES** section 
  needs to be filled

* **Dates**:
  ::

    Ystart,Mstart = 2013,1 # Starting month
    Yend,Mend = 2013,3    # Ending month 

Representing the starting and ending month of you simulation

* **Product informations**:
  ::

    product_id = 'cmems_mod_glo_phy_my_0.083deg_P1M-m'
    variables = ["zos","uo","vo","thetao","so"]
    multi_files = False

The desired product type is specified with ``product_id`` along the variables you 
want to download. In case of large dataset, ``multi_files`` option allows to 
put each variable in its own netcdf.

List of all the product is available on `Copernicus website <https://data.marine.copernicus.eu/products>`_  
but here a list of common used datasets:

.. list-table::

  * - Dai9ly Reanalysis
    - ``cmems_mod_glo_phy_my_0.083deg_P1D-m``
  * - Monthly Reanalysis
    - ``cmems_mod_glo_phy_my_0.083deg_P1M-m``
  * - Daily Analysis and Forecast 
    - ``cmems_mod_glo_phy_anfc_0.083deg_P1D-m``

* **CROCO informations**:
  ::

    croco_dir = '../'
    croco_grd = 'croco_grd.nc'

* **Output informations**:
  ::

    output_dir = './GLORYS_DATA/'
    output_prefix = 'glo12-reana-daily'

After completing all fields, lauch the script::

    python doawnload_glorys_data_copernicus_cli.py

Login and password will be asked and data downloaded on put inside ``output_dir``.


Download river discharge
^^^^^^^^^^^^^^^^^^^^^^^^

``download_glofas_river.py`` makes possible to download data for rivers discharges.
Data are stored on `Climate Data Store website <https://cds.climate.copernicus.eu/cdsapp#!/home>`_ 
It uses `CDS API <for data request <https://cds.climate.copernicus.eu/api-how-to#>`_ for data request. 

.. note:: 

    Follow instructions on `CDS API <https://cds.climate.copernicus.eu/api-how-to#>`_ 
    to create the file ``$HOME/.cdsapirc`` This file should contains your login and 
    parssword informations

River discharge can be found `here <https://cds.climate.copernicus.eu/cdsapp#!/dataset/cems-glofas-historical?tab=overview>`_ 
and is retreive from the Global Flood Awareness System (GloFAS).

.. note::

    User also need to have CDO to convert grib to netcdf format

To use ``download_glofas_river.py``, **USER CHANGES** section 
needs to be filled:

* **Dates**:
  ::

    Ystart,Mstart = 2013,1 # Starting month
    Yend,Mend = 2013,3    # Ending month 

Representing the starting and ending month of you simulation

* **Product informations**:
  ::

    product_id = 'cems-glofas-historical'

Other dataset of interest can be found on the website, look around and just put 
the desired ``product_id``

 **CROCO informations**:
  ::

    croco_dir = '../'
    croco_grd = 'croco_grd.nc'

* **Output informations**:
  ::

    convert2netcdf = True
    output_dir = './DATA_RIVER/'
    output_name = 'cems_glofas'

If ``convert2netcdf`` is set to ``False`` data will be in grib format and no 
conversion is done.

To use this script use the command::

    python download_glofas_river.py


