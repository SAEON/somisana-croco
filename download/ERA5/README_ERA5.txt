This readme and scripts are copied from the official croco_tools-v1.3.1
The original download scripts prompts you to download data into a model configuration dir
Our workflow is to rather downloading the ERA5 data into the DATASETS_CROCOTOOLS directory,
located in the root directory of this repo. This is so that
multiple croco configurations can make use of the same downloaded files

The easiest approach may be to just copy all the files in this directory  into a DATASETS_CROCOTOOLS/ERA5/ directory
and edit the era5_crocotools_param.py file appropriately

DATASETS_CROCOTOOLS is in .gitignore so nothing in there will be tracked in the repo

everything below is unchanged from the original readme
----------------------------------------------------------

Here is a set of matlab and python routines allowing to :

1- Download native ERA-5 hourly atmospheric data into netcdf format (Run/DATA/ConfigName_ERA5_native/ directory) 

2- Convert the raw data into a format useable by the ONLINE_BULK capability of croco ocean model. Run/DATA/ConfigName_ERA5/ directory) 
	- apply correct unit transformation
	- rename the variable with "CROCO" compatible names

3- ( Optional : If not using ONLINE_BULK capabilities) : create from these converted data at the ERA5 resolution, new CROCO "frc" and/or "blk" files (similar to CFSR ones) interpolated onto the croco grid

#####################
--> Pre-requisite : First the user need to install the ERA5 python API : https://cds.climate.copernicus.eu/api-how-to

--> Edit the era5_crocotools_param.py parameter file 

--> Then to download the ERA5 data (step 1 above)
    ./ERA5_request.py

--> Then to convert the ERA5 data with unit and name into a "CROCO online bulk" compatible format (unit and names) (step 2 above):
    ./ERA5_convert.py

--> If the user want to interpolate the data onto the croco model grid in order to crerate frc/blk netcdf file (Optional step 3 above) 

    -->  adapt the crocotools_param.m section :
    %--------------------------------------------------
    %  Options for make_ERA5	
    %--------------------------------------------------

    --> Finally, executing make_ERA5 under matlab should work fine.
