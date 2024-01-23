GGF comment - Here, I'm only doing option 1 below i.e. just the download of raw files
I do the processessing steps in the configuration directory 


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
