This readme and scripts are based on the official croco_tools-v1.3.1

The scripts are adapted to be unrelated to a particular croco configuration i.e. the file structure is not prescribed

Also, the protocol for downloading ERA5 was changed in Sep 2024 and these scripts take the changes into account (this was not the case in the official croco_tools at the time of writing)

Here is what you would need to do to download ERA5 data

--> Pre-requisite : First the user needs to install the ERA5 python API : https://cds.climate.copernicus.eu/how-to-api

--> Edit the era5_crocotools_param.py parameter file 

--> Then to download the ERA5 data (step 1 above)
    ./ERA5_request.py

--> Then to convert the ERA5 data with unit and name into a "CROCO online bulk" compatible format (unit and names):
    ./ERA5_convert.py

