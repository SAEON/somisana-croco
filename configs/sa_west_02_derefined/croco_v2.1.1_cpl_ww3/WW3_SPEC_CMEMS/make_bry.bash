#!/bin/bash

# with the wavespectra environment activated!

python /home/$USER/code/somisana-ww3/cli.py make_bry_cmems_monthy \
	--input_dir /home/gfearon/code/somisana-croco/DATASETS_CROCOTOOLS/CMEMS_WAV/eez/ \
	--lon_file ../WW3_GRID/lon.dat \
	--lat_file ../WW3_GRID/lat.dat \
	--mask_file ../WW3_GRID/mask.dat \
	--output_dir $PWD \
	--month_start 2010-01 \
	--month_end 2010-02 

