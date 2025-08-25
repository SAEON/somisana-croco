#!/bin/bash

# we need a way of specifying the suffix of the file
python /home/gfearon/code/somisana-croco/cli.py make_ini_inter \
	--input_dir /home/gfearon/code/somisana-croco/DATASETS_CROCOTOOLS/GLORYS/gulf \
	--output_dir /home/gfearon/code/somisana-croco/configs/test_02/croco_v1.3.1/GLORYS.1 \
	--month_start 2015-01 \
	--Yorig 1993

