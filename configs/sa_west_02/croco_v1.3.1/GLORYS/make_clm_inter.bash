#!/bin/bash

python /home/gfearon/code/somisana-croco/cli.py make_clm_inter \
	--input_dir /home/gfearon/code/somisana-croco/DATASETS_CROCOTOOLS/GLORYS/swcape \
	--output_dir $PWD \
	--month_start 2013-08 \
	--month_end 2013-08 \
	--Yorig 1993

