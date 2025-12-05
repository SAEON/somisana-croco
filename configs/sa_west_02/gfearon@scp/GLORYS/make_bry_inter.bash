#!/bin/bash

python /home/gfearon/code/somisana-croco/cli.py make_bry_inter \
	--input_dir /home/gfearon/code/somisana-croco/DATASETS_CROCOTOOLS/GLORYS/swcape \
	--output_dir $PWD \
	--month_start 1993-02 \
	--month_end 2019-12 \
	--Yorig 1993

