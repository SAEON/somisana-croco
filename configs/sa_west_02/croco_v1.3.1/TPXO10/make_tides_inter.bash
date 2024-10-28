#!/bin/bash

python /home/gfearon/code/somisana-croco/cli.py make_tides_inter \
	--input_dir /home/gfearon/code/somisana-croco/DATASETS_CROCOTOOLS/TPXO10/ \
	--output_dir /home/gfearon/code/somisana-croco/configs/sa_west_02/croco_v1.3.1/TPXO10 \
	--month_start 2009-01 \
	--month_end 2009-12 \
	--Yorig 1993

