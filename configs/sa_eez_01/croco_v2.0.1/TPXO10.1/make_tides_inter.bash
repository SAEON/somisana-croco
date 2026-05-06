#!/bin/bash

python /home/gfearon/code/somisana-croco/cli.py make_tides_inter \
	--input_dir /home/gfearon/code/somisana-croco/DATASETS_CROCOTOOLS/TPXO10/ \
	--output_dir $PWD \
	--month_start 2008-01 \
	--month_end 2013-12 \
	--Yorig 1993

