#!/bin/bash

python /home/gfearon/somisana-croco/cli.py make_clm_inter \
	--input_dir /mnt/lustre/groups/ERTH1103/data/GLORYS/eez \
	--output_dir $PWD \
	--month_start 2008-01 \
	--month_end 2013-12 \
	--Yorig 1993

