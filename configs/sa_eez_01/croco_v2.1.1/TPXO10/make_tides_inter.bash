#!/bin/bash

python /home/gfearon/somisana-croco/cli.py make_tides_inter \
	--input_dir /mnt/lustre/groups/ERTH1103/data/TPXO10/ \
	--output_dir $PWD \
	--month_start 1993-03 \
	--month_end 2025-11 \
	--Yorig 1993
