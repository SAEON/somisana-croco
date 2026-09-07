#!/bin/bash

python /home/g.rautenbach/Projects/somisana-croco/cli.py make_tides_inter \
	--input_dir /home/g.rautenbach/Data/TPXO10/ \
	--output_dir $PWD \
	--month_start 2011-02 \
	--month_end 2011-03 \
	--Yorig 1993
