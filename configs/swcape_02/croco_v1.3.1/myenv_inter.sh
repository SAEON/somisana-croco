#!/bin/bash

# Set the environment variables for an inter-annual run
# just putting in the variables which we may want to make configurable
# we can always add more here if we want

# source code
export SOURCE=/home/gfearon/code/croco-v1.3.1/OCEAN/
#
# dir with compile options
export EXENAME=C01
#
# MPI partitioning
export MPI_NUM_X=2
export MPI_NUM_Y=8
export MPI_NUM_PROCS=$(($MPI_NUM_X * $MPI_NUM_Y))
#
# boundary and atmospheric forcing
export ATMOS_BULK=ERA5
export ATMOS_FRC=NA
export OGCM=GLORYS
#
# runtime input
export INNAME=I01
