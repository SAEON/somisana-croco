#!/bin/bash

# Set the environment variables for an inter-annual run
# just putting in the variables which we may want to make configurable
# we can always add more here if we want

# source code
export SOURCE=/home/$USER/code/croco-v1.3.1/OCEAN/
#
# dir with compile options
export EXENAME=C03
#
# MPI partitioning
export MPI_NUM_X=5
export MPI_NUM_Y=3
export MPI_NUM_PROCS=$(($MPI_NUM_X * $MPI_NUM_Y))
#
# boundary and atmospheric forcing
export ATMOS_BULK=ERA5
export ATMOS_FRC=NA
export OGCM=GLORYS
#
# tidal forcing
export TIDE_FRC=TIDETPXO
#
# runtime input
export INNAME=I02
