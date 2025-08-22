#!/bin/bash

# Set the environment variables for an inter-annual run
# just putting in the variables which we may want to make configurable
# we can always add more here if we want

# source code
export SOURCE=/home/gfearon/croco-v2.0.1/OCEAN/

# dir with compile options
export EXENAME=C09

# MPI partitioning
export MPI_NUM_X=30
export MPI_NUM_Y=4
export MPI_NUM_PROCS=$(($MPI_NUM_X * $MPI_NUM_Y))

# boundary and atmospheric forcing
export BULK_FILES=0
export ATMOS_BULK=ERA5
#
export FORCING_FILES=0
export ATMOS_FRC=NA
#
export CLIMATOLOGY_FILES=0
export BOUNDARY_FILES=1
export OGCM=GLORYS
#
export TIDE_FILES=1
export TIDE_FRC=TPXO10

# runtime input
export INNAME=I01
