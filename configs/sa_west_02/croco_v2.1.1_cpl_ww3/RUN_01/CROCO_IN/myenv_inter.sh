#!/bin/bash

# Set the environment variables for an inter-annual run
# just putting in the variables which we may want to make configurable
# we can always add more here if we want

# source code
export CROCO_SRC=/home/gfearon/croco-v2.1.1/OCEAN/

# MPI partitioning
export CROCO_MPI_NUM_X=4
export CROCO_MPI_NUM_Y=14
export CROCO_MPI_NUM_PROCS=$(($CROCO_MPI_NUM_X * $CROCO_MPI_NUM_Y))

# boundary and atmospheric forcing
export BULK_FILES=0
export ATMOS_BULK=ERA5
#
export FORCING_FILES=0
export ATMOS_FRC=NA
#
export CLIMATOLOGY_FILES=1
export BOUNDARY_FILES=0
export OGCM=GLORYS
#
export TIDE_FILES=1
export TIDE_FRC=TPXO10

# Model time step [seconds]
#
CROCO_DT0=60
CROCO_NFAST=60
# Number of hours for output files
NH_AVG=24
NH_HIS=24 # not used in the existing config - see croco_inter.in
NH_AVGSURF=24 # not used in the existing config - see croco_inter.in 
NH_HISSURF=1
NH_STA=1 # not used in the existing config - undef STATION in cppdefs.h
# Convert to seconds
NS_AVG=$((NH_AVG*3600))
NS_HIS=$((NH_HIS*3600))
NS_AVGSURF=$((NH_AVGSURF*3600))
NS_HISSURF=$((NH_HISSURF*3600))
NS_STA=$((NH_STA*3600)) 
#
# Number total of grid levels (1: No child grid)
#
NLEVEL=1
