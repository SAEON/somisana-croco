#!/bin/bash

# Set the environment variables for an inter-annual run
# just putting in the variables which we may want to make configurable
# we can always add more here if we want

# source code
export WW3_EXE_DIR=/home/gfearon/WW3/model/exe_ow/

# directories corresponding to this configuration
export WW3_SPEC_DIR_NAME=WW3_SPEC_CMEMS # name of the directory in the directory level above this one
export WW3_INPUTDIR_SRF=/mnt/lustre/groups/ERTH1103/data/ERA5/eez_for_croco

# MPI settings
export WW3_MPI_NUM_PROCS=64

