#!/bin/bash

# Set the environment variables for a forecast run
# just putting in the variables which we may want to make configurable
# A lot of variables in myenv_inter.sh aren't here as they are already defined as part of the operational workflow

# provide the path to the croco source code
export SOURCE=/home/$USER/code/croco-v1.3.1/OCEAN/

# MPI partitioning
export MPI_NUM_X=5
export MPI_NUM_Y=18
export MPI_NUM_PROCS=$(($MPI_NUM_X * $MPI_NUM_Y))
