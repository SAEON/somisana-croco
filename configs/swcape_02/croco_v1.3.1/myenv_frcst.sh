#!/bin/bash

# Set the environment variables for a forecast run
# just putting in the variables which we may want to make configurable
# we can always add more here if we want

# source code
# this is the path inside our docker image used to run the model
export SOURCE=/croco-v1.3.1/OCEAN/
#
# MPI partitioning
export MPI_NUM_X=3
export MPI_NUM_Y=10
export MPI_NUM_PROCS=$(($MPI_NUM_X * $MPI_NUM_Y))
