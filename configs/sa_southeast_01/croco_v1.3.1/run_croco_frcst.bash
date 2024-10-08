#!/bin/bash
#
#unalias cp
#unalias mv
#limit coredumpsize unlimited
CP=/bin/cp
MV=/bin/mv
LN=/bin/ln
MK=/bin/mkdir
########################################################
#  Define files and run parameters
########################################################

# This script is much simpler than run_croco_inter.bash because all the workflow for operational system gets handled elsewhere (see .github/workflows in the root of this repo)
# But you can run it locally if the scratch dir is already set up with all the necessary input files

if [ -z "$1" ]; then
  echo "Error: Missing argument. Please provide directory name where the run will get done e.g. C01_I99_MERCATOR_GFS"
  exit 1
fi

RUNNAME=$1
source myenv_frcst.sh

RUNDIR=`pwd`
SCRATCHDIR=$RUNDIR/$RUNNAME/scratch
OUTPUTDIR=$RUNDIR/$RUNNAME/output

cd $SCRATCHDIR

#  COMPUTE
#
date
mpirun -np $MPI_NUM_PROCS ./croco croco.in > croco.out
date
#
# Test if the run finised properly
echo "Test croco.out"
status=`tail -2 croco.out | grep DONE | wc -l`
if [[ $status == 1 ]]; then
  echo "All good"
  echo "XXXXXXXX"
else
  echo
  echo "Warning: run not finished properly"
  echo
  tail -20 croco.out
  echo
  exit 1
fi
#
#  Archive
#
$MV -f croco_his.nc ${OUTPUTDIR}/croco_his.nc
$MV -f croco_rst.nc ${OUTPUTDIR}/croco_rst.nc
$MV -f croco_avg.nc ${OUTPUTDIR}/croco_avg.nc
$MV -f croco_his_surf.nc ${OUTPUTDIR}/croco_his_surf.nc
$MV -f croco_avg_surf.nc ${OUTPUTDIR}/croco_avg_surf.nc
# $MV -f croco_sta.nc ${OUTPUTDIR}/croco_sta.nc
#
#############################################################
