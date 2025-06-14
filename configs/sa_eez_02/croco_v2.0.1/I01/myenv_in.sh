#!/bin/bash

# environment variables related to the runtime input file(s)
# these would typically be in the croco run script, but we're putting
# them here as they relate directly to what gets written to 
# the runtime input (i.e. *.in) file
# 
#
# Model time step [seconds]
#
DT=180
DTFAST=60
#      Number of hours for averages files
#      24 is 1 days
#      48 is 2 days
#      72 is 3 days
#      96 is 4 days
#      120 is 5 days
NH_AVG=24 #6
NH_HIS=24 #1
NH_AVGSURF=1 #6
NH_HISSURF=24 #1
NH_STA=1
#
# Time refinement coefficient (factor to apply to time-step at each child level)
#
T_REF=3
#
# Number total of grid levels (1: No child grid)
#
NLEVEL=1
#
# these are the time extents for this run
#NY_START=1993
#NY_END=2019
#NM_START=1
#NM_END=12
NY_START=2009
NY_END=2013
NM_START=1
NM_END=12
#
# Set month format at 1 or 2 digits (for input and output files): "%01d" = 1 digit/ "%02d" = 2 digit
MTH_FORMAT="%02d"
#
# Number of year that are considered to be part of the spin-up (i.e. 365 days per year)
NY_SPIN=0
#
#  Restart file - RSTFLAG=0 --> No Restart
#		  RSTFLAG=1 --> Restart
#
RSTFLAG=0
#
#  Time Schedule  -  TIME_SCHED=0 --> yearly files
#                    TIME_SCHED=1 --> monthly files
#
TIME_SCHED=1
#
