#!/bin/bash

# environment variables related to the runtime input file(s)
# these would typically be in the croco run script, but we're putting
# them here as they relate directly to what gets written to 
# the runtime input file(s)
# 
#
# Model time step [seconds]
#
DT=60
DTFAST=60
#      Number of hours for averages files
#      24 is 1 days
#      48 is 2 days
#      72 is 3 days
#      96 is 4 days
#      120 is 5 days
NH_AVG=1 #6
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
NY_START=2011
NY_END=2011
NM_START=2
NM_END=4
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
RSTFLAG=1
#
#  Time Schedule  -  TIME_SCHED=0 --> yearly files
#                    TIME_SCHED=1 --> monthly files
#
TIME_SCHED=1
#
