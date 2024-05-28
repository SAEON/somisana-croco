#!/bin/bash

# environment variables related to the runtime input file(s)
# these would typically be in the croco run script, but we're putting
# them here as they relate directly to what gets written to 
# the runtime input (i.e. *.in) file
# 
#
# Model time step [seconds]
#
DT=60
DTFAST=60
#
# Output frequency (in hours)
#
# full domain output
NH_AVG=1
NH_HIS=24
# surface output
NH_AVGSURF=1
NH_HISSURF=24
# station output
NH_STA=1
# restart file
# (writing at 6 hourly intervals gives us flexibility for restarting from a previous run which might have been run say 18 hours ago)
NH_RST=6
#
# Time refinement coefficient (factor to apply to time-step at each child level)
#
T_REF=3
#
# Number total of grid levels (1: No child grid)
#
NLEVEL=1
#
# We don't define if the run is from a restart file or not here - this is handled in the operational workflow
# 
