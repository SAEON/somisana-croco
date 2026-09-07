#!/bin/bash

# This script is designed to run on the chpc
# it submits a series of run_croco_inter.pbs job scripts, each covering MONTHS_PER_JOB
# months, chained with afterok dependencies. Jobs are chunked so each one fits inside
# the walltime requested in run_croco_inter.pbs.

# Months per PBS job.
# A month takes ~50 min on 64 procs; we budget 1 hr/month as a safety factor against
# the 48 hr walltime in run_croco_inter.pbs (48 hr is the CHPC queue maximum), so
# 44 months uses 44 of those 48 hrs: ~4 hrs slack at the budgeted rate, ~11 hrs at
# the measured rate.
# 44 is the sweet spot for this run window. 45, 46 and 47 all need the same 9 jobs but
# give away slack for nothing; 48 would fill the walltime exactly and get killed; 43
# drops to 10 jobs. So 44 is the fewest jobs with the most headroom.
# If you change the walltime, change this to match.
MONTHS_PER_JOB=44

CURRENT_DATE="1993-03-01"
END_DATE="2025-11-01"
RSTFLAG=0 # set to 1 if you want to restart from a previous months output
prev_jobid=""

# we need to give the absolute path to the script, even though we are running this bash script from the same directory!
SCRIPT="/home/gfearon/lustre/somisana-croco/configs/sa_eez_01/croco_v2.1.1/run_croco_inter.pbs"

while [ "$(date -d "$CURRENT_DATE" +%Y%m)" -le "$(date -d "$END_DATE" +%Y%m)" ]; do

    # add MONTHS_PER_JOB-1 months to get the end month covered by this job
    END_DATE_JOB=$(date -d "$CURRENT_DATE +$((MONTHS_PER_JOB - 1)) month" +%Y-%m-01)

    # the run window is not necessarily a whole number of MONTHS_PER_JOB jobs, so clamp
    # the last job to END_DATE - otherwise we'd run past the end of the forcing
    if [ "$(date -d "$END_DATE_JOB" +%Y%m)" -gt "$(date -d "$END_DATE" +%Y%m)" ]; then
        END_DATE_JOB=$(date -d "$END_DATE" +%Y-%m-01)
    fi
    
    if [ -z "$prev_jobid" ]; then
        # First job (no dependency)
        jobid=$(qsub -v START_DATE=$CURRENT_DATE,END_DATE=$END_DATE_JOB,RSTFLAG=$RSTFLAG $SCRIPT)
    else
        # include dependency of previous job
        jobid=$(qsub -W depend=afterok:$prev_jobid -v START_DATE=$CURRENT_DATE,END_DATE=$END_DATE_JOB,RSTFLAG=$RSTFLAG $SCRIPT)
    fi

    echo "Submitted run starting at month $CURRENT_DATE as job $jobid"
    echo
    prev_jobid=$(echo $jobid | cut -d. -f1)  # Extract numeric job ID    

    # Advance by one month
    CURRENT_DATE=$(date -d "$END_DATE_JOB +1 month" +%Y-%m-01)

    # Always set restart to 1 after the first job
    RSTFLAG=1

done

