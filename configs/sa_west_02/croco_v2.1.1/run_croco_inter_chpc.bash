#!/bin/bash

# This script is designed to run on the chpc
# it submits a series of run_croco_inter.pbs job scripts, each of which is for a few months 
# this is needed to ensure we stay under the 48 hr job length for each job submission

CURRENT_DATE="2008-01-01"
END_DATE="2013-12-01"
RSTFLAG=0 # set to 1 if you want to restart from a previous months output
prev_jobid=""

# we need to give the absolute path to the script, even though we are running this bash script from the same directory!
SCRIPT="/home/gfearon/lustre/somisana-croco/configs/sa_west_02/croco_v2.1.1/run_croco_inter.pbs"

while [ "$(date -d "$CURRENT_DATE" +%Y%m)" -le "$(date -d "$END_DATE" +%Y%m)" ]; do

    # using jobs of 4 month increments, so add 3 months to get the end dates for this job
    END_DATE_JOB=$(date -d "$CURRENT_DATE +3 month" +%Y-%m-01)
    
    NY_START=$(date -d "$CURRENT_DATE" +%Y)
    NM_START=$(date -d "$CURRENT_DATE" +%-m)
    NY_END=$(date -d "$END_DATE_JOB" +%Y)
    NM_END=$(date -d "$END_DATE_JOB" +%-m)
    
    if [ -z "$prev_jobid" ]; then
        # First job (no dependency)
        jobid=$(qsub -v NY_START=$NY_START,NY_END=$NY_END,NM_START=$NM_START,NM_END=$NM_END,RSTFLAG=$RSTFLAG $SCRIPT)
    else
        # include dependency of previous job
        jobid=$(qsub -W depend=afterok:$prev_jobid -v NY_START=$NY_START,NY_END=$NY_END,NM_START=$NM_START,NM_END=$NM_END,RSTFLAG=$RSTFLAG $SCRIPT)
    fi

    echo "Submitted run starting at month $CURRENT_DATE as job $jobid"
    echo
    prev_jobid=$(echo $jobid | cut -d. -f1)  # Extract numeric job ID    

    # Advance by one month
    CURRENT_DATE=$(date -d "$END_DATE_JOB +1 month" +%Y-%m-01)

    # Always set restart to 1 after the first month
    RSTFLAG=1

done

