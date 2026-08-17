#!/bin/bash

# This script is designed to run on the chpc
# it submits a series of run_inter.pbs job scripts, each of which is for a few months 
# this is needed to ensure we stay under the 48 hr job length for each job submission

CURRENT_DATE="2008-01-01"
END_DATE="2013-12-31"
RSTFLAG=0 # set to 1 if you want to restart from a previous months output
prev_jobid=""

# which run directory to integrate (RUN_01 = full wave coupling,
# RUN_02 = conservative wave effects only - see RUN_02/README.md)
RUN_NAME=${RUN_NAME:-RUN_02}

CONFIG_DIR="/home/gfearon/lustre/somisana-croco/configs/sa_west_02/croco_v2.1.1_cpl_ww3"

# we need to give the absolute path to the script, even though we are running this bash script from the same directory!
SCRIPT="$CONFIG_DIR/run_inter.pbs"

# run_inter.pbs has static #PBS -o/-e directives that cannot interpolate
# $RUN_NAME, so override them here to keep each run's logs separate
PBS_LOGS="-o $CONFIG_DIR/$RUN_NAME/stdout -e $CONFIG_DIR/$RUN_NAME/stderr -N $RUN_NAME"

while [ "$(date -d "$CURRENT_DATE" +%Y%m)" -le "$(date -d "$END_DATE" +%Y%m)" ]; do

    # using jobs of 4 month increments, so add 3 months to get the end dates for this job
    END_DATE_JOB=$(date -d "$CURRENT_DATE +3 month" +%Y-%m-01)
    
    NY_START=$(date -d "$CURRENT_DATE" +%Y)
    NM_START=$(date -d "$CURRENT_DATE" +%-m)
    NY_END=$(date -d "$END_DATE_JOB" +%Y)
    NM_END=$(date -d "$END_DATE_JOB" +%-m)

    if [ -z "$prev_jobid" ]; then
        # First job (no dependency)
        jobid=$(qsub $PBS_LOGS -v START_DATE=$CURRENT_DATE,END_DATE=$END_DATE_JOB,RSTFLAG=$RSTFLAG,RUN_NAME=$RUN_NAME $SCRIPT)
    else
        # include dependency of previous job
        jobid=$(qsub $PBS_LOGS -W depend=afterok:$prev_jobid -v START_DATE=$CURRENT_DATE,END_DATE=$END_DATE_JOB,RSTFLAG=$RSTFLAG,RUN_NAME=$RUN_NAME $SCRIPT)
    fi

    echo "Submitted $RUN_NAME starting at month $CURRENT_DATE as job $jobid"
    echo
    prev_jobid=$(echo $jobid | cut -d. -f1)  # Extract numeric job ID    

    # Advance by one month
    CURRENT_DATE=$(date -d "$END_DATE_JOB +1 month" +%Y-%m-01)

    # Always set restart to 1 after the first job
    RSTFLAG=1

done

