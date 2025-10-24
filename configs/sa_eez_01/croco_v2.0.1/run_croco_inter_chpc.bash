#!/bin/bash

CURRENT_DATE="2008-01-01"
END_DATE="2013-12-31"
RSTFLAG=0 # set to 1 if you want to restart from a previous months output
prev_jobid=""
   
# I want to submit command line inputs to run_croco_inter.pbs for NY, NM and RSTFLAG
# but when I do this the #PBS directives are not read!
# So I am defining them here as inputs to qsub
# This is a bit hacky, and could cause confusion when editing the #PBS directives in run_croco_inter.pbs and not seeing any change in behaviour
# but this is all I could come up with at the time 
QSUB_OPTIONS="-l select=4:ncpus=16:mpiprocs=16 -P ERTH1103 -q normal -l walltime=24:00:00"
QSUB_OPTIONS+=" -o /home/gfearon/run_croco_inter_lengau/stdout -e /home/gfearon/run_croco_inter_lengau/stderr"

while [ "$(date -d "$CURRENT_DATE" +%Y%m)" -le "$(date -d "$END_DATE" +%Y%m)" ]; do

    # using jobs of 4 month increments, so add 3 months to get the end dates for this job
    END_DATE_JOB=$(date -d "$CURRENT_DATE +3 month" +%Y-%m-01)
    
    NY_START=$(date -d "$CURRENT_DATE" +%Y)
    NM_START=$(date -d "$CURRENT_DATE" +%-m)
    NY_END=$(date -d "$END_DATE_JOB" +%Y)
    NM_END=$(date -d "$END_DATE_JOB" +%-m)
    
    if [ -z "$prev_jobid" ]; then
        # First job (no dependency)
        jobid=$(qsub $QSUB_OPTIONS -- run_croco_inter.pbs $NY_START $NY_END $NM_START $NM_END $RSTFLAG)
    else
        # include dependency of previous job
        jobid=$(qsub -W depend=afterok:$prev_jobid $QSUB_OPTIONS -- run_croco_inter.pbs $NY_START $NY_END $NM_START $NM_END $RSTFLAG)
    fi

    echo "Submitted run starting at month $CURRENT_DATE as job $jobid"
    echo
    prev_jobid=$(echo $jobid | cut -d. -f1)  # Extract numeric job ID    

    # Advance by one month
    CURRENT_DATE=$(date -d "$END_DATE_JOB +1 month" +%Y-%m-01)

    # Always set restart to 1 after the first month
    RSTFLAG=1

done

