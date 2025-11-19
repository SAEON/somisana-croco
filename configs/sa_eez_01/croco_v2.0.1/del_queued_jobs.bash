#!/bin/bash
for job in $(qstat -u gfearon | awk 'NR>5 {print $1}'); do
    qdel "$job"
done
