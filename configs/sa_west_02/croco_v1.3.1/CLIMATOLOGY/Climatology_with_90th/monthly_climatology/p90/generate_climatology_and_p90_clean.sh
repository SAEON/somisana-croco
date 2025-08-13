#!/bin/bash

# === CONFIGURATION ===
#INPUT_DIR="/media/nc.memela/Leviathan/MODELS/somisana-croco/configs/SA_West/croco_v1.3.1/C01_I01_GLORYS_ERA5/output"
INPUT_DIR="/mnt/ocims-somisana/sa-west/v1.0/hindcasts/GLORYS-ERA5/raw"
OUTDIR="$(pwd)"
LOGFILE="$OUTDIR/p90_generation.log"
mkdir -p "$OUTDIR"

# Clear log
echo "🌍 Starting p90 climatology generation..." > "$LOGFILE"

# === PROCESS EACH MONTH ===
for MONTH in $(seq -w 1 12); do
    echo "📆 Processing month: $MONTH" | tee -a "$LOGFILE"

    # Gather files for this month across all years (1993–2023)
    FILE_LIST=""
    for YEAR in $(seq 1993 2019); do
        FILE="${INPUT_DIR}/croco_avg_Y${YEAR}M${MONTH}.nc"
        if [ -f "$FILE" ]; then
            FILE_LIST="$FILE_LIST $FILE"
        else
            echo "⚠️ Skipping missing: $FILE" | tee -a "$LOGFILE"
        fi
    done

    # Stack into one monthly timeseries
    STACKED="${OUTDIR}/temp_timeseries_month${MONTH}.nc"
    cdo mergetime $FILE_LIST $STACKED 2>>"$LOGFILE"

    # Compute min and max (for p90 bounds)
    TIMMIN="${OUTDIR}/temp_min_month${MONTH}.nc"
    TIMMAX="${OUTDIR}/temp_max_month${MONTH}.nc"
    cdo timmin $STACKED $TIMMIN 2>>"$LOGFILE"
    cdo timmax $STACKED $TIMMAX 2>>"$LOGFILE"

    # Compute climatology mean and 90th percentile
    CLIM_MEAN="${OUTDIR}/temp_clim_month${MONTH}.nc"
    P90_FILE="${OUTDIR}/temp_p90_month${MONTH}.nc"
    cdo timmean $STACKED $CLIM_MEAN 2>>"$LOGFILE"
    cdo timpctl,90 $STACKED $TIMMIN $TIMMAX $P90_FILE 2>>"$LOGFILE"

    echo "✅ Done: $CLIM_MEAN and $P90_FILE" | tee -a "$LOGFILE"

    # Delete large intermediate file
    echo "🧹 Removing $STACKED to save space..." | tee -a "$LOGFILE"
    rm -f "$STACKED"
done

echo "🎉 All done. Output files saved in $OUTDIR" | tee -a "$LOGFILE"

