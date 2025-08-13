#!/bin/bash
set -e

# === CONFIGURATION ===
INPUT_DIR="/media/nc.memela/Leviathan/MODELS/somisana-croco/configs/SA_West/croco_v1.3.1/C01_I01_GLORYS_ERA5/output"
OUTDIR="$(pwd)"
LOGFILE="$OUTDIR/p90_generation.log"

echo "🔁 Starting monthly stack and 90th percentile computation: $(date)" > "$LOGFILE"

for MONTH in $(seq -w 1 12); do
    echo "" | tee -a "$LOGFILE"
    echo "📅 [INFO] Processing month $MONTH..." | tee -a "$LOGFILE"

    STACKED_FILE="$OUTDIR/stacked_M${MONTH}.nc"
    P90_FILE="$OUTDIR/monthly_p90_M${MONTH}.nc"

    # Skip if already done
    if [ -f "$P90_FILE" ]; then
        echo "✅ [SKIP] $P90_FILE already exists." | tee -a "$LOGFILE"
        continue
    fi

    # Stack CROCO monthly files (e.g., croco_avg_Y????M01.nc)
    echo "📦 [INFO] Stacking files for month $MONTH..." | tee -a "$LOGFILE"
    MONTHLY_FILES=($(ls "$INPUT_DIR"/croco_avg_Y????M${MONTH}.nc 2>/dev/null))

    if [ ${#MONTHLY_FILES[@]} -eq 0 ]; then
        echo "⚠️  [WARN] No input files found for month $MONTH. Skipping." | tee -a "$LOGFILE"
        continue
    fi

    cdo mergetime "${MONTHLY_FILES[@]}" "$STACKED_FILE" >> "$LOGFILE" 2>&1

    # Compute 90th percentile
    echo "📊 [INFO] Computing 90th percentile for month $MONTH..." | tee -a "$LOGFILE"
    cdo -P 4 timpctl,90 "$STACKED_FILE" "$P90_OUT" >> "$LOGFILE" 2>&1


    echo "✅ [DONE] Month $MONTH complete." | tee -a "$LOGFILE"
done

echo "" | tee -a "$LOGFILE"
echo "🎉 All months processed. See $LOGFILE for details." | tee -a "$LOGFILE"

