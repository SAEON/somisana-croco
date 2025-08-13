#!/bin/bash
set -e

# === CONFIGURATION ===
INPUT_DIR="/media/nc.memela/Leviathan/MODELS/somisana-croco/configs/SA_West/croco_v1.3.1/C01_I01_GLORYS_ERA5/output"
OUTDIR="$(pwd)"
LOGFILE="$OUTDIR/simplified_p90_generation.log"

echo "🔁 Starting monthly stack and 90th percentile computation: $(date)" > "$LOGFILE"

# === PROCESS EACH MONTH ===
for MONTH in $(seq -w 1 12); do
    echo "" | tee -a "$LOGFILE"
    echo "📅 [INFO] Processing month $MONTH..." | tee -a "$LOGFILE"

    STACKED_FILE="$OUTDIR/stacked_M${MONTH}.nc"
    P90_OUT="$OUTDIR/monthly_p90_M${MONTH}.nc"

    # Skip if already done
    if [ -f "$P90_OUT" ]; then
        echo "✅ [INFO] $P90_OUT already exists. Skipping." | tee -a "$LOGFILE"
        continue
    fi

    # Create stacked file if not exists
    if [ ! -f "$STACKED_FILE" ]; then
        echo "📦 [INFO] Stacking files for month $MONTH..." | tee -a "$LOGFILE"
        FILELIST=$(find "$INPUT_DIR" -type f -name "croco_avg_Y????M${MONTH}.nc" | sort)
        if [ -z "$FILELIST" ]; then
            echo "⚠️  [WARN] No files found for month $MONTH in $INPUT_DIR" | tee -a "$LOGFILE"
            continue
        fi
        cdo -O mergetime $FILELIST "$STACKED_FILE" >> "$LOGFILE" 2>&1
    else
        echo "📦 [INFO] Using existing $STACKED_FILE" | tee -a "$LOGFILE"
    fi

    # Compute 90th percentile
    echo "📊 [INFO] Computing 90th percentile for month $MONTH..." | tee -a "$LOGFILE"
    cdo -P 4 timpctl,90 "$STACKED_FILE" "$P90_OUT" >> "$LOGFILE" 2>&1
done

echo "" | tee -a "$LOGFILE"
echo "✅ Finished all months: $(date)" | tee -a "$LOGFILE"

