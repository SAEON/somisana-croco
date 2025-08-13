#!/bin/bash

# === Setup ===
INPUT_DIR="/media/nc.memela/Leviathan/MODELS/somisana-croco/configs/SA_West/croco_v1.3.1/C01_I01_GLORYS_ERA5/output/"
OUT_DIR="/home/nc.memela/Projects/somisana-croco/configs/sa_west_02/croco_v1.3.1/CLIMATOLOGY/Climatology_with_90th"
mkdir -p "$OUT_DIR"

LOGFILE="${OUT_DIR}/climatology_generation.log"
echo "Climatology generation started at $(date)" > "$LOGFILE"

export OMP_NUM_THREADS=4
VAR="temp"

# === Process each month ===
for MONTH in $(seq -w 1 12); do
    echo "[INFO] Processing month ${MONTH}..." | tee -a "$LOGFILE"
    
    # === Mean climatology for temp ===
    MEAN_OUT="${OUT_DIR}/monthly_mean_M${MONTH}.nc"
    if [[ ! -f "$MEAN_OUT" ]]; then
        FILES=$(ls ${INPUT_DIR}/croco_avg_Y*M${MONTH}.nc 2>/dev/null)
        if [[ -z "$FILES" ]]; then
            echo "[WARN] No files found for month ${MONTH}. Skipping mean." | tee -a "$LOGFILE"
        else
            echo "[INFO] Computing monthly mean for ${MONTH}..." | tee -a "$LOGFILE"
            cdo -P 4 -selname,${VAR} -ensmean $FILES "$MEAN_OUT" >> "$LOGFILE" 2>&1
        fi
    else
        echo "[SKIP] Mean output already exists: $MEAN_OUT" | tee -a "$LOGFILE"
    fi

    # === 90th percentile climatology for temp ===
    STACKED="${OUT_DIR}/stacked_M${MONTH}.nc"
    P90_OUT="${OUT_DIR}/monthly_p90_M${MONTH}.nc"
    
    if [[ ! -f "$P90_OUT" ]]; then
        echo "[INFO] Stacking time series for percentile computation..." | tee -a "$LOGFILE"
        cdo -P 4 -selname,${VAR} -mergetime ${INPUT_DIR}/croco_avg_Y*M${MONTH}.nc "$STACKED" >> "$LOGFILE" 2>&1

        echo "[INFO] Computing 90th percentile for ${MONTH}..." | tee -a "$LOGFILE"
        cdo -P 4 -timselpercentile,30,90 "$STACKED" "$P90_OUT" >> "$LOGFILE" 2>&1

        rm -f "$STACKED"
    else
        echo "[SKIP] 90th percentile output already exists: $P90_OUT" | tee -a "$LOGFILE"
    fi
done

# === Combine monthly outputs ===
echo "[INFO] Concatenating all mean climatologies..." | tee -a "$LOGFILE"
cdo -O -P 4 mergetime ${OUT_DIR}/monthly_mean_M*.nc ${OUT_DIR}/monthly_mean_clim.nc

echo "[INFO] Concatenating all 90th percentile climatologies..." | tee -a "$LOGFILE"
cdo -O -P 4 mergetime ${OUT_DIR}/monthly_p90_M*.nc ${OUT_DIR}/monthly_p90_clim.nc

echo "✅ All climatologies generated and saved in $OUT_DIR" | tee -a "$LOGFILE"

