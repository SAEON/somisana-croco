#!/bin/bash

# --- Configuration ---
# Path to the python CLI script
CLI_PATH="/home/philile/somisana-croco/cli.py"

# Input files
FORECAST_FILE="/mnt/ocims-somisana/public-facing/sa-west/v1.0/forecasts/latest/MERCATOR-GFS/croco_avg.nc"
CAT_FILE="/mnt/ocims-somisana/public-facing/sa-west/v1.0/forecasts/latest/MERCATOR-GFS/forecast_mhw_mcs_categories.nc"
CLIM_FILE="/mnt/ocims-somisana/public-facing/sa-west/v1.0/hindcasts/GLORYS-ERA5/climatology/day_of_year_climatology.nc"

# Output directory for the images and GIFs
OUT_DIR="/mnt/ocims-somisana/sa-west/v1.0/forecasts/latest/MERCATOR-GFS"

# Define forecast dates
START_DATE=$(date -d "-5 days" +%Y-%m-%d)
END_DATE=$(date -d "+5 days" +%Y-%m-%d)

# Reference year for CROCO time
YORIG=2000

# --- Execution ---
echo "Generating operational MHW/MCS plots..."

python3 "$CLI_PATH" plot_mhw_forecast \
    --forecast_file "$FORECAST_FILE" \
    --cat_file "$CAT_FILE" \
    --clim_file "$CLIM_FILE" \
    --out_dir "$OUT_DIR" \
    --start_date "$START_DATE" \
    --end_date "$END_DATE" \
    --Yorig "$YORIG"

