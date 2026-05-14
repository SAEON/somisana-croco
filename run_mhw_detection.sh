#!/bin/bash

# Path to the python CLI script in my repo
CLI_PATH="/home/philile/somisana-croco/cli.py"

# Path to the latest forecast temperature data
TEMP_FILE="/mnt/ocims-somisana/public-facing/sa-west/v1.0/forecasts/latest/MERCATOR-GFS/croco_avg.nc"

# Path to the pre-built 4D climatology file
CLIM_FILE="/mnt/ocims-somisana/public-facing/sa-west/v1.0/hindcasts/GLORYS-ERA5/climatology/day_of_year_climatology.nc"

# Output path for the detected categories
OUTPUT_DIR="/mnt/ocims-somisana/sa-west/v1.0/forecasts/latest/MERCATOR-GFS"
OUTPUT_FILE="$OUTPUT_DIR/forecast_mhw_mcs_categories.nc"

# Reference year for CROCO time (seconds since Yorig-01-01)
YORIG=2000

# Rows processed at once (batch_size)
BATCH_SIZE=5

# --- Execution ---
echo "Starting Marine Heatwave and Marine coldspell detection for forecast data..."

python3 "$CLI_PATH" detect_mhw_forecast \
    --temp_file "$TEMP_FILE" \
    --clim_file "$CLIM_FILE" \
    --fname_out "$OUTPUT_FILE" \
    --temp_var "temp" \
    --Yorig "$YORIG" \
    --batch_size "$BATCH_SIZE"

