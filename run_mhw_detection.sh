#!/bin/bash

# --- Configuration ---
# Path to the python CLI script in your repo
CLI_PATH="/home/philile/somisana-croco/cli.py"

# Path to the latest forecast temperature data
TEMP_FILE="/mnt/ocims-somisana/public-facing/sa-west/v1.0/forecasts/latest/MERCATOR-GFS/croco_avg.nc"

# Path to the pre-built 4D climatology file
CLIM_FILE="/mnt/ocims-somisana/public-facing/sa-west/v1.0/hindcasts/GLORYS-ERA5/products/mhw/Climatology_4D_Unified.nc"

# Output path for the detected categories
# Creating a dedicated output folder in your home directory is recommended
OUTPUT_DIR="/home/philile/mhw_outputs"
mkdir -p "$OUTPUT_DIR"
OUTPUT_FILE="$OUTPUT_DIR/forecast_mhw_categories.nc"

# Reference year for CROCO time (seconds since Yorig-01-01)
YORIG=2000

# Rows processed at once (batch_size)
# Using 5 as per the CLI default
BATCH_SIZE=5

# --- Execution ---
echo "Starting Marine Heatwave detection for forecast data..."

python3 "$CLI_PATH" detect_mhw_forecast \
    --temp_file "$TEMP_FILE" \
    --clim_file "$CLIM_FILE" \
    --fname_out "$OUTPUT_FILE" \
    --temp_var "temp" \
    --Yorig "$YORIG" \
    --batch_size "$BATCH_SIZE"

if [ $? -eq 0 ]; then
    echo "Detection completed successfully."
    echo "Output: $OUTPUT_FILE"
else
    echo "Error: MHW detection failed."
    exit 1
fi