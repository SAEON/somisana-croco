#!/bin/bash

# === Output filenames ===
CLIM_OUT="climatology_12months.nc"
P90_OUT="p90_12months.nc"

# === Input pattern ===
CLIM_PATTERN="temp_clim_month*.nc"
P90_PATTERN="temp_p90_month*.nc"

echo "📦 Concatenating climatology files into $CLIM_OUT..."
cdo mergetime $CLIM_PATTERN $CLIM_OUT

echo "📦 Concatenating p90 files into $P90_OUT..."
cdo mergetime $P90_PATTERN $P90_OUT

echo "✅ Done. Output files: $CLIM_OUT and $P90_OUT"

