# #!/usr/bin/env python3
# import xarray as xr
# import matplotlib.pyplot as plt
# import numpy as np
# import os

# # === CONFIGURATION ===
# data_dir = os.getcwd()
# lon_point = 18.0    # ⚠️ Your desired longitude
# lat_point = -34.0   # ⚠️ Your desired latitude
# temp_var = "temp"   # Name of temperature variable
# depth_dim = "s_rho" # Depth dimension name in your CROCO files

# # === STORAGE ===
# months = []
# clim_values = []
# p90_values = []

# # === LOOP THROUGH MONTHS ===
# for m in range(1, 13):
#     month_str = f"{m:02d}"
#     clim_file = f"{data_dir}/temp_clim_month{month_str}.nc"
#     p90_file = f"{data_dir}/temp_p90_month{month_str}.nc"

#     print(f"📦 Loading Month {month_str}: {os.path.basename(p90_file)}")

#     try:
#         # Open datasets
#         ds_clim = xr.open_dataset(clim_file)
#         ds_p90 = xr.open_dataset(p90_file)

#         # Promote lon/lat to coordinates for both datasets
#         if 'lon_rho' in ds_clim.variables and 'lat_rho' in ds_clim.variables:
#             ds_clim = ds_clim.set_coords(['lon_rho', 'lat_rho'])
#         if 'lon_rho' in ds_p90.variables and 'lat_rho' in ds_p90.variables:
#             ds_p90 = ds_p90.set_coords(['lon_rho', 'lat_rho'])


#         # Extract values at the surface (last index in s_rho)
#         val_clim = ds_clim[temp_var].sel(lon_rho=lon_point, lat_rho=lat_point, method='nearest').isel({depth_dim: -1}).values.item()
#         val_p90 = ds_p90[temp_var].sel(lon_rho=lon_point, lat_rho=lat_point, method='nearest').isel({depth_dim: -1}).values.item()

#         months.append(m)
#         clim_values.append(val_clim)
#         p90_values.append(val_p90)

#     except Exception as e:
#         print(f"❌ Skipping month {month_str}: {e}")

# # === PLOT ===
# plt.figure(figsize=(10, 5))
# plt.plot(months, clim_values, marker='o', label='Climatology (mean)')
# plt.plot(months, p90_values, marker='s', label='90th Percentile')
# plt.title(f'Surface Temperature Climatology vs P90\n({lat_point}°, {lon_point}°)')
# plt.xlabel('Month')
# plt.ylabel('Temperature (°C)')
# plt.xticks(months)
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
# plt.savefig("climatology_vs_p90_surface_plot.png", dpi=300)
# plt.show()

#!/usr/bin/env python3
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import os

# === User configuration ===
lon_idx = 100   # Index in xi_rho direction
lat_idx = 300   # Index in eta_rho direction
depth_dim = "s_rho"
temp_var = "temp"

# === Setup ===
clim_template = "temp_clim_month{:02d}.nc"
p90_template = "temp_p90_month{:02d}.nc"
months = range(1, 13)

clim_vals = []
p90_vals = []

print("📈 Extracting temperature at grid index (eta, xi):", lat_idx, lon_idx)

for month in months:
    fname_clim = clim_template.format(month)
    fname_p90 = p90_template.format(month)

    if not (os.path.exists(fname_clim) and os.path.exists(fname_p90)):
        print(f"⚠️  Skipping month {month:02d}: file(s) missing")
        continue

    print(f"📦 Processing Month {month:02d}")

    try:
        ds_clim = xr.open_dataset(fname_clim)
        ds_p90 = xr.open_dataset(fname_p90)

        val_clim = ds_clim[temp_var].isel(
            xi_rho=lon_idx, eta_rho=lat_idx, **{depth_dim: -1}
        ).values.item()

        val_p90 = ds_p90[temp_var].isel(
            xi_rho=lon_idx, eta_rho=lat_idx, **{depth_dim: -1}
        ).values.item()

        clim_vals.append(val_clim)
        p90_vals.append(val_p90)

    except Exception as e:
        print(f"❌ Failed to extract month {month:02d}: {e}")
        continue

# === Plot ===
if clim_vals and p90_vals:
    months_label = np.arange(1, len(clim_vals)+1)

    plt.figure(figsize=(10, 5))
    plt.plot(months_label, clim_vals, label="Climatology", marker='o')
    plt.plot(months_label, p90_vals, label="90th Percentile", marker='s')
    plt.title("Surface Temperature Comparison (Monthly)")
    plt.xlabel("Month")
    plt.ylabel("Temperature (°C)")
    plt.grid(True)
    plt.xticks(ticks=months_label)
    plt.legend()
    plt.tight_layout()
    plt.savefig("clim_vs_p90_surface_comparison.png", dpi=150)
    plt.show()
else:
    print("⚠️ No valid data to plot.")

