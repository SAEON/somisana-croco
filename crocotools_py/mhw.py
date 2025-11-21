
def run_mhw_detection(
    fname,
    var_str="temp",
    mode="grid",             # "grid", "area_mean", "nearest"
    time=slice(None),        # same idea as get_var(): slice of datetimes or slice(None)
    level=-1,                # s_rho index: -1=surface, 0=bottom; None=>full water column (grid only)
    Yorig=1993,              # only used as a fallback if units are unhelpful
    pctile=90,               # MHW percentile
    clim_start_year=None,    # climatology baseline years; if None -> full span
    clim_end_year=None,
    target_lon=None,         # used when mode="nearest"
    target_lat=None,
    grdname=None,            # grid file for nearest-point search; if None uses fname
    Bottom=None,             # passed to find_nearest_point if you want depth-aware nearest
    write_climatology=False, # write seas/thresh-only file
    nc_out=None,             # main output filename
    clim_out=None,           # optional seas/thresh filename
    engine="netcdf4",
):
    """
    Lean marine heatwave detection for CROCO output.

    This version writes ONLY daily diagnostics:

      - seas(time)
      - thresh(time)
      - anom_relSeas(time)
      - anom_relThresh(time)
      - is_mhw(time)       (0/1)
      - category_code(time) (0–4)

    It deliberately does NOT:
      - save raw temperature in the output file
      - save per-event tables from mhw_out

    Parameters
    ----------
    fname : str, list, or xr.Dataset
        CROCO output file name / glob pattern / list of files, or open Dataset.
    var_str : str
        Variable name to analyse (e.g. "temp").
    mode : {"grid", "area_mean", "nearest"}
        "grid"      – full MHW fields on the CROCO grid
        "area_mean" – spatial mean time series
        "nearest"   – single time series at nearest rho point to (target_lon, target_lat)
    time : slice or datetime-like
        Time selection (e.g. slice("1993-01-01","1994-12-31")) or slice(None).
    level : int or None
        s_rho index. -1=surface, 0=bottom.
        If None and mode="grid", detection is done on the full water column
        (per s_rho, eta_rho, xi_rho). For 1-D modes, level must be an integer.
    Yorig : int, optional
        Only used as a fallback if time units do not contain "since".
    pctile : float
        Percentile used for threshold (usually 90).
    clim_start_year, clim_end_year : int or None
        Climatology baseline years. If None, the full span in the data is used.
    target_lon, target_lat : float
        Required if mode="nearest".
    grdname : str or xr.Dataset, optional
        Grid file for find_nearest_point(); if None, uses fname.
    Bottom : float, optional
        Passed through to find_nearest_point() if you want a minimum depth.
    write_climatology : bool
        If True, also write a small seas+thresh-only NetCDF.
    nc_out : str, optional
        Main output filename. If None, a name is constructed.
    clim_out : str, optional
        Climatology-only filename. If None and write_climatology=True, a name is constructed.
    engine : str
        NetCDF backend used in to_netcdf().
    """

    def _summarize_mhw_output(da, mhw_out, clim, clim_period_label="(full span)"):
        """Console-only summary, no effect on output content."""
        t0 = pd.to_datetime(da.time.values[0]).date()
        t1 = pd.to_datetime(da.time.values[-1]).date()
        years_local = pd.to_datetime(da.time.values).year
        N = int(mhw_out["n_events"])

        print("\n===== MHW detection summary =====")
        print(f"Time span: {t0} → {t1}  (years {years_local.min()}–{years_local.max()})")
        print(f"Climatology period used: {clim_period_label}")
        print(f"Events detected: {N}")

        if N == 0:
            return

        starts = np.array(mhw_out["date_start"])
        ends   = np.array(mhw_out["date_end"])
        peaks  = np.array(mhw_out["date_peak"])
        dur    = np.array(mhw_out["duration"], dtype=int)
        Imax   = np.array(mhw_out["intensity_max"], dtype=float)
        Imean  = np.array(mhw_out["intensity_mean"], dtype=float)
        cat    = np.array(mhw_out["category"])

        k = min(5, N)
        idx_top = np.argsort(-Imax)[:k]
        print("\nTop events by peak intensity:")
        for rank, ii in enumerate(idx_top, 1):
            print(
                f" {rank}. {starts[ii]} → {ends[ii]} (peak {peaks[ii]}); "
                f"dur={dur[ii]} d; I_max={Imax[ii]:.2f} °C; "
                f"I_mean={Imean[ii]:.2f} °C; cat={cat[ii]}"
            )

        imax_dur = int(np.argmax(dur))
        print("\nLongest event:")
        print(
            f" - {starts[imax_dur]} → {ends[imax_dur]} "
            f"(dur={dur[imax_dur]} d; I_max={Imax[imax_dur]:.2f} °C; "
            f"cat={cat[imax_dur]})"
        )

        delta_local = np.nanmean(np.asarray(clim["thresh"]) - np.asarray(clim["seas"]))
        print(f"\nMean (threshold - seasonal) = {delta_local:.3f} °C (should be > 0)")

    # -------------------------------
    # open data (no handle_time here)
    # -------------------------------
    print("")
    print(f"  Running run_mhw_detection() for {var_str} in mode='{mode}'")

    if isinstance(fname, xr.Dataset):
        ds = fname.copy()

    elif isinstance(fname, (list, tuple)):
        ds = xr.open_mfdataset(
            fname,
            compat="override",
            decode_times=False,
            data_vars="minimal",
            coords="minimal",
            engine=engine,
        )

    elif isinstance(fname, str):
        # string can be a single file or a glob pattern
        ds = xr.open_mfdataset(
            fname,
            compat="override",
            decode_times=False,
            data_vars="minimal",
            coords="minimal",
            engine=engine,
        )
    else:
        raise TypeError("fname must be a Dataset, list of files, or file pattern string.")

    if var_str not in ds:
        ds.close()
        raise KeyError(f"Variable '{var_str}' not found in dataset.")

    # -------------------------------
    # decode time from CF units
    # -------------------------------
    time_da = ds["time"]
    time_vals = time_da.values.astype("float64")
    units = time_da.attrs.get("units", "")

    if "since" in units:
        unit_part, _, ref_str = units.partition(" since ")
        unit_part = unit_part.strip().lower()
        ref_ts = pd.to_datetime(ref_str.strip())

        if "second" in unit_part:
            t_pd_full = ref_ts + pd.to_timedelta(time_vals, unit="s")
        elif "day" in unit_part:
            t_pd_full = ref_ts + pd.to_timedelta(time_vals, unit="D")
        else:
            ds.close()
            raise ValueError(f"Unsupported time unit in {units!r}")
    else:
        # fallback: assume "days since Yorig-01-01"
        if Yorig is None:
            ds.close()
            raise ValueError(f"Time units missing 'since' and Yorig is None: {units!r}")
        ref_ts = pd.Timestamp(f"{int(Yorig)}-01-01T00:00:00")
        t_pd_full = ref_ts + pd.to_timedelta(time_vals, unit="D")

    # apply decoded datetimes as coordinate
    ds = ds.assign_coords(time=("time", t_pd_full.to_numpy(dtype="datetime64[ns]")))

    # optional time selection from caller
    if not (isinstance(time, slice) and time == slice(None)):
        ds = ds.sel(time=time)

    da_full = ds[var_str]

    if "s_rho" not in da_full.dims:
        ds.close()
        raise ValueError("run_mhw_detection currently expects var_str with s_rho dimension.")

    # -------------------------------
    # decide depth selection label
    # -------------------------------
    if level is None:
        depth_label = "full_water_column"
    else:
        level = int(level)
        if level == -1:
            depth_label = "surface"
        elif level == 0:
            depth_label = "bottom"
        else:
            depth_label = f"s_rho{level}"

    # -------------------------------
    # pick level / full column
    # -------------------------------
    if level is None:
        if mode != "grid":
            ds.close()
            raise ValueError("level=None (full water column) only supported for mode='grid'.")
        surf = da_full
    else:
        surf = da_full.isel(s_rho=level)  # (time, eta_rho, xi_rho)

    # daily resampling
    surf = surf.resample(time="1D").mean()
    surf = surf.chunk({"time": -1})

    # years from resampled time axis
    t_pd = pd.to_datetime(surf["time"].values)
    years = t_pd.year
    print("  Time span:", years.min(), "→", years.max())

    # -------------------------------
    # climatology years (always valid)
    # -------------------------------
    if clim_start_year is None:
        clim_start_use = int(years.min())
    else:
        clim_start_use = int(clim_start_year)

    if clim_end_year is None:
        clim_end_use = int(years.max())
    else:
        clim_end_use = int(clim_end_year)

    # if user-provided outside the range, fall back to full span
    if (
        (clim_start_use < years.min())
        or (clim_end_use > years.max())
        or (clim_start_use > clim_end_use)
    ):
        clim_start_use = int(years.min())
        clim_end_use   = int(years.max())

    clim_year_str = f"{clim_start_use}_{clim_end_use}"

    # ordinals for mhw.detect (days as integers)
    t_ord = np.array([ts.to_pydatetime().toordinal() for ts in t_pd], dtype=int)

    attrs_obj = CROCO_Attrs()

    if nc_out is None:
        nc_out = f"mhw_output_{clim_year_str}_{depth_label}.nc"

    # -------------------------------
    # GRID MODE
    # -------------------------------
    if mode == "grid":
        if "mask_rho" in ds:
            surf_masked = surf.where(ds["mask_rho"] > 0)
        else:
            surf_masked = surf

        def _detect_1d(temp_1d, tord, cstart, cend):
            """
            1-D helper for apply_ufunc in grid mode.

            Returns seas, thresh, is_mhw (no raw temp, no event table).
            """
            x = np.asarray(temp_1d, dtype=float)
            L = x.size

            if L == 0 or np.all(np.isnan(x)):
                seas_local = np.full(L, np.nan, dtype="float32")
                thr_local  = np.full(L, np.nan, dtype="float32")
                flag_local = np.zeros(L, dtype="int8")
                return seas_local, thr_local, flag_local

            if np.any(np.isnan(x)):
                ii = np.flatnonzero(~np.isnan(x))
                if ii.size < 5:
                    seas_local = np.full(L, np.nan, dtype="float32")
                    thr_local  = np.full(L, np.nan, dtype="float32")
                    flag_local = np.zeros(L, dtype="int8")
                    return seas_local, thr_local, flag_local
                x = np.interp(np.arange(L), ii, x[ii])

            clim_period = [int(cstart), int(cend)]

            mhw_out_local, clim_local = mhw.detect(
                tord,
                x,
                climatologyPeriod=clim_period,
                pctile=pctile,
                windowHalfWidth=5,
                smoothPercentile=True,
                smoothPercentileWidth=31,
                minDuration=5,
                joinAcrossGaps=True,
                maxGap=2,
                maxPadLength=False,
                coldSpells=False,
                alternateClimatology=False,
                Ly=False,
            )

            seas_local = np.asarray(clim_local["seas"],   dtype="float32")
            thr_local  = np.asarray(clim_local["thresh"], dtype="float32")

            flag_local = np.zeros(L, dtype="int8")
            N_local = int(mhw_out_local["n_events"])
            for ev in range(N_local):
                i0 = int(mhw_out_local["index_start"][ev])
                i1 = int(mhw_out_local["index_end"][ev])
                flag_local[i0:i1+1] = 1

            return seas_local, thr_local, flag_local

        seas, thresh, is_mhw = xr.apply_ufunc(
            _detect_1d,
            surf_masked,
            input_core_dims=[["time"]],
            output_core_dims=[["time"], ["time"], ["time"]],
            output_dtypes=[np.float32, np.float32, np.int8],
            vectorize=True,
            dask="parallelized",
            kwargs={"tord": t_ord, "cstart": clim_start_use, "cend": clim_end_use},
        )

        # use temp internally but do NOT write it
        temp = surf_masked.astype("float32")

        anom_relSeas   = (temp - seas).astype("float32")
        anom_relThresh = (temp - thresh).astype("float32")

        delta = xr.where(thresh > seas, thresh - seas, 0.0)

        category_code = xr.zeros_like(temp, dtype=np.int8)
        category_code = xr.where(temp >= thresh,                   1, category_code)
        category_code = xr.where(temp >= (thresh + 1.0 * delta),   2, category_code)
        category_code = xr.where(temp >= (thresh + 2.0 * delta),   3, category_code)
        category_code = xr.where(temp >= (thresh + 3.0 * delta),   4, category_code)

        coords = {dim: surf_masked[dim] for dim in surf_masked.dims}

        ds_out = xr.Dataset(
            data_vars=dict(
                seas=seas,
                thresh=thresh,
                anom_relSeas=anom_relSeas,
                anom_relThresh=anom_relThresh,
                is_mhw=is_mhw,
                category_code=category_code,
            ),
            coords=coords,
            attrs=dict(
                title="Gridded MHW daily fields on CROCO grid",
                method="Hobday et al. (2016) via marineHeatWaves.detect()",
                percentile=str(pctile),
                windowHalfWidth="5 days",
                smoothPercentileWidth="31 days",
                minDuration="5 days",
                joinAcrossGaps="true",
                maxGap="2 days",
                climatologyPeriod=f"{clim_start_use}-{clim_end_use}",
                depth_selection=depth_label,
                source="somisana-croco",
            ),
        )

        # attach static fields if available
        for vname in ("lon_rho", "lat_rho", "mask_rho"):
            if vname in ds:
                ds_out = ds_out.assign_coords({vname: ds[vname]})

        # apply CROCO_Attrs to known vars
        for name in ("seas", "thresh", "anom_relSeas", "anom_relThresh", "is_mhw", "category_code"):
            if name in ds_out:
                ds_out[name] = change_attrs(attrs_obj, ds_out[name], name)

        encoding = {name: {"zlib": True, "complevel": 4} for name in ds_out.data_vars}
        ds_out.to_netcdf(nc_out, engine=engine, encoding=encoding)

        print(f"\n  [MHW grid] NetCDF written: {nc_out}")
        print("  Sizes:", dict(ds_out.sizes))
        print("  Vars:", ", ".join(sorted(list(ds_out.data_vars))))

        if write_climatology:
            if clim_out is None:
                clim_out = f"mhw_climatology_{pctile}th_{clim_year_str}_{depth_label}.nc"
            clim_vars = ds_out[["seas", "thresh"]].astype("float32")
            try:
                encoding_clim = {
                    vname: {"zlib": True, "complevel": 4}
                    for vname in clim_vars.data_vars
                }
                clim_vars.to_netcdf(clim_out, encoding=encoding_clim)
            except Exception as e:
                print(f"  compressed write of climatology failed ({e}); retrying without compression.")
                clim_vars.to_netcdf(clim_out)
            print(f"  Climatology file written: {clim_out}")

        ds.close()
        return ds_out

    # -------------------------------
    # 1-D MODES: area_mean / nearest
    # -------------------------------
    if level is None:
        ds.close()
        raise ValueError("1-D modes require a single s_rho level; set level (e.g. -1 for surface).")

    da2d = surf  # (time, eta_rho, xi_rho)

    if mode == "area_mean":
        if "mask_rho" in ds:
            w = xr.where(ds["mask_rho"] > 0, 1.0, np.nan)
            da1d = (da2d * w).mean(dim=("eta_rho", "xi_rho"), skipna=True)
        else:
            da1d = da2d.mean(dim=("eta_rho", "xi_rho"), skipna=True)
        use_area_mean = True
        lon_used, lat_used = None, None

    elif mode == "nearest":
        da2d = surf

        if grdname is None:
            if isinstance(fname, str):
                grdname_local = fname
            else:
                grdname_local = fname[0]
        else:
            grdname_local = grdname

        j, i = find_nearest_point(grdname_local, target_lon, target_lat, Bottom)
        da1d = da2d.isel(eta_rho=j, xi_rho=i)

        use_area_mean = False
        lon_used = float(da2d.lon_rho.values[j, i]) if "lon_rho" in da2d.coords else None
        lat_used = float(da2d.lat_rho.values[j, i]) if "lat_rho" in da2d.coords else None

    else:
        ds.close()
        raise ValueError("mode must be 'grid', 'area_mean', or 'nearest'.")

    # 1-D series
    t_pd_1d  = pd.to_datetime(da1d["time"].values)
    t_ord_1d = np.array([ts.to_pydatetime().toordinal() for ts in t_pd_1d], dtype=int)
    temp_1d  = da1d.load().values.astype(float)

    mhw_out, clim = mhw.detect(
        t_ord_1d,
        temp_1d,
        climatologyPeriod=[clim_start_use, clim_end_use],
        pctile=pctile,
        windowHalfWidth=5,
        smoothPercentile=True,
        smoothPercentileWidth=31,
        minDuration=5,
        joinAcrossGaps=True,
        maxGap=2,
        maxPadLength=False,
        coldSpells=False,
        alternateClimatology=False,
        Ly=False,
    )

    clim_label = f"{clim_start_use}-{clim_end_use}"

    _summarize_mhw_output(da1d, mhw_out, clim, clim_period_label=clim_label)

    # --- build lean 1D dataset (time only; no events, no raw temp) ---
    time_vals_1d = da1d["time"].values.astype("datetime64[ns]")
    T = time_vals_1d.size
    N = int(mhw_out["n_events"])

    seas_1d   = np.asarray(clim["seas"],   dtype=np.float32)
    thresh_1d = np.asarray(clim["thresh"], dtype=np.float32)
    temp32_1d = temp_1d.astype("float32")

    anom_relSeas_1d   = temp32_1d - seas_1d
    anom_relThresh_1d = temp32_1d - thresh_1d

    is_mhw_1d = np.zeros(T, dtype="int8")
    for ev in range(N):
        i0 = int(mhw_out["index_start"][ev])
        i1 = int(mhw_out["index_end"][ev])
        is_mhw_1d[i0:i1+1] = 1

    delta_1d = np.where(thresh_1d > seas_1d, thresh_1d - seas_1d, 0.0).astype("float32")
    category_code_1d = np.zeros(T, dtype="int8")
    category_code_1d = np.where(temp32_1d >= thresh_1d,                   1, category_code_1d)
    category_code_1d = np.where(temp32_1d >= (thresh_1d + 1.0 * delta_1d), 2, category_code_1d)
    category_code_1d = np.where(temp32_1d >= (thresh_1d + 2.0 * delta_1d), 3, category_code_1d)
    category_code_1d = np.where(temp32_1d >= (thresh_1d + 3.0 * delta_1d), 4, category_code_1d)

    ds_out = xr.Dataset(
        data_vars=dict(
            seas=("time", seas_1d),
            thresh=("time", thresh_1d),
            anom_relSeas=("time", anom_relSeas_1d),
            anom_relThresh=("time", anom_relThresh_1d),
            is_mhw=("time", is_mhw_1d),
            category_code=("time", category_code_1d),
        ),
        coords=dict(
            time=("time", time_vals_1d),
        ),
        attrs=dict(
            title="1-D MHW daily diagnostics",
            method="Hobday et al. (2016) via marineHeatWaves.detect()",
            percentile=str(pctile),
            windowHalfWidth="5 days",
            smoothPercentileWidth="31 days",
            minDuration="5 days",
            joinAcrossGaps="true",
            maxGap="2 days",
            climatologyPeriod=clim_label,
            input_reduction="area_mean" if use_area_mean else "nearest_point",
            depth_selection=depth_label,
            source="somisana-croco",
        ),
    )

    for name in ("seas", "thresh", "anom_relSeas", "anom_relThresh", "is_mhw", "category_code"):
        if name in ds_out:
            ds_out[name] = change_attrs(attrs_obj, ds_out[name], name)

    if not use_area_mean and (lon_used is not None) and (lat_used is not None):
        ds_out = ds_out.assign_attrs(point_lon=float(lon_used), point_lat=float(lat_used))

    encoding = {name: {"zlib": True, "complevel": 4} for name in ds_out.data_vars}
    ds_out.to_netcdf(nc_out, engine=engine, encoding=encoding)

    print(f"\n  [MHW 1-D] NetCDF written: {nc_out}")
    print("  Variables written:", ", ".join(sorted(list(ds_out.data_vars))))
    print("  Dims:", dict(ds_out.dims))

    if write_climatology:
        if clim_out is None:
            clim_out = f"mhw_climatology_{pctile}th_{clim_year_str}_{depth_label}.nc"

        clim_vars = ds_out[["seas", "thresh"]].astype("float32")
        clim_vars.attrs.update(
            dict(
                title="Seasonal climatology and MHW percentile threshold",
                method="Hobday et al. (2016) via marineHeatWaves.detect()",
                climatologyPeriod=clim_label,
                source="somisana-croco",
            )
        )

        try:
            encoding_clim = {
                vname: {"zlib": True, "complevel": 4}
                for vname in clim_vars.data_vars
            }
            clim_vars.to_netcdf(clim_out, encoding=encoding_clim)
        except Exception as e:
            print(f"  compressed write of climatology failed ({e}); retrying without compression.")
            clim_vars.to_netcdf(clim_out)

        print(f"  Climatology file written: {clim_out}")

    ds.close()
    return ds_out

