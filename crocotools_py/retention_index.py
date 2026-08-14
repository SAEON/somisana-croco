"""
Coastal retention index — box-based, surface, hourly release / daily output.

Pipeline:
  1. build_coast_boxes_km   -> physical (km) sized boxes walked along your
                                coast_order path, replacing the fixed-degree
                                boxes used in marineheatwaves.py's flag map.
  2. seed_opendrift_boxes   -> seeds particles hourly into each box using
                                OpenDrift's origin_marker so each particle
                                keeps a record of which box it started in.
  3. compute_retention      -> post-processes the OpenDrift output file into
                                a per-box, per-time retention fraction using
                                "currently inside the origin box" (particles
                                may leave and return), with beaching tracked
                                as a separate category by default.

Requires: pyproj, shapely, opendrift  (pip install pyproj shapely opendrift --break-system-packages)
"""

import numpy as np
import xarray as xr
from pyproj import Geod
from shapely.geometry import Point, Polygon

GEOD = Geod(ellps="WGS84")


# 1. Box geometry in physical units

def _densify_path(targets, coast_order, spacing_km=1.0):
    """Walk the coast_order polyline and return a densely sampled
    (lon, lat, cumulative_km) path using geodesic interpolation."""
    dense_lons, dense_lats = [], []
    for k in range(len(coast_order) - 1):
        lon0, lat0 = targets[coast_order[k]]
        lon1, lat1 = targets[coast_order[k + 1]]
        _, _, dist_m = GEOD.inv(lon0, lat0, lon1, lat1)
        npts = max(1, int(dist_m / 1000.0 / spacing_km))
        if k == 0:
            dense_lons.append(lon0)
            dense_lats.append(lat0)
        if npts > 0:
            intermediate = GEOD.npts(lon0, lat0, lon1, lat1, npts)
            for lo, la in intermediate:
                dense_lons.append(lo)
                dense_lats.append(la)
        dense_lons.append(lon1)
        dense_lats.append(lat1)

    dense_lons = np.array(dense_lons)
    dense_lats = np.array(dense_lats)

    seg_km = np.zeros(len(dense_lons))
    for i in range(1, len(dense_lons)):
        _, _, d = GEOD.inv(dense_lons[i - 1], dense_lats[i - 1], dense_lons[i], dense_lats[i])
        seg_km[i] = d / 1000.0
    cum_km = np.cumsum(seg_km)
    return dense_lons, dense_lats, cum_km


def _square_polygon_km(center_lon, center_lat, box_km):
    """Build an approximately-square polygon (geodesic) of side box_km,
    centered on (center_lon, center_lat)."""
    half = (box_km / 2.0) * 1000.0  # meters
    corners = []
    for az in (315, 45, 135, 225):  # NW, NE, SE, SW
        lon_c, lat_c, _ = GEOD.fwd(center_lon, center_lat, az, half * np.sqrt(2))
        corners.append((lon_c, lat_c))
    return Polygon(corners)


def build_coast_boxes_km(targets, coast_order, box_km=25.0, offshore_km=-5.0,
                          step_km=25.0, mask_check_fn=None):
    """
    Replaces the fixed-degree box construction in marineheatwaves.py's
    plot_flag_map with physically-sized (km) boxes.

    Parameters
    ----------
    targets : dict[str, (lon, lat)]      -- same as TARGETS in your script
    coast_order : list[str]              -- same as coast_order
    box_km : float                       -- box side length, km
    offshore_km : float                  -- perpendicular offset from coast,
                                             negative = seaward (matches your
                                             existing sign convention)
    step_km : float                      -- spacing between box centers, km
    mask_check_fn : callable(lon, lat) -> bool, optional
                                          -- return True if (lon, lat) is wet
                                             water on your CROCO grid. If
                                             given, boxes centered on land
                                             are dropped and a warning is
                                             printed (you'll still want to
                                             filter individual seed points
                                             at seeding time — see step 2).

    Returns
    -------
    list of dicts: {'id', 'nearest_site', 'center_lon', 'center_lat', 'polygon'}
    """
    dense_lons, dense_lats, cum_km = _densify_path(targets, coast_order)

    boxes = []
    box_id = 0
    for bd in np.arange(0, cum_km[-1], step_km):
        cx = np.interp(bd, cum_km, dense_lons)
        cy = np.interp(bd, cum_km, dense_lats)

        idx = max(1, min(np.searchsorted(cum_km, bd), len(cum_km) - 1))
        az_fwd, _, _ = GEOD.inv(dense_lons[idx - 1], dense_lats[idx - 1],
                                 dense_lons[idx], dense_lats[idx])
        az_perp = az_fwd + 90.0  # perpendicular to local coast direction

        off_lon, off_lat, _ = GEOD.fwd(cx, cy, az_perp, offshore_km * 1000.0)

        if mask_check_fn is not None and not mask_check_fn(off_lon, off_lat):
            print(f"  [skip] box at {bd:.0f} km ({off_lon:.3f}, {off_lat:.3f}) centers on land")
            continue

        nearest_site = min(targets, key=lambda s: GEOD.inv(off_lon, off_lat, *targets[s])[2])

        boxes.append({
            "id": box_id,
            "nearest_site": nearest_site,
            "center_lon": off_lon,
            "center_lat": off_lat,
            "polygon": _square_polygon_km(off_lon, off_lat, box_km),
        })
        box_id += 1

    print(f"Built {len(boxes)} boxes, {box_km} km square, {step_km} km spacing")
    return boxes


# 2. OpenDrift seeding template

def seed_opendrift_boxes(boxes, croco_file, start_time, end_time,
                          particles_per_box_per_hour=20,
                          horizontal_diffusivity=1.0,
                          outfile="retention_run.nc"):
    """
    Template for the seeding + run step. Fill in your actual CROCO file
    path and adjust config to match your earlier decisions (surface-only,
    stranding action explicit, RK4 advection).
    """
    from opendrift.readers import reader_ROMS_native
    from opendrift.models.oceandrift import OceanDrift
    import pandas as pd

    o = OceanDrift(loglevel=20)
    reader = reader_ROMS_native.Reader(croco_file)
    o.add_reader(reader)

    # --- config matching the surface-only, first-pass decisions made earlier ---
    o.set_config('drift:vertical_mixing', False)
    o.set_config('drift:advection_scheme', 'runge-kutta4')
    o.set_config('drift:horizontal_diffusivity', horizontal_diffusivity)
    o.set_config('general:coastline_action', 'stranding')  # explicit, not 'previous'
    o.set_config('drift:current_uncertainty', 0)

    release_times = pd.date_range(start_time, end_time, freq='1H')

    for box in boxes:
        lons, lats = box["polygon"].exterior.xy
        o.seed_within_polygon(
            lons=list(lons), lats=list(lats),
            number=particles_per_box_per_hour * len(release_times),
            time=list(release_times),
            z=0,  # surface
            origin_marker=box["id"],
        )

    o.run(end_time=pd.Timestamp(end_time), time_step=3600,
          time_step_output=86400,  # daily output, as decided
          outfile=outfile)

    return outfile


# 3. Retention calculation

def compute_retention(outfile, boxes, count_beached_as_retained=False):
    """
    Reads the OpenDrift output trajectory file and computes, per box and
    per output time, the retention fraction using "currently inside the
    origin box at time t" (particles may leave and return).

    Denominator per box = number of that box's particles that were
    successfully placed in water at t=0 (excludes any seeded directly on
    land, which strand immediately and would otherwise deflate the index
    artificially).

    Returns an xarray.Dataset: dims (box, time), variables:
      - retained_frac   : fraction currently inside origin box, status active
      - stranded_frac    : fraction beached, reported separately unless
                            count_beached_as_retained=True (then folded into
                            retained_frac instead)
      - n_seeded         : particles successfully seeded in water, per box
    """
    ds = xr.open_dataset(outfile)

    lon = ds['lon'].values          # (trajectory, time)
    lat = ds['lat'].values
    status = ds['status'].values    # 0=active, 1=stranded (check your OpenDrift version's status_dict)
    origin = ds['origin_marker'].values[:, 0]  # constant per trajectory
    times = ds['time'].values
    n_time = len(times)
    n_boxes = len(boxes)

    retained_frac = np.full((n_boxes, n_time), np.nan)
    stranded_frac = np.full((n_boxes, n_time), np.nan)
    n_seeded = np.zeros(n_boxes, dtype=int)

    for b, box in enumerate(boxes):
        poly = box["polygon"]
        traj_mask = origin == box["id"]

        # denominator: particles from this box valid (not stranded) at t=0
        valid0 = traj_mask & (status[:, 0] == 0)
        n_seeded[b] = valid0.sum()
        if n_seeded[b] == 0:
            continue

        for t in range(n_time):
            sel = np.where(valid0)[0]
            lo_t, la_t, st_t = lon[sel, t], lat[sel, t], status[sel, t]

            active = st_t == 0
            stranded = st_t == 1

            inside_active = np.array([
                poly.contains(Point(lo, la)) for lo, la, a in zip(lo_t, la_t, active) if a
            ])
            n_inside_active = inside_active.sum() if len(inside_active) else 0

            if count_beached_as_retained:
                inside_stranded = np.array([
                    poly.contains(Point(lo, la)) for lo, la, s in zip(lo_t, la_t, stranded) if s
                ])
                n_inside_stranded = inside_stranded.sum() if len(inside_stranded) else 0
                retained_frac[b, t] = (n_inside_active + n_inside_stranded) / n_seeded[b]
            else:
                retained_frac[b, t] = n_inside_active / n_seeded[b]
                stranded_frac[b, t] = stranded.sum() / n_seeded[b]

    out = xr.Dataset(
        {
            "retained_frac": (("box", "time"), retained_frac),
            "stranded_frac": (("box", "time"), stranded_frac),
            "n_seeded": (("box",), n_seeded),
        },
        coords={
            "box": [b["id"] for b in boxes],
            "nearest_site": ("box", [b["nearest_site"] for b in boxes]),
            "time": times,
        },
    )
    return out


if __name__ == "__main__":
    # import TARGETS and coast_order directly from marineheatwaves.py
    from marineheatwaves import TARGETS, coast_order

    boxes = build_coast_boxes_km(
        targets=TARGETS,
        coast_order=coast_order,
        box_km=25.0,
        offshore_km=-5.0,
        step_km=25.0,
    )

    outfile = seed_opendrift_boxes(
        boxes,
        croco_file="/home/pmvula/PRODUCTS/FORECAST/sa_west_02/croco_v1.3.1/MERCATOR-GFS/RETENTION/data/croco_avg.nc",
        start_time="2026-08-08",
        end_time="2026-08-16",
    )

    retention = compute_retention(outfile, boxes, count_beached_as_retained=False)
    retention.to_netcdf("/home/pmvula/PRODUCTS/FORECAST/sa_west_02/croco_v1.3.1/MERCATOR-GFS/RETENTION/data/retention_by_box.nc")
    print(retention)
