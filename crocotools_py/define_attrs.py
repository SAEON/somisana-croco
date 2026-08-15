'''
CF-compliant attributes for the files this repo writes.

Variable attributes
-------------------
ATTRS: dict of variable name -> VarAttrs for scalars and grid variables.
VECTOR_ATTRS: dict of variable name -> (grid-aligned, eastnorth-rotated,
              section-rotated) VarAttrs triple for vector components. The
              section-rotated form is "across-section" for x-side keys
              (u, ubar, sustr, ...) and "along-section" for y-side keys
              (v, vbar, svstr, ...).

apply_attrs(): applies CF attributes to a DataArray, with a warning if
               the variable is not in the registry.

This registry is the CF baseline: the three fields every variable must have.
It deliberately does not hold the richer, run-specific attributes some
variables carry (the MHW/MCS category's flag_values and min_duration_days, the
stratification's reference_pressure_dbar and target/deep depths, ...), because
those depend on the settings a particular run was made with and so cannot live
in a static dict. Those are written by whatever computes the variable and are
carried through by postprocess.get_var(), which applies this registry on top of
them rather than in place of them.

Global attributes
-----------------
global_attrs(): builds the global attribute set for an output file, shared by
                products.py and regridding.py so that every file this repo
                writes describes itself the same way. See also coverage_attrs()
                and history_line().
'''

from dataclasses import dataclass
from datetime import datetime, timezone

import numpy as np
import pandas as pd

@dataclass
class VarAttrs:
    long_name: str
    units: str
    standard_name: str

# Scalar and grid variables (same regardless of rotation)
ATTRS = {
    'xi_rho':    VarAttrs('x-dimension of the grid', '1', 'x_grid_index'),
    'eta_rho':   VarAttrs('y-dimension of the grid', '1', 'y_grid_index'),
    'lon_rho':   VarAttrs('Longitude', 'degrees_east', 'longitude'),
    'lat_rho':   VarAttrs('Latitude', 'degrees_north', 'latitude'),
    'h':         VarAttrs('Depth of the sea floor', 'm', 'sea_floor_depth'),
    'mask':      VarAttrs('Land-sea mask (1=water, 0=land)', '1', 'land_binary_mask'),
    'depth':     VarAttrs('Depth', 'm', 'depth'),
    'zeta':      VarAttrs('Sea Surface Elevation', 'm', 'sea_surface_elevation'),
    'zeta_anom': VarAttrs('Sea Surface Elevation Anomaly', 'm', 'sea_surface_elevation_anomaly'),
    'temp':      VarAttrs('Sea Water Temperature', 'degC', 'sea_water_temperature'),
    'temp_anom': VarAttrs('Sea Water Temperature Anomaly', 'degC', 'sea_water_temperature_anomaly'),
    'salt':      VarAttrs('Sea Water Salinity', '1', 'sea_water_salinity'),
    'salt_anom': VarAttrs('Sea Water Salinity Anomaly', '1', 'sea_water_salinity_anomaly'),
    'w':         VarAttrs('Upward seawater velocity', 'm s-1', 'averaged vertical momentum component'),
    # Keep in step with marineheatwaves.category_attrs(), which is what actually
    # writes this variable's metadata - this entry only sets the three fields below.
    'category':  VarAttrs('MHW/MCS combined event category', '1', 'status_flag'),
    'sst_front': VarAttrs('Sea Surface Temperature Horizontal Front Magnitude', 'degC km-1', 'sea_surface_temperature_gradient_magnitude'),
    'density_deep':         VarAttrs('Potential density (referenced to 0 dbar) at the deep reference depth, or at the seabed where shallower', 'kg m-3', 'sea_water_potential_density'),
    'density_target_depth': VarAttrs('Potential density (referenced to 0 dbar) at the target depth', 'kg m-3', 'sea_water_potential_density'),
    'stratification':       VarAttrs('Water column stratification (potential density at the deep reference depth minus potential density at the target depth, both referenced to 0 dbar)', 'kg m-3', 'sea_water_potential_density_difference'),
}

# Vector variables: (grid-aligned attrs, eastnorth-rotated attrs, section-rotated attrs)
VECTOR_ATTRS = {
    'u':         (VarAttrs('Sea water velocity in x direction', 'm s-1', 'sea_water_x_velocity'),
                  VarAttrs('Eastward component of baroclinic velocity', 'm s-1', 'baroclinic_eastward_sea_water_velocity'),
                  VarAttrs('Across-section component of baroclinic velocity', 'm s-1', 'baroclinic_across_section_sea_water_velocity')),
    'u_anom':    (VarAttrs('Sea water velocity in x direction Anomaly', 'm s-1', 'sea_water_x_velocity_anomaly'),
                  VarAttrs('Eastward component of baroclinic velocity anomaly', 'm s-1', 'baroclinic_eastward_sea_water_velocity_anomaly'),
                  VarAttrs('Across-section component of baroclinic velocity anomaly', 'm s-1', 'baroclinic_across_section_sea_water_velocity_anomaly')),
    'v':         (VarAttrs('Sea water velocity in y direction', 'm s-1', 'sea_water_y_velocity'),
                  VarAttrs('Northward component of baroclinic velocity', 'm s-1', 'baroclinic_northward_sea_water_velocity'),
                  VarAttrs('Along-section component of baroclinic velocity', 'm s-1', 'baroclinic_along_section_sea_water_velocity')),
    'v_anom':    (VarAttrs('Sea water velocity in y direction Anomaly', 'm s-1', 'sea_water_y_velocity_anom'),
                  VarAttrs('Northward component of baroclinic velocity anomaly', 'm s-1', 'baroclinic_northward_sea_water_velocity_anomaly'),
                  VarAttrs('Along-section component of baroclinic velocity anomaly', 'm s-1', 'baroclinic_along_section_sea_water_velocity_anomaly')),
    'ubar':      (VarAttrs('Barotropic velocity of sea water in x direction', 'm s-1', 'barotropic_sea_water_x_velocity'),
                  VarAttrs('Eastward component of barotropic velocity', 'm s-1', 'barotropic_eastward_sea_water_velocity'),
                  VarAttrs('Across-section component of barotropic velocity', 'm s-1', 'barotropic_across_section_sea_water_velocity')),
    'vbar':      (VarAttrs('Barotropic velocity of sea water in y direction', 'm s-1', 'barotropic_sea_water_y_velocity'),
                  VarAttrs('Northward component of barotropic velocity', 'm s-1', 'barotropic_northward_sea_water_velocity'),
                  VarAttrs('Along-section component of barotropic velocity', 'm s-1', 'barotropic_along_section_sea_water_velocity')),
    'sustr':     (VarAttrs('Wind stress on sea surface in x direction', 'N m-2', 'surface_downward_x_stress'),
                  VarAttrs('Eastward component of surface stress', 'N m-2', 'surface_eastward_stress'),
                  VarAttrs('Across-section component of surface stress', 'N m-2', 'surface_across_section_stress')),
    'svstr':     (VarAttrs('Wind stress on sea surface in y direction', 'N m-2', 'surface_downward_y_stress'),
                  VarAttrs('Northward component of surface stress', 'N m-2', 'surface_northward_stress'),
                  VarAttrs('Along-section component of surface stress', 'N m-2', 'surface_along_section_stress')),
    'bustr':     (VarAttrs('Stress due to sea water on sea floor in x direction', 'N m-2', 'stress_due_to_sea_water_on_sea_floor_in_x_direction'),
                  VarAttrs('Eastward component of bottom stress', 'N m-2', 'bottom_eastward_stress'),
                  VarAttrs('Across-section component of bottom stress', 'N m-2', 'bottom_across_section_stress')),
    'bvstr':     (VarAttrs('Stress due to sea water on sea floor in y direction', 'N m-2', 'stress_due_to_sea_water_on_sea_floor_in_y_direction'),
                  VarAttrs('Northward component of bottom stress', 'N m-2', 'bottom_northward_stress'),
                  VarAttrs('Along-section component of bottom stress', 'N m-2', 'bottom_along_section_stress')),
}

def apply_attrs(da, var_str, rotated=False, section=False):
    """
    Apply CF-compliant attributes to a DataArray.

    Parameters
    ----------
    da       : xarray DataArray to update
    var_str  : variable name to look up in the registry
    rotated  : if True, use east/north attrs for vector variables.
               Ignored when section=True.
    section  : if True, use section-rotated (across/along) attrs for vector
               variables. Takes precedence over `rotated`.

    Returns the DataArray (modified in place).
    """
    if var_str in ATTRS:
        meta = ATTRS[var_str]
    elif var_str in VECTOR_ATTRS:
        if section:
            meta = VECTOR_ATTRS[var_str][2]
        elif rotated:
            meta = VECTOR_ATTRS[var_str][1]
        else:
            meta = VECTOR_ATTRS[var_str][0]
    else:
        print(f"WARNING: no CF compliant attributes defined for '{var_str}' - metadata not updated to be CF compliant")
        return da
    da.attrs['long_name'] = meta.long_name
    da.attrs['units'] = meta.units
    da.attrs['standard_name'] = meta.standard_name
    return da


# Global attributes
#
# Everything below builds the global attribute set of an output file. It lives
# here rather than in products.py or regridding.py because both need it and
# neither should import the other.

INSTITUTION = ('South African Environmental Observation Network (SAEON), '
               'a facility of the National Research Foundation (NRF)')

REFERENCES = ('Project: Sustainable Ocean Modelling Initiative: a South AfricaN '
              'Approach (SOMISANA; https://somisana.ac.za/); '
              'Tools: https://github.com/SAEON/somisana-croco')

CONVENTIONS = 'CF-1.8, ACDD-1.3'

# Global attributes carried forward from the file a step read its data from.
#
# This is a whitelist, not a blanket copy, because a step can equally be run on
# raw CROCO output - and a raw CROCO header carries SRCS, CPP-options, the
# rst/his/avg/grd filenames and the whole timestep and mixing configuration,
# none of which describes the file being written. A raw input simply has none
# of the keys below, so it needs no special-casing: nothing is carried and the
# output gets the descriptive set built from scratch.
#
# 'history' is handled separately (it is appended to, not replaced), and the
# geospatial/temporal coverage is not carried at all - coverage_attrs() derives
# it from the output dataset, which stays correct even when a step changes the
# grid, as tier 3 regridding does.
CARRIED_ATTRS = ['domain']


def history_line(action):
    """
    One CF style 'history' entry: a UTC timestamp and what was done.

    Each step that writes to a file adds a line, so the file records what was
    done to it and when, in the order it happened.
    """
    stamp = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S')
    return f'{stamp} UTC: {action}'


def extend_history(existing, action):
    """Append a history line to whatever history the input file already had."""
    return '\n'.join(filter(None, [str(existing or '').strip(),
                                   history_line(action)]))


def coverage_attrs(ds):
    """
    The ACDD geospatial and temporal coverage attributes of ds.

    Derived from the dataset being written rather than carried over from its
    input, so that a step which changes the grid or the time axis cannot leave
    stale coverage behind - tier 3 regridding builds a new regular grid, and
    carrying the curvilinear extents across would misreport it.

    Handles both the curvilinear grids (lon_rho/lat_rho) and the rectilinear
    ones tier 3 produces (longitude/latitude).
    """
    attrs = {}

    for lon_name, lat_name in (('lon_rho', 'lat_rho'), ('longitude', 'latitude')):
        if lon_name in ds.variables and lat_name in ds.variables:
            lon, lat = np.asarray(ds[lon_name]), np.asarray(ds[lat_name])
            attrs.update({'geospatial_lon_min': float(np.nanmin(lon)),
                          'geospatial_lon_max': float(np.nanmax(lon)),
                          'geospatial_lat_min': float(np.nanmin(lat)),
                          'geospatial_lat_max': float(np.nanmax(lat)),
                          'geospatial_lon_units': 'degrees_east',
                          'geospatial_lat_units': 'degrees_north'})
            break

    if 'time' in ds.variables and ds['time'].size:
        times = pd.to_datetime(np.asarray(ds['time'].values))
        attrs.update({'time_coverage_start': times.min().strftime('%Y-%m-%dT%H:%M:%SZ'),
                      'time_coverage_end': times.max().strftime('%Y-%m-%dT%H:%M:%SZ')})

    attrs['cdm_data_type'] = 'Grid'
    return attrs


def global_attrs(title, summary, source, ds, action, src_attrs=None,
                 doi_link=None, extra=None):
    """
    Build the global attributes for a file written by this repo.

    Parameters
    ----------
    title     : short name for what this file is
    summary   : one paragraph describing the file. File level only - anything
                specific to a single variable belongs on that variable, not
                here, so that it cannot go stale when the variable changes
    source    : the file (or files) this one was built from
    ds        : the dataset about to be written, read for its coverage
    action    : what this step did, appended to 'history'
    src_attrs : global attributes of the input file, if it had any. Only the
                keys in CARRIED_ATTRS are taken, plus 'history', which is
                extended rather than replaced - that chain is what carries the
                provenance of a products file through the regridding tiers
    doi_link  : bare DOI, written as a full https://doi.org/ URL
    extra     : any further attributes to add, e.g. the vertical coordinate
                parameters a CROCO format file has to keep

    Returns a dict ready to assign to Dataset.attrs.
    """
    src_attrs = dict(src_attrs or {})

    attrs = {'title': title,
             'summary': summary,
             'institution': INSTITUTION,
             'source': source,
             'references': REFERENCES,
             'Conventions': CONVENTIONS}

    for key in CARRIED_ATTRS:
        if str(src_attrs.get(key, '')).strip():
            attrs[key] = src_attrs[key]

    attrs.update(coverage_attrs(ds))

    if extra:
        attrs.update(extra)
    if doi_link is not None:
        attrs['doi'] = f'https://doi.org/{doi_link}'

    attrs['date_created'] = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    attrs['history'] = extend_history(src_attrs.get('history'), action)
    return attrs
