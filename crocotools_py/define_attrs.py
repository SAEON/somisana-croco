'''
CF-compliant attributes for CROCO output variables.

ATTRS: dict of variable name -> VarAttrs for scalars and grid variables.
VECTOR_ATTRS: dict of variable name -> (grid-aligned VarAttrs, rotated VarAttrs)
              for vector components that have different metadata depending on
              whether they have been rotated to east/north components.

apply_attrs(): applies CF attributes to a DataArray, with a warning if
               the variable is not in the registry.
'''

from dataclasses import dataclass

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
}

# Vector variables: (grid-aligned attrs, rotated/east-north attrs)
VECTOR_ATTRS = {
    'u':         (VarAttrs('Sea water velocity in x direction', 'm s-1', 'sea_water_x_velocity'),
                  VarAttrs('Eastward component of baroclinic velocity', 'm s-1', 'baroclinic_eastward_sea_water_velocity')),
    'u_anom':    (VarAttrs('Sea water velocity in x direction Anomaly', 'm s-1', 'sea_water_x_velocity_anomaly'),
                  VarAttrs('Eastward component of baroclinic velocity anomaly', 'm s-1', 'baroclinic_eastward_sea_water_velocity_anomaly')),
    'v':         (VarAttrs('Sea water velocity in y direction', 'm s-1', 'sea_water_y_velocity'),
                  VarAttrs('Northward component of baroclinic velocity', 'm s-1', 'baroclinic_northward_sea_water_velocity')),
    'v_anom':    (VarAttrs('Sea water velocity in y direction Anomaly', 'm s-1', 'sea_water_y_velocity_anom'),
                  VarAttrs('Northward component of baroclinic velocity anomaly', 'm s-1', 'baroclinic_northward_sea_water_velocity_anomaly')),
    'ubar':      (VarAttrs('Barotropic velocity of sea water in x direction', 'm s-1', 'barotropic_sea_water_x_velocity'),
                  VarAttrs('Eastward component of barotropic velocity', 'm s-1', 'barotropic_eastward_sea_water_velocity')),
    'vbar':      (VarAttrs('Barotropic velocity of sea water in y direction', 'm s-1', 'barotropic_sea_water_y_velocity'),
                  VarAttrs('Northward component of barotropic velocity', 'm s-1', 'barotropic_northward_sea_water_velocity')),
    'sustr':     (VarAttrs('Wind stress on sea surface in x direction', 'N m-2', 'surface_downward_x_stress'),
                  VarAttrs('Eastward component of surface stress', 'N m-2', 'surface_eastward_stress')),
    'svstr':     (VarAttrs('Wind stress on sea surface in y direction', 'N m-2', 'surface_downward_y_stress'),
                  VarAttrs('Northward component of surface stress', 'N m-2', 'surface_northward_stress')),
    'bustr':     (VarAttrs('Stress due to sea water on sea floor in x direction', 'N m-2', 'stress_due_to_sea_water_on_sea_floor_in_x_direction'),
                  VarAttrs('Eastward component of bottom stress', 'N m-2', 'bottom_eastward_stress')),
    'bvstr':     (VarAttrs('Stress due to sea water on sea floor in y direction', 'N m-2', 'stress_due_to_sea_water_on_sea_floor_in_y_direction'),
                  VarAttrs('Northward component of bottom stress', 'N m-2', 'bottom_northward_stress')),
}

def apply_attrs(da, var_str, rotated=False):
    """
    Apply CF-compliant attributes to a DataArray.

    Parameters
    ----------
    da       : xarray DataArray to update
    var_str  : variable name to look up in the registry
    rotated  : if True, use east/north attrs for vector variables

    Returns the DataArray (modified in place).
    """
    if var_str in ATTRS:
        meta = ATTRS[var_str]
    elif var_str in VECTOR_ATTRS:
        meta = VECTOR_ATTRS[var_str][1 if rotated else 0]
    else:
        print(f"WARNING: no CF compliant attributes defined for '{var_str}' - metadata not updated to be CF compliant")
        return da
    da.attrs['long_name'] = meta.long_name
    da.attrs['units'] = meta.units
    da.attrs['standard_name'] = meta.standard_name
    return da
