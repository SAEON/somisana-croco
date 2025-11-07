'''
Here we define attributes of CROCO variables to make them CF-Compliant
We have different classes for variables as read from the output file(s)
and then abother for one with rotated vectors i.e. vectors which represent East/North components
'''

class VariableMetadata:
    def __init__(self, long_name, units, standard_name):
        self.long_name = long_name
        self.units = units
        self.standard_name = standard_name

    def __repr__(self):
        return f"VariableMetadata(long_name='{self.long_name}', units='{self.units}', standard_name='{self.standard_name}')"

class CROCO_Attrs_RotatedVectors:
    def __init__(self):
        self.xi_rho = VariableMetadata('x-dimension of the grid', '1', 'x_grid_index')
        self.eta_rho= VariableMetadata('y-dimension of the grid', '1', 'y_grid_index')
        self.lon_rho= VariableMetadata('Longitude', 'degrees_east', 'longitude')
        self.lat_rho= VariableMetadata('Latitude', 'degrees_north', 'latitude')
        self.zeta   = VariableMetadata('Sea Surface Elevation', 'm', 'sea_surface_elevation')
        self.zeta_anom   = VariableMetadata('Sea Surface Elevation Anomaly', 'm', 'sea_surface_elevation_anomaly')
        self.temp   = VariableMetadata('Sea Water Temperature', 'degC', 'sea_water_temperature')
        self.temp_anom   = VariableMetadata('Sea Water Temperature Anomaly', 'degC', 'sea_water_temperature_anomaly')
        self.salt   = VariableMetadata('Sea Water Salinity', '1', 'sea_water_salinity')
        self.salt_anom   = VariableMetadata('Sea Water Salinity Anomaly', '1', 'sea_water_salinity_anomaly')
        self.h      = VariableMetadata('Depth of the sea floor', 'm', 'sea_floor_depth')
        self.mask   = VariableMetadata('Land-sea mask (1=water, 0=land)', '1', 'land_binary_mask')
        self.depth  = VariableMetadata('Depth', 'm', 'depth')
        self.u      = VariableMetadata('Eastward component of baroclinic velocity', 'm s-1', 'baroclinic_eastward_sea_water_velocity')
        self.u_anom      = VariableMetadata('Eastward component of baroclinic velocity anomaly', 'm s-1', 'baroclinic_eastward_sea_water_velocity_anomaly')
        self.sustr  = VariableMetadata('Eastward component of surface stress', 'N m-2', 'surface_eastward_stress')
        self.bustr  = VariableMetadata('Eastward component of bottom stress', 'N m-2', 'bottom_eastward_stress')
        self.ubar   = VariableMetadata('Eastward component of barotropic velocity', 'm s-1', 'barotropic_eastward_sea_water_velocity')
        self.v      = VariableMetadata('Northward component of baroclinic velocity', 'm s-1', 'baroclinic_northward_sea_water_velocity')
        self.v_anom      = VariableMetadata('Northward component of baroclinic velocity anomaly', 'm s-1', 'baroclinic_northward_sea_water_velocity_anomaly')
        self.svstr  = VariableMetadata('Northward component of surface stress', 'N m-2', 'surface_northward_stress')
        self.bvstr  = VariableMetadata('Northward component of bottom stress', 'N m-2', 'bottom_northward_stress')
        self.vbar   = VariableMetadata('Northward component of barotropic velocity', 'm s-1', 'barotropic_northward_sea_water_velocity')
        self.w      = VariableMetadata('Upward seawater velocity', 'm s-1', 'averaged vertical momentum component')

class CROCO_Attrs:
    def __init__(self):
        self.xi_rho = VariableMetadata('x-dimension of the grid', '1', 'x_grid_index')
        self.eta_rho= VariableMetadata('y-dimension of the grid', '1', 'y_grid_index')
        self.lon_rho= VariableMetadata('Longitude', 'degrees_east', 'longitude')
        self.lat_rho= VariableMetadata('Latitude', 'degrees_north', 'latitude')       
        self.zeta   = VariableMetadata('Sea Surface Elevation', 'm', 'sea_surface_elevation')
        self.zeta_anom   = VariableMetadata('Sea Surface Elevation Anomaly', 'm', 'sea_surface_elevation_anomaly')
        self.temp   = VariableMetadata('Sea Water Temperature', 'degC', 'sea_water_temperature')
        self.temp_anom   = VariableMetadata('Sea Water Temperature Anomaly', 'degC', 'sea_water_temperature_anomaly')
        self.salt   = VariableMetadata('Sea Water Salinity', '1', 'sea_water_salinity')
        self.salt_anom   = VariableMetadata('Sea Water Salinity Anomaly', '1', 'sea_water_salinity_anomaly')
        self.h      = VariableMetadata('Depth of the sea floor', 'm', 'sea_floor_depth')
        self.mask   = VariableMetadata('Land-sea mask (1=water, 0=land)', '1', 'land_binary_mask')
        self.depth  = VariableMetadata('Depth', 'm', 'depth')
        self.u      = VariableMetadata('Sea water velocity in x direction', 'm s-1', 'sea_water_x_velocity')
        self.u_anom      = VariableMetadata('Sea water velocity in x direction Anomaly', 'm s-1', 'sea_water_x_velocity_anomaly')
        self.sustr  = VariableMetadata('Wind stress on sea surface in x direction', 'N m-2', 'surface_downward_x_stress')
        self.bustr  = VariableMetadata('Stress due to sea water on sea floor in x direction', 'N m-2', 'stress_due_to_sea_water_on_sea_floor_in_x_direction')
        self.ubar   = VariableMetadata('Barotropic velocity of sea water in x direction', 'm s-1', 'barotropic_sea_water_x_velocity')  
        self.v      = VariableMetadata('Sea water velocity in y direction', 'm s-1', 'sea_water_y_velocity')
        self.v_anom      = VariableMetadata('Sea water velocity in y direction Anomaly', 'm s-1', 'sea_water_y_velocity_anom')
        self.svstr  = VariableMetadata('Wind stress on sea surface in y direction', 'N m-2', 'surface_downward_y_stress')
        self.bvstr  = VariableMetadata('Stress due to sea water on sea floor in y direction', 'N m-2', 'stress_due_to_sea_water_on_sea_floor_in_y_direction')
        self.vbar   = VariableMetadata('Barotropic velocity of sea water in y direction', 'm s-1', 'barotropic_sea_water_y_velocity')
        self.w      = VariableMetadata('Upward seawater velocity', 'm s-1', 'averaged vertical momentum component')
        # ------------------------------------------------------------------
        # MHW daily fields (for run_mhw_detection outputs)
        # ------------------------------------------------------------------
        self.seas = VariableMetadata(
            "Seasonal climatology of sea water temperature",
            "degC",
            "sea_water_temperature_climatology",
        )
        self.thresh = VariableMetadata(
            "Marine heatwave threshold (e.g. 90th percentile)",
            "degC",
            "sea_water_temperature_threshold",
        )
        self.anom_relSeas = VariableMetadata(
            "Temperature anomaly relative to seasonal climatology",
            "degC",
            "sea_water_temperature_anomaly_relative_to_climatology",
        )
        self.anom_relThresh = VariableMetadata(
            "Temperature anomaly relative to MHW threshold",
            "degC",
            "sea_water_temperature_anomaly_relative_to_threshold",
        )
        self.is_mhw = VariableMetadata(
            "Marine heatwave day flag (1=in event, 0=not in event)",
            "1",
            "marine_heatwave_day_flag",
        )
        self.category_code = VariableMetadata(
            "Marine heatwave category code (1=Moderate,2=Strong,3=Severe,4=Extreme)",
            "1",
            "marine_heatwave_category_code",
        )