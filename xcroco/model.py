"""
Define dictionnaries to change the model names by standart names
"""
class Model:
    
    def __init__(self, name, ):
        self.name = name

        if name == "croco_xios":
            self.rename_vars = {
                # surface wind stress
                "sustr"        : "xtau_sfc_u",        # x-wind stress component
                "svstr"        : "ytau_sfc_v",        # y-wind stress component
                "wstr"         : "tau_sfc",           # wind stress module
                # surface fluxes
                "radsw"        : "sw_sfc_down",       # surface_downward_shortwave_flux
                "shflx_rlw"    : "lw_sfc_down",       # surface_downward_longwave_flux
                "shflux"       : "hf_net_sfc_down",   # surface_downward_shortwave_flux
                "shflx_lat"    : "lh_sfc_down",       # surface_downward_latent_heat_flux
                "shflx_sen"    : "sh_sfc_down",       # surface_downward_sensible_heat_flux
                "swflux"       : "ep_flux_sfc_down",  # surface_downward_freshwater_flux (E-P)
                # ocean surface
                "zeta"         : "z_sfc",             # sea_surface_height
                "ssh"          : "z_sfc",             # sea_surface_height
                # currents
                "u"            : "xcur",              # x-current component
                "v"            : "ycur",              # y-current component
                "w"            : "zcur",              # z-current component
                "omega"        : "zcur_sc",           # z-current component in sigma coordinate
                "ubar"         : "xcur_btrope",       # x-barotropic current component
                "vbar"         : "ycur_btrope",       # y-barotropic current component 
                # MLD
                "hbls"         : "mld_turb",          # ocean_turbulent_mixed_layer_depth
                "hbl"          : "mld_turb",          # ocean_turbulent_mixed_layer_depth
                # tracers et al.
                "temp"         : "temp",              # ocean potential temperature
                "salt"         : "salt",              # ocean salinity
                # other diag variables
                "Akt"          : "zdiff_tra_w",       # ocean_vertical_tracer_diffusivity_at_w_points
                "Akv"          : "zdiff_mtm_w",       # ocean_vertical_momentum_diffusivity
                # bottom layer
                "bostr"        : "tau_bot",           # ocean bottom stress
                "hbbl"         : "bbl",               # bottom boundary layer
                # runoff
                "Qbar"         : "runoff",            # river runoff
                "temp_src"     : "runoff_temp",       # runoff temperature
                "salt_src"     : "runoff_salt",       # runoff salinity
                "runoff_name"  : "runoff_name",       # runoff_name
                # grid et al.
                "h"            : "h",                 # terrain height or bathy (positive upward, static)
                "nav_lat_rho"  : "lat",               # latitude at mass points
                "lat_rho"      : "lat",               # latitude at mass points
                "nav_lon_rho"  : "lon",               # longitude at mass points
                "lon_rho"      : "lon",               # longitude at mass points
                "nav_lat_u"    : "lat_u",             # latitude at u points
                "lat_u"        : "lat_u",             # latitude at u points
                "nav_lon_u"    : "lon_u",             # longitude at u points
                "lon_u"        : "lon_u",             # longitude at u points
                "nav_lat_v"    : "lat_v",             # latitude at v points
                "lat_v"        : "lat_v",             # latitude at v points
                "nav_lon_v"    : "lon_v",             # longitude at v points
                "lon_v"        : "lon_v",             # longitude at v points
                "nav_lon"      : "lon",               # longitude at mass points
                "nav_lat"      : "lat",               # latitude at mass points
                "mask_rho"     : "mask",              # land mask
                "f"            : "f",                 # Coriolis parameter (1/seconds) at RHO-points
                "VertCoordType": "vtransform" ,       # new vertical S-coordinates
                # time
                "time_counter" : "t",                 # time
                # dimensions
                "xi_rho"       : "x",                 # x dimension at mass points
                "eta_rho"      : "y",                 # y dimension at mass points
                "x_rho"        : "x",                 # x dimension at mass points
                "y_rho"        : "y",                 # x dimension at mass points
                "xi_u"         : "x_u",               # x dimension at u points
                "eta_u"        : "y",                 # y dimension at u points
                "x_u"          : "x_u",               # x dimension at u points
                "y_u"          : "y",                 # y dimension at u points
                "xi_v"         : "x",                 # x dimension at v points
                "eta_v"        : "y_v",               # y dimension at v points
                "x_v"          : "x",                 # x dimension at v points
                "y_v"          : "y_v",               # y dimension at v points
                "x_w"          : "x",                 # x dimension at mass points
                "y_w"          : "y",                 # y dimension at mass points
                "s_rho"        : "s",                 # sigma dimension at rho level
                "s_w"          : "s_w",               # sigma dimension at w level
            }
            
        elif name == "croco_native":
            self.rename_vars = {
                # surface wind stress
                "sustr"        : "xtau_sfc_u",        # x-wind stress component
                "svstr"        : "ytau_sfc_v",        # y-wind stress component
                "wstr"         : "tau_sfc",           # wind stress module
                # surface fluxes
                "radsw"        : "sw_sfc_down",       # surface_downward_shortwave_flux
                "shflx_rlw"    : "lw_sfc_down",       # surface_downward_longwave_flux
                "shflux"       : "hf_net_sfc_down",   # surface_downward_shortwave_flux
                "shflx_lat"    : "lh_sfc_down",       # surface_downward_latent_heat_flux
                "shflx_sen"    : "sh_sfc_down",       # surface_downward_sensible_heat_flux
                "swflux"       : "ep_flux_sfc_down",  # surface_downward_freshwater_flux (E-P)
                # ocean surface
                "zeta"         : "z_sfc",             # sea_surface_height
                "ssh"          : "z_sfc",             # sea_surface_height
                # currents
                "u"            : "xcur",              # x-current component
                "v"            : "ycur",              # y-current component
                "w"            : "zcur",              # z-current component
                "omega"        : "zcur_sc",           # z-current component in sigma coordinate
                "ubar"         : "xcur_btrope",       # x-barotropic current component
                "vbar"         : "ycur_btrope",       # y-barotropic current component 
                # MLD
                "hbls"         : "mld_turb",          # ocean_turbulent_mixed_layer_depth
                "hbl"          : "mld_turb",          # ocean_turbulent_mixed_layer_depth
                # tracers et al.
                "temp"         : "temp",              # ocean potential temperature
                "salt"         : "salt",              # ocean salinity
                # other diag variables
                "Akt"          : "zdiff_tra_w",       # ocean_vertical_tracer_diffusivity_at_w_points
                "Akv"          : "zdiff_mtm_w",       # ocean_vertical_momentum_diffusivity
                # bottom layer
                "bostr"        : "tau_bot",           # ocean bottom stress
                "hbbl"         : "bbl",               # bottom boundary layer
                # grid et al.
                "h"            : "h",                 # terrain height or bathy (positive upward, static)
                "lat_rho"      : "lat",               # latitude at mass points
                "lon_rho"      : "lon",               # longitude at mass points
                "lat_u"        : "lat_u",             # latitude at u points
                "lon_u"        : "lon_u",             # longitude at u points
                "lat_v"        : "lat_v",             # latitude at v points
                "lon_v"        : "lon_v",             # longitude at v points
                "mask_rho"     : "mask",              # land mask
                "f"            : "f",                 # Coriolis parameter (1/seconds) at RHO-points
                "VertCoordType": "vtransform" ,       # new vertical S-coordinates
                # time
                "time"         : "t",                 # time
                # dimensions
                "xi_rho"       : "x",                 # x dimension at mass points
                "eta_rho"      : "y",                 # y dimension at mass points
                "x_rho"        : "x",                 # x dimension at mass points
                "y_rho"        : "y",                 # x dimension at mass points
                "xi_u"         : "x_u",               # x dimension at u points
                "eta_u"        : "y",                 # y dimension at u points
                "x_u"          : "x_u",               # x dimension at u points
                "y_u"          : "y",                 # y dimension at u points
                "xi_v"         : "x",                 # x dimension at v points
                "eta_v"        : "y_v",               # y dimension at v points
                "x_v"          : "x",                 # x dimension at v points
                "y_v"          : "y_v",               # y dimension at v points
                "x_w"          : "x",                 # x dimension at mass points
                "y_w"          : "y",                 # y dimension at mass points
                "s_rho"        : "s",                 # sigma dimension at rho level
                "s_w"          : "s_w",               # sigma dimension at w level
            }
            
        else:
            print('Model name not defined. See module model.py')
            

        self.dims_var = {
            'sc_r' : ['s'],
            'sc_w' : ['s_w'],
        }