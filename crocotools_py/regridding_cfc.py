import crocotools_py.postprocess as post
import numpy as np
import xarray as xr
import os,sys
from datetime import datetime,timedelta
from dask import delayed
import dask.array as da
from scipy.interpolate import griddata
import matplotlib.path as mplPath
from glob import glob
import netCDF4 as nc

def writing_cfc_nc(file_out, Vars, Dims, _epoch, info_dir, regrid_ver):
    """
    This function takes inputs from the outputs of the CROCO (Coastal and Regional Ocean COmmunity) model and saves it as a CF-Comipliant netCDF file.

    INPUT:
    file_out    : Directory and name of the netCDF file to save the outputs. 
    Vars        : Dictionary of the variables to save. 
    Dims        : Dictionary of the dimensions to save.
    _epoch      : Reference datetime in datetime format from when the model was run (i.e. datetime.datetime(YYYY,MM,DD,HH,MM).
    info_dir    : Directory which stores information regarding the model run - params.py file
    regrid_ver  : Regrid version which determines which regridding function was used. Valid for eith: 1, 2 or 3 which corrisponds to regrid1_cf_compliant, regrid2_cf_compliant and regrid3_cf_compliant, respectively.

    OUTPUT:
    Saves a CF-Compliant netCDF4, file_out, which can be used for either research, publication or visualisation. The user is recommended to read the regrid1_cf_compliant,  regrid2_cf_compliant and regrid3_cf_compliant functions to know exactly what they are looking for. 
    """

    sys.path.append(info_dir)
    import croco_cf_compliance as complaince

    print('')
    print('Creating: ', file_out)
    print('')

    # read variables
    zeta  = Vars['zeta']
    temp  = Vars['temp']
    salt  = Vars['salt']
    u     = Vars['u']
    v     = Vars['v']      

    # read dimensions
    time    = Dims['time']
    lon_rho = Dims['lon_rho']
    lat_rho = Dims['lat_rho']

    if (np.ndim(lon_rho) == 1) & (np.ndim(lat_rho) == 1):
        lat_min, lat_max = float(np.min(lat_rho[:].data)), float(np.max(lat_rho[:].data))
        lon_min, lon_max = float(np.min(lon_rho[:].data)), float(np.max(lon_rho[:].data))
    elif (np.ndim(lon_rho) == 2) & (np.ndim(lat_rho) == 2):
        lat_min, lat_max = float(np.min(lat_rho[:,0].data)), float(np.max(lat_rho[:,0].data))
        lon_min, lon_max = float(np.min(lon_rho[0,:].data)), float(np.max(lon_rho[0,:].data))
    else:
        print('Error: incorrect number of dimensions.')
        sys.exit()

    
    if regrid_ver == 1:
        eta = Dims['eta_rho']
        xi = Dims['xi_rho']
        s_rho = Dims['s_rho']
        depth = Vars['depth']
        h = Vars['h']
    elif regrid_ver == 2:
        eta = Dims['eta_rho']
        xi = Dims['xi_rho']
        depth = Dims['depth']
        h = Vars['h']
    elif regrid_ver == 3:
        depth = Dims['depth']
    else:
        print("Error: Regrid version does not exist. Available versions are 1, 2 & 3.")
        sys.exit
        
    # Read dimension sizes
    Nt,Nz,Ny,Nx = np.shape(temp)

    # Write netCDF file 
    ncf = nc.Dataset(file_out, "w", format="NETCDF4",clobber=True)

    ncf.title = complaince.title
    ncf.institution = complaince.institution
    ncf.source = complaince.source
    creation_time = datetime.now()
    ncf.history = "Created " + str(datetime.strftime(creation_time,"%Y-%m-%d %H:%M:%S"))
    ncf.references = complaince.references
    ncf.comment = complaince.comment
    ncf.Conventions = complaince.conventions
    ncf.summary = complaince.summary   
    ncf.keywords = complaince.keywords
    ncf.license = complaince.license
    ncf.license_url = complaince.license_url
    ncf.date_created = complaince.date_created
    ncf.creator_name = complaince.creator_name
    ncf.creator_email = complaince.creator_email
    ncf.project = complaince.project
    ncf.publisher_name = complaince.publisher_name
    ncf.publisher_email = complaince.publisher_email
    ncf.geospatial_lat_min = lat_min
    ncf.geospatial_lat_max = lat_max
    ncf.geospatial_lat_units = complaince.geospatial_lat_units
    ncf.geospatial_lon_min = lon_min
    ncf.geospatial_lon_max = lon_max
    ncf.geospatial_lon_units = complaince.geospatial_lon_units
    ncf.geospatial_vertical_min = float(np.min(depth.data))
    ncf.geospatial_vertical_max = float(np.max(depth.data))
    ncf.geospatial_vertical_units = complaince.geospatial_vertical_units
    ncf.geospatial_vertical_positive = complaince.geospatial_vertical_positive
    ncf.time_coverage_start = complaince.time_coverage_start
    ncf.time_coverage_end = complaince.time_coverage_end
    ncf.date_modified = complaince.date_modified
    ncf.date_metadata_modified = complaince.date_metadata_modified
    ncf.instrument = complaince.instrument
    ncf.lineage = complaince.lineage
    
    # Write dimensions
    dim_time = ncf.createDimension("time", Nt)
    
    if regrid_ver == 1:
        dim_srho = ncf.createDimension("sigma_level", Nz)
    elif (regrid_ver == 2) | (regrid_ver == 3):
        dim_zlev = ncf.createDimension("z", Nz)
    else:
        print("Error: Regrid version does not exist. Available versions are 1, 2 & 3.")
        sys.exit
        
    dim_eta = ncf.createDimension("y", Ny) 
    dim_xi  = ncf.createDimension("x", Nx)

    # write variables
    nc_time  = ncf.createVariable("time",'f8',("time",))
    nc_time.standard_name = "time"
    nc_time.long_name = "Time"
    nc_time.units = "seconds since " + str(_epoch)
    nc_time.calendar = "standard"
    nc_time.valid_min = time[0] 
    nc_time.valid_max = time[-1] 
    nc_time[:] = time[:]

    if regrid_ver == 1:
        nc_srho  = ncf.createVariable("ocean_sigma_coordinate",float,("sigma_level"))
        nc_srho.standard_name = "ocean_sigma_coordinate"
        nc_srho.long_name = "Ocean Sigma Coordinate"
        nc_srho.units = "1"
        nc_srho.valid_min = float(np.min(s_rho.data))
        nc_srho.valid_max = float(np.max(s_rho.data))
        nc_srho.axis = "Z"
        nc_srho[:] = s_rho.astype("float")
    elif (regrid_ver == 2) | (regrid_ver == 3):
        nc_zlev  = ncf.createVariable("depth",float,("z"))
        nc_zlev.standard_name = "depth"
        nc_zlev.long_name = "Water Depth"
        nc_zlev.units = "m"
        nc_zlev.valid_min = float(np.min(depth.data))
        nc_zlev.valid_max = float(np.max(depth.data))
        nc_zlev.axis = "Z"
    else:
        print("Error: Regrid version does not exist. Available versions are 1, 2 & 3.")
        sys.exit

    if (regrid_ver == 1) | (regrid_ver == 2):
        nc_eta = ncf.createVariable("y_index",float,("y",))
        nc_eta.standard_name = "y_index"
        nc_eta.long_name = "Y index"
        nc_eta.units = "1"
        nc_eta.valid_min = float(np.min(eta.data))
        nc_eta.valid_max = float(np.max(eta.data))
        nc_eta.axis	= "Y"
        nc_eta[:] = eta.astype("float")
    
        nc_xi = ncf.createVariable("x_index",float,("x",))
        nc_xi.standard_name = "x_index"
        nc_xi.long_name = "X index"
        nc_xi.units = "1"
        nc_xi.valid_min = float(np.min(xi.data))
        nc_xi.valid_max = float(np.max(xi.data))
        nc_xi.axis	= "X"
        nc_xi[:] = xi.astype("float")

    else:
        pass

    if (np.ndim(lon_rho) == 1) & (np.ndim(lat_rho) == 1):
        nc_lon = ncf.createVariable("longitude",float,("x",))
        nc_lon.standard_name = "longitude"
        nc_lon.long_name = "Longitude"
        nc_lon.units = "degrees"
        nc_lon.valid_min = float(np.nanmin(lon_rho[:].data))
        nc_lon.valid_max = float(np.nanmax(lon_rho[:].data))
        nc_lon.coordinates	= "X"
        nc_lon[:] = lon_rho.astype("float")
        
        nc_lat = ncf.createVariable("latitude",float,("y",))
        nc_lat.standard_name = "latitude"
        nc_lat.long_name = "Latitude"
        nc_lat.units = "degrees"
        nc_lat.valid_min = float(np.nanmin(lat_rho[:].data))
        nc_lat.valid_max = float(np.nanmax(lat_rho[:].data))
        nc_lat.coordinates	= "Y"
        nc_lat[:] = lat_rho.astype("float")
        
    elif (np.ndim(lon_rho) == 2) | (np.ndim(lat_rho) == 2):
        nc_lon = ncf.createVariable("longitude",float,("y","x",))
        nc_lon.standard_name = "longitude"
        nc_lon.long_name = "Longitude"
        nc_lon.units = "degrees"
        nc_lon.valid_min = float(np.nanmin(lon_rho[0,:].data))
        nc_lon.valid_max = float(np.nanmax(lon_rho[0,:].data))
        nc_lon.coordinates	= "Y X"
        nc_lon[:] = lon_rho.astype("float")
        
        nc_lat = ncf.createVariable("latitude",float,("y","x",))
        nc_lat.standard_name = "latitude"
        nc_lat.long_name = "Latitude"
        nc_lat.units = "degrees"
        nc_lat.valid_min = float(np.nanmin(lat_rho[:,0].data))
        nc_lat.valid_max = float(np.nanmax(lat_rho[:,0].data))
        nc_lat.coordinates	= "Y X"
        nc_lat[:] = lat_rho.astype("float")
        
        nc_h  = ncf.createVariable("sea_floor_depth",float,("y","x",))
        nc_h.standard_name = "sea_floor_depth"
        nc_h.long_name = "Depth of the sea floor"
        nc_h.units = "m"
        nc_h.valid_min = float(np.nanmin(h.data)) #(RECOMMENDED)	
        nc_h.valid_max = float(np.nanmax(h.data)) #(RECOMMENDED)	
        nc_h.coordinates = "Y X" #(RECOMMENDED)
        nc_h[:] = h.astype("float")
        
    else:
        print('Error: Incorrect number of dimensions')
        sys.exit()

    nc_zeta  = ncf.createVariable("sea_surface_elevation",float,("time","y","x",))
    nc_zeta.standard_name = "sea_surface_elevation"
    nc_zeta.long_name = "Sea Surface Elevation"
    nc_zeta.units = "m"
    nc_zeta.valid_min = float(np.nanmin(zeta.data)) #(RECOMMENDED)	
    nc_zeta.valid_max = float(np.nanmax(zeta.data)) #(RECOMMENDED)	
    nc_zeta.coordinates = "time y x" #(RECOMMENDED)
    nc_zeta[:] = zeta.astype("float")

    if regrid_ver == 1:
        nc_temp  = ncf.createVariable("sea_water_temperature",float,("time","sigma_level","y","x",))
        nc_temp.coordinates = "time sigma_level y x" 
        
    elif (regrid_ver == 2) | (regrid_ver == 3):
        nc_temp  = ncf.createVariable("sea_water_temperature",float,("time","z","y","x",))
        nc_temp.coordinates = "time z y x"
        
    else:
        print("Error: Regrid version does not exist. Available versions are 1, 2 & 3.")
        sys.exit
        
    nc_temp.standard_name = "sea_water_temperature"
    nc_temp.long_name = "Sea Water Temperature"
    nc_temp.units = "degC"
    nc_temp.valid_min = float(np.nanmin(temp)) #(RECOMMENDED)	
    nc_temp.valid_max = float(np.nanmax(temp)) #(RECOMMENDED)	
    nc_temp[:] = temp.astype("float")

    if regrid_ver == 1:
        nc_salt   = ncf.createVariable("sea_water_salinity",float,("time","sigma_level","y","x",))
        nc_salt.coordinates = "time sigma_level y x"
        
    elif (regrid_ver == 2) | (regrid_ver == 3):
        nc_salt   = ncf.createVariable("sea_water_salinity",float,("time","z","y","x",))
        nc_salt.coordinates = "time z y x"
        
    else:
        print("Error: Regrid version does not exist. Available versions are 1, 2 & 3.")
        sys.exit

    nc_salt.standard_name = "sea_water_salinity"
    nc_salt.long_name = "Sea Water Salinity"
    nc_salt.units = "1e-3"
    nc_salt.valid_min = float(np.nanmin(salt)) #(RECOMMENDED)	
    nc_salt.valid_max = float(np.nanmax(salt)) #(RECOMMENDED)	
    nc_salt[:] = salt.astype("float")

    if regrid_ver == 1:
        nc_u   = ncf.createVariable("eastward_sea_water_velocity",float,("time","sigma_level","y","x",))
        nc_u.coordinates = "time sigma_level y x"
    elif (regrid_ver == 2) | (regrid_ver == 3):
        nc_u   = ncf.createVariable("eastward_sea_water_velocity",float,("time","z","y","x",))
        nc_u.coordinates = "time z y x"
    else:
        print("Error: Regrid version does not exist. Available versions are 1, 2 & 3.")
        sys.exit
    nc_u.standard_name = "eastward_sea_water_velocity"
    nc_u.long_name = "Eastward Sea Water Velocity"
    nc_u.units = 'm s-1'
    nc_u.valid_min = float(np.nanmin(u)) #(RECOMMENDED)	
    nc_u.valid_max = float(np.nanmax(u)) #(RECOMMENDED)	
    nc_u[:] = u.astype("float")

    if regrid_ver == 1:
        nc_v = ncf.createVariable("northward_sea_water_velocity",float,("time","sigma_level","y","x",))
        nc_v.coordinates = "time sigma_level y x"
    elif (regrid_ver == 2) | (regrid_ver == 3):
        nc_v = ncf.createVariable("northward_sea_water_velocity",float,("time","z","y","x",))
        nc_v.coordinates = "time z y x"
    else:
        print("Error: Regrid version does not exist. Available versions are 1, 2 & 3.")
        sys.exitz
    nc_v.standard_name = "northward_sea_water_velocity"
    nc_v.long_name = "Northward Sea Water Velocity"
    nc_v.units = 'm s-1'
    nc_v.valid_min = float(np.nanmin(v)) #(RECOMMENDED)	
    nc_v.valid_max = float(np.nanmax(v)) #(RECOMMENDED)	
    nc_v[:] = v.astype("float")

    if regrid_ver == 1:
        nc_depth = ncf.createVariable("depth",float,("time","sigma_level","y","x",))
        nc_depth.standard_name = "depth"
        nc_depth.long_name = "Water Depth"
        nc_depth.units = "m"
        nc_depth.valid_min = float(np.nanmin(depth))
        nc_depth.valid_max = float(np.nanmax(depth))
        nc_depth.coordinates = "time sigma_level y x" #(RECOMMENDED)
        nc_depth[:] = depth.astype("float")
    else:
        pass
    
    ncf.close()   

def regrid1_cf_compliant(fname_in,info_dir,out_dir=None):
    """
    Function to regrid the croco variable (sea-surface height, temperature, salinity, u, v) onto the rho points of the grid cell. 
    This function also computes the depths of the sigma layers.
    Finally, it saves the variables in a netCDF file that is CF-Complient. 

    INPUTS:
    fname_in : Path to the croco file/s that is of interest to be regrided. Can be a specific file, or a list of files or one can use a 
    wildcard (*) to include all the files in a specific directory.  
    out_dir  : Path to the directory to save outputs.
    info_dir : Path to croco directory with model info (i.e. /somisana-croco/configs/sa_west_02/croco_v1.3.1/MERCATOR/)
    
    OUTPUTS:
    Saves a netCDF file in the out_dir that is CF-compliant.
    """
    sys.path.append(info_dir)
    import crocotools_param as param
    Yorig,Morig,Dorig = param.Yorig,param.Morig,param.Dorig
    _epoch = datetime( int(Yorig), int(Morig), int(Dorig) )
    
    if type(fname_in) == str:
        if fname_in.find('*') > 0:
            fname_in = glob(fname_in)
        elif fname_in.contains(','):
            fname_list = [x for x in fname_in.split(',')]
            fname_in = fname_list
        else:
            fname_in = [fname_in]
    else:
        print("fname_in: " , fname_in)
        print('Error: unkown input format. Input variable fname_in needs to be str or list.')
        sys.exit()
    
    for file in fname_in:
        print('')
        print('Opening: ', file)
        print('')
        
        print("Extracting the model output variables we need")


        # using the get_var() function so masking is included
        h = post.get_var(file, "h")
        temp = post.get_var(file, "temp", ref_date=_epoch)
        salt = post.get_var(file, "salt", ref_date=_epoch)
        zeta = post.get_var(file, "zeta", ref_date=_epoch)
        u, v = post.get_uv(file, ref_date=_epoch)

        # get the depth levels of the sigma layers
        print("Computing depth of sigma levels")
        depth = post.get_depths(post.get_ds(file)) # get_depths now takes an xarray dataset as input

        ds = xr.open_dataset(file)
        time_steps = ds.time.values
        ds.close()

        # Create a dictionary with dimensions
        Dims = {

            "time": time_steps.astype('timedelta64[s]').astype(int),
            "eta_rho": temp.eta_rho.values,
            "xi_rho": temp.xi_rho.values,
            "lon_rho": temp.lon_rho.values,
            "lat_rho": temp.lat_rho.values,
            "s_rho": temp.s_rho.values

        }

        Vars = {

            'h':h.values,
            'zeta':zeta.values,
            'temp':temp.values,
            'salt':salt.values,
            'u':u.values,
            'v':v.values,
            'depth':depth,

        }

        # Ensure the directory for the specified output exists
        if out_dir is None:
            out_dir = os.path.join(os.path.dirname(file), 'regrid_tier_1/',)
            os.makedirs(os.path.dirname(out_dir), exist_ok=True)
        else:
            os.makedirs(os.path.dirname(out_dir), exist_ok=True)

        fname_out = os.path.abspath(os.path.join(os.path.dirname(out_dir), os.path.basename(file[:-3]) + '_t1.nc'))

        if os.path.exists(fname_out):
            os.remove(fname_out)
        else:
            pass

        writing_cfc_nc(fname_out, Vars, Dims, _epoch, info_dir, regrid_ver=1)

    print('Done')
    print('')

def regrid2_cf_compliant(fname_in,info_dir,depths=None,out_dir=None):
    print(' depths : ', depths)
    """
    Function that takes the output of regrid-tier1 as input and
    regrids the sigma levels to constant z levels, including the surface and bottom layers
    output variables are the same as tier 1, only depths is now a dimension with the user specified values.

    INPUTS:
    fname_in : Path to the croco file/s that is of interest to be regrided. Can be a specific file, or a list of files or one can use a
    wildcard (*) to include all the files in a specific directory.
    out_dir  : Path to the directory to save outputs.
    info_dir : Path to croco directory with model info (i.e. /somisana-croco/configs/sa_west_02/croco_v1.3.1/MERCATOR/)
    depths : Depths of the horizontal slices.

    OUTPUTS:
    Saves a netCDF file in the out_dir that is CF-compliant.
    """
    sys.path.append(info_dir)
    import crocotools_param as param
    Yorig,Morig,Dorig = param.Yorig,param.Morig,param.Dorig
    _epoch = datetime(int(Yorig),int(Morig),int(Dorig))

    if depths == None:
        depths = [0,-5,-10,-20,-50,-75,-100,-200,-500,-1000,-99999]
    else:
        pass

    if type(fname_in) == str:
        if fname_in.find('*') < 0:
            fname_in = [fname_in]
        elif fname_in.find('*') > 0:
            fname_in = glob(fname_in)
        else:
            print('Error: unkown input format. Input variable fname_in needs to be str or list.')
            sys.exit()
    elif type(fname_in) == list:
        fname_in=fname_in
    else:
        print('Error: unkown input format. Input variable fname_in needs to be str or list.')
        sys.exit()

    print("Extracting the tier 1 re-gridded model output")

    for file in fname_in:
        print('')
        print('Opening: ', file)
        print('')

        ds = xr.open_dataset(file)

        depth_in=ds.depth.values
        temp_in=ds.sea_water_temperature.values
        salt_in=ds.sea_water_salinity.values
        u_in=ds.eastward_sea_water_velocity.values
        v_in=ds.northward_sea_water_velocity.values
        zeta_in=ds.sea_surface_elevation.values
        h_in=ds.sea_floor_depth.values
        lon_rho = ds.longitude.values
        lat_rho = ds.latitude.values

        # Dimensions
        s_rho=ds.sigma_level.values
        xi_rho=ds.x.values
        eta_rho=ds.y.values
        time_in=ds.time.values
        time_out = [(ts.astype('datetime64[s]').astype(int) - np.datetime64(_epoch).astype('datetime64[s]').astype(int)) for ts in time_in]
        ds.close()

        T,N,M,L=np.shape(depth_in)

        # set up the output arrays
        depth_out = np.array(depths).astype(float)
        temp_out=np.zeros((T,len(depth_out),M,L))
        salt_out=temp_out.copy()
        u_out=temp_out.copy()
        v_out=temp_out.copy()

        print("Doing the vertical interpolation to the constant depth levels")

        for d, depth in enumerate(depth_out):
            if depth==0: # surface layer
                temp_out[:,d,:,:]=temp_in[:,N-1,:,:]
                salt_out[:,d,:,:]=salt_in[:,N-1,:,:]
                u_out[:,d,:,:]=u_in[:,N-1,:,:]
                v_out[:,d,:,:]=v_in[:,N-1,:,:]
            elif depth==-99999: # bottom layer
                temp_out[:,d,:,:]=temp_in[:,0,:,:]
                salt_out[:,d,:,:]=salt_in[:,0,:,:]
                u_out[:,d,:,:]=u_in[:,0,:,:]
                v_out[:,d,:,:]=v_in[:,0,:,:]
            else:
                print("Depth = ", depth, " m")
                for t in np.arange(T):
                    temp_out[t,d,:,:]=post.hlev(temp_in[t,::], depth_in[t,::], depth)
                    salt_out[t,d,:,:]=post.hlev(salt_in[t,::], depth_in[t,::], depth)
                    u_out[t,d,:,:]=post.hlev(u_in[t,::], depth_in[t,::], depth)
                    v_out[t,d,:,:]=post.hlev(v_in[t,::], depth_in[t,::], depth)
        
        # Create a dictionary with dimensions
        Dims = {

            "time" : time_out,
            "eta_rho" : eta_rho,
            "xi_rho" : xi_rho,
            "lon_rho" : lon_rho,
            "lat_rho" : lat_rho,
            "depth" : depth,  

        }
        
        Vars = {
            
            "h" : h_in,
            "zeta" : zeta_in,
            "temp" : temp_out,
            "salt" : salt_out,
            "u" : u_out,
            "v" : v_out,
        }

        # Ensure the directory for the specified output exists
        if out_dir is None:
            out_dir = os.path.join(os.path.dirname(file), '..', 'regrid_tier_2/')
            os.makedirs(os.path.dirname(out_dir), exist_ok=True)
        else:
            os.makedirs(os.path.dirname(out_dir), exist_ok=True)
        
        fname_out = os.path.abspath(os.path.join(os.path.dirname(out_dir), os.path.basename(file[:-6])+'_t2.nc'))
        
        if os.path.exists(fname_out):
            os.remove(fname_out)
        else:
            pass
        
        writing_cfc_nc(fname_out, Vars, Dims, _epoch, info_dir,regrid_ver=2)

    print('Done')
    print('')

def regrid3_cf_compliant(fname_in,info_dir,spacing=None,out_dir=None):
    """
    Function that takes the output of regrid-tier2 as input and
    regrids the curvilinear grid to a rectilinear grid for visulaisation.
    output variables are sea-surface height, temperature, salinity, u and v.

    INPUTS:
    fname_in : Path to the croco file/s that is of interest to be regrided. Can be a specific file, or a list of files or one can use a 
    wildcard (*) to include all the files in a specific directory.  
    out_dir  : Path to the directory to save outputs.
    info_dir : Path to croco directory with model info (i.e. /somisana-croco/configs/sa_west_02/croco_v1.3.1/MERCATOR/)
    spacing  : spacing of the rectilinear grid (float).
    
    OUTPUTS:
    Saves a netCDF file in the out_dir that is CF-compliant.
    """

    sys.path.append(info_dir)
    import crocotools_param as param
    Yorig,Morig,Dorig = param.Yorig,param.Morig,param.Dorig
    _epoch = datetime(int(Yorig),int(Morig),int(Dorig))
    
    if spacing == None: 
        spacing = 0.01
    else:
        pass

    if type(fname_in) == str:
        if fname_in.find('*') < 0:
            fname_in = [fname_in]
        elif fname_in.find('*') > 0:
            fname_in = glob(fname_in)
        else:
            print('Error: unkown input format. Input variable fname_in needs to be str or list.')
            sys.exit()
    elif type(fname_in) == list:
        fname_in=fname_in
    else:
        print('Error: unkown input format. Input variable fname_in needs to be str or list.')
        sys.exit()

    for file in fname_in:
        print('')
        print('Opening: ', file)
        print('')

        ds = xr.open_dataset(file,decode_times=False)
        zeta_in = ds.sea_surface_elevation.values
        temp_in = ds.sea_water_temperature.values
        salt_in = ds.sea_water_salinity.values
        u_in = ds.eastward_sea_water_velocity.values
        v_in = ds.northward_sea_water_velocity.values
        time_in = ds.time.values.astype(int)
        depth = ds.depth.values
        lon_rho = ds.longitude.values
        lat_rho = ds.latitude.values
        Nt, Nz, Ny, Nx = np.shape(temp_in)
        lon_rho_1d = np.ravel(lon_rho)
        lat_rho_1d = np.ravel(lat_rho)

        # input for griddata function later
        lonlat_input = np.array([lon_rho_1d, lat_rho_1d]).T

        print("Generating the regular horizontal output grid")
        # get the model boundary polygon
        lon_boundary = np.hstack(
            (lon_rho[0:, 0], lon_rho[-1, 1:-1], lon_rho[-1::-1, -1], lon_rho[0, -2::-1])
        )
        lat_boundary = np.hstack(
            (lat_rho[0:, 0], lat_rho[-1, 1:-1], lat_rho[-1::-1, -1], lat_rho[0, -2::-1])
        )
    
        # find the corners of the output regular grid (just big enough to cover the model domain)
        lon_min = np.floor(np.min(lon_boundary) / spacing) * spacing
        lon_max = np.ceil(np.max(lon_boundary) / spacing) * spacing
        lat_min = np.floor(np.min(lat_boundary) / spacing) * spacing
        lat_max = np.ceil(np.max(lat_boundary) / spacing) * spacing
    
        # generate the regular grid
        Nlon = int(np.rint((lon_max - lon_min) / spacing)) + 1
        Nlat = int(np.rint((lat_max - lat_min) / spacing)) + 1
        lon_out = np.linspace(lon_min, lon_max, num=Nlon, endpoint=True)
        lat_out = np.linspace(lat_min, lat_max, num=Nlat, endpoint=True)
        lon_out_grd, lat_out_grd = np.meshgrid(lon_out, lat_out)
    
        # get a mask for the output grid which tells us which points are inside the CROCO model grid
        poly_boundary = mplPath.Path(np.array([lon_boundary, lat_boundary]).T)
        mask_out = np.zeros_like(lon_out_grd)
        for y in np.arange(Nlat):
            for x in np.arange(Nlon):
                if poly_boundary.contains_point((lon_out_grd[y, x], lat_out_grd[y, x])):
                    mask_out[y, x] = 1
                else:
                    mask_out[y, x] = np.nan

        print('')
        print("Interpolating the model output onto the regular horizontal output grid")
    
        @delayed
        def compute_2d_chunk(t, variable, method="nearest"):
            return (
                griddata(
                    lonlat_input,
                    np.ravel(variable[t, ::]),
                    (lon_out_grd, lat_out_grd),
                    method,
                )
                * mask_out
            )
    
        @delayed
        def compute_3d_chunk(t, variable, n, method="nearest"):
            return (
                griddata(
                    lonlat_input,
                    np.ravel(variable[t, n, ::]),
                    (lon_out_grd, lat_out_grd),
                    method,
                )
                * mask_out
            )
    
        # Separate lists for each time step
        zeta_out = []
        temp_out_time = []
        salt_out_time = []
        u_out_time = []
        v_out_time = []
        for t in np.arange(Nt):
            # Lists for each depth level
            temp_out_depth = []
            salt_out_depth = []
            u_out_depth = []
            v_out_depth = []
    
            zeta_out.append(
                da.from_delayed(
                    compute_2d_chunk(t, zeta_in),
                    shape=(Nlat, Nlon),
                    dtype=float,
                )
            )
            for n in np.arange(Nz):
                #print(
                #    f"Timestep {str(t+1).zfill(3)}/{Nt}. Depth {str(np.round(depth[n]).astype(int)).zfill(5)}m (lvl {str(n+1).zfill(2)}/{Nz})."
                #)
                temp_out_depth.append(
                    da.from_delayed(
                        compute_3d_chunk(t, temp_in, n),
                        shape=(Nlat, Nlon),
                        dtype=float,
                    )
                )
                salt_out_depth.append(
                    da.from_delayed(
                        compute_3d_chunk(t, salt_in, n),
                        shape=(Nlat, Nlon),
                        dtype=float,
                    )
                )
                u_out_depth.append(
                    da.from_delayed(
                        compute_3d_chunk(t, u_in, n),
                        shape=(Nlat, Nlon),
                        dtype=float,
                    )
                )
                v_out_depth.append(
                    da.from_delayed(
                        compute_3d_chunk(t, v_in, n),
                        shape=(Nlat, Nlon),
                        dtype=float,
                    )
                )
    
    
            # Stack the depth dimension and append to the time list
            temp_out_time.append(da.stack(temp_out_depth, axis=0))
            salt_out_time.append(da.stack(salt_out_depth, axis=0))
            u_out_time.append(da.stack(u_out_depth, axis=0))
            v_out_time.append(da.stack(v_out_depth, axis=0))
    
        # Stack the time dimension
        zeta_out = da.stack(zeta_out, axis=0)
        temp_out = da.stack(temp_out_time, axis=0)
        salt_out = da.stack(salt_out_time, axis=0)
        u_out = da.stack(u_out_time, axis=0)
        v_out = da.stack(v_out_time, axis=0)
    
        # Create a dictionary with dimensions
        Dims = {
            
            "time": time_in,
            "lon_rho": lon_out,
            "lat_rho": lat_out,
            "depth": ds.depth.values
            
        }
        
        Vars = {
        
            'zeta':np.array(zeta_out),
            'temp':np.array(temp_out),
            'salt':np.array(salt_out),
            'u':np.array(u_out),
            'v':np.array(v_out)
            
        }
    
        # Ensure the directory for the specified output exists
        if out_dir is None:
            out_dir = os.path.join(os.path.dirname(file), '..', 'regrid_tier_3/')
            os.makedirs(os.path.dirname(out_dir), exist_ok=True)
        else:
            os.makedirs(os.path.dirname(out_dir), exist_ok=True)
        
        fname_out = os.path.abspath(os.path.join(os.path.dirname(out_dir), os.path.basename(file[:-6])+'_t3.nc'))

        if os.path.exists(fname_out):
            os.remove(fname_out)
        else:
            pass
        
        writing_cfc_nc(fname_out, Vars, Dims, _epoch, info_dir, regrid_ver=3)

    print('Done')
    print('')

if __name__ == "__main__":
    #fname_in = '/home/g.rautenbach/Data/models/sa_southeast/croco_avg_Y2009M09.nc'
    fname_in = '/home/g.rautenbach/Data/models/sa_southeast/regrid_tier_1/croco_avg_Y2009M09_t1.nc'
    info_dir = '/home/g.rautenbach/Projects/somisana-croco/configs/sa_west_02/croco_v1.3.1/MERCATOR/'
    depths = [0,-5,-10]
    spacing='0.01'
    #regrid1_cf_compliant(fname_in,info_dir,out_dir=None)
    #regrid2_cf_compliant(fname_in,info_dir,depths,out_dir=None)
    regrid3_cf_compliant(fname_in,info_dir,spacing,out_dir=None)
